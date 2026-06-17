library(shiny)
library(shinythemes)
library(tidymodels)
library(themis)
library(protr)
library(xgboost)
library(ranger)
library(kernlab)
library(kknn)
library(tidyverse)
library(Biostrings)
library(Matrix)
library(svglite)
library(DT)
library(Peptides)
library(workflowsets)
library(plotly)
library(htmlwidgets)

# --- 1. CONFIGURACIÓN ---
options(shiny.maxRequestSize = 800 * 1024^2)
###models_bundle <- tryCatch({ readRDS("models/Ensemble_4m.rds") }, error = function(e) NULL)
BLAST_DB <- "C:/Users/jllavin/Desktop/SelCys_data/BLASTdb/uniprot_db" 
BLAST_PATH <- "C:/ncbi-blast-2.17.0+/bin/blastp.exe"

# --- 2. FUNCIONES DE EXTRACCIÓN Y BLAST ---
clean_id_func <- function(id_str) {
  id_str <- as.character(id_str)
  if (grepl("^sp\\|", id_str)) return(strsplit(id_str, "\\|")[[1]][2])
  return(strsplit(id_str, "[ ;|]")[[1]][1])
}

safe_extractAAC <- function(seq) {
  tryCatch({
    res <- extractAAC(seq)
    if (is.null(res) || length(res) != 20 || any(is.na(res))) return(NULL)
    return(res)
  }, error = function(e) NULL)
}

safe_extractCTDC <- function(seq) {
  tryCatch({
    res <- extractCTDC(seq)
    if (is.null(res) || any(is.na(res))) return(NULL)
    return(res)
  }, error = function(e) NULL)
}

extract_features_logic <- function(seqs) {
  if (length(seqs) == 0) return(NULL)
  seqs <- toupper(gsub("[ \n\r\t]", "", seqs))
  seqs <- gsub("U", "C", seqs)
  results_list <- vector("list", length(seqs))
  for (i in seq_along(seqs)) {
    aac <- safe_extractAAC(seqs[i]); ctdc <- safe_extractCTDC(seqs[i])
    if (!is.null(aac) && !is.null(ctdc)) {
      aac_df <- as.data.frame(t(aac)); names(aac_df) <- paste0("AAC_", names(aac_df))
      ctdc_df <- as.data.frame(t(ctdc)); names(ctdc_df) <- paste0("GLC_", names(ctdc_df))
      results_list[[i]] <- cbind(aac_df, ctdc_df, Sequence = seqs[i])
    }
  }
  results_list <- results_list[!sapply(results_list, is.null)]
  if (length(results_list) == 0) return(NULL)
  bind_rows(results_list)
}

run_blast_evidence <- function(fasta_path, evalue_cut) {
  evalue_str <- paste0("1e", evalue_cut)
  cmd <- sprintf('"%s" -query "%s" -db "%s" -outfmt "6 qseqid stitle pident" -max_target_seqs 1 -evalue %s -num_threads 4', 
                 BLAST_PATH, fasta_path, BLAST_DB, evalue_str)
  res <- tryCatch(system(cmd, intern = TRUE), error = function(e) NULL)
  if(is.null(res) || length(res) == 0) return(NULL)
  read.table(text = res, sep = "\t", quote = "", fill = TRUE, 
             col.names = c("qseqid", "Description", "pident"))
}

# --- 3. UI COMPLETA ---
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("SelCys: SelenoCysteine proteins detector via ML + BLAST"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta_in", "Upload FASTA", accept = c(".fasta", ".fa")),
      # NUEVO: Selector de modelo
      selectInput("model_type", "Choose model:",
                  choices = c("General (Eukaryotes/Prokaryotes)" = "General", 
                              "Bacteria Only" = "Bacterial")),
      
      checkboxInput("enable_blast", "Enable BLAST analysis (BLAST)", value = TRUE),
      sliderInput("evalue_cut", "E-value cut-off (log10):", min = -10, max = -1, value = -5, step = 1, pre = "1e"),
      helpText("Notice: BLAST analysis can take several minutes."),
      actionButton("predict_btn", "Run Analysis", class = "btn-primary btn-block"),
      hr(),
      downloadButton("download_tsv", "Download Results (TSV)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Prediction Table", DTOutput("results_table")),
        tabPanel("3D Display Space", 
                 plotlyOutput("pca_3d_plot", height = "600px"),
                 hr(),
                 downloadButton("down_pca_3d_html", "Export 3D PCA (HTML)"))
      )
    )
  )
)

# --- 4. SERVER COMPLETO ---
server <- function(input, output) {
  
  current_model <- reactive({
    filename <- if(input$model_type == "General") "SelCysGeneral.rds" else "SelCysBacterial.rds"
    path <- file.path("models", filename)
    
    if(file.exists(path)) {
      return(readRDS(path))
    } else {
      showNotification(paste("Error: No se encuentra el archivo", filename), type = "error")
      return(NULL)
    }
  })
  
  pca_plot_obj <- reactiveVal(NULL)
  
  results_data <- eventReactive(input$predict_btn, {
    req(input$fasta_in, current_model()) 
    
    raw <- readAAStringSet(input$fasta_in$datapath)
    ids_limpios <- sapply(names(raw), clean_id_func)
    names(raw) <- ids_limpios
    
    df_features <- extract_features_logic(as.character(raw))
    if(is.null(df_features)) stop("Error: No se pudieron extraer descriptores.")
    
    # --- PASO CRÍTICO: PREPARACIÓN PARA EL MODELO ---
    # Necesitamos un dataframe que incluya las columnas que el modelo "exige" (ID, Sequence)
    # pero asegurando que los tipos sean compatibles con el modelo.
    df_para_modelo <- df_features
    df_para_modelo$ID <- ids_limpios[1:nrow(df_features)]
    # La columna Sequence ya está en df_features gracias a extract_features_logic
    
    # Predicción: Pasamos el dataframe completo. 
    # Tidymodels usará la receta interna del archivo .rds para ignorar ID y Sequence 
    # durante el cálculo (siempre que tengan roles de "id" en el workflow).
    preds_list <- map(current_model(), ~predict(.x, df_para_modelo, type = "prob"))
    prob_selcys <- rowMeans(do.call(cbind, map(preds_list, ~ .x$.pred_SelCys)))
    
    # Construcción final uniendo todo para la tabla
    final_df <- data.frame(
      ID = df_para_modelo$ID,
      Veredict = ifelse(prob_selcys > 0.5, "Selenoprotein", "Standard"),
      Description = "No evaluation",
      pident = NA,
      Confidence = round(pmax(prob_selcys, 1-prob_selcys) * 100, 2),
      Sequence = df_para_modelo$Sequence,
      stringsAsFactors = FALSE
    )
    
    # BLAST y MERGE directo
    # 4. Análisis BLAST con merge seguro
    if(input$enable_blast) {
      positives <- final_df %>% filter(Veredict == "Selenoprotein")
      if(nrow(positives) > 0) {
        temp_fasta <- tempfile(fileext = ".fasta")
        writeXStringSet(AAStringSet(setNames(positives$Sequence, positives$ID)), temp_fasta)
        blast_res <- run_blast_evidence(temp_fasta, input$evalue_cut)
        
        if(!is.null(blast_res) && nrow(blast_res) > 0) {
          # 1. Nos aseguramos de tener solo un hit por ID
          blast_res_unique <- blast_res %>%
            group_by(qseqid) %>%
            slice(1) %>% 
            ungroup() %>%
            rename(blast_pident = pident, blast_desc = Description) # Renombramos para evitar conflictos
          
          # 2. Join limpio usando los nombres renombrados
          final_df <- final_df %>% 
            left_join(blast_res_unique, by = c("ID" = "qseqid")) %>%
            mutate(
              pident = ifelse(!is.na(blast_pident), blast_pident, pident),
              Description = ifelse(!is.na(blast_desc), blast_desc, Description)
            ) %>% 
            select(ID, Veredict, Description, pident, Confidence, Sequence)
        }
      }
    }
    
    feat_num <- df_features %>% select(-Sequence)
    pca_res <- prcomp(feat_num, scale. = TRUE)
    final_df$PC1 <- pca_res$x[,1]; final_df$PC2 <- pca_res$x[,2]; final_df$PC3 <- pca_res$x[,3]
    return(final_df)
  })
  
  output$results_table <- renderDT({ 
    datatable(results_data() %>% select(-PC1, -PC2, -PC3), options=list(scrollX=TRUE)) 
  })
  
  output$pca_3d_plot <- renderPlotly({
    req(results_data())
    p <- plot_ly(results_data(), x=~PC1, y=~PC2, z=~PC3, color=~Veredict, 
                 text=~ID, hoverinfo="text+x+y+z", type='scatter3d', mode='markers')
    pca_plot_obj(p); p
  })
  
  output$down_pca_3d_html <- downloadHandler(
    filename = "PCA_3D.html", content = function(file) { saveWidget(pca_plot_obj(), file) }
  )
  
  output$download_tsv <- downloadHandler(
    filename = "Results.tsv", content = function(file) { write.table(results_data() %>% select(-PC1, -PC2, -PC3), file, sep="\t", row.names=FALSE) }
  )
}

shinyApp(ui, server)