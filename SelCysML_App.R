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
library(processx)

# --- 1. CONFIGURATION ---
options(shiny.maxRequestSize = 800 * 1024^2)
BLAST_CONFIG <- list(
  EXE  = "C:/ncbi-blast-2.17.0+/bin/blastp.exe", 
  BSEP = "C:/Users/jllavin/Desktop/SelCys_data/BLASTdb/selenoprots_BSepDB",
  UNI  = "C:/Users/jllavin/Desktop/SelCys_data/BLASTdb/uniprot_db"
)

# --- 2. EXTRACTION AND BLAST FUNCTIONS ---
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

run_blast_safe <- function(fasta_path, db_path, evalue_cut) {
  evalue_str <- paste0("1e", evalue_cut)
  out_file <- tempfile(fileext = ".tsv")
  
  tryCatch({
    processx::run(
      command = BLAST_CONFIG$EXE,
      args = c("-query", fasta_path, 
               "-db", db_path, 
               "-outfmt", "6 qseqid stitle pident", 
               "-max_target_seqs", "1", 
               "-evalue", evalue_str),
      stdout = out_file,
      error_on_status = FALSE
    )
  }, error = function(e) return(NULL))
  
  if (!file.exists(out_file) || file.size(out_file) == 0) return(NULL)
  
  read.table(out_file, sep = "\t", quote = "", fill = TRUE, 
             col.names = c("qseqid", "Description", "pident"))
}

# --- 3. UI ---
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("SelCys: SelenoCysteine protein detector via ML + BLAST"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta_in", "Upload FASTA", accept = c(".fasta", ".fa")),
      selectInput("model_type", "Choose model:", 
                  choices = c("General (Eukaryotes/Prokaryotes)" = "General", 
                              "Bacteria Only" = "Bacterial")),
      checkboxInput("enable_blast", "Enable BLAST analysis", value = TRUE),
      sliderInput("evalue_cut", "E-value cut-off (log10):", min = -10, max = -1, value = -5, step = 1, pre = "1e"),
      helpText("Note: BLAST analysis may take several minutes."),
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

# --- 4. SERVER ---
server <- function(input, output) {
  current_model <- reactive({
    filename <- if(input$model_type == "General") "SelCysGeneral.rds" else "SelCysBacterial.rds"
    path <- file.path("models", filename)
    if(file.exists(path)) return(readRDS(path)) else { 
      showNotification(paste("Error: File not found", filename), type = "error"); return(NULL) 
    }
  })
  
  pca_plot_obj <- reactiveVal(NULL)
  
  results_data <- eventReactive(input$predict_btn, {
    req(input$fasta_in, current_model()) 
    raw <- readAAStringSet(input$fasta_in$datapath)
    ids_limpios <- sapply(names(raw), clean_id_func)
    names(raw) <- ids_limpios
    df_features <- extract_features_logic(as.character(raw))
    if(is.null(df_features)) stop("Error: Could not extract descriptors.")
    
    df_para_modelo <- df_features
    df_para_modelo$ID <- ids_limpios[1:nrow(df_features)]
    
    preds_list <- map(current_model(), ~predict(.x, df_para_modelo, type = "prob"))
    prob_selcys <- rowMeans(do.call(cbind, map(preds_list, ~ .x$.pred_SelCys)))
    
    final_df <- data.frame(
      ID = df_para_modelo$ID,
      Veredict = ifelse(prob_selcys > 0.5, "Selenoprotein", "Standard"),
      Pred_test = as.character("No evaluation"),
      Description = as.character("No evaluation"),
      pident = as.numeric(NA),
      Confidence = round(pmax(prob_selcys, 1-prob_selcys) * 100, 2),
      Sequence = df_para_modelo$Sequence,
      stringsAsFactors = FALSE
    )
    
    if(input$enable_blast) {
      positives <- final_df %>% filter(Veredict == "Selenoprotein")
      if(nrow(positives) > 0) {
        temp_fasta <- tempfile(fileext = ".fasta")
        writeXStringSet(AAStringSet(setNames(positives$Sequence, positives$ID)), temp_fasta)
        
        res_bsep <- run_blast_safe(temp_fasta, BLAST_CONFIG$BSEP, input$evalue_cut)
        if(!is.null(res_bsep)) {
          for(i in seq_len(nrow(res_bsep))) {
            idx <- which(final_df$ID == res_bsep$qseqid[i])
            if(length(idx) > 0) final_df$Pred_test[idx] <- as.character(res_bsep$Description[i])
          }
        }
        res_uni <- run_blast_safe(temp_fasta, BLAST_CONFIG$UNI, input$evalue_cut)
        if(!is.null(res_uni)) {
          for(i in seq_len(nrow(res_uni))) {
            idx <- which(final_df$ID == res_uni$qseqid[i])
            if(length(idx) > 0) {
              final_df$Description[idx] <- as.character(res_uni$Description[i])
              final_df$pident[idx] <- res_uni$pident[i]
            }
          }
        }
        unlink(temp_fasta)
      }
    }
    
    final_df <- final_df %>%
      mutate(
        Evidence_Status = case_when(
          Veredict == "Selenoprotein" & Pred_test != "No evaluation" ~ "Validated (ML + BSepDB)",
          Veredict == "Selenoprotein" & Pred_test == "No evaluation" ~ "Alert: Potential FP (ML only)",
          Veredict == "Standard" & Pred_test != "No evaluation"    ~ "Alert: Potential FN (BLAST+)",
          TRUE ~ "Standard (Concordant)"
        ),
        Confidence_Adjusted = case_when(
          Veredict == "Selenoprotein" & Pred_test != "No evaluation" ~ 99.9,
          Veredict == "Selenoprotein" & Pred_test == "No evaluation" ~ round(Confidence * 0.3, 1),
          TRUE ~ Confidence
        )
      ) %>%
      select(-Sequence, everything(), Sequence) # Moving Sequence to the end
    
    feat_num <- df_features %>% select(-Sequence)
    pca_res <- prcomp(feat_num, scale. = TRUE)
    final_df$PC1 <- pca_res$x[,1]; final_df$PC2 <- pca_res$x[,2]; final_df$PC3 <- pca_res$x[,3]
    return(final_df)
  })
  
  output$results_table <- renderDT({ 
    data_res <- results_data()
    vista <- data_res[, !(names(data_res) %in% c("PC1", "PC2", "PC3"))]
    datatable(vista, options=list(scrollX=TRUE)) 
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
    filename = "Results.tsv", content = function(file) { 
      data_out <- results_data()[, !(names(results_data()) %in% c("PC1", "PC2", "PC3"))]
      write.table(data_out, file, sep="\t", row.names=FALSE) 
    }
  )
}

shinyApp(ui, server)