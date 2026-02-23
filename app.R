# ---------------------------
# Load libraries
# ---------------------------
library(shiny)
library(dplyr)
install.packages("rsconnect")

# ---------------------------
# Load data
# ---------------------------
# Use file.choose() for now, if you want to select manually
expr_mat <- read.csv(file.choose(), row.names = 1, check.names = FALSE)
meta_df  <- read.csv(file.choose(), stringsAsFactors = FALSE)
umap_df  <- read.csv(file.choose(), stringsAsFactors = FALSE)

# ---------------------------
# Helper functions
# ---------------------------

# Scale 0–100 for plotting
scale_0_100 <- function(x){
  rng <- range(x, na.rm = TRUE)
  if(rng[1] == rng[2]) return(rep(50, length(x)))
  (x - rng[1]) / (rng[2] - rng[1]) * 100
}

# Compute gene statistics for a cell type
compute_gene_stats <- function(expr_mat, meta_df, target_ct){
  
  in_cells  <- meta_df$cell_id[meta_df$cell_type == target_ct]
  out_cells <- meta_df$cell_id[meta_df$cell_type != target_ct]
  
  xin  <- expr_mat[in_cells, , drop = FALSE]
  xout <- expr_mat[out_cells, , drop = FALSE]
  
  det_in  <- colMeans(xin > 0)
  det_out <- colMeans(xout > 0)
  
  mean_in  <- colMeans(xin)
  mean_out <- colMeans(xout)
  
  diff <- mean_in - mean_out
  
  data.frame(
    gene = colnames(expr_mat),
    det_in = det_in,
    det_out = det_out,
    mean_in = mean_in,
    mean_out = mean_out,
    diff = diff,
    stringsAsFactors = FALSE
  )
}

# Pick marker gene
pick_marker_gene <- function(gene_stats_df){
  gene_stats_df %>%
    arrange(desc(diff), desc(det_in)) %>%
    slice(1) %>%
    pull(gene)
}

# ---------------------------
# Shiny UI
# ---------------------------
ui <- fluidPage(
  
  titlePanel("Single-Cell RNA-seq Marker Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "cell_type",
        "Select Cell Type:",
        choices = unique(meta_df$cell_type)
      )
    ),
    
    mainPanel(
      plotOutput("umap_plot"),
      verbatimTextOutput("marker_info"),
      tableOutput("gene_table")
    )
  )
)

# ---------------------------
# Shiny Server
# ---------------------------
server <- function(input, output){
  
  # Reactive gene stats
  gene_stats <- reactive({
    compute_gene_stats(expr_mat, meta_df, input$cell_type)
  })
  
  # Reactive marker gene
  marker_gene <- reactive({
    pick_marker_gene(gene_stats())
  })
  
  # Prepare merged UMAP data with marker expression
  umap_data <- reactive({
    merged <- meta_df %>%
      left_join(umap_df, by = "cell_id")
    merged$marker_expr <- expr_mat[merged$cell_id, marker_gene()]
    merged$marker_scaled <- scale_0_100(merged$marker_expr)
    merged
  })
  
  # UMAP Plot
  output$umap_plot <- renderPlot({
    
    df <- umap_data()
    
    # Color by marker expression
    cols <- colorRampPalette(c("lightgrey", "blue", "red"))(101)
    point_cols <- cols[round(df$marker_scaled) + 1]
    
    plot(df$UMAP_1, df$UMAP_2,
         col = point_cols,
         pch = 16,
         main = paste("UMAP - Marker:", marker_gene()),
         xlab = "UMAP 1",
         ylab = "UMAP 2")
    
    legend("topright", legend = c("Low","High"),
           fill = c("blue","red"))
  })
  
  # Marker info text
  output$marker_info <- renderText({
    gs <- gene_stats()
    mg <- marker_gene()
    diff_value <- gs$diff[gs$gene == mg]
    paste0(
      "Selected Cell Type: ", input$cell_type, "\n",
      "Marker Gene: ", mg, "\n",
      "Specificity Score (diff): ", round(diff_value, 4)
    )
  })
  
  # Gene statistics table
  output$gene_table <- renderTable({
    gene_stats()
  })
}

# ---------------------------
# Run App
# ---------------------------
shinyApp(ui = ui, server = server)
library(rsconnect)
rsconnect::setAccountInfo(name='boufahjanarmin777',
                          token='3BA897526AAC60279E79B43A31B432D6',
                          secret='<SECRET>')
rsconnect::deployApp()
getwd()
list.files()