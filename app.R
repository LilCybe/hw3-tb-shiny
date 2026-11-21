library(shiny)
library(ggplot2)

TBnanostring <- readRDS("TBnanostring.rds")
# Keep TB status out of the dimension reduction and clustering.
genes <- as.matrix(TBnanostring[, -1])

# Compute the three reductions once when the app starts.
pca_fit <- prcomp(genes, center = TRUE, scale. = FALSE)
set.seed(2026)
umap_fit <- uwot::umap(
  genes, n_neighbors = 15, min_dist = 0.1, n_components = 2,
  n_threads = 1, n_sgd_threads = 1, verbose = FALSE
)
# Subtract each gene's minimum to make every NMF input nonnegative.
genes_positive <- sweep(genes, 2, apply(genes, 2, min), "-")
nmf_fit <- NMF::nmf(
  genes_positive, rank = 2, method = "lee", nrun = 10, seed = 2026,
  .options = list(parallel = FALSE)
)

# Each matrix has one patient per row and two plotting coordinates.
reductions <- list(
  PCA = pca_fit$x[, 1:2],
  UMAP = umap_fit,
  NMF = NMF::basis(nmf_fit)
)
colnames(reductions$UMAP) <- c("UMAP 1", "UMAP 2")
colnames(reductions$NMF) <- c("Component 1", "Component 2")

ui <- fluidPage(
  titlePanel("TB gene expression: dimension reduction and k-means"),
  sidebarLayout(
    sidebarPanel(
      selectInput("method", "Dimension reduction", c("PCA", "UMAP", "NMF"),
                  selected = "PCA", selectize = FALSE),
      numericInput("clusters", "Number of clusters (k)", 2, min = 1, max = 9,
                   step = 1),
      helpText("179 patients and 107 genes: 79 TB and 100 LTBI."),
      helpText("Color shows cluster; shape shows TB status. Black crosses mark cluster centers."),
      width = 3
    ),
    mainPanel(
      plotOutput("cluster_plot", height = "420px"),
      fluidRow(
        column(7, h4("Clusters vs. TB status"), tableOutput("status_table")),
        column(5, h4("Adjusted Rand index"), textOutput("ari"))
      ),
      width = 9
    )
  )
)

server <- function(input, output, session) {
  # Switching methods selects stored coordinates without refitting the reduction.
  selected_data <- reactive(reductions[[input$method]])

  clusters <- reactive({
    req(input$clusters %in% 1:9)
    # Use 25 starting partitions and keep the one with the smallest within-cluster sum of squares.
    set.seed(2026)
    kmeans(selected_data(), centers = input$clusters, nstart = 25)
  })

  output$cluster_plot <- renderPlot({
    coordinates <- selected_data()
    fit <- clusters()
    plot_data <- data.frame(
      Dimension1 = coordinates[, 1], Dimension2 = coordinates[, 2],
      Cluster = factor(fit$cluster), TB_Status = TBnanostring$TB_Status
    )
    centers <- data.frame(Dimension1 = fit$centers[, 1],
                          Dimension2 = fit$centers[, 2])
    ggplot(plot_data, aes(Dimension1, Dimension2, color = Cluster,
                         shape = TB_Status)) +
      geom_point(size = 2.7, alpha = 0.85) +
      # Reuse the same fit for points and centers so they describe the same partition.
      geom_point(data = centers, aes(Dimension1, Dimension2),
                 inherit.aes = FALSE, shape = 4, size = 6, stroke = 1.6) +
      scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73",
                                    "#CC79A7", "#E69F00", "#56B4E9",
                                    "#A6761D", "#666666", "#332288")) +
      scale_shape_manual(values = c(TB = 16, LTBI = 17)) +
      labs(title = paste(input$method, "with", input$clusters, "clusters"),
           x = colnames(coordinates)[1], y = colnames(coordinates)[2]) +
      theme_minimal(base_size = 14)
  })

  output$status_table <- renderTable({
    # TB status is used only to assess the unsupervised clusters.
    as.data.frame.matrix(table(Cluster = clusters()$cluster,
                               TB_Status = TBnanostring$TB_Status))
  }, rownames = TRUE, striped = TRUE)

  output$ari <- renderText({
    # ARI compares pairwise group membership and adjusts for chance agreement.
    sprintf("%.4f", mclust::adjustedRandIndex(clusters()$cluster,
                                             TBnanostring$TB_Status))
  })
}

shinyApp(ui = ui, server = server)
