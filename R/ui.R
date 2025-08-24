# R/ui.R
#' @keywords internal
ui <- function() {
  shiny::fluidPage(
    titlePanel("CITEscope"), #title
    sidebarLayout(
      # Sidebar contains: File Upload, Action Button, Upload Status, Download
      sidebarPanel(
        width = 3,
        # Data Upload and Progress Status
        h4("Please select the upload type"),
        radioButtons(
          inputId = "uploadMode",
          label   = "Please select the upload type：",
          choices = c("ZIP File (.zip)" = "zip",
                      "Upload Files Separately"   = "files"),
          selected = "zip"
        ),

        # zip file
        conditionalPanel(
          condition = "input.uploadMode == 'zip'",
          fileInput("zipFile", "upload 10x data (.zip)", accept = ".zip")
        ),

        # upload files separately
        conditionalPanel(
          condition = "input.uploadMode == 'files'",
          fileInput(
            "mtxFile", "upload matrix.mtx(.gz)：", accept = c(".mtx",".gz"), multiple = FALSE
          ),
          fileInput(
            "bcFile",  "upload barcodes.tsv(.gz)：", accept = c(".tsv",".gz"), multiple = FALSE
          ),
          fileInput(
            "ftFile",  "upload features.tsv(.gz)：", accept = c(".tsv",".gz"), multiple = FALSE
          ),
          hr(),
          h5("Or upload one expression matrix file TSV："),
          fileInput(
            "tsvFile", "Upload a single .tsv(.gz)：", accept = c(".tsv",".gz"), multiple = FALSE
          )
        ),
        actionButton("btnPreprocess","Preprocessing"), br(), textOutput("status_preprocess"), hr(),
        actionButton("btnWNN","WNN Clustering + UMAP"), br(), textOutput("status_wnn"), hr(),
        actionButton("btnSingleR","SingleR Automatic Annotation"), br(), textOutput("status_singleR"), hr(),
        actionButton("btnInfection","Infection Status Classification"), br(), textOutput("status_infection"), hr(),
        actionButton("btnModuleScore","Module Scoring"), br(), textOutput("status_module"), hr(),
        actionButton("btnTI", "Trajectory Inference"), br(), textOutput("status_ti"), hr(),
        actionButton("btnDE","Differential Expression Analysis"), br(), textOutput("status_de")
      ),

      mainPanel(
        tabsetPanel(
          id = "mainTabs",
          selected = "Preview",

          # Preview label
          tabPanel("Preview",
                   h4("Meta Data (Top 10 Rows)"),
                   tableOutput("tblPreview"),
                   br(),
                   downloadButton("dlPreviewMeta", "download full metadata CSV")
          ),

          # WNN Clustering label
          tabPanel("WNN Clustering",
                   fluidRow(
                     column(8, h4("WNN UMAP + Clustering"), plotOutput("plot1")),
                     column(4, downloadButton("dlPlot1", "download"))
                   ),
                   fluidRow(
                     column(8, h4("SCT.weight on UMAP"), plotOutput("plot2")),
                     column(4, downloadButton("dlPlot2", "download"))
                   ),
                   fluidRow(
                     column(8, h4("ADT.weight on UMAP"), plotOutput("plot3")),
                     column(4, downloadButton("dlPlot3", "download"))
                   )
          ),


          # SingleR Automatic Annotation
          tabPanel("SingleR Automatic Annotation",
                   fluidRow(
                     column(6,
                            selectInput(
                              "species", "Species",
                              choices = c("Mouse", "Human"),
                              selected = "Mouse"
                            )
                     ),
                     column(6,
                            uiOutput("refset_ui")  # dynamic reference set
                     )
                   ),
                   hr(),
                   fluidRow(
                     column(5,
                            selectizeInput(
                              inputId = "markerSelect",
                              label = "Choose or Input Core Markers (Maximum 5)：",
                              choices = NULL,
                              multiple = TRUE,
                              options = list(
                                maxItems = 5,
                                create = TRUE,
                                placeholder = 'Choose or Input Core Markers')
                            ),
                            br(),
                            helpText("Default candidate markers in current data："),
                            DT::dataTableOutput("tblValidMarkers")
                     ),
                     column(7,
                            helpText("All genes in current dataset (searchable & filterable)："),
                            DT::dataTableOutput("tblAllGenes")
                     )
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("SingleR Cell Type Annotation"), plotOutput("plotSingleR")),
                     column(4, downloadButton("dlPlotSingleR", "download Annotation plot"))
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("ADT Marker UMAP"), plotOutput("plotADTmarkers")),
                     column(4, downloadButton("dlPlotADTmarkers", "download ADT plot"))
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("SCT Top Features"), plotOutput("plotSCTmarkers")),
                     column(4, downloadButton("dlPlotSCTmarkers", "download SCT plot"))
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("Custom Marker UMAP"), plotOutput("plotMarkers")),
                     column(4, downloadButton("dlPlotMarkers", "download Marker UMAP plot"))
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("Custom Marker Expression Distribution"), plotOutput("plotMarkerVln")),
                     column(4, downloadButton("dlPlotMarkerVln", "download Marker Violin plot"))
                   )
          ),

          # Infection Status Classification
          tabPanel(
            title = "Infection Status Classification",
            fluidRow(
              column(6,h4("NP Expression (Violin Plot)"),
                     plotOutput("violinNP"),
                     br(),
                     downloadButton("dlPlotViolin", "download Violin")
              ),
              column(6,h4("Infection Proportion (Bar Plot)"),
                     plotOutput("barInfection"),
                     br(),
                     downloadButton("dlPlotBar", "download Bar")
              )
            ),
            hr(),
            fluidRow(
              column(8, h4("UMAP Visualization of Infection Status"), plotOutput("plotInfection")),
              column(4, downloadButton("dlPlotInfection", "download"))
            )
          ),

          # Module Scoring
          tabPanel("Module Scoring",
                   fluidRow(
                     column(6,
                            selectInput(
                              inputId = "isgList",
                              label   = "Select ISG List：",
                              choices = c(
                                "Hallmark IFN-α+IFN-γ" = "H",
                                "MSigDB C7 IFN Gene Sets" = "C7",
                                "Cross-Species Core ISGs (Human/Mouse/Cow/Sheep)" = "Core",
                                "Viral Sensor/Signaling"  = "PRR",
                                " Upload Custom List (CSV File)" = "Custom"
                              ),
                              selected = "H"
                            ),
                            fileInput(
                              inputId = "customISG",
                              label = "Upload Custom ISG List (Optional)：",
                              accept = ".csv"
                            ),
                            helpText(
                              "if you need more ISG：visit",
                              tags$a(href="https://www.interferome.org/", target="_blank","Interferome v2.0"),
                              "Search online and download gene sets for the desired species/treatment conditions."
                            )
                     )
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("Merge Modules Score UMAP"), plotOutput("plotModuleCombined")),
                     column(4, downloadButton("dlPlotModuleCombined", "download"))
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("ISG Module Score"), plotOutput("plotISG")),
                     column(4, downloadButton("dlPlotISG", "download ISG plot"))
                   ),
                   fluidRow(
                     column(8, h4("Exh Module Score"), plotOutput("plotExh")),
                     column(4, downloadButton("dlPlotExh", "download Exh plot"))
                   ),
                   fluidRow(
                     column(8, h4("Cytokine/Inflammatory Module Score (UMAP)"), plotOutput("plotCytokine")),
                     column(4, downloadButton("dlPlotCytokine", "download Cytokine plot"))
                   ),
                   hr(),
                   fluidRow(
                     column(6,
                            selectInput(
                              "ms_groupby", "Group by:",
                              choices  = c("infection_status", "seurat_clusters"),
                              selected = "infection_status"
                            )
                     ),
                     column(6,
                            selectizeInput(
                              "ms_modules", "Select modules to display:",
                              choices  = c("ISG_Score", "Exh_Score", "Cytokine_Score"),
                              selected = c("ISG_Score", "Exh_Score", "Cytokine_Score"),
                              multiple = TRUE,
                              options  = list(maxItems = 3)
                            )
                     )
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("DotPlot: Group × Module"), plotOutput("plotMSDot", height = "420px")),
                     column(4, br(), downloadButton("dlMSDot", "Download DotPlot"))
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("Heatmap: Module Scores (Group Means)"), plotOutput("plotMSHeat", height = "420px")),
                     column(4, br(), downloadButton("dlMSHeat", "Download Heatmap"))
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("RidgePlot: Score Distributions"), plotOutput("plotMSRidge", height = "520px")),
                     column(4, br(), downloadButton("dlMSRidge", "Download RidgePlot"))
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("Pseudobulk Bar: Mean ± SEM"), plotOutput("plotMSBar", height = "400px")),
                     column(4, br(), downloadButton("dlMSBar", "Download Bar Plot"))
                   ),
                   hr(),
                   fluidRow(
                     column(8, h4("ISG_Score vs NP-ADT (Correlation)"), plotOutput("plotISGvsNP", height = "400px")),
                     column(4, br(), downloadButton("dlISGvsNP", "Download Scatter Plot"))
                   )
          ),

          # Trajectory Inference (Monocle3)
          tabPanel("Trajectory Inference",
                   fluidRow(
                     column(5,
                            h4("Settings"),
                            radioButtons(
                              "ti_root_mode", "Root selection strategy",
                              choices = c("Manual cluster" = "manual",
                                          "Highest ADT marker (quantile)" = "adt"),
                              selected = "adt"
                            ),
                            uiOutput("ti_root_cluster_ui"),
                            uiOutput("ti_adt_feature_ui"),
                            conditionalPanel(
                              "input.ti_root_mode == 'adt'",
                              sliderInput("ti_adt_quantile", "Root quantile (ADT):",
                                          min = 0.5, max = 0.99, value = 0.9, step = 0.01
                              )
                            ),
                            helpText("RNA determines the trajectory. ADT is used as a phenotypic overlay for annotation.")
                     ),
                     column(7,
                            h4("UMAP + Learned Graph (pseudotime)"),
                            plotOutput("plotTI_cells", height = "460px"),
                            downloadButton("dlTIcells", "download UMAP+Graph")
                     )
                   ),
                   hr(),
                   fluidRow(
                     column(6,
                            h4("UMAP colored by ADT"),
                            plotOutput("plotTI_ADT", height = "420px"),
                            downloadButton("dlTIadt", "download ADT-colored UMAP")
                     ),
                     column(6,
                            h4("ADT vs Pseudotime"),
                            plotOutput("plotTI_ADTvsPT", height = "420px"),
                            downloadButton("dlTIcorr", "download ADT~pseudotime")
                     )
                   ),
                   hr(),
                   fluidRow(
                     column(
                       8,
                       h4("Pseudotime by Group (Ridge)"),
                       plotOutput("plotTI_ptRidge", height = "420px")
                     ),
                     column(
                       4,
                       br(),
                       downloadButton("dlTI_ptRidge", "download")
                     )
                   )
          ),

          # Differential Expression Analysis
          tabPanel("Differential Expression Analysis",
                   h4("Differential Gene Results"), DT::dataTableOutput("deTable"), hr(),
                   fluidRow(column(6,plotOutput("plotMA")), column(6,plotOutput("plotVolcano"))), hr(),
                   fluidRow(column(12,h4("Top10 DEGs Heatmap"),plotOutput("plotHeatmap", height = "400px"))), hr(),
                   fluidRow(column(12,h4("GO/KEGG enrichment"),plotOutput("plotGOKEGG", height = "700px"))), hr(),
                   downloadButton("dlDECSV","download DEG CSV"),
                   downloadButton("dlGOCSV","download GO+KEGG CSV"),
                   downloadButton("dlPlotMA","download MA plot"),
                   downloadButton("dlPlotVolcano","download volcano plot"),
                   downloadButton("dlPlotHeatmap","download heatmap"),
                   downloadButton("dlPlotGOKEGG","download enrichment plot")
          )
        )
      )
    )
  )
}
