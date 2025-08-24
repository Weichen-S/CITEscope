# R/server.R
#' @keywords internal
server <- function(input, output, session) {
  # Reset plan and stop clusters
  shiny::onStop(function() {
    try(future::plan(future::sequential), silent = TRUE)
    try(future:::ClusterRegistry("stop"), silent = TRUE)
  })

    ## Global reactiveValues & States
    rvPlots <- reactiveValues(
      p1 = NULL, p2 = NULL, p3 = NULL,
      hist = NULL, inf = NULL,
      sr = NULL, adt = NULL, sct = NULL,
      mod = NULL, isg = NULL, ifn = NULL, exh = NULL, cyto = NULL,
      ma = NULL, volcano = NULL, heatmap = NULL, go = NULL
    )

    rvSeurat <- reactiveVal(NULL)
    deResults <- reactiveVal(NULL)
    maData <- reactiveVal(NULL)
    volcanoData <- reactiveVal(NULL)
    heatmapData <- reactiveVal(NULL)
    enrichResults <- reactiveVal(NULL)
    rvCDS <- reactiveVal(NULL)
    deAssay <- reactiveVal(NULL)
    status_pre <- reactiveVal("Awaiting Preprocessing...")
    status_wnn <- reactiveVal("Awaiting WNN...")
    status_inf <- reactiveVal("Awaiting Infection Classification...")
    status_sr <- reactiveVal("Awaiting SingleR...")
    status_module <- reactiveVal("Awaiting Module Scoring...")
    status_de <- reactiveVal("Awaiting DE Analysis...")
    status_ti <- reactiveVal("Awaiting TI...")

    # output
    output$status_preprocess <- renderText(status_pre())
    output$status_wnn <- renderText(status_wnn())
    output$status_infection  <- renderText(status_inf())
    output$status_singleR <- renderText(status_sr())
    output$status_module <- renderText(status_module())
    output$status_de <- renderText(status_de())
    output$status_ti <- renderText(status_ti())


    ### Preprocessing
    observeEvent(input$btnPreprocess, {
      status_pre("Preprocessing…")

      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::sequential)

      withProgress(message = "Preprocessing in Progress", value = 0, {

        # Create and Clear Temporary Directory
        incProgress(0.1, detail = "Create Temporary Directory…")
        tmpDir <- file.path(tempdir(), "upload_data")
        if (dir.exists(tmpDir)) unlink(tmpDir, recursive = TRUE)
        dir.create(tmpDir, recursive = TRUE)

        # unzip or Copy Input Files
        incProgress(0.1, detail = "Unzip/Copy Raw Files…")
        if (input$uploadMode == "zip") {
          req(input$zipFile)
          unzip(input$zipFile$datapath, exdir = tmpDir)

          # automatically locate directory containing 10X files
          dirs <- c(tmpDir, list.dirs(tmpDir, recursive = FALSE))
          dataDir <- dirs[sapply(dirs, function(d)
            all(c("matrix.mtx", "features.tsv", "barcodes.tsv") %in% basename(list.files(d))))][1]
        } else {
          req(input$mtxFile, input$bcFile, input$ftFile)
          file.copy(input$mtxFile$datapath, file.path(tmpDir, "matrix.mtx"))
          file.copy(input$bcFile$datapath, file.path(tmpDir, "barcodes.tsv"))
          file.copy(input$ftFile$datapath, file.path(tmpDir, "features.tsv"))
          dataDir <- tmpDir
        }

        # read matrix & features
        incProgress(0.2, detail = "read matrix.mtx & features.tsv…")
        feature_df <- read.delim(
          file.path(dataDir, "features.tsv"),
          header = FALSE, stringsAsFactors = FALSE
        )
        full_mat <- ReadMtx(
          mtx = file.path(dataDir, "matrix.mtx"),
          cells = file.path(dataDir, "barcodes.tsv"),
          features = file.path(dataDir, "features.tsv"),
          feature.column = 2
        )

        # separate RNA and ADT by 3rd column
        incProgress(0.2, detail = "sperate RNA & ADT data…")
        adt_idx <- which(feature_df[,3] == "Antibody Capture")
        rna_idx <- setdiff(seq_len(nrow(full_mat)), adt_idx)
        rna_counts <- full_mat[rna_idx, , drop = FALSE]
        adt_counts <- full_mat[adt_idx, , drop = FALSE]

        # create Seurat object & normalize RNA
        incProgress(0.2, detail = "create Seurat object & SCTransform (RNA)…")
        seu <- CreateSeuratObject(
          counts = rna_counts,
          project = "upload",
          min.cells = 3,
          min.features= 200
        )
        seu <- SCTransform(
          object  = seu,
          assay   = "RNA",
          method  = "glmGamPoi",
          verbose = FALSE
        )

        # add ADT and perform CLR + Scaling
        incProgress(0.1, detail = "add & normalize ADT (CLR + Scale)…")
        common_cells <- intersect(colnames(seu), colnames(adt_counts))
        seu[["ADT"]] <- CreateAssayObject(counts = adt_counts[, common_cells, drop = FALSE])
        seu <- NormalizeData(seu, assay = "ADT", normalization.method = "CLR", verbose = FALSE)
        seu <- ScaleData(seu, assay = "ADT", verbose = FALSE)

        # design variable features for ADT assay
        incProgress(0.1, detail = "ADT: Identify variable features…")
        seu <- FindVariableFeatures(
          object = seu,
          assay = "ADT",
          selection.method = "vst",
          verbose = FALSE
        )

        incProgress(0.1, detail = "Preprocessing Completed！")
        rvSeurat(seu)
      })

      # Render Preview Table
      output$tblPreview <- renderTable(
        head(rvSeurat()@meta.data, 10),
        rownames = TRUE
      )
      status_pre("Preprocessing Completed！")
    })



    ##  WNN Clustering + UMAP
    observeEvent(input$btnWNN, {
      req(rvSeurat())
      status_wnn("Running WNN Clustering + UMAP…")
      withProgress(message = " WNN Clustering in Progress...",
                   value = 0, {
                     seu <- rvSeurat()

                     # RNA PCA
                     incProgress(0.2, detail = "RNA PCA…")
                     seu <- RunPCA(seu,
                                   assay = "SCT",
                                   reduction.name = "pca",
                                   npcs = 30,
                                   verbose = FALSE)

                     # ADT PCA
                     incProgress(0.4, detail = "ADT PCA…")
                     seu <- RunPCA(seu,
                                   assay = "ADT",
                                   reduction.name = "apca",
                                   npcs = 30,
                                   verbose = FALSE)

                     # Dynamically retrieve the number of available PCs
                     pca_n <- ncol( Embeddings(seu, "pca" ) )
                     apca_n <- ncol( Embeddings(seu, "apca") )
                     dims_list <- list(
                       pca = seq_len(pca_n),
                       apca = seq_len(apca_n)
                     )

                     # Calculate Multimodal Nearest Neighbors
                     incProgress(0.6, detail = "Calculate Multimodal Nearest Neighbors…")
                     seu <- FindMultiModalNeighbors(
                       object = seu,
                       reduction.list = c("pca","apca"),
                       dims.list = dims_list,
                       modality.weight.name = c("SCT.weight","ADT.weight"),
                       verbose = FALSE
                     )

                     # UMAP + Clustering
                     incProgress(0.8, detail = "UMAP + Clustering…")
                     seu <- RunUMAP(seu,
                                    nn.name = "weighted.nn",
                                    reduction = "wnn",
                                    reduction.name = "wnn.umap",
                                    umap.method = "uwot",
                                    verbose = FALSE)
                     seu <- FindClusters(seu,
                                         graph.name = "wsnn",
                                         resolution = 0.5,
                                         verbose = FALSE)

                     incProgress(1, detail = "finish")
                     rvSeurat(seu)
                   })

      # ggplot
      output$plot1 <- renderPlot({
        p <- DimPlot(rvSeurat(), reduction = "wnn.umap", label = TRUE) +
          ggtitle("WNN UMAP + Clusters")
        rvPlots$p1 <- p
        p
      })
      output$plot2 <- renderPlot({
        p <- FeaturePlot(rvSeurat(),
                         features = "SCT.weight",
                         reduction = "wnn.umap",
                         pt.size = 0.5) +
          ggtitle("SCT.weight on WNN UMAP")
        rvPlots$p2 <- p
        p
      })
      output$plot3 <- renderPlot({
        p <- FeaturePlot(rvSeurat(),
                         features = "ADT.weight",
                         reduction = "wnn.umap",
                         pt.size = 0.5) +
          ggtitle("ADT.weight on WNN UMAP")
        rvPlots$p3 <- p
        p
      })

      status_wnn("WNN Clustering Completed")
    })


    ###  SingleR Automatic Annotation
    default_markers <- c(
      "CD14","CD68","CSF1R","CD163",
      "CD3D","CD3E","CD4","CD8A",
      "CD19","CD79A","MS4A1",
      "NCAM1","NCR1","NKG7","GNLY",
      "EPCAM","KRT8","KRT18",
      "PDGFRA","PDGFRB","COL1A1",
      "PECAM1","VWF","CLDN5"
    )

    #  species dynamic reference set
    observeEvent(input$species, {
      if (is.null(input$species)) return(NULL)
      if (input$species == "Human") {
        refs <- c(
          "Human Primary Cell Atlas" = "HumanPrimaryCellAtlasData",
          "Blueprint/ENCODE (Blood/Immune)" = "BlueprintEncodeData",
          "Monaco Immune" = "MonacoImmuneData",
          "DICE (Immune)" = "DatabaseImmuneCellExpressionData",
          "Novershtern Hematopoietic" = "NovershternHematopoieticData"
        )
      } else {
        refs <- c(
          "ImmGen (Mouse immune)" = "ImmGenData",
          "MouseRNAseq (general)" = "MouseRNAseqData"
        )
      }
      output$refset_ui <- renderUI({
        selectInput("refset", "Reference", choices = refs, selected = refs[1])
      })
    }, ignoreInit = FALSE)

    # load function
    .get_ref <- function(species, refset) {
      if (species == "Human") {
        switch(refset,
               "HumanPrimaryCellAtlasData"  = celldex::HumanPrimaryCellAtlasData(),
               "BlueprintEncodeData" = celldex::BlueprintEncodeData(),
               "MonacoImmuneData" = celldex::MonacoImmuneData(),
               "DatabaseImmuneCellExpressionData" = celldex::DatabaseImmuneCellExpressionData(),
               "NovershternHematopoieticData" = celldex::NovershternHematopoieticData()
        )
      } else {
        switch(refset,
               "ImmGenData" = celldex::ImmGenData(),
               "MouseRNAseqData" = celldex::MouseRNAseqData()
        )
      }
    }
    ## Update Marker Dropdown & Validation Table
    observeEvent(rvSeurat(), {
      seu <- rvSeurat()
      avail <- intersect(default_markers, rownames(seu[["SCT"]]))

      # core Marker 4 lists ——
      output$tblValidMarkers <- DT::renderDataTable({
        n  <- length(default_markers)
        half <- n / 2
        pres <- default_markers %in% avail

        # Build Four Columns：Marker1 / Present1 / Marker2 / Present2
        df_core <- data.frame(
          Marker1 = default_markers[1:half],
          Present1 = pres[1:half],
          Marker2 = default_markers[(half+1):n],
          Present2 = pres[(half+1):n],
          stringsAsFactors = FALSE
        )

        DT::datatable(
          df_core,
          rownames = FALSE,
          options  = list(dom = 't', pageLength = half)
        )
      })

      # Full Gene Dynamics Table
      output$tblAllGenes <- DT::renderDataTable({
        df_all <- data.frame(
          Feature = rownames(seu[["SCT"]]),
          stringsAsFactors = FALSE
        )
        DT::datatable(
          df_all,
          rownames = FALSE,
          options  = list(pageLength = 10, lengthMenu = c(10, 50, 100))
        )
      })

      # update dropdown
      if (length(avail) == 0) {
        showNotification("No preset core markers found in the current dataset.",
                         type = "warning")
        updateSelectizeInput(
          session, "markerSelect",
          choices = rownames(seu[["SCT"]]),
          selected = NULL,
          server = TRUE,
          options = list(create = TRUE, placeholder = "Please enter or select a gene")
        )
      } else {
        updateSelectizeInput(
          session, "markerSelect",
          choices  = avail,
          selected = head(avail, 2),
          server = TRUE,
          options  = list(create = TRUE)
        )
      }
    })

    ### SingleR Annotation
    observeEvent(input$btnSingleR, {
      req(rvSeurat())
      status_sr("Running SingleR Annotation...")
      withProgress(message = "SingleR Annotation in Progress...", value = 0, {
        incProgress(0.3, detail = "Loading Reference Data...")
        seu <- rvSeurat()

        # Load reference set    default / mouse
        species_sel <- if (is.null(input$species)) "Mouse" else input$species
        refset_sel  <- if (is.null(input$refset))  (if (species_sel == "Human") "HumanPrimaryCellAtlasData" else "MouseRNAseqData") else input$refset
        ref <- .get_ref(species_sel, refset_sel)

        incProgress(0.6, detail = "Running SingleR()...")
        test_mat <- GetAssayData(seu, assay = "SCT", layer = "data")
        pred     <- SingleR(test = test_mat, ref = ref, labels = ref$label.main)

        seu$SingleR_label <- pred$labels
        incProgress(1, detail = "finished")
        rvSeurat(seu)
      })

      # ggplot
      seu <- rvSeurat()
      p_sr <- DimPlot(seu, reduction="wnn.umap", group.by="SingleR_label", repel=TRUE) +
        ggtitle("SingleR Annotation") + theme_minimal(base_size=14) +
        theme(plot.title=element_text(hjust=0.5))

      DefaultAssay(seu) <- "ADT"
      p_adt <- FeaturePlot(seu, features=rownames(seu[["ADT"]]), reduction="wnn.umap",
                           pt.size=1.5, ncol=4) +
        ggtitle("ADT Marker UMAP") + theme_minimal(base_size=14) +
        theme(plot.title=element_text(hjust=0.5))

      DefaultAssay(seu) <- "SCT"
      seu <- FindVariableFeatures(seu, assay="SCT", selection.method="vst", verbose=FALSE)
      top6 <- head(VariableFeatures(seu), 6)
      p_sct <- FeaturePlot(seu, features=top6, reduction="wnn.umap",
                           pt.size=1.5, ncol=3) +
        ggtitle("SCT Top Features") +
        theme_minimal(base_size=14) +
        theme(plot.title=element_text(hjust=0.5))

      rvPlots$sr <- p_sr;  output$plotSingleR <- renderPlot(p_sr, height=400)
      rvPlots$adt <- p_adt; output$plotADTmarkers <- renderPlot(p_adt, height=400)
      rvPlots$sct <- p_sct; output$plotSCTmarkers <- renderPlot(p_sct, height=400)

      status_sr("SingleR Annotation Completed")
    })

    # Custom Marker Selection & Visualization
    observeEvent(input$markerSelect, {
      req(rvSeurat())
      seu <- rvSeurat()
      sel <- input$markerSelect
      present <- intersect(sel, rownames(seu))
      validate(need(length(present) > 0, "Please select at least one valid marker！"))

      # Plot List for a Single Gene
      plots_list <- FeaturePlot(
        object = seu,
        features = present,
        reduction = "wnn.umap",
        pt.size = 1.5,
        combine = FALSE
      )

      # Add individual titles to each subplot
      plots_list <- lapply(seq_along(plots_list), function(i) {
        plots_list[[i]] +
          ggtitle(present[i]) +
          theme_minimal(base_size = 14) +
          theme(
            plot.title = element_text(hjust = 0.5)
          )
      })

      p_mark <- patchwork::wrap_plots(plots_list, ncol = min(length(plots_list), 3))

      plots_vln <- VlnPlot(
        object = seu,
        features = present,
        group.by = "seurat_clusters",
        pt.size = 0,
        combine = FALSE
      )
      plots_vln <- lapply(seq_along(plots_vln), function(i) {
        plots_vln[[i]] +
          ggtitle(present[i]) +
          theme_minimal(base_size = 14) +
          theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1)
          )
      })
      p_vln <- patchwork::wrap_plots(plots_vln, ncol = 1)

      rvPlots$markers <- p_mark
      rvPlots$markerVln <- p_vln
      output$plotMarkers <- renderPlot(p_mark, height = 400)
      output$plotMarkerVln <- renderPlot(p_vln, height = 800)
    })



    ###  Infection Status Classification
    observeEvent(input$btnInfection, {
      req(rvSeurat())
      status_inf("Running Infection Status Classification…")

      withProgress(message = "Infection Status Classification in Progress...",
                   value = 0, {
                     seu <- rvSeurat()

                     # check ADT assay
                     incProgress(0.1, detail = "check ADT assay…")
                     if (!"ADT" %in% names(seu@assays)) {
                       showNotification("ADT assay not detected, classification unavailable！",
                                        type = "error")
                       return(NULL)
                     }

                     # find NP marker
                     incProgress(0.2, detail = "find NP marker…")
                     adt_feats <- rownames(seu[["ADT"]])
                     np_feats <- grep("^NP", adt_feats, value = TRUE)
                     np_feat <- np_feats[1]

                     # extract and classify
                     incProgress(0.4, detail = "Extract and Classify…")
                     adt_data <- GetAssayData(seu, assay="ADT", layer="data")
                     np_vals <- adt_data[np_feat, ]
                     thresh <- median(np_vals) + 2 * sd(np_vals)
                     status_vec <- factor(ifelse(np_vals > thresh, "infected", "uninfected"),
                                          levels = c("uninfected","infected"))
                     seu <- AddMetaData(seu, metadata=status_vec, col.name="infection_status")
                     rvSeurat(seu)

                     # generate Plot
                     incProgress(0.8, detail = "generate Plot…")
                     # violin plot
                     violin_plot <- VlnPlot(
                       object = seu,
                       assay = "ADT",
                       features = np_feat,
                       group.by = "infection_status",
                       pt.size = 0,
                       combine = FALSE
                     )[[1]] + ggtitle(paste0(np_feat, " by infection status"))

                     # bar plot
                     bar_df <- as.data.frame(table(status = status_vec))
                     bar_plot <- ggplot(bar_df, aes(x=status, y=Freq, fill=status)) +
                       geom_col() +
                       geom_text(aes(label=Freq), vjust=-0.5) +
                       theme_minimal() +
                       labs(x="status", y="cells number",
                            title="Proportion of Cells by Infection Status") +
                       guides(fill="none")

                     # UMAP
                     inf_plot <- DimPlot(
                       seu,
                       reduction = "wnn.umap",
                       group.by  = "infection_status"
                     ) + ggtitle("Infection status UMAP")

                     incProgress(1, detail = "completed")
                   })

      rvPlots$hist <- violin_plot
      rvPlots$bar <- bar_plot
      rvPlots$inf <- inf_plot

      # Render
      output$violinNP <- renderPlot({ violin_plot })
      output$barInfection <- renderPlot({ bar_plot })
      output$plotInfection<- renderPlot({ inf_plot })
      status_inf("Infection Status Classification Completed")
    })




    ### Module Scoring Analysis ----
    observeEvent(input$btnModuleScore, {
      req(rvSeurat())
      status_module("Running Module Scoring ...")

      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::sequential)

      # helper function
      get_ms_cols <- function(seu, selected = NULL) {
        all_cols <- grep("_Score$", colnames(seu@meta.data), value = TRUE)
        if (is.null(selected) || length(selected) == 0) return(all_cols)

        # match the complete column name directly
        direct <- intersect(selected, all_cols)

        norm <- function(x) tolower(gsub("[^A-Za-z0-9]+", "", x)) # Remove non-alphanumeric characters and convert to lowercase.

        base_cols <- sub("_Score$", "", all_cols) # ISG/Exh/Cytokine
        sel_base <- sub("_Score$", "", selected) # ISG/ISG Score

        by_base <- all_cols[norm(base_cols) %in% norm(sel_base)]

        # tips: ignore Aa
        if (length(c(direct, by_base)) == 0) {
          cand <- unlist(lapply(sel_base, function(s)
            grep(paste0("^", gsub("[^A-Za-z0-9]+", "", s), "$"),
                 gsub("[^A-Za-z0-9]+", "", base_cols), ignore.case = TRUE, value = FALSE)))
          by_grep <- all_cols[unique(cand)]
        } else {
          by_grep <- character(0)
        }

        unique(c(direct, by_base, by_grep))
      }


      get_grp_col <- function(seu) {
        if ("group" %in% colnames(seu@meta.data)) return("group")
        if ("seurat_clusters" %in% colnames(seu@meta.data)) return("seurat_clusters")
        NULL
      }

      withProgress(message = "Module Scoring in Progress...", value = 0, {

        incProgress(0.2, detail = "Preparing data...")
        seu <- rvSeurat()

        #  group
        if (!"group" %in% colnames(seu@meta.data)) {
          if ("seurat_clusters" %in% colnames(seu@meta.data)) {
            seu$group <- as.character(seu$seurat_clusters)
          } else {
            seu$group <- "sample"
          }
          showNotification("default grouping column：group", type = "warning")
        }

        # define gene sets
        incProgress(0.3, detail = "Loading gene sets...")
        m_df <- msigdbr(species = "Mus musculus", category = "H")
        c7_df <- msigdbr(species = "Mus musculus", category = "C7")

        modules <- list(
          ISG = switch(input$isgList,
                       "H" = unique(c(
                         m_df %>% dplyr::filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>% dplyr::pull(gene_symbol),
                         m_df %>% dplyr::filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% dplyr::pull(gene_symbol)
                       )),
                       "C7" = c7_df %>% dplyr::filter(grepl("INTERFERON", gs_name)) %>% dplyr::pull(gene_symbol) %>% unique(),
                       "Core" = unique(c(
                         msigdbr(species = "Homo sapiens", category = "H") %>%
                           dplyr::filter(gs_name %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE")) %>%
                           dplyr::pull(gene_symbol),
                         m_df %>% dplyr::filter(gs_name %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE")) %>%
                           dplyr::pull(gene_symbol)
                       )),
                       "PRR" = c("DDX58","IFIH1","MAVS","TBK1","IRF3","IRF7"),
                       "Custom" = {
                         if (!is.null(input$customISG)) {
                           read.csv(input$customISG$datapath, header = TRUE)[,1] |> as.character()
                         } else character(0)
                       }
          ),
          Exh = c7_df %>% dplyr::filter(grepl("EXHAUST", gs_name)) %>% dplyr::pull(gene_symbol) %>% unique(),
          Cytokine = m_df %>% dplyr::filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% dplyr::pull(gene_symbol) %>% unique()
        )

        # calculate module score
        map_genes_case_insensitive <- function(query_genes, target_rownames) {
          lut <- setNames(target_rownames, toupper(target_rownames))
          mapped <- lut[unique(toupper(query_genes))]
          unique(stats::na.omit(mapped))
        }


        incProgress(0.5, detail = "Calculating scores...")

        # ensure scoring is performed on RNA
        # has been normalized
        if (!"RNA" %in% names(seu@assays)){
          showNotification("The RNA assay is not present in the object; cannot calculate the module score.",
                           type = "error")
          rvSeurat(seu); return(NULL)
        }
        prev_assay <- DefaultAssay(seu)
        on.exit({ DefaultAssay(seu) <- prev_assay }, add = TRUE)
        DefaultAssay(seu) <- "RNA"

        rna_data <- GetAssayData(seu, assay = "RNA", slot = "data")
        if (inherits(rna_data, "dgCMatrix")) {
          if (Matrix::nnzero(rna_data) == 0) {
            seu <- NormalizeData(seu, assay = "RNA", normalization.method = "LogNormalize")
          }
        } else {
          if (all(rna_data == 0)) {
            seu <- NormalizeData(seu, assay = "RNA", normalization.method = "LogNormalize")
          }
        }


        rna_genes <- rownames(seu[["RNA"]])
        scored_any <- FALSE

        # check if the data exists
        for (mod in names(modules)) {
          genes_present <- map_genes_case_insensitive(modules[[mod]], rna_genes)
          if (length(genes_present) < 5) {
            showNotification(paste0("[", mod, "] Fewer than 5 valid genes; skipped."), type = "warning")
            next
          }

          tryCatch({
            seu <- AddModuleScore(
              seu,
              features = list(genes_present),
              assay = "RNA",
              name = paste0(mod, "_"),
              ctrl = 25,
              nbin = 24,
              seed = 1
            )

            score_col <- grep(paste0("^", mod, "_\\d+$"), colnames(seu@meta.data), value = TRUE)
            if (length(score_col) == 1) {
              seu@meta.data[[paste0(mod, "_Score")]] <- seu@meta.data[[score_col]]
              seu@meta.data[[score_col]] <- NULL
              scored_any <- TRUE
            } else {
              showNotification(paste0("[", mod, "] Scoring output column not found, rename step skipped."),
                               type = "warning")
            }
          }, error = function(e) {
            showNotification(paste0("[", mod, "] AddModuleScore failed：", conditionMessage(e)),
                             type = "error")
          })
        }
        if (!scored_any) {
          showNotification("No module scores could be calculated. Please check if the gene names match the RNA gene names in the object.",
                           type = "error")
          rvSeurat(seu); return(NULL)
        }
        rvSeurat(seu)

        incProgress(0.8, detail = "Generating visualizations...")


        # 2. combined plot
        output$plotModuleCombined <- renderPlot({
          p_list <- list(
            FeaturePlot(rvSeurat(), "ISG_Score", reduction = "wnn.umap", cols = viridis::viridis(10)) +
              theme(legend.position = "none"),
            FeaturePlot(rvSeurat(), "Exh_Score", reduction = "wnn.umap", cols = viridis::viridis(10)) +
              theme(legend.position = "none"),
            FeaturePlot(rvSeurat(), "Cytokine_Score", reduction = "wnn.umap", cols = viridis::viridis(10))
          )
          p <- patchwork::wrap_plots(p_list, ncol = 2)
          rvPlots$mod <- p
          p
        }, height = 400)

        # ISG
        output$plotISG <- renderPlot({
          p <- FeaturePlot(rvSeurat(), features = "ISG_Score", reduction = "wnn.umap",
                           cols = viridis::viridis(10), pt.size = 0.5) +
            ggtitle("ISG Module Score") + theme_minimal()
          rvPlots$isg <- p
          p
        }, height = 400)

        # Exh
        output$plotExh <- renderPlot({
          p <- FeaturePlot(rvSeurat(), features = "Exh_Score", reduction = "wnn.umap",
                           cols = viridis::viridis(10), pt.size = 0.5) +
            ggtitle("Exhaustion Module Score") + theme_minimal()
          rvPlots$exh <- p
          p
        }, height = 400)

        # Cytokine
        output$plotCytokine <- renderPlot({
          p <- FeaturePlot(rvSeurat(), features = "Cytokine_Score", reduction = "wnn.umap",
                           cols = viridis::viridis(10), pt.size = 0.5) +
            ggtitle("Inflammatory Module Score") + theme_minimal()
          rvPlots$cyto <- p
          p
        }, height = 400)

        # heatmap(group mean)
        output$plotMSHeat <- renderPlot({
          seu <- rvSeurat()
          grp <- if ("group" %in% colnames(seu@meta.data)) "group" else "seurat_clusters"
          feats <- get_ms_cols(seu, input$ms_modules)
          if (length(feats) == 0) feats <- grep("_Score$", colnames(seu@meta.data), value = TRUE)

          mean_scores <- seu@meta.data |>
            dplyr::group_by(.data[[grp]]) |>
            dplyr::summarise(dplyr::across(dplyr::all_of(feats), mean, na.rm = TRUE)) |>
            as.data.frame()

          mat <- as.matrix(mean_scores[,-1, drop = FALSE])
          rownames(mat) <- mean_scores[[1]]

          rvPlots$msheat_mat <- mat  # matrix

          pheatmap::pheatmap(mat, color = viridis::viridis(100),
                             cluster_rows = TRUE, cluster_cols = TRUE,
                             scale = "column", main = "Module Scores Heatmap",
                             angle_col = 45, display_numbers = TRUE, number_format = "%.2f")
        }, height = 600)

        # dot plot（Group × Module）
        output$plotMSDot <- renderPlot({
          feats <- get_ms_cols(rvSeurat()); grp <- get_grp_col(rvSeurat())
          validate(need(length(feats) > 0, "No module scores available"),
                   need(!is.null(grp), "Grouping column not found."))
          p <- Seurat::DotPlot(rvSeurat(), features = feats, group.by = grp) +
            scale_color_viridis_c() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
          rvPlots$msdot <- p  #save rvPlots
          p
        })

        # RidgePlot
        output$plotMSRidge <- renderPlot({
          feats <- get_ms_cols(seu); grp <- get_grp_col(seu)
          validate(need(length(feats) > 0, "No module scores available"),
                   need(!is.null(grp), "Grouping column not found."))
          p <- Seurat::RidgePlot(seu, features = feats, group.by = grp)
          rvPlots$msridge <- p      # save
          p
        }, height = 600)

        # Bar (Mean ± SEM)
        output$plotMSBar <- renderPlot({
          feats <- get_ms_cols(seu); grp <- get_grp_col(seu)
          validate(need(!is.null(grp), "Grouping column not found."),
                   need(length(feats) > 0, "No module scores available"))

          df <- seu@meta.data[, c(grp, feats), drop = FALSE] %>%
            tidyr::pivot_longer(cols = dplyr::all_of(feats),
                                names_to = "module",
                                values_to = "score") %>%
            dplyr::group_by(.data[[grp]], module) %>%
            dplyr::summarise(mean = mean(score, na.rm = TRUE),
                             sem  = stats::sd(score, na.rm = TRUE) / sqrt(dplyr::n()),
                             .groups = "drop")

          p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[grp]], y = mean, fill = module)) +
            ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.7) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sem, ymax = mean + sem),
                                   position = ggplot2::position_dodge(0.8), width = 0.2) +
            ggplot2::labs(title = "Module Scores by Group", x = "", y = "Score (Mean ± SEM)") +
            ggplot2::theme_minimal()
          rvPlots$msbar <- p  # save
          p
        }, height = 500)

        # ISG-NP correlation plot
        output$plotISGvsNP <- renderPlot({
          req(rvSeurat())
          seu <- rvSeurat()

          validate(
            need("ADT" %in% names(seu@assays), "ADT assay no found"),
            need(any(grepl("^NP", rownames(seu[["ADT"]]))), "No features starting with 'NP' were found in ADT征"),
            need("ISG_Score" %in% colnames(seu@meta.data), "metadata lose ISG_Score colume")
          )

          # ADT: features x cells -> cells x features
          adt_mat <- t(as.matrix(GetAssayData(seu, assay = "ADT", slot = "data")))
          np_feature <- grep("^NP", colnames(adt_mat), value = TRUE)[1]
          validate(need(!is.na(np_feature), "No NP-related ADT features found."))

          # align by intersection of cell names
          cells <- intersect(rownames(seu@meta.data), rownames(adt_mat))
          validate(need(length(cells) > 0, "No cells available for alignment."))

          grp <- get_grp_col(seu)
          GroupVec <- if (!is.null(grp)) seu@meta.data[cells, grp, drop = TRUE] else factor("All Cells")

          plot_data <- data.frame(
            ISG_Score = seu@meta.data[cells, "ISG_Score", drop = TRUE],
            NP_ADT = adt_mat[cells, np_feature],
            Group = GroupVec,
            row.names = cells
          )

          cor_value <- round(stats::cor(plot_data$ISG_Score, plot_data$NP_ADT,
                                        use = "pairwise.complete.obs"), 3)

          ggplot(plot_data, aes(x = NP_ADT, y = ISG_Score, color = Group)) +
            geom_point(alpha = 0.6, size = 2) +
            geom_smooth(method = "lm", se = TRUE, color = "darkred") +
            scale_color_viridis_d() +
            labs(
              x = paste0(np_feature, " (ADT Expression)"),
              y = "ISG Module Score",
              title = paste("ISG-NP Correlation (r =", cor_value, ")"),
              subtitle = "Linear regression fit with 95% confidence interval"
            ) +
            theme_minimal(base_size = 14) +
            theme(legend.position = "bottom",
                  plot.title = element_text(face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5))
        }, height = 500)

        # download
        output$dlISGvsNP <- downloadHandler(
          filename = function() paste0("ISG_NP_Correlation_", Sys.Date(), ".png"),
          content = function(file) {
            req(rvSeurat())
            seu <- rvSeurat()

            adt_mat <- t(as.matrix(GetAssayData(seu, assay = "ADT", slot = "data")))
            np_feature <- grep("^NP", colnames(adt_mat), value = TRUE)[1]
            cells <- intersect(rownames(seu@meta.data), rownames(adt_mat))
            grp <- get_grp_col(seu)
            GroupVec <- if (!is.null(grp)) seu@meta.data[cells, grp, drop = TRUE] else factor("All Cells")

            plot_data <- data.frame(
              ISG_Score = seu@meta.data[cells, "ISG_Score", drop = TRUE],
              NP_ADT = adt_mat[cells, np_feature],
              Group = GroupVec,
              row.names = cells
            )

            p <- ggplot(plot_data, aes(x = NP_ADT, y = ISG_Score, color = Group)) +
              geom_point(alpha = 0.6, size = 2) +
              geom_smooth(method = "lm", se = TRUE, color = "darkred") +
              labs(x = paste0(np_feature, " Expression"), y = "ISG Module Score")

            ggsave(file, plot = p, width = 8, height = 6, dpi = 300, bg = "white")
          }
        )
      })
    })



    ### Trajectory Inference

    # build a monocle3 cell_data_set (cds) from a Seurat object
    .make_cds_from_seurat <- function(seu) {
      expr <- NULL
      # Prefer counts if available; otherwise fall back to data

      if ("RNA" %in% names(seu@assays)) {
        expr <- tryCatch(GetAssayData(seu, assay = "RNA", layer = "counts"),
                         error = function(e) NULL)
      }
      if (is.null(expr) && "SCT" %in% names(seu@assays)) {
        expr <- tryCatch(GetAssayData(seu, assay = "SCT", layer = "counts"),
                         error = function(e) NULL)
      }
      if (is.null(expr)) {
        assay_use <- if ("RNA" %in% names(seu@assays)) "RNA" else "SCT"
        expr <- GetAssayData(seu, assay = assay_use, layer = "data")
      }

      # Prefer SeuratWrappers

      if (requireNamespace("SeuratWrappers", quietly = TRUE)) {
        cds <- SeuratWrappers::as.cell_data_set(seu)
        return(cds)
      } else {
        gene_meta <- data.frame(gene_short_name = rownames(expr), row.names = rownames(expr))
        cell_meta <- seu@meta.data
        cds <- monocle3::new_cell_data_set(expr,
                                           cell_metadata = cell_meta,
                                           gene_metadata = gene_meta)
        return(cds)
      }
    }

    # largest cluster label among seurat_clusters

    .largest_cluster <- function(seu) {
      cl <- as.character(seu$seurat_clusters)
      names(sort(table(cl), decreasing = TRUE))[1]
    }

    # UI: manual root cluster (default to largest cluster)

    output$ti_root_cluster_ui <- renderUI({
      req(rvSeurat())
      seu <- rvSeurat()
      cl_vec <- as.character(seu$seurat_clusters)
      choices <- sort(unique(cl_vec))
      default_cl <- .largest_cluster(seu)

      selectInput(
        "ti_root_cluster",
        "Root cluster (manual):",
        choices  = choices,
        selected = default_cl
      )
    })

    # ADT feature dropdown for coloring & ADT-root mode

    output$ti_adt_feature_ui <- renderUI({
      req(rvSeurat())
      seu <- rvSeurat()
      if (!"ADT" %in% names(seu@assays)) {
        return(helpText("No ADT assay detected. TI can run without ADT coloring/root."))
      } else {
        np_default <- grep("^NP", rownames(seu[["ADT"]]), value = TRUE)[1]
        selectInput("ti_adt_feature", "ADT feature to overlay / select root:",
                    choices  = rownames(seu[["ADT"]]),
                    selected = ifelse(length(np_default) > 0, np_default, rownames(seu[["ADT"]])[1]))
      }
    })

    #  pseudotime ridge by group

    output$plotTI_ptRidge <- renderPlot({
      cds <- req(rvCDS())
      seu <- req(rvSeurat())

      # Ensure cells are aligned

      cells <- intersect(colnames(cds), rownames(seu@meta.data))
      validate(need(length(cells) > 0, "No overlapping cells between cds and Seurat."))

      pt_vec <- tryCatch(monocle3::pseudotime(cds), error = function(e) rep(NA_real_, ncol(cds)))
      pt_vec <- as.numeric(pt_vec[cells])

      grp <- if ("infection_status" %in% colnames(seu@meta.data)) {
        seu@meta.data[cells, "infection_status", drop = TRUE]
      } else {
        seu@meta.data[cells, "seurat_clusters",  drop = TRUE]
      }

      df <- data.frame(Pseudotime = pt_vec, Group = grp)
      df <- df[is.finite(df$Pseudotime) & !is.na(df$Group), , drop = FALSE]

      p <- ggplot2::ggplot(df, ggplot2::aes(Pseudotime, y = Group, fill = Group)) +
        ggridges::geom_density_ridges(alpha = 0.8, scale = 1.1, rel_min_height = 0.01) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(title = "Pseudotime by Group (Ridge)", y = "")

      rvPlots$ti_ptRidge <- p
      p
    })

    # Main TI pipeline

    observeEvent(input$btnTI, {
      req(rvSeurat())
      status_ti("Running Monocle3 TI…")
      withProgress(message = "Monocle3 Trajectory Inference…", value = 0, {
        incProgress(0.1, detail = "Build cell_data_set from Seurat")
        seu <- rvSeurat()
        cds <- .make_cds_from_seurat(seu)

        adt_col_name <- NULL
        if ("ADT" %in% names(seu@assays) && !is.null(input$ti_adt_feature)) {
          adt_feats <- rownames(seu[["ADT"]])
          selected_adt_feature <- input$ti_adt_feature

          if (selected_adt_feature %in% adt_feats) {

            adt_mat <- GetAssayData(seu, assay = "ADT", layer = "data")


            common_cells <- intersect(colnames(cds), colnames(adt_mat))
            if (length(common_cells) > 0) {

              adt_col_name <- paste0("ADT_", selected_adt_feature)


              cd <- SummarizedExperiment::colData(cds)
              adt_values <- rep(NA, ncol(cds))
              names(adt_values) <- colnames(cds)
              adt_values[common_cells] <- as.numeric(adt_mat[selected_adt_feature, common_cells])
              cd[[adt_col_name]] <- adt_values
              SummarizedExperiment::colData(cds) <- cd
            }
          }
        }

        incProgress(0.3, detail = "preprocess → reduce_dimension → cluster_cells")
        cds <- monocle3::preprocess_cds(cds, num_dim = 30, method = "PCA")
        cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
        cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")

        # Select root cells
        incProgress(0.6, detail = "learn_graph & order_cells")
        root_cells <- character(0)

        if (identical(input$ti_root_mode, "manual")) {
          # Manual mode
          cl_vec <- as.character(seu$seurat_clusters)
          default_cl <- .largest_cluster(seu)
          chosen_cl  <- if (!is.null(input$ti_root_cluster) && nzchar(input$ti_root_cluster)) {
            input$ti_root_cluster
          } else {
            default_cl
          }
          cand <- rownames(seu@meta.data)[cl_vec == chosen_cl]
          root_cells <- intersect(cand, colnames(cds))

        } else if (identical(input$ti_root_mode, "adt") &&
                   !is.null(adt_col_name) &&
                   adt_col_name %in% colnames(SummarizedExperiment::colData(cds))) {


          vals <- SummarizedExperiment::colData(cds)[[adt_col_name]]
          # remove NA
          valid_vals <- vals[!is.na(vals)]
          if (length(valid_vals) > 0) {
            q <- tryCatch(stats::quantile(valid_vals, probs = input$ti_adt_quantile, na.rm = TRUE),
                          error = function(e) NA_real_)
            if (is.finite(q)) {
              root_cells <- names(vals)[vals >= q & !is.na(vals)]
            }
            if (length(root_cells) < 50) {
              # select ADT highest cell
              valid_indices <- which(!is.na(vals))
              ord <- order(vals[valid_indices], decreasing = TRUE)
              root_cells <- names(vals)(valid_indices[head(ord, min(200, length(valid_indices)))])
            }
          }
        }

        if (length(root_cells) == 0 && "seurat_clusters" %in% colnames(seu@meta.data)) {
          cl_vec <- as.character(seu@meta.data[colnames(cds), "seurat_clusters", drop = TRUE])
          top_cl <- names(sort(table(cl_vec), decreasing = TRUE))[1]
          cand   <- names(cl_vec)[cl_vec == top_cl]
          root_cells <- head(cand, 200)
        }

        if (length(root_cells) == 0) {
          root_cells <- head(colnames(cds), min(200, ncol(cds)))
        }

        cds <- monocle3::learn_graph(cds, use_partition = TRUE)
        cds <- monocle3::order_cells(cds, root_cells = root_cells)

        rvCDS(cds)

        # UMAP + learned graph colored by pseudotime
        p_cells <- monocle3::plot_cells(
          cds,
          color_cells_by = "pseudotime",
          label_groups_by_cluster = FALSE,
          label_leaves = TRUE,
          label_branch_points = TRUE,
          graph_label_size = 2
        )
        output$plotTI_cells <- renderPlot(p_cells)
        rvPlots$ti_cells <- p_cells

        # UMAP colored by ADT feature
        if (!is.null(adt_col_name) && adt_col_name %in% colnames(SummarizedExperiment::colData(cds))) {
          # check
          adt_values <- SummarizedExperiment::colData(cds)[[adt_col_name]]
          if (any(!is.na(adt_values))) {
            p_adt <- monocle3::plot_cells(
              cds,
              color_cells_by = adt_col_name,
              label_groups_by_cluster = FALSE,
              label_leaves = FALSE,
              label_branch_points = FALSE,
              graph_label_size = 2
            )
          } else {
            p_adt <- ggplot2::ggplot() + ggplot2::theme_void() +
              ggplot2::ggtitle("No valid ADT values available for coloring.")
          }
        } else {
          p_adt <- ggplot2::ggplot() + ggplot2::theme_void() +
            ggplot2::ggtitle("No ADT features available for coloring.")
        }
        output$plotTI_ADT <- renderPlot(p_adt)
        rvPlots$ti_adt <- p_adt

        # ADT vs Pseudotime with GAM smooth
        if (!is.null(adt_col_name) &&
            adt_col_name %in% colnames(SummarizedExperiment::colData(cds))) {

          adt_values <- SummarizedExperiment::colData(cds)[[adt_col_name]]
          pseudotime_vals <- monocle3::pseudotime(cds, reduction_method = "UMAP")


          common_cells <- intersect(names(adt_values), names(pseudotime_vals))
          if (length(common_cells) > 10) {
            df <- data.frame(
              pseudotime = pseudotime_vals[common_cells],
              ADT = adt_values[common_cells],
              row.names = common_cells
            )
            df <- df[is.finite(df$pseudotime) & is.finite(df$ADT), ]

            if (nrow(df) > 10) {
              p_corr <- ggplot2::ggplot(df, ggplot2::aes(pseudotime, ADT)) +
                ggplot2::geom_point(alpha = 0.3, size = 0.6) +
                ggplot2::geom_smooth(method = "gam", formula = y ~ mgcv::s(x, k = 5), se = TRUE) +
                ggplot2::theme_minimal(base_size = 13) +
                ggplot2::labs(x = "Pseudotime", y = gsub("^ADT_", "", adt_col_name),
                              title = "ADT vs Pseudotime (GAM smooth)")
            } else {
              p_corr <- ggplot2::ggplot() + ggplot2::theme_void() +
                ggplot2::ggtitle("Insufficient data for correlation analysis.")
            }
          } else {
            p_corr <- ggplot2::ggplot() + ggplot2::theme_void() +
              ggplot2::ggtitle("Not enough overlapping cells for correlation analysis.")
          }
        } else {
          p_corr <- ggplot2::ggplot() + ggplot2::theme_void() +
            ggplot2::ggtitle("No ADT features available for correlation analysis.")
        }
        output$plotTI_ADTvsPT <- renderPlot(p_corr)
        rvPlots$ti_corr <- p_corr

        incProgress(1, detail = "Completed")
        status_ti("Trajectory Inference Completed")
      })
    })





    ### Differential Expression Analysis ###
    observeEvent(input$btnDE, {
      req(rvSeurat())
      status_de("Running Differential Expression Analysis...")

      #close future
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::sequential)

      #max
      old_max <- getOption("future.globals.maxSize")
      on.exit(options(future.globals.maxSize = old_max), add = TRUE)
      options(future.globals.maxSize = 64 * 1024^3)  # 64 GB

      withProgress(message = "Differential Expression Analysis in Progress...",
                   value = 0, {
                     # Is there an existing infection grouping
                     incProgress(0.05, detail = "check infection status…")
                     seu <- rvSeurat()
                     validate(
                       need("infection_status" %in% colnames(seu@meta.data),
                            "Please run the Infection Status Classification module first")
                     )

                     # Select the assay and slot for DE.
                     # Record assay for heatmap usage.
                     incProgress(0.15, detail = "Select assay & slot…")
                     assays <- names(seu@assays)
                     assayUse <- if ("SCT" %in% assays) "SCT" else "RNA"
                     slotUse <- "data"
                     DefaultAssay(seu) <- assayUse
                     Idents(seu) <- seu$infection_status
                     deAssay(assayUse)

                     # defferential genes
                     incProgress(0.25, detail = "Step 1/5: FindMarkers()…")
                     de_res <- FindMarkers(
                       object = seu,
                       ident.1 = "infected",
                       ident.2 = "uninfected",
                       assay = assayUse,
                       slot  = slotUse,
                       logfc.threshold = 0.25,
                       min.pct = 0.1
                     )
                     de_df <- as.data.frame(de_res) %>% tibble::rownames_to_column("gene")
                     deResults(de_df)

                     # MA data
                     incProgress(0.45, detail = "Step 2/5: construct MA data…")
                     avg_expr <- AverageExpression(
                       object = seu,
                       assays = assayUse,
                       slot = slotUse,
                       group.by = "infection_status"
                     )[[assayUse]]
                     ma_df <- avg_expr %>%
                       as.data.frame() %>%
                       tibble::rownames_to_column("gene") %>%
                       dplyr::rename(infected = infected, uninfected = uninfected) %>%
                       dplyr::mutate(
                         A = log2((infected + uninfected) / 2 + 1),
                         M = log2(infected + 1) - log2(uninfected + 1)
                       )
                     maData(ma_df)

                     # volcano data
                     incProgress(0.65, detail = "Step 3/5: construct Volcano plot data…")
                     volcano_df <- de_df %>%
                       dplyr::mutate(
                         negLog10P = -log10(p_val_adj + 1e-300),
                         sig = (p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
                       )
                     volcanoData(volcano_df)

                     # heatmap top10 genes
                     incProgress(0.80, detail = "Step 4/5: prepare Top10 genes for Heatmap…")

                     genes_in_assay <- rownames(seu[[assayUse]])
                     top10_all <- de_df %>% dplyr::arrange(p_val_adj) %>% dplyr::pull(gene)
                     top10 <- intersect(top10_all, genes_in_assay) %>% head(10)
                     if (length(top10) == 0) {
                       showNotification("Top DE genes not present in selected assay; heatmap will be empty.",
                                        type = "error")
                     } else {

                       seu <- tryCatch(
                         {
                           ScaleData(seu, assay = assayUse, features = top10, verbose = FALSE)
                         },
                         error = function(e) seu
                       )
                     }
                     heatmapData(top10)
                     rvSeurat(seu)

                     # GO/KEGG enrichment
                     incProgress(0.90, detail = "Step 5/5: GO/KEGG enrichment…")

                     sig_genes <- intersect(
                       de_df %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>% dplyr::pull(gene),
                       rownames(seu)
                     )

                     if (length(sig_genes) > 0) {
                       entrez_ids <- AnnotationDbi::mapIds(
                         org.Mm.eg.db,
                         keys = sig_genes,
                         column = "ENTREZID",
                         keytype  = "SYMBOL",
                         multiVals = "first"
                       ) %>% stats::na.omit()
                       ego_bp <- clusterProfiler::enrichGO(
                         gene  = entrez_ids,
                         OrgDb = org.Mm.eg.db,
                         keyType = "ENTREZID",
                         ont  = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         readable  = TRUE
                       )
                       ekegg <- clusterProfiler::enrichKEGG(
                         gene  = entrez_ids,
                         organism  = "mmu",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05
                       ) %>% clusterProfiler::setReadable(OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
                     } else {
                       ego_bp <- NULL; ekegg <- NULL
                     }
                     enrichResults(list(go = ego_bp, kegg = ekegg))

                     incProgress(1, detail = "finished")
                     status_de("Differential Expression Analysis Completed")
                   })
    })

    # result render
    # DE table
    output$deTable <- DT::renderDataTable({
      req(deResults())
      deResults()
    })

    # MA Plot
    output$plotMA <- renderPlot({
      df <- req(maData())
      top50  <- head(deResults()$gene, 50)
      p <- ggplot(df, aes(x = A, y = M)) +
        geom_point(alpha = 0.3, size = 1) +
        geom_point(data = df[df$gene %in% top50, ], color = "red", size = 1.5) +
        theme_minimal() + labs(title = "MA Plot",
                               x = "Average log2(Expression+1)",
                               y = "Log2 Fold Change")
      rvPlots$ma <- p; p
    })

    # volcano
    output$plotVolcano <- renderPlot({
      df <- req(volcanoData())
      p <- ggplot(df, aes(x = avg_log2FC, y = negLog10P, color = sig)) +
        geom_point(alpha = 0.5) +
        scale_color_manual(values = c("grey70", "red")) +
        theme_minimal() +
        labs(title = "Volcano Plot",
             x = "avg_log2FC",
             y = "-log10(adj.p-value)")
      rvPlots$volcano <- p; p
    })

    # Top10 DEG heatmap
    output$plotHeatmap <- renderPlot({
      genes <- req(heatmapData())
      seu   <- rvSeurat()
      assayHM <- if (!is.null(deAssay())) deAssay() else DefaultAssay(seu)
      feats   <- intersect(genes, rownames(seu[[assayHM]]))
      validate(need(length(feats) > 0,
                    "No requested features found in the selected assay."))

      p <- DoHeatmap(
        object = seu,
        features = feats,
        group.by = "infection_status",
        assay = assayHM,
        slot  = "data"
      ) + NoLegend() + ggtitle("Top10 DEGs Heatmap")

      rvPlots$heatmap <- p
      p
    })

    # GO+KEGG enrichment
    output$plotGOKEGG <- renderPlot({
      enr <- req(enrichResults())
      if (is.null(enr$go) || is.null(enr$kegg)) {
        plot.new(); text(0.5, 0.5, "No Significant DEG Enrichment Results", cex = 1.2)
      } else {
        p <- dotplot(enr$go, showCategory = 20) | dotplot(enr$kegg, showCategory = 20)
        rvPlots$go <- p; p
      }
    })



    # download handlers
    output$dlPreviewMeta <- downloadHandler(
      filename = function() "metadata_full.csv",
      content = function(file) {
        req(rvSeurat())
        write.csv(rvSeurat()@meta.data, file, row.names = TRUE)
      }
    )
    output$dlPlot1 <- downloadHandler(
      filename = function() "WNN_UMAP_Clusters.png",
      content = function(file) {
        ggsave(file, plot = rvPlots$p1, width = 6, height = 5)
      }
    )
    output$dlPlot2 <- downloadHandler(
      filename = function() "SCT_weight_UMAP.png",
      content = function(file) {
        ggsave(file, plot = rvPlots$p2, width = 6, height = 5)
      }
    )
    output$dlPlot3 <- downloadHandler(
      filename = function() "ADT_weight_UMAP.png",
      content = function(file) {
        ggsave(file, plot = rvPlots$p3, width = 6, height = 5)
      }
    )

    ## Download Infection Status Classification Plot ##
    output$dlPlotViolin <- downloadHandler(
      filename = function() "NP_violin.png",
      content = function(file) ggsave(file, plot = rvPlots$hist, width = 6, height = 5)
    )
    output$dlPlotBar <- downloadHandler(
      filename = function() "Infection_bar.png",
      content = function(file) ggsave(file, plot = rvPlots$bar, width = 6, height = 5)
    )
    output$dlPlotInfection <- downloadHandler(
      filename = function() "Infection_Status_UMAP.png",
      content = function(file) ggsave(file, plot = rvPlots$inf, width = 6, height = 5)
    )

    ## Download SingleR Annotation Plot ##
    output$dlPlotSingleR <- downloadHandler(
      filename = function() "SingleR_Annotation_UMAP.png",
      content = function(file) {
        ggsave(file, plot = rvPlots$sr, width = 6, height = 5)
      }
    )
    output$dlPlotADTmarkers <- downloadHandler(
      filename = function() "ADT_Markers_UMAP.png",
      content = function(file) {
        ggsave(file, plot = rvPlots$adt, width = 8, height = 6)
      }
    )
    output$dlPlotSCTmarkers <- downloadHandler(
      filename = function() "SCT_VariableFeatures_UMAP.png",
      content = function(file) {
        ggsave(file, plot = rvPlots$sct, width = 8, height = 6)
      }
    )

    output$dlPlotMarkers <- downloadHandler(
      filename = function() "Custom_Marker_UMAP.png",
      content = function(file) ggsave(file, plot = rvPlots$markers,
                                      width = 10, height = 7, dpi = 300, bg = "white")
    )
    output$dlPlotMarkerVln <- downloadHandler(
      filename = function() "Custom_Marker_Violin.png",
      content = function(file) ggsave(file, plot = rvPlots$markerVln,
                                      width = 8, height = 10, dpi = 300, bg = "white")
    )

    ## Download Module Scoring Plot ##
    output$dlPlotModuleCombined <- downloadHandler(
      filename = function() "ModuleScores_Combined_UMAP.png",
      content = function(file) {
        ggsave(file, plot = rvPlots$mod, width = 12, height = 8)
      }
    )
    output$dlPlotISG <- downloadHandler(
      filename = function() "ISG_ModuleScore_UMAP.png",
      content = function(file) {
        ggsave(file, plot = rvPlots$isg, width = 6, height = 5)
      }
    )

    output$dlPlotExh <- downloadHandler(
      filename = function() "Exhaustion_ModuleScore_UMAP.png",
      content = function(file) {
        ggsave(file, plot = rvPlots$exh, width = 6, height = 5)
      }
    )
    output$dlPlotCytokine <- downloadHandler(
      filename = function() "Cytokine_ModuleScore_UMAP.png",
      content = function(file) {
        ggsave(file, plot = rvPlots$cyto, width = 6, height = 5)
      }
    )

    output$dlMSRidge <- downloadHandler(
      filename = function() "MS_RidgePlot.png",
      content = function(file) ggsave(file, plot = rvPlots$msridge,
                                      width = 8, height = 6, dpi = 300, bg = "white")
    )
    output$dlMSBar <- downloadHandler(
      filename = function() "MS_BarPlot.png",
      content = function(file) ggsave(file, plot = rvPlots$msbar,
                                      width = 8, height = 6, dpi = 300, bg = "white")
    )


    output$dlDECSV <- downloadHandler(
      filename = function() paste0("DEG.csv"),
      content = function(f) write.csv(deResults(), f, row.names=FALSE)
    )
    output$dlGOCSV <- downloadHandler(
      filename = function() "GO_KEGG.csv",
      content = function(f) {
        enr <- enrichResults()
        go_df <- if (!is.null(enr$go)) cbind(source="GO" , as.data.frame(enr$go))   else NULL
        kegg_df <- if (!is.null(enr$kegg)) cbind(source="KEGG", as.data.frame(enr$kegg)) else NULL
        out <- dplyr::bind_rows(go_df, kegg_df)
        if (is.null(out) || nrow(out)==0) out <- data.frame(message="No significant terms")
        write.csv(out, f, row.names = FALSE)
      }
    )
    output$dlTI_ptRidge <- downloadHandler(
      filename = function() "TI_Pseudotime_Ridge.png",
      content = function(file) ggsave(file, plot = rvPlots$ti_ptRidge,
                                      width = 8, height = 6, dpi = 300, bg = "white")
    )

    output$dlPlotMA <- downloadHandler(
      filename = function() "DE_MA_plot.png",
      content = function(file) ggsave(file, plot = rvPlots$ma,
                                      width = 7, height = 5, dpi = 300, bg = "white")
    )
    output$dlPlotVolcano <- downloadHandler(
      filename = function() "DE_Volcano_plot.png",
      content = function(file) ggsave(file, plot = rvPlots$volcano,
                                      width = 7, height = 5, dpi = 300, bg = "white")
    )
    output$dlPlotHeatmap <- downloadHandler(
      filename = function() "DE_Top10_Heatmap.png",
      content = function(file) ggsave(file, plot = rvPlots$heatmap,
                                      width = 8, height = 5, dpi = 300, bg = "white")
    )
    output$dlPlotGOKEGG <- downloadHandler(
      filename = function() "DE_Enrichment.png",
      content = function(file) ggsave(file, plot = rvPlots$go,
                                      width = 10, height = 6, dpi = 300, bg = "white")
    )

    #
    output$dlMSDot <- downloadHandler(
      filename = function() "MS_DotPlot.png",
      content = function(file) {
        req(rvPlots$msdot)
        ggsave(file, plot = rvPlots$msdot, width = 8, height = 6, dpi = 300, bg = "white")
      }
    )
    output$dlMSHeat <- downloadHandler(
      filename = function() "MS_Heatmap.png",
      content = function(file) {
        req(rvPlots$msheat_mat)
        png(file, width = 1800, height = 1400, res = 300)
        pheatmap::pheatmap(rvPlots$msheat_mat, color = viridis::viridis(100),
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           scale = "column", main = "Module Scores Heatmap",
                           angle_col = 45, display_numbers = TRUE, number_format = "%.2f")
        dev.off()
      }
    )


    # Download TI Plot
    output$dlTIcells <- downloadHandler(
      filename = function() "TI_UMAP_Pseudotime_Graph.png",
      content = function(file) ggsave(file, plot = rvPlots$ti_cells,
                                      width = 7, height = 6, dpi = 300)
    )
    output$dlTIadt <- downloadHandler(
      filename = function() "TI_UMAP_ADT.png",
      content = function(file) ggsave(file, plot = rvPlots$ti_adt,
                                      width = 7, height = 6, dpi = 300)
    )
    output$dlTIcorr <- downloadHandler(
      filename = function() "ADT_vs_Pseudotime.png",
      content = function(file) ggsave(file, plot = rvPlots$ti_corr,
                                      width = 7, height = 6, dpi = 300)
    )
}
