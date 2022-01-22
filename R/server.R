#' tinySEV.server
#'
#' @param objects A named list of (paths to)
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} objects.
#' @param uploadMaxSize The maximum upload size. Set to zero to disable upload.
#' @param maxPlot The maximum number of features to allow for plotting heatmaps.
#' @param genelists An optional named list of genes/features which will be 
#'   flagged in the gene tab.
#'
#' @return A shiny server function.
#' @export
#' @import shiny ggplot2 SummarizedExperiment SEtools waiter
#' @importFrom plotly ggplotly renderPlotly event_data
#' @importFrom shinydashboard updateTabItems
#' @importFrom shinyjs showElement hideElement
#' @importFrom DT datatable renderDT
#' @importFrom ComplexHeatmap draw
#' @importFrom S4Vectors metadata
tinySEV.server <- function(objects=NULL, uploadMaxSize=50*1024^2, maxPlot=500,
                           genelists=list(), genelistTab=length(genelists)>0){
  options(shiny.maxRequestSize=uploadMaxSize)
  genelists <- genelists[lengths(genelists)>0]

  if(!is.null(objects) && is.null(names(objects))){
    names(objects) <- paste("Object ",seq_along(objects))
    x <- sapply(objects, FUN=function(x){
      if(is.character(x))
        return(gsub("\\.SE\\.rds$|\\.rds$", "", basename(x), ignore.case=TRUE))
      return(NULL)
    })
    x[which(is.null(x))] <- names(objects)[which(is.null(x))]
    names(objects) <- make.unique(x, sep=" ")
  }

  grepGene <- function(x,g){
    if(!is.character(x)){
      g <- grepGene(row.names(x), g)
      return(x[g,drop=FALSE])
    }
    if(all(g %in% x)) return(g)
    g <- paste0("^",g,"\\.|^",g,"$|\\.",g,"$")
    g <- lapply(g,FUN=function(i) grep(i, x, value=TRUE, ignore.case=TRUE))
    return(unique(unlist(g)))
  }

  getDef <- function(se,var){
    if(length(var)>1){
      y <- unlist(lapply(var, FUN=function(x) getDef(se,x)))
      y <- y[!sapply(y,is.null)]
      if(length(y)==0) return(NULL)
      return(y)
    }
    if(is.null(se@metadata$default_view[[var]])) return(NULL)
    se@metadata$default_view[[var]]
  }

  function(input, output, session) {

    output$uploadMenu <- renderUI({
      if(!(uploadMaxSize>0)) return(NULL)
      menuSubItem("Upload object", tabName="tab_fileinput")
    })
    output$menu_genelist <- renderUI({
      if(genelistTab && length(genelists)==0) return(NULL)
      menuItem("Gene lists", tabName="tab_genelists")
    })
    output$maxGenes <- renderText(maxPlot)

    SEinit <- function(x){
      if(is.null(assayNames(x)))
        assayNames(x) <- paste0("assay",1:length(assays(x)))
      if(ncol(rowData(x))==0) rowData(x)$name <- row.names(x)
      updateSelectizeInput(session, "gene_input", selected=isolate(selGene()),
                           choices=sort(unique(row.names(x))), server=TRUE)
      colvars <- colnames(colData(x))
      updateSelectInput(session, "assay_input", choices=assayNames(x),
                        selected=getDef(x, "assay"))
      updateSelectInput(session, "assay_input2", choices=assayNames(x),
                        selected=getDef(x, "assay"))
      updateSelectizeInput(session, "hm_order", choices=colvars)
      updateSelectizeInput(session, "hm_anno", choices=colvars,
                           selected=getDef(x, c("groupvar","colvar")))
      updateSelectizeInput(session, "hm_gaps", choices=colvars,
                           selected=getDef(x, "gridvar"))
      updateSelectInput(session, "select_groupvar", choices=colvars,
                        selected=getDef(x, "groupvar"))
      updateSelectInput(session, "select_colorvar", choices=colvars,
                        selected=getDef(x, "colvar"))
      updateSelectInput(session, "select_gridvars", choices=colvars,
                        selected=getDef(x, "gridvar"))
      if(!is.null(getDef(x, "assay")))
        updateCheckboxInput(session, "hm_scale",
                            value=!grepl("FC$",getDef(x, "assay")))
      x
    }

    SEs <- reactiveValues()
    for(nn in names(objects)) SEs[[nn]] <- objects[[nn]]
    updateSelectInput(session, "object", choices=names(objects))
    if(length(genelists)>0)
      updateSelectInput(session, "genelist_input", choices=names(genelists))

    SE <- reactive({
      if(is.null(input$object) || input$object=="" ||
         is.null(SEs[[input$object]])) return(NULL)
      if(is.character(SEs[[input$object]]))
        SEs[[input$object]] <- readRDS(SEs[[input$object]])
      SEinit(SEs[[input$object]])
    })
    
    observeEvent(input$file, {
      tryCatch({
        if(!is.null(input$file)){
          x <- readRDS(input$file$datapath)
          if(is(x,"SummarizedExperiment")){
            SEname <- gsub("\\.SE\\.rds$|\\.rds$","",
                           basename(input$file$name), ignore.case=TRUE)
            SEs[[SEname]] <- x
            updateSelectInput(session, "object", selected=SEname,
                              choices=union(names(objects), names(SEs)))
          }else{
            stop("The object is not a SummarizedExperiment!")
          }
        }}, error=function(e){
          showModal(modalDialog(easyClose=TRUE, title="Error with upload",
            "The file was not recognized. Are you sure that it is a R .rds file?",
            tags$pre(e)))
        })
    })


    output$SEout <- renderPrint({
      if(is.null(SE())) return(NULL)
      print(SE())
      md <- metadata(SE())
      for(f in c("title","name","source","description")){}
      if(length(md <- md[intersect(c("title","name","source","description"),
                                   names(md))])>0){
        cat("
Object metadata:
")
        for(f in names(md)) cat(f, ": ", md[[f]], "\n")
      }
    })
    
    output$SEout2 <- renderText({
      if(is.null(SE())) return("No object loaded")
      paste("A SummarizedExperiment with ", ncol(SE()), "samples, ",
            length(DEAs()), " differential expression analyses and ",
            length(EAs()), " enrichment analyses.")
    })
    output$fileout <- renderPrint({
      if(is.null(SE())) return(NULL)
      print(SE())
    })

    output$features <- renderDT({
      if(is.null(SE())) return(NULL)
      RD <- rowData(SE())
      RD <- RD[,unlist(sapply(RD, is.vector)),drop=FALSE]
      datatable( as.data.frame(RD), filter="top", class="compact",
                 options=list( pageLength=30, dom = "fltBip" ),
                 extensions=c("ColReorder") )
    }, server = TRUE)

    output$samples <- renderDT({
      if(is.null(SE())) return(NULL)
      datatable( as.data.frame(colData(SE())), filter="top", class="compact",
                 options=list( pageLength=30, dom = "fltBip" ),
                 extensions=c("ColReorder") )
    })

    ############
    ### BEGIN DEAs

    DEAs <- reactive({
      if(is.null(SE())) return(list())
      RD <- rowData(SE())
      deas <- grep("^DEA\\.", colnames(RD), value=TRUE)
      deas <- deas[unlist(sapply(deas, FUN=function(x) is.data.frame(RD[[x]]) ||
                            is(RD[[x]], "DFrame") ))]
      if(length(deas)==0) return(list())
      lapply(setNames(deas, gsub("^DEA\\.","",deas)), FUN=function(x) RD[[x]])
    })
    
    output$menu_DEA <- renderUI({
      if(is.null(DEAs()) || length(DEAs())==0)
        return(menuItem("DEA results", tabName="tab_dea",
                    badgeLabel="N/A", badgeColor="red"))
      menuItem("DEA results", tabName="tab_dea", badgeLabel=length(DEAs()),
               badgeColor="aqua")
    })

    output$dea_input <- renderUI({
      hideElement("dea_box")
      if(length(DEAs())==0)
        return(tags$h5(icon("exclamation-triangle"),
                       "The object contains no DEA."))
      showElement("dea_box")
      selectInput("dea", label="Differential expression analysis",
                  choices=names(DEAs()))
    })

    DEA <- reactive({
      if(is.null(input$dea) || input$dea=="" || is.null(DEAs()[[input$dea]]))
        return(NULL)
      dea <- DEAs()[[input$dea]]
      dea[order(dea$FDR, dea$PValue),]
    })

    observeEvent(input$dea_geneFilt, {
      if(length(g <- input$dea_table_rows_all)>0){
        g <- row.names(DEA())[head(g,maxPlot)]
        updateTextAreaInput(session, "input_genes", value=paste(g, collapse=", "))
        updateTabItems(session, "main_tabs", selected="tab_heatmap")
      }
    })
    observeEvent(input$dea_geneSel, {
      if(length(g <- input$dea_table_rows_selected)>0){
        g <- row.names(DEA())[head(g,maxPlot)]
        updateTextAreaInput(session, "input_genes", value=paste(g, collapse=", "))
        updateTabItems(session, "main_tabs", selected="tab_heatmap")
      }
    })

    output$dea_overview <- renderText({
      if(is.null(dea <- DEA())) return(NULL)
      paste("A differential expression analysis across", sum(!is.na(dea$FDR)),
      "features, ", sum(dea$FDR<0.05, na.rm=TRUE), "of which are at FDR<0.05.")
    })

    output$dea_pvalues <- renderPlot({
      if(is.null(dea <- DEA())) return(NULL)
      hist(dea$PValue, xlab="Unadjusted p-values", main="")
    })

    output$dea_table <- renderDT({
      if(is.null(dea <- DEA())) return(NULL)
      datatable( dround(as.data.frame(dea), digits=3, roundGreaterThan1=TRUE),
                 filter="top", class="compact", extensions=c("ColReorder"),
                 options=list( pageLength=30, dom = "fltBip" ) )
    })

    output$dea_volcano <- renderPlotly({
      if(is.null(dea <- DEA())) return(NULL)
      dea <- head(dea, 2000)
      if(is.null(dea$logFC))
        dea$logFC <- apply(
          as.matrix(dea[,grep("logFC|log2FC",colnames(dea)),drop=FALSE]), 1,
          FUN=function(x) x[which.max(abs(x))])
      if(is.null(dea$MeanExpr)){
        if(!is.null(dea$logCPM)) dea$MeanExpr <- dea$logCPM
        if(!is.null(dea$baseMean)) dea$MeanExpr <- dea$baseMean
      }
      dea$feature <- row.names(dea)
      if(is.null(dea$MeanExpr)){
        p <- ggplot(dea, aes(logFC, -log10(FDR), Feature=feature, FDR=FDR,
                             PValue=PValue))
      }else{
        p <- ggplot(dea, aes(logFC, -log10(FDR), Feature=feature, FDR=FDR,
                             PValue=PValue, MeanExpr=MeanExpr, colour=MeanExpr))
      }
      p <- p + geom_vline(xintercept=0, linetype="dashed") + geom_point() +
        theme_classic()
      plotlyObs$resume()
      ggplotly(p, tooltip=c("Feature","logFC","PValue","FDR"), source="volcano")
    })

    plotlyObs <- observeEvent(event_data("plotly_click", "volcano",
                                         priority="event"), suspended=TRUE, {
      if(is.null(dea <- DEA())) return(NULL)
      event <- event_data("plotly_click", "volcano")
      if(!is.list(event) || is.null(event$pointNumber) && event$pointNumber>=0)
        return(NULL)
      g <- row.names(dea)[as.integer(event$pointNumber+1)]
      selGene(g)
      updateTabItems(session, "main_tabs", "tab_gene")
    })

    ### END DEAs
    ############
    ### BEGIN EAs

    EAs <- reactive({
      if(is.null(SE()) || is.null(eas <- metadata(SE())$EA) || !is.list(eas) ||
         length(eas)==0) return(list())
      if(is.list(eas[[1]]) && !is.data.frame(eas[[1]]) &&
         !is(eas[[1]], "DFrame")) eas <- unlist(eas, recursive=FALSE)
      eas
    })
    
    output$menu_Enrichments <- renderUI({
      if(length(EAs())==0) return(menuItem("Enrichments", tabName="tab_ea", 
                                           badgeLabel="N/A", badgeColor="red"))
      menuItem("Enrichments", tabName="tab_ea", badgeLabel=length(EAs()),
               badgeColor="aqua")
    })

    output$ea_input <- renderUI({
      hideElement("ea_box")
      if(length(EAs())==0)
        return(tags$h5(icon("exclamation-triangle"),
                       "The object contains no enrichment analysis"))
      eas <- EAs()
      eas <- setNames(names(eas), gsub("DEA\\.","",names(eas)))
      showElement("ea_box")
      selectInput("ea", label="Enrichment analysis", choices=eas)
    })

    EA <- reactive({
      if(is.null(input$ea) || input$ea=="" || is.null(ea <- EAs()[[input$ea]]))
        return(NULL)
      ea
    })

    output$ea_table <- renderDT({
      if(is.null(ea <- EA())) return(NULL)
      if(!isTRUE(input$ea_showGenes)) ea$genes <- ea$leadingEdge <- NULL
      datatable( dround(as.data.frame(ea), digits=3, roundGreaterThan1=TRUE),
                 filter="top", class="compact", extensions=c("ColReorder"),
                 options=list( pageLength=30, dom = "fltBip", colReorder=TRUE ) )
    })

    observeEvent(input$ea_geneSel, {
      if(is.null(ea <- EA()) || is.null(ea$genes)) return(NULL)
      if(length(RN <- input$ea_table_rows_selected)>0){
        g <- head(unique(unlist(.getWordsFromString(ea[RN,"genes"]))), maxPlot)
        updateTextAreaInput(session, "input_genes", value=paste(g, collapse=", "))
        updateTabItems(session, "main_tabs", selected="tab_heatmap")
      }
    })

    ### END EAs
    ############
    ### BEGIN HEATMAP

    selGenes <- reactive({
      g <- gsub(",|\n|\r|;,"," ",input$input_genes)
      g <- strsplit(g," ",fixed=T)[[1]]
      g <- unique(g[which(g!="")])
      grepGene(row.names(SE()), g)
    })


    output$heatmap <- renderPlot({
      validate( need(!is.null(SE()), "No SummarizedExperiment loaded.") )
      g <- selGenes()
      validate( need(length(g)>0, paste("No gene selected in the 'Genes' tab,",
                      "or the selected genes are not found in the object.")) )
      if(length(g)>2 && input$hm_clusterRow){
        srow <- 1:ncol(SE())
      }else{
        srow <- NULL
      }
      se <- SE()[g,]
      o <- input$hm_order
      if(is.null(o)) o <- c()
      for(f in rev(o)) se <- se[,order(colData(se)[[f]])]
      breaks <- NULL
      if(grepl("logFC|log2FC|scaledLFC|zscore", input$assay_input2,
               ignore.case=TRUE) || input$hm_scale)
        breaks <- 1-(input$hm_breaks/100)
      draw(sechm(se, g, do.scale=input$hm_scale,
           assayName=input$assay_input2,
           sortRowsOn=srow, anno_columns=input$hm_anno, gaps_at=input$hm_gaps,
           cluster_cols=input$hm_clusterCol, cluster_rows=FALSE,
           breaks=breaks), merge_legends=TRUE)

    })

    ### END HEATMAP
    ############




    ############
    ### START GENE TAB
    
    selGene <- reactiveVal()
    observeEvent(input$gene_input, {
      selGene(input$gene_input)
    })
    
    output$gene_plot <- renderPlot({
      d <- tryCatch(meltSE(SE(), selGene()),
                    error=function(x) NULL)
      validate( need(!is.null(d) && nrow(d)>0,
                  "No gene selected. Select one in the dropdown list above.
You can type the first few letters in the box and select from the matching suggestions."))
      gr <- input$select_groupvar
      if(input$asfactor) d[[gr]] <- factor(d[[gr]])
      if(input$select_plotpoints){
        p <- ggplot(d, aes_string(input$select_groupvar, input$assay_input,
                                  colour=input$select_colorvar))
      }else{
        p <- ggplot(d, aes_string(input$select_groupvar, input$assay_input,
                                  fill=input$select_colorvar))
      }
      if(input$plottype_input=="violin plot"){
        p <- p + geom_violin()
      }else{
        p <- p + geom_boxplot(outlier.shape = NA)
      }
      if(input$select_plotpoints)
        p <- p + geom_point(position = position_jitterdodge())
      p <- p + theme_classic() + ggtitle(input$gene_input) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


      if(!is.null(input$select_gridvars)){
        form <- as.formula(paste0("~",paste(input$select_gridvars, collapse="+")))
        if(input$select_freeaxis){
          p <- p + facet_wrap(form, scales="free_y")
        }else{
          p <- p + facet_wrap(form)
        }
      }
      p + ggtitle(selGene())
    })
    
    output$gene_inList <- renderText({
      if(is.null(g <- selGene()) || g=="" || length(genelists)==0)
        return("")
      if(grepl(".+\\.[a-zA-Z].+", g)) g <- gsub("^[^.]+\\.","",g)
      g <- strsplit(g,"/")[[1]]
      x <- which(sapply(genelists, FUN=function(x){ any(g %in% x) }))
      if(length(x)==0)
        return("This gene is included in none of the registered genelists.")
      out <- "This gene is included in the following list(s):\n"
      paste(out, names(genelists)[x], collapse="\n")
    })
    

    ### END GENE TAB
    ############
    ### START GENELISTS TAB
    
    output$genelist_size <- renderText({
      if(is.null(input$genelist_input) || 
         is.null(gl <- genelists[[input$genelist_input]])) return(NULL)
      length(gl)
    })
    output$genelist_out <- renderText({
      if(is.null(input$genelist_input) || 
         is.null(gl <- genelists[[input$genelist_input]])) return(NULL)
      paste(gl, collapse=", ")
    })
    observeEvent(input$btn_importGenelist, {
      if(!is.null(input$genelist_input) &&  
         !is.null(g <- genelists[[input$genelist_input]]))
      g <- head(g,maxPlot)
      updateTextAreaInput(session, "input_genes", value=paste(g, collapse=", "))
      updateTabItems(session, "main_tabs", selected="tab_heatmap")
    })
    
    
    ### END GENELISTS TAB
    ############
    

    observeEvent(input$quickStart, showModal(.getHelp("general")))
    observeEvent(input$help_SE, showModal(.getHelp("SE")))
    observeEvent(input$help_gassay, showModal(.getHelp("assay")))
    observeEvent(input$help_ggroup, showModal(.getHelp("group")))
    observeEvent(input$help_ggrid, showModal(.getHelp("grid")))
    observeEvent(input$help_gfreeaxes, showModal(.getHelp("grid")))
    observeEvent(input$help_hmassay, showModal(.getHelp("assay")))
    observeEvent(input$help_hmscale, showModal(.getHelp("scale")))
    observeEvent(input$help_hmtrim, showModal(.getHelp("scaletrim")))
    observeEvent(input$help_genelists, showModal(.getHelp("genelists")))
    
    
    waiter_hide()
  }
}
