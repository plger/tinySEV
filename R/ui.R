#' A shiny SummarizedExperiment viewer
#'
#' @name tinySEV
#' @rdname tinySEV
#' @aliases tinySEV
#' @examples
#' ui <- tinySEV.iu()
#' server <- tinySEV.server(objects=c(SE1="path/to/SE1.rds",
#'                                    SE2="path/to/SE2.rds"))
NULL

#' tinySEV.ui
#'
#' @param title The title of the app (displayed in the header)
#' @param waiterContent Optional content of the loading mask; should be a
#' `tagList`, NULL to use default, or FALSE to disable the waiter
#' @param about Optional content of the introduction page (NULL to disable
#' intro page)
#'
#' @return a shiny UI
#' @export
#' @import shiny shinydashboard shinyjqui waiter
#' @importFrom shinycssloaders withSpinner
#' @importFrom plotly plotlyOutput
#' @importFrom shinyjs useShinyjs
#' @importFrom DT DTOutput
tinySEV.ui <- function(title="tinySEV", waiterContent=NULL, about=NULL){

  if(is.null(waiterContent) || isTRUE(waiterContent))
    waiterContent <- tagList(
      tags$h3("Please wait while the application is initialized..."), spin_1())
  if(isFALSE(waiterContent)){
    waiterContent <- NULL
  }else{
    waiterContent <- waiter_show_on_load(html=waiterContent)
  }
  aboutMenu <- NULL
  if(!is.null(about)) aboutMenu <- menuItem("About", tabName="tab_about")

  shinyUI( dashboardPage(
    dashboardHeader(title=title),
    dashboardSidebar(
      sidebarMenu(id="main_tabs", aboutMenu,
        .modify_stop_propagation(
          menuItem("Object", startExpanded=TRUE,
             selectInput("object", label=NULL, choices=c()),
             menuSubItem("Overview", tabName="tab_object"),
             menuSubItem("Samples", tabName="tab_samples"),
             menuSubItem("Features", tabName="tab_features"),
             menuItemOutput("uploadMenu"))),
        menuItem("DEA results", tabName="tab_dea"),
        menuItem("Enrichments", tabName="tab_ea"),
        menuItem("Plot gene", tabName="tab_gene"),
        .modify_stop_propagation(menuItem("Heatmap", startExpanded=TRUE,
          menuSubItem("Genes", tabName="tab_hm_genes"),
          menuSubItem("Heatmap", tabName="tab_heatmap")
        ))
      )
    ),
    dashboardBody(
      tags$head(tags$style(HTML("
        .sidebar-menu li.treeview, .sidebar-menu li.treeview:hover a{
        	background-color: #2c3b41;
        }
      "))),
      use_waiter(), useShinyjs(), waiterContent,
      tabItems(
        tabItem("tab_object",
                box(width=12, title="Object overview",
                    tags$p(withSpinner(textOutput("SEout2"))),
                    withSpinner(verbatimTextOutput("SEout")))),
        tabItem("tab_fileinput",
          box(width=7,
              fileInput("file", "Choose SE .rds file", multiple=FALSE,
                        accept=c(".rds",".RDS")),
              withSpinner(verbatimTextOutput("fileout")) )
        ),
        tabItem("tab_samples",
                box(width=12, withSpinner(DTOutput("samples")))),
        tabItem("tab_features",
                box(width=12, withSpinner(DTOutput("features")))),
  	    tabItem("tab_dea", uiOutput("dea_input"),
  	       tabBox(id="dea_box", width=12,
  	         tabPanel("Overview", textOutput("dea_overview"),
  	           withSpinner(plotOutput("dea_pvalues", width="400px", height="300px"))),
  	         tabPanel("Table", withSpinner(DTOutput("dea_table")),
  	                  actionButton("dea_geneFilt", "Transfer filtered genes (max 500) to the heatmap"),
  	                  actionButton("dea_geneSel", "Transfer selected rows (max 500) to the heatmap")),
  	         tabPanel("Volcano", withSpinner(plotlyOutput("dea_volcano")),
  	                  tags$p("You may click on a gene to view it."))
  	       )
  	     ),
  	     tabItem("tab_ea",
           box(width=12, id="ea_box", fluidRow(
                 column(8, uiOutput("ea_input")),
                 column(4, checkboxInput("ea_showGenes", "Show genes", value=FALSE))
               ),
               tags$div(style="width: 100%; overflow-x: scroll; font-size: 80%",
                        withSpinner(DTOutput("ea_table"))),
               actionButton("ea_geneSel", "Transfer genes from selected rows (max 500) to the heatmap"))),
        tabItem("tab_gene",
          column(4, selectizeInput("gene_input", "Select Gene", choices=c(), multiple=FALSE)),
          box(width=8, title="Options", collapsible=TRUE, collapsed=TRUE,
            selectInput("assay_input", "Assay", choices=c(), multiple=FALSE),
            column(6,
              selectInput("plottype_input", "Type of Plot",
                          choices=c("violin plot","box plot"), multiple=FALSE),
              checkboxInput('select_plotpoints','Plot Points', value=TRUE),
              checkboxInput('select_logaxis','Logarithmic Axis', value=FALSE)),
           column(6,
              selectInput("select_groupvar", "Group by", choices=c(), multiple=FALSE),
              checkboxInput('asfactor','As factor', value=TRUE),
              selectInput("select_colorvar", "Color by", choices=c(), multiple=FALSE),
              selectizeInput("select_gridvars", "Grid by", choices=c(),
                             options = list(maxItems = 2), multiple=TRUE),
              checkboxInput('select_freeaxis','Free Axis', value=TRUE))
          ),
          box(width=12, withSpinner(shinyjqui::jqui_resizable(plotOutput("gene_plot"))))
        ),
        tabItem("tab_hm_genes", box(width=12, title="Select genes to plot",
          textAreaInput('input_genes','Genes to plot', width="90%", rows=10,
            placeholder="Enter genes symbols separated by commas, spaces, or line breaks..."))),
        tabItem("tab_heatmap",
             box( width=12, title="Heatmap parameters", collapsible=TRUE,
                  column(4, selectInput("assay_input2", "Assay", choices=c(), multiple=FALSE),
                         checkboxInput('hm_scale', 'Scale rows'),
                         sliderInput("hm_breaks", "Color scale trim",
                                     min=98, max=100, value=99.5, step=0.25)
                         ),
                  column(4, selectizeInput('hm_anno', "Column annotation", choices=c(), multiple=T),
                         selectizeInput('hm_gaps', "Gaps at", choices=c(), multiple=T)
                  ),
                  column(4, selectizeInput('hm_order', "Column ordering", choices=c(), multiple=T),
                         checkboxInput('hm_clusterCol','Cluster columns', value=F),
                         checkboxInput('hm_clusterRow','Sort rows', value=T)
                  )
             ),
             box(width=12, title="Heatmap",
                 withSpinner(shinyjqui::jqui_resizable(plotOutput("heatmap", height=600))))
        ),
        tabItem("tab_about", about)
    ), tags$div(style="clear: both;")
  )))
}
