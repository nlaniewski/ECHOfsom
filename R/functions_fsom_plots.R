#' Pheatmap of FlowSOM Clusters/Nodes
#'
#' @param dat.mat numeric matrix
#' @param color.type c("divergent","sequential"); heatmap colorscale
#' @param break.vec a numeric vector for generating color breaks
#' @param ... additional arguments for pheatmap::pheatmap
#'
#' @return a pheatmap
#' @export
#'
pheat <- function(dat.mat, color.type = c("divergent","sequential"),break.vec = seq(0, 1, by = 0.05),...){
  if(is.null(rownames(dat.mat))){
    rownames(dat.mat) <- 1:nrow(dat.mat)
  }
  if(any(grepl('node|cluster',colnames(dat.mat)))){
    dat.mat <- dat.mat[,!colnames(dat.mat) %in% c("node","cluster")]
  }
  dat.mat <- dat.mat[,order(colnames(dat.mat))]
  color.type <- match.arg(color.type)
  color.type <- switch(color.type,
                       divergent = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7,name ="RdYlBu")))(length(break.vec)),
                       sequential = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name ="Greens"))(length(break.vec))
  )
  pheatmap::pheatmap(dat.mat,
                     color = color.type,
                     breaks = break.vec,
                     #silent = F,
                     ...
  )
}

#' Cluster medians heatmap: plotly heatmap with hclust dendrogram ordering
#'
#' @param cluster.medians cluster median values; returned from 'generate.cluster.medians()'
#' @param break.vec a numeric vector for generating color breaks
#'
#' @return Plotly heatmap
#' @export
#'
plheat <- function(cluster.medians,break.vec = seq(0, 1, by = 0.05)){
  #prepare cluster data
  dat <- as.matrix(cluster.medians[,!colnames(cluster.medians) %in% 'cluster'])
  c.order <- rev(unlist(stats::as.dendrogram(stats::hclust(stats::dist(dat)))))
  m.order <- unlist(stats::as.dendrogram(stats::hclust(stats::dist(t(dat)))))
  dat <- dat[c.order,m.order]
  rownames(dat) <- c.order
  #plotly heatmap
  plotly.heatmap <- plotly::plot_ly(x=colnames(dat),
                                    y=rownames(dat),
                                    z=dat,
                                    zauto = FALSE,
                                    zmin = min(break.vec),
                                    zmax = max(break.vec),
                                    type = "heatmap",
                                    colors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name ="Greens"))(length(break.vec)),
                                    source = 'cluster.heatmap'
  )
  plotly.heatmap <- plotly::layout(plotly.heatmap,
                                   yaxis = list(tickfont = list(size = 10),
                                                type = "category")
  )
  return(plotly.heatmap)
}

#' Plot a single ECHO-specific boxplot
#'
#' @param echo.melted.mdatframe a 'reshape2::melt()' data.frame; returned from 'echo.melt()'
#' @param variable.value cluster/node value (#) to plot
#' @param x.var x-axis variable to plot; factor
#' @param drop.col.var a paired character 'c(factor,level)'; drops a factor level
#'
#' @return a ggplot boxplot
#' @export
#'
echo.boxplot.single <- function(echo.melted.mdatframe,variable.value,x.var = 'visit',drop.col.var=NULL){
  ##
  if(any(grepl('cluster',colnames(echo.melted.mdatframe)))){
    i <- 'cluster'
  }else if(any(grepl('node',colnames(echo.melted.mdatframe)))){
    i <- 'node'
  }else{
    stop("Expect a column named 'cluster' or 'node'; is this a melted frame (returned from 'echo.melt')?")
  }
  ##
  dat <- droplevels(echo.melted.mdatframe[echo.melted.mdatframe[,i]==variable.value,])
  title.name = paste(paste0(stringr::str_to_title(i),":"),unique(dat[,i]))
  ##
  y.name <- as.name(grep("% of",colnames(dat),value = T))
  ##
  if(!is.null(drop.col.var)&length(drop.col.var)==2){
    dat <- droplevels(dat[dat[,drop.col.var[1]]!=drop.col.var[2],])
  }

  n.samples <- table(dat[,x.var])
  caption.name <- paste(paste0(names(n.samples),":"),n.samples, collapse = "   ")

  p <- ggplot2::ggplot(dat) +
    ggplot2::aes_string(x=x.var,y=y.name,color=x.var,fill=x.var) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(1),outlier.shape = NA,alpha = 0.2,lwd = 0.05) +
    ggplot2::geom_jitter(size = 1, width = 0.2, alpha = 0.6) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 10),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 10),
                   #strip.text = element_text(size = 8),
                   strip.text.x = ggplot2::element_text(size = 8, margin = ggplot2::margin(b = 1, t = 1)),
                   strip.background = ggplot2::element_rect(size = 0.1),
                   legend.title = ggplot2::element_text(size = 10),
                   legend.text = ggplot2::element_text(size = 10),
                   panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 0.1)) +
    #facet_wrap(~cluster, scales = "free_y") +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(title = title.name,caption = caption.name)
  ##
  if(x.var=='visit'){
    colors.group.visit <- c("V4" = '#999999',
                            "V6" = '#E69F00',
                            "V7" = '#56B4E9',
                            "Adult" = '#009E73'
    )
    p <- p +
      ggplot2::scale_color_manual(name = "Blood Draw",
                                  labels = c("Birth", "Month 6", "Month 12", "Adult"),
                                  values = colors.group.visit) +
      ggplot2::scale_fill_manual(values = colors.group.visit) +
      ggplot2::scale_x_discrete(labels = c("V4" = "Birth", "V6" = "Month 6", "V7" = "Month 12"))
  }else if(x.var=='condition'){
    p <- p +
      ggplot2::scale_color_manual(name = paste("Stimulation","Condition",sep = "\n"),
                         labels = c("SEB", "Unstim."),
                         values = c("red","blue")) +
      ggplot2::scale_fill_manual(values = c("red","blue"))
  }
  return(p)
}

#' Interactive Shiny; various plots of FlowSOM results
#'
#' @param fsom.somnambulated.rds.path a file.path to exisiting fsom; returned from 'fsom.somnambulation()'
#' @param fsom.somnambulated object returned from 'fsom.somnambulation.mod()'
#' @param markers character vector of length 2; marker pairs for initial plots
#'
#' @return launches an interactive window ('shiny::shinyApp()')
#' @export
#'
somnambulate <- function(fsom.somnambulated.rds.path=NULL,fsom.somnambulated=NULL,markers = NULL){

  ##ggplot function
  gg.func.bivariate <- function(dat,...){
    ggplot2::ggplot(dat,ggplot2::aes_string(...)) +
      ggplot2::geom_hex(bins = 200) +
      viridis::scale_fill_viridis(option = "plasma", limits = c(0, 50), oob = scales::squish)
  }
  ##choose file
  if(is.null(fsom.somnambulated.rds.path)&is.null(fsom.somnambulated)){
    f.path <- file.choose()
    message("loading chosen fsom.rds...")
    fsom.somnambulated <- readRDS(f.path)
  }
  ##provide path
  if(!is.null(fsom.somnambulated.rds.path)&is.null(fsom.somnambulated)){
    fsom.somnambulated <- readRDS(fsom.somnambulated.rds.path)
  }
  ##provide environment object
  if(is.null(fsom.somnambulated.rds.path)&!is.null(fsom.somnambulated)){
    message("Using environment object...")
  }
  ##check class
  if(!inherits(fsom.somnambulated,"FlowSOM Somnambulated")){
    stop("Need the return object from 'fsom.somnambulation()'")
  }

  ##
  c.names <- colnames(fsom.somnambulated$dat.all)
  dat.all.max.rows <- nrow(fsom.somnambulated$dat.all)

  if(!is.null(markers)){
    m1 <- markers[1];m2 <- markers[2]
  }else{
    m1 <- m2 <- NULL
  }

  ####

  ##ui
  ui <- shiny::fluidPage(
    #useShinyjs(),
    # Application title
    shiny::titlePanel("SOMnambulate: Sleep Walking Through a High-Dimensional Dream",
                      windowTitle = "SOMnambulate"),

    shiny::tabsetPanel(

      shiny::tabPanel("Cyto Plots",
                      shiny::sidebarLayout(
                        shiny::sidebarPanel(
                          #
                          shiny::selectInput(inputId = "marker1",
                                             label = "Marker (x):",
                                             choices = c.names,
                                             selected = m1),
                          shiny::selectInput(inputId = "marker2",
                                             label = "Marker (y):",
                                             choices = c.names,
                                             selected = m2),
                          shiny::selectInput(inputId = "mdat1",
                                             label = "Meta-data variable:",
                                             choices = c('visit','condition','batch'),
                                             selected = 'visit'),
                          shiny::selectInput(inputId = "nc",
                                             label = "Node or Cluster:",
                                             choices = c('node','cluster'),
                                             selected = 'cluster'),
                          shiny::selectInput(inputId = "nc.val",
                                             label = "Cluster/Node #:",
                                             choices = ""),
                          shiny::numericInput(inputId = "rowsamp",
                                              label = "# of 'All Events' to display:",
                                              value = 1E5,
                                              min = 1E5,
                                              max = dat.all.max.rows,
                                              step = 1E5),
                          #
                          width = 2,
                          id = 'sidebar'
                        ),

                        shiny::mainPanel(
                          shiny::fluidRow(
                            shiny::plotOutput("plots_2_1")
                          ),
                          shiny::fluidRow(
                            if(!is.null(fsom.somnambulated$mdats)){
                              shiny::column(5, shiny::plotOutput("plots_3"))
                            },
                            #column(6, plotOutput("plots_heat"))
                            shiny::column(7, plotly::plotlyOutput("plotly_heat"))
                          ),
                          width = 10
                        )
                      )
      ),

      shiny::tabPanel("Boxplots",
                      shiny::fillPage(
                        shiny::plotOutput("boxplots_faceted", height = "80vh")
                      )
      ),

      shiny::tabPanel("Heatmap - Full",
                      shiny::fillPage(
                        shiny::plotOutput("pheat", height = "80vh")
                      )
      )
    )
  )

  ##server
  server = function(input, output) {
    ##
    shiny::observe({
      x <- input$nc
      choices = seq(fsom.somnambulated$nc.vals[x])
      shiny::updateSelectInput(inputId = 'nc.val',
                               label = paste(stringr::str_to_title(x),"#:"),
                               choices = choices)
    })
    ##
    row.index <- shiny::reactive({
      set.seed(20040501)
      i <- sample(1:dat.all.max.rows,input$rowsamp)
    })
    ##
    plot1 <- shiny::reactive({
      gg.func.bivariate(data.frame(fsom.somnambulated$dat.all[row.index(),]),x = input$marker1,y = input$marker2) +
        ggplot2::xlim(fsom.somnambulated$xy.lims[,input$marker1]) +
        ggplot2::ylim(fsom.somnambulated$xy.lims[,input$marker2]) +
        ggplot2::labs(title = "All Events",
                      subtitle = paste(length(row.index()), "of", dat.all.max.rows, "displayed"))
    })
    plot2 <- shiny::reactive({
      shiny::validate(shiny::need(input$nc.val,"populating cluster/node values based on 'updateSelectInput'"))
      if(input$nc=='node'){
        title.sub <- fsom.somnambulated$node.titles[[input$nc.val]]
        dat <- data.frame(fsom.somnambulated$cluster.mats[[fsom.somnambulated$metaclustering[as.numeric(input$nc.val)]]])
        dat <- dat[dat$node==input$nc.val,]
        p <- gg.func.bivariate(dat,x = input$marker1,y = input$marker2)
      }else if(input$nc=='cluster'){
        title.sub <- fsom.somnambulated$cluster.titles[[input$nc.val]]
        p <- gg.func.bivariate(data.frame(fsom.somnambulated$cluster.mats[[input$nc.val]]),x = input$marker1,y = input$marker2)
      }
      p <- p +
        ggplot2::xlim(fsom.somnambulated$xy.lims[,input$marker1]) +
        ggplot2::ylim(fsom.somnambulated$xy.lims[,input$marker2]) +
        ggplot2::labs(title = paste(paste0(stringr::str_to_title(input$nc),":"),input$nc.val),
                      subtitle = title.sub) +
        ggplot2::guides(fill = "none")
      return(p)
    })
    if(!is.null(fsom.somnambulated$mdats)){
      plot3 <- shiny::reactive({
        shiny::validate(shiny::need(input$nc.val,"populating cluster/node values based on 'updateSelectInput'"))
        i <- input$nc
        echo.boxplot.single(echo.melted.mdatframe = fsom.somnambulated$mdats[[i]],
                            variable.value = input$nc.val,
                            x.var = input$mdat1)
      })
    }
    # boxplots.faceted <- reactive({
    #   boxplot.list[[input$mdat1]]
    # })

    output$plots_2_1 <- shiny::renderPlot({
      do.call(ggpubr::ggarrange,list(plot2(),plot1()))
    })
    output$plots_3 <- shiny::renderPlot({
      plot3()
    })

    output$plotly_heat <- plotly::renderPlotly(fsom.somnambulated$heatmaps$pl.heat)

    output$pheat <- shiny::renderPlot({fsom.somnambulated$heatmaps$pheat})

    # output$boxplots_faceted <- renderPlot({
    #   boxplots.faceted()
    # })

    ##plotly heatmap click data
    clicks <- shiny::reactiveValues(dat = data.frame(marker1 = NA, marker2 = NA))

    click <- shiny::reactive({
      plotly::event_data("plotly_click", priority = 'event', source = 'cluster.heatmap')
    })
    #
    shiny::observeEvent(eventExpr = click(),{
      if(is.na(clicks$dat$marker1)&is.na(clicks$dat$marker2)){
        clicks$dat$marker1 <- click()$x
      }else if(is.na(clicks$dat$marker2)){
        clicks$dat$marker2 <- click()$x
        shiny::updateSelectInput(inputId = 'marker1',
                                 selected = clicks$dat$marker1)
        shiny::updateSelectInput(inputId = 'marker2',
                                 selected = clicks$dat$marker2)
        shiny::updateSelectInput(inputId = 'nc.val',
                                 selected = stringr::str_extract(click()$y,"[0-9]+"))
        clicks$dat$marker1 <- NA
        clicks$dat$marker2 <- NA
      }
    })
  }
  # Run app ----
  shiny::shinyApp(ui, server)
}
##
##

# echo.boxplots <- function(mdat = fsom$mdat,cluster.counts=fsom$cluster.counts,factor.col.and.level.to.drop=NULL,...){
#   ##
#   dat.melt <- echo.melt(...)
#   ##
#   y.name <- grep("% of", colnames(dat.melt),value = T)
#   ##
#   colors.group.visit <- c("V4" = '#999999',
#                           "V6" = '#E69F00',
#                           "V7" = '#56B4E9',
#                           "Adult" = '#009E73'
#   )
#   ##
#   if(!is.null(factor.col.and.level.to.drop)){
#     colors.group.visit <- colors.group.visit[!names(colors.group.visit) %in% factor.col.and.level.to.drop]
#     p <- ggplot(data = droplevels(dat.melt[dat.melt[,factor.col.and.level.to.drop[1]] != factor.col.and.level.to.drop[2], ]))
#   }else{
#     p <- ggplot(data = dat.melt)
#   }
#   ##
#   p.list <- sapply(c('visit','condition','batch'),function(i){
#     i.f <- as.formula(paste0('~',i))
#     p <- p +
#       aes_(x = i.f, y = as.name(y.name), color = i.f, fill = i.f) +
#       geom_boxplot(position = position_dodge(1), outlier.shape = NA, alpha = 0.2, lwd = 0.05) +
#       geom_jitter(size = 1, width = 0.2, alpha = 0.3) +
#       facet_wrap(~cluster, scales = "free_y")
#     if(i=='visit'){
#       p <- p +
#         scale_color_manual(name = "Blood Draw",
#                            labels = c("Birth", "Month 6", "Month 12", "Adult"),
#                            values = colors.group.visit) +
#         scale_fill_manual(values = colors.group.visit) +
#         scale_x_discrete(labels = c("V4" = "Birth", "V6" = "Month 6", "V7" = "Month 12"))
#     }
#     p <- p +
#       theme_bw() +
#       theme(axis.text.x = element_blank(),
#             axis.text.y = element_text(size = 10),
#             axis.title.x = element_blank(),
#             axis.title.y = element_text(size = 10),
#             #strip.text = element_text(size = 8),
#             strip.text.x = element_text(size = 8, margin = margin(b = 1, t = 1)),
#             strip.background = element_rect(size = 0.1),
#             legend.title = element_text(size = 10),
#             legend.text = element_text(size = 10),
#             panel.border = element_rect(fill = NA, color = "black", size = 0.1)) +
#       #facet_wrap(~cluster, scales = "free_y") +
#       guides(fill = "none")
#     return(p)
#   },simplify = F)
#   ##
# }
