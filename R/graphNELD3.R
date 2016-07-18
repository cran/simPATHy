# graphNEL with D3 JavaScript
#
# Plot a graphNEL object, with an associated covariance matrix, using the D3 JavaScript force directed newtwork.
# @param graph A graphNEL object.
# @param S The sample covariance matrix. If supplied, it specifies the colour of the edges between nodes. If \code{NULL} the default colour is \code{#ccc}.
# @param colNode A character string specifying the colour of the nodes. The colour node is common for all the nodes.
# @param limColEdges A vector of two elements with the upper and lower bounds of the S values.
# @param edgeClick Activate edge animation with click.
# @param width,height Height and width of the plot in pixels. Must be a number, which will be coerced to a string and have 'px' append.
# @references https://d3js.org (Micheal Bostock).
# @references http://www.htmlwidgets.org (Ramnath Vaidyanathan, Kenton Russell, and RStudio).
# @references https://christophergandrud.github.io/networkD3/ (Christopher Gandrud, JJ Allaire, & Kent Russell)
# @seealso \code{\link[simPATHy]{plotGraphNELD3}},  \code{\link[simPATHy]{easyLookShiny}}, \code{\link[simPATHy]{fitSgraph}}
graphNELD3 <- function(graph,S=NULL,colNode = "#c0c0c0",limColEdges=NULL,edgeClick=FALSE, width = NULL, height = NULL) {

  legend<-FALSE
  dataLegend<-data.frame()

  type<-NULL
  if (!(class(graph) == "graphNEL"))
    stop("graph argument is not a graphNEL object")
  else {
    if (gRbase::is.DG.graphNEL(graph)) type <- "dag"
    if (gRbase::is.UG.graphNEL(graph)) type <- "ug"
    if (is.null(type)) stop("graph argument is not a valid object")
  }

  nodes<-data.frame(name=graph::nodes(graph),stringsAsFactors = FALSE)
  Elist<-gRbase::edgeList(graph)
  E<-sapply(Elist,function(e,n=nodes$name) match(e,n)-1 )

  # option: colour edges
  if( is.null(limColEdges) | is.null(S) ) colEdges<-rep("#ccc",length(Elist))
  else{
      # check colEdges
      # sto rappresentando correlazioni/covarianze
      if(length(limColEdges)>2) stop("length of limColEdges must be 2") #scrivi meglio
      if(limColEdges[1]>0 & limColEdges[2]<0) stop("limColEdges parameter not valid") # scrivi meglio

      values<-sapply(Elist,function(e,s=S) s[e[1],e[2]] )
      colEdges<-giveColors(values,lim = limColEdges)

      legend<-TRUE

      ns<-4

      step1<-abs(limColEdges[1]/ns)
      step2<-abs(limColEdges[2]/ns)
      dataLegend<-data.frame(value=c(limColEdges[1]+(step1* (0:(ns-1)) ),0,0+( (1:ns)*step2 )))

  }


  edges<-data.frame(source=E[1,],target=E[2,],colour=colEdges,stringsAsFactors = FALSE)

  options<-list(type=type,colNode=colNode,legend=legend,edgeClick=edgeClick)

  x<-list(nodes=nodes,edges=edges,dataLegend=dataLegend,options=options)
  # create widget
  htmlwidgets::createWidget(
    name = 'graphNELD3',
    x= x,
    width = width,
    height = height,
    htmlwidgets::sizingPolicy(padding = 10, browser.fill = TRUE),
    package = 'simPATHy'
  )
}

#' Shiny bindings for plotGraphNELD3
#'
#' Output and render functions for using plotGraphNELD3 within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId Output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a graphNELD3
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name graphNELD3-shiny
#' @export
graphNELD3Output <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'graphNELD3', width, height, package = 'simPATHy')
}

#' @rdname graphNELD3-shiny
#' @export
renderGraphNELD3 <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, graphNELD3Output, env, quoted = TRUE)
}

