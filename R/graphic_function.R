# Give colors
#
# Given a value it returns ..
# @param val a vectors of values
# @param lim a vector of two elements
# @param digits integer indicating the number of decimal places to be used.
# @param outCol the colors assigned to the out of bounds values.
# @return Formatted colours
giveColors<-function(val,lim,digits=4,outCol="#c0c0c0"){
  #library(grDevices)
  val<-round(val,digits)

  CRPal_neg<-grDevices::colorRampPalette(c("red4","white"))
  CRPal_pos<-grDevices::colorRampPalette(c("white","royalblue4"))

  seq_neg<-seq(lim[1],0,by=1/10^digits)
  seq_pos<-seq(0,lim[2],by=1/10^digits)

  COL<-c(CRPal_neg(length(seq_neg)),CRPal_pos(length(seq_pos)-1))
  seq<-round(c(seq_neg,seq_pos[-1]),digits)
  colours<-COL[match(val,seq)]
  colours[is.na(colours)]<-outCol

  colours
}


#' Plot correlation or partial correlation matrix
#'
#' Plot a correlation or partial correlation matrix with the possibility to emphasize the graphical structure.
#' @param S1,S2 Sample covariance matrix. If \code{S2} supplied, the difference between the two corresponding correlation or partial correlation matrices is plotted.
#' @param type Character string specifying which matrix is to be plotted. Either \code{cor} for correlation matrix, or \code{pcor} for partial correlation matrix.
#' @param graph A graphNEL object.
#' @param path A list of edges in edgesList format (see \pkg{gRbase}).
#' @param main The main title.
#' @param colLim Numeric vector of length two specifying the lower and upper bound of the color range (see Details).
#' @param legendColor Logical value indicating whether the color legend should be added to the plot.
#' @details
#' If the \code{graph} is supplied, the zero elements of the adjacency matrix are represented as shaded squares, whereas non-zero elements are represented as squares with grey borderline. \cr \cr
#' Admissible values for \code{colLim} are contained in the interval \code{[-1,1]} when \code{S2=NULL}, otherwise the admissible interval is \code{[-2,2]}.
#' When an element is outside of the \code{colLim} interval, it is colored gray.
#' @return Correlation or partial correlation matrix plot.
#' @examples
#' if( require(gRbase) & require(graph)){
#'  graph <- gRbase::ug(~a:b, ~a:c, ~c:d, ~b:d, ~b:c)
#'
#'  S <- matrix(c(2,  0.8,0.5,-0.3,
#'               0.8,1.5,0.6,-0.7,
#'               0.5,0.6,1,  0.7,
#'               -0.3,-0.7,0.7,3), ncol=4,nrow=4)
#'  colnames(S) <- rownames(S) <- graph::nodes(graph)
#'
#'  # Plot the correlation matrix of S
#'  plotCorGraph(S)
#'
#'
#'  S<-fitSgraph(graph = graph,S = S)
#'  # Change the color range
#'  plotCorGraph(S, colLim=c(-0.5,0.5))
#'
#'  # Visualize the adjacency matrix
#'  plotCorGraph(S, type="cor", graph = graph)
#'
#'  # Show the partial correlation matrix
#'  plotCorGraph(S, type="pcor", graph = graph)
#'
#'  # Plot the difference between two matrices
#'  S2 <- S
#'  # Change the element c~a
#'  S2["a","c"] <- S2["c","a"]<- -0.1
#'  plotCorGraph(S1=S, S2=S2)
#'  plotCorGraph(S1=S, S2=S2, type="pcor")
#'
#'  S2<-fitSgraph(graph = graph,S = S2)
#'  # Highlight the graphical structure
#'  plotCorGraph(S1=S, S2=S2, type="pcor",graph = graph)
#'  # Highlight the element c~a
#'  plotCorGraph(S1=S, S2=S2, type="pcor",graph = graph,path = list(c("a","c")))
#'
#' }
#' @export
plotCorGraph<-function(S1,type="cor",S2=NULL,graph=NULL,path=NULL,main="",colLim=c(-1,1),legendColor=TRUE){

  a<-0.5
  digits<-4

  ## check0: colLim
  if(is.null(S2)){
    if( !(colLim[1]<= -0.1  & colLim[1]>=-1) ) stop("colLim[1] must be between -0.1 and -1")
    if( !(colLim[2]<= 1  & colLim[2]>=0.1) ) stop("colLim[2] must be between 0.1 and 1")
  } else {
    if( !(colLim[1]<= -0.1  & colLim[1]>=-2) ) stop("colLim[1] must be between -0.1 and -2")
    if( !(colLim[2]<= 2  & colLim[2]>=0.1) ) stop("colLim[2] must be between 0.1 and 2")
  }

  ## check1: matrix
  checkMatrix(S1,name=deparse(substitute(S1)))
  if(!is.null(S2)){
    checkMatrix(S2,name=deparse(substitute(S2)))
    ## check1.1: compatibility S1 e S2
    if( !all(dim(S1)==dim(S2)) | !all(colnames(S1)==colnames(S2))  ) stop("S1 and S2 have different sizes or different elements")
    diff<-TRUE
  } else diff<-FALSE

  ## check2: type
  if( !(type %in% c("cor","pcor")) ) stop("Invalid parameter type: must be cor or pcor")

  ## check3: graph
  if(!is.null(graph)){
    Gtype<-NULL
    if(!(class(graph)=="graphNEL")) stop("graph argument is not a graphNEL object")
    else {
      if(gRbase::is.DG(graph)) Gtype<-"dag"
      if(gRbase::is.UG(graph)) Gtype<-"ug"
    }
    if(is.null(Gtype)) stop("graph argument is not a valid object")
  }

  if(type=="cor"){
    # correlation case
    if(!diff) S<- riscala(S1)
    else S<- riscala(S1) - riscala(S2)
    Text<-"Correlation Matrix"

  } else {
    # partial correlation case
    if(!diff) S<- riscala(solve(S1))
    else S<- riscala(solve(S1)) - riscala(solve(S2))
    Text<-"Partial Correlation Matrix"

    indDiag<-!(diag(ncol(S))==1)
    S[indDiag]<- -S[indDiag]

  }
  S<-round(S,digits = digits)
  name<-colnames(S)
  n<-ncol(S)


  ## check4: path
  if(!is.null(path)){
    if(!(class(path)=="list")) stop("Invalid class path")
    if( !all(unique(unlist(path)) %in% name) ) stop("All path elements must be in S matrix")
  }


  ## Plot0: param plot
  graphics::par(pty="s",mar=c(1,3,3,4))
  graphics::plot(1,xlim=c(a,nrow(S)+a+1.5),ylim=c(a-1.5,nrow(S)+a),type="n",axes=FALSE,xlab="",ylab="",main=main)
  graphics::axis(2, line=-1.5,at = (1:nrow(S)),labels=name[nrow(S):1],tick=FALSE,las=2,cex.axis=0.6)
  graphics::axis(3, line=-1.5,at = 1:nrow(S), labels=name[1:nrow(S)],tick=FALSE,las=2,cex.axis=0.6)


  ## Plot0: plot S
  indRect<-which(!S==0,arr.ind=TRUE)

  if(nrow(indRect)>0){
    colRect<-giveColors(val=S[indRect],lim=colLim,digits=digits)
    invisible(sapply(1:nrow(indRect),function(i,iR=indRect,cR=colRect)  plotRect(coord =as.vector(iR[i,]) ,col =cR[i] ,a =a ,n =n) ))
  }

  ## Plot1: graph structure
  # se si fornisce graph
  # devo distinguere dag e ug
  if(!is.null(graph)){

    if(Gtype=="dag"){
      # directed graph
      #if(type=="cor") A<- gRbase::as.adjMAT(graph) + t(gRbase::as.adjMAT(graph))
      #else A<- gRbase::as.adjMAT(gRbase::moralize(graph))
      A<- gRbase::as.adjMAT(gRbase::moralize(graph))

    } else A<- gRbase::as.adjMAT(graph) # undirected graph

    diag(A)<-1
    # metto il bordo
    indBord<-which(!A==0,arr.ind=TRUE)
    indDens<-which(A==0,arr.ind=TRUE)

    if(nrow(indDens)>0) invisible(sapply(1:nrow(indDens),function(i,iR=indDens)  plotRect(coord =as.vector(iR[i,]) ,col="white",density = 30,a =a ,n =n) ))
    if(nrow(indBord)>0) invisible(sapply(1:nrow(indBord),function(i,iR=indBord)  plotRect(coord =as.vector(iR[i,]) ,col=NA,border = "gray50",a =0.4 ,n =n) ))
    # density a quelli vincolati
  }

  ## Plot2: path
  if(!is.null(path)){

    indPath<-matrix(t(sapply(path,function(e,N=name) match(e, N) )),ncol=2)
    indPath<-rbind(indPath,cbind(indPath[,2],indPath[,1]))
    invisible(sapply(1:nrow(indPath),function(i,iR=indPath)  plotRect(coord =as.numeric(iR[i,]) ,col=NA,border = "green",a =0.3 ,n =n) ))

  }

  if(legendColor){
    addScaleColor(n =n-0.5,min =colLim[1] ,max =colLim[2] ,a =a )
  }
  graphics::rect(xleft=a,ybottom=a,xright=nrow(S)+a,ytop=nrow(S)+a,border="gray50")
  graphics::mtext(text = Text,side = 1,line=-1.5,col="gray50",cex=0.7)

}



# Add Legend
#
# Add a Legend Color
# @param n number of elements
# @param min min value
# @param max max value
# @param a dist parameter
# @return add a legend to the plot
addScaleColor<-function(n,min=-1,max=1,a=0.5){

  seq_neg<-seq(min,0,length.out = 51)
  seq_pos<-seq(0,max,length.out = 51)

  s<-c(seq_neg,seq_pos[-1])
  s<-round(s,4)
  ns<-length(s)

  ytop<- n
  ybottom<- 1

  passo<-c(0,cumsum(rep((ytop-ybottom)/ns,ns-1)))

  vec_ytop<-ytop-passo
  vec_ybottom<- vec_ytop-passo[2]

  Xleft<-rep(n+1.5,ns)
  Xrigth<-rep(n+1.5+a,ns)
  col<-giveColors(val = s,lim =c(min,max) )

  invisible(sapply(1:ns,function(i) graphics::rect(xleft = Xleft[i],ybottom =vec_ybottom[i],xright =Xrigth[i]  ,ytop = vec_ytop[i] ,border = col[i],col=col[i]) ))
  graphics::rect(xleft = Xleft[1],ybottom =vec_ybottom[ns],xright =Xrigth[1]  ,ytop = vec_ytop[1] ,border = "gray")

  t<-seq(from=1,to=ns,length.out = 11)

  graphics::text(x = (Xrigth+0.5)[t] , y = ((vec_ybottom+vec_ytop)/2)[t], labels = round(s[t],2),cex=a,col="gray50")
}

# Plot a rectangle
#
# Plot a single rectangle
# @param coord position
# @param col colour
# @param n number of variable
# @param border color borde
# @param a dist parameter
# @param density number density
# @return a rect formatted element
plotRect<-function(coord,col,a,n,border="white",density=NULL){
  x<-coord[2]
  y<- (n+1)-coord[1]
  graphics::rect(xleft=x-a,ybottom=y-a,xright=x+a,ytop=y+a,col=col,border=border,density=density)
}

# Check matrix
#
# Check format matrix
# @param M squared matrix
# @param name name error
# @return TRUE or FALSE
checkMatrix<-function(M,name){

  #if( !(class(M)=="matrix") ) stop(paste0(name," must be class matrix"))
  if( !( sum(class(M)=="matrix")>0  ) ) stop(paste0(name," must be class matrix"))
  if(!(ncol(M)==nrow(M)) ) stop(paste0("the matrix ",name," are not square"))
  if(!all(rownames(M)==colnames(M)) ) stop(paste0(name,": column and row elements must be equal"))

}


#' Choose a path in a graph from an interactive shiny app
#'
#' Choose a path in a graph from an interactive shiny app with the rigth format for simPATHy function.
#' @param graph A graphNEL object.
#' @return Selected path with the rigth format for simPATHy function.
#' @seealso \code{\link[simPATHy]{simPATHy}}
#' @examples
#'  if(require(gRbase)){
#'   graph <- gRbase::dag(~c:a, ~c:b, ~d:c, ~e:d)
#'
#'   # Launch the interactive plot
#'   # path <- getPathShiny(graph)
#'  }
#' @export
getPathShiny<-function(graph){

  appGetPath<-shiny::shinyApp(ui=shiny::fluidPage(

    #### Header ####
    shiny::fluidRow(
      shiny::column(12,
             shiny::p("getPath",shiny::span("simPATHy",style="font-weight:100; color: #66b2ff"),
               style = "text-align: center; font-size: 30px")
      )
    ),

    #### Buttons ####
    shiny::fluidRow(
      shiny::column(12,
             shiny::p(shiny::actionButton(inputId = "get", label = "Get path",width ="100px" ),
               shiny::actionButton(inputId = "exit",label =  "Exit", width ="100px"),
               style="text-align: center")
      )

    ),

    #### Plot ####
    shiny::fluidRow(
      shiny::column(12,
             graphNELD3Output(outputId ="plot")
      )
    )

  ),server =function(input,output){


    shiny::observe({
      #### Exit button behaviour ####
      if(input$exit>0) shiny::stopApp(NULL)

      #### Get button behaviour ####
      if(input$get>0) {
        value<-as.numeric(input$mydata)
        if(length(value)==0){ path<-NULL  }
        else{
          value<-sapply(value,function(v,n=paste(graph::nodes(graph))) n[v])
          path<-list()
          n<-length(input$mydata)/2
          for(i in 1:n) path<-append(path,list(c(value[i],value[n+i])))
        }

        shiny::stopApp( path  )
      }

    })


    #### Plot ####
    output$plot <- renderGraphNELD3({
      graphNELD3(graph = graph,edgeClick = TRUE,colNode = "#99ccff")
    })

  })

  path<-shiny::runApp(appGetPath)
  return(path)
}



#' Dynamic plot of a graph
#'
#' Dynamic plot of a graphNEL object with the possibility to emphasize the strength of relations between nodes, represented by either a pairwise correlation or a partial correlation coefficient. \cr \cr
#' The interactive graph is an implementation of the javascript D3.js package (force-layout) for undirected and directed graphNEL objects (see references).
#' @param graph A graphNEL object.
#' @param S1,S2 Sample covariance matrix. If \code{S1} is supplied edges between nodes are colored in accordance with pairwise correlation or partial correlation coefficients. If \code{S2} supplied, the difference between the two corresponding correlation or partial correlation matrices is plotted.
#' @param type Character string specifying which matrix is to be used. Either cor for correlation matrix, or pcor for partial correlation matrix.
#' @param colLim Numeric vector of length two specifying the lower and upper bound of the color range (see Details).
#' @param legendColor ogical value indicating whether the color legend should be added to the plot.
#' @param colNode A character string specifying the colour of the nodes. The colour node is common for all nodes.
#' @return Dynamic plot of a graphNEL object.
#' @details Admissible values for \code{colLim} are contained in the interval \code{[-1,1]} when \code{S2=NULL}, otherwise the admissible interval is \code{[-2,2]}. When an element is outside of the colLim interval, it is colored gray and represented as a dashed link.
#' @references \url{https://d3js.org} (Micheal Bostock).
#' @references \url{http://www.htmlwidgets.org} (Ramnath Vaidyanathan, Kenton Russell, and RStudio).
#' @references \url{https://christophergandrud.github.io/networkD3/} (Christopher Gandrud, JJ Allaire, & Kent Russell)
#' @examples
#' if(require(gRbase) & require(graph)){
#'   graph <- gRbase::ug(~a:b, ~a:c, ~c:d, ~b:d, ~b:c)
#'   # Plot a graphNEL
#'   plotGraphNELD3(graph)
#'
#'   # Plot a graphNEL coloring edges in correspondance with pairwise correlation coefficients
#'   S <- matrix(c(2,  0.8,0.5,-0.3,
#'                 0.8,1.5,0.6,-0.7,
#'                 0.5,0.6,1,  0.7,
#'                 -0.3,-0.7,0.7,3), ncol=4,nrow=4)
#'   colnames(S) <- rownames(S) <- graph::nodes(graph)
#'   plotGraphNELD3(graph, S1=S)
#'
#'   # Plot a graphNEL coloring edges in correspondance with partial correlation coefficients
#'   plotGraphNELD3(graph, S1=S, type="pcor")
#'
#'   # Change the color range
#'   plotGraphNELD3(graph, S1=S, type="cor", colLim=c(-0.7,0.8))
#'   # Change nodes color
#'   plotGraphNELD3(graph, S1=S, type="cor", colNode = "pink")
#'
#'   # Plot the difference between two graphical models
#'   S2 <- S
#'   S2[1,3] <- S2[3,1]<- -0.1
#'   plotGraphNELD3(graph,S1=S, S2=S2)
#'
#' }
#' @export
plotGraphNELD3<-function(graph,type="cor",S1=NULL,S2=NULL,colLim=c(-1,1),legendColor=TRUE,colNode = "#c0c0c0"){


  ## check0: colLim
  if(is.null(S2)){
    if( !(colLim[1]<= -0.1  & colLim[1]>=-1) ) stop("colLim[1] must be between -0.1 and -1")
    if( !(colLim[2]<= 1  & colLim[2]>=0.1) ) stop("colLim[2] must be between 0.1 and 1")
  } else {
    if( !(colLim[1]<= -0.1  & colLim[1]>=-2) ) stop("colLim[1] must be between -0.1 and -2")
    if( !(colLim[2]<= 2  & colLim[2]>=0.1) ) stop("colLim[2] must be between 0.1 and 2")
  }


  ## check0: graph
  Gtype<-NULL
  if (!(class(graph) == "graphNEL"))
    stop("graph argument is not a graphNEL object")
  else {
    if (gRbase::is.DG(graph)) Gtype <- "dag"
    if (gRbase::is.UG(graph)) Gtype <- "ug"
    if (is.null(Gtype)) stop("graph argument is not a valid object")
  }


  S<-NULL
  ## check2: matrix
  if(!is.null(S1)){

    #colNode = "#c0c0c0"colNode = "#99ccff"
    ## check1: type
    if( !(type %in% c("cor","pcor")) ) stop("Invalid parameter type: must be cor or pcor")
    checkMatrix(S1,name=deparse(substitute(S1)))

    if(!is.null(S2)){
      # S1 - S2
      checkMatrix(S2,name=deparse(substitute(S2)))
      if(type=="cor") S<- riscala(S1) - riscala(S2)
      else S<- riscala(solve(S1)) - riscala(solve(S2))

    } else{
      # only S1
      if(type=="cor") S<- riscala(S1)
      else S<- riscala(solve(S1))
    }
  }

  if(type=="pcor"){
    indDiag<-!(diag(ncol(S))==1)
    S[indDiag]<- -S[indDiag]
  }


  graphNELD3(graph = graph,S = S,limColEdges = colLim,colNode = colNode)

}



#' Visual dysregulation summary
#'
#' A Shiny application for visual summary of dysregulation.
#' @param resObj The output of simPATHy function
#' @param graph The graphNEL object given to the simPATHy function to obtain resObj.
#' @param heightGraph,heightMatrix The height of the graph and correlation matrix plots in pixels. Must be a number, which will be coerced to a string and have 'px' append.
#' @seealso \code{\link[simPATHy]{simPATHy}}, \code{\link[simPATHy]{plotGraphNELD3}}, \code{\link[simPATHy]{plotCorGraph}}, \code{\link[simPATHy]{easyLookDys}}
#' @return Interactive plots for exploring the output of simPATHy.
#' @export
easyLookShiny<-function(resObj,graph,heightGraph=NULL,heightMatrix=NULL){

  #if(!(nrow %in% c(2,3)) ) stop("nrow must be 2 or 3")
  tab<-easyLookDys(resObj)
  tab$strength<-round(tab$strength,2)
  tab[,4:7]<-round(tab[,4:7],4)

  Hm<-NULL; Hg<-NULL
  ## due o tre righe
  row2<-TRUE
  if( !is.null(heightGraph) | !is.null(heightMatrix) ){

    row2<-FALSE

    if(!is.null(heightGraph)) Hg<-paste0(heightGraph,"px") else Hg<-"400px"
    if(!is.null(heightMatrix)) Hm<-paste0(heightMatrix,"px") else Hm<-"400px"
  }



  ## Sidebar
  sidebar <- shinydashboard::dashboardSidebar(

    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("Legend colors",
                               shiny::conditionalPanel(condition = "input.matrix==3",
                                 shiny::sliderInput(inputId ="rpos3",label="Positive range",ticks=FALSE,min = 0.1,max =2 ,value =1 ,step =0.10 ),
                                 shiny::sliderInput(inputId ="rneg3",label="Negative range",ticks=FALSE,min = 0.1,max =2 ,value =1 ,pre="-",step =0.10 )
                               ),
                               shiny::conditionalPanel(condition = "input.matrix==1 | input.matrix==2",
                                  shiny::sliderInput(inputId ="rpos",label="Positive range",ticks=FALSE,min = 0.1,max =1 ,value =1 ,step =0.10 ),
                                  shiny::sliderInput(inputId ="rneg",label="Negative range",ticks=FALSE,min = 0.1,max =1 ,value =1 ,pre="-",step =0.10 )
                               ),
                               shiny::br(),
                               icon = shiny::icon("sliders")
      ),
      shinydashboard::menuItem("Matrix",
                               shiny::radioButtons(inputId = "matrix",label =NULL ,choices =list("S1"=1,"S2"=2,"S1-S2"=3),selected = 1 ),
                               shiny::br(),
                               icon = shiny::icon("th")
      ),
      shinydashboard::menuItem("Type",
                               shiny::radioButtons(inputId = "type",label ="" ,choices =list("Correlation"='cor',"Partial Correlation"='pcor'),selected = 'cor' ),
                               shiny::br(),
                               icon = shiny::icon("gears")
      )
    ),
    shiny::hr(),
    shiny::p(shiny::actionButton(inputId = "exit",label =  "Exit", width ="100px"),style = "text-align: center")
  )

  ## Body
  if( row2 ){
    body<- shinydashboard::dashboardBody(

      shiny::fluidRow(
        shinydashboard::box(graphNELD3Output(outputId = "Graph"),
            title="Graph",status = "primary",solidHeader = TRUE,collapsible = TRUE,width = 6),
        shinydashboard::box(shiny::plotOutput(outputId = "Cor"),
            title="Matrix",status = "primary",solidHeader = TRUE,collapsible = TRUE,width = 6)
      ),
      shiny::fluidRow(

        shinydashboard::box(title="Info",
                            shiny::dataTableOutput(outputId = "DysInfo"),
                            status = "primary",solidHeader = TRUE,collapsible = TRUE,width = 12,collapsed = TRUE)

      )
    )
  } else {

    body<- shinydashboard::dashboardBody(

      shiny::fluidRow(
        shinydashboard::box(graphNELD3Output(outputId = "Graph",height=Hg),
                            title="Graph",status = "primary",solidHeader = TRUE,collapsible = TRUE,width = 12)
        ),
      shiny::fluidRow(
        shinydashboard::box(shiny::plotOutput(outputId = "Cor",height=Hm),
                            title="Matrix",status = "primary",solidHeader = TRUE,collapsible = TRUE,width = 12,collapsed = TRUE)
      ),
      shiny::fluidRow(

        shinydashboard::box(title="Info",
                            shiny::dataTableOutput(outputId = "DysInfo"),
                            status = "primary",solidHeader = TRUE,collapsible = TRUE,width = 12,collapsed = TRUE)

      )
    )

  }

  ## server function
  server<-function(input,output){

    shiny::observe({
      if(input$exit>0) shiny::stopApp()
    })

    output$DysInfo<-shiny::renderDataTable({
      tab
    })

    output$Cor<-shiny::renderPlot({

      if(input$matrix==1) { S1<-resObj$S1; S2<-NULL; lim<-c(-input$rneg,input$rpos)}
      if(input$matrix==2) { S1<-resObj$S2; S2<-NULL; lim<-c(-input$rneg,input$rpos)}
      if(input$matrix==3) { S1<-resObj$S1; S2<-resObj$S2; lim<-c(-input$rneg3,input$rpos3)}
      plotCorGraph(S1 = S1,S2 = S2,type = input$type,colLim = lim ,graph = graph,path = resObj$path)
    })

    output$Graph<-renderGraphNELD3({
      if(input$matrix==1) { S1<-resObj$S1; S2<-NULL; lim<-c(-input$rneg,input$rpos) }
      if(input$matrix==2) { S1<-resObj$S2; S2<-NULL; lim<-c(-input$rneg,input$rpos) }
      if(input$matrix==3) { S1<-resObj$S1; S2<-resObj$S2;  lim<-c(-input$rneg3,input$rpos3) }

      plotGraphNELD3(graph,S1 = S1,S2 = S2,type = input$type,colLim = lim)
    })
  }
  ## call shinyApp
  shiny::shinyApp(ui=shinydashboard::dashboardPage(shinydashboard::dashboardHeader(title = "simPATHy"),sidebar,body), server)
}


