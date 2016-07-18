#' Chimera data
#'
#' A matrix containing the expression values of 3405 genes deriving from Affimetrix single
#' channel technology, consisting of 41 observations from one experimental condition
#' (absence of BCR/ABL gene arrangment, class 1), and 37 observations from another
#' experimental condition (presence of BCR/ABL gene arrangment, class 2).
#' @docType data
#' @keywords datasets
#' @name chimera
#' @usage chimera
#' @format A matrix with 8405 genes (rows) and 78 samples (columns).
#' @source Sabina Chiaretti, Xiaochun Li, Robert Gentleman, Antonella Vitale, Marco Vignetti, Franco Mandelli, Jerome Ritz, and Robin Foa Gene expression profile of adult T-cell acute lymphocytic leukemia identifies distinct subsets of patients with different response to therapy and survival. Blood, 1 April 2004, Vol. 103, No. 7.
#' @examples data(chimera)
NULL

#' Find one path in a graph
#'
#' Find one shortest path in the graph between two given nodes.
#' @param graph A directed or undirected graph represented as a graphNEL object.
#' @param from,to The nodes (character node id) giving the first and the last nodes of the path to be calculated. If \code{NULL} then the \code{from} and \code{to} nodes are randomly choosen.
#' @return A list of edges in edgesList format (see \pkg{gRbase}).
#' @seealso \code{\link[igraph]{get.all.shortest.paths}}
#' @export
generatePath<-function(graph,from=NULL,to=NULL){

  # check graph
  if(!(class(graph)=="graphNEL")) stop("argument graph is not a graphNEL object")
  else {
    if(gRbase::is.DG.graphNEL(graph)) type<-"dag"
    if(gRbase::is.UG.graphNEL(graph)) type<-"ug"

    if(is.null(type)) stop("argument graph is not an admitted object")
  }
  nodes<-graph::nodes(graph)

  cond<-sum(is.null(from),is.null(to))
  if(cond==1) stop("from and to arguments must both be NULL or character")

  #check from & to
  if( sum(is.null(from),is.null(to))==0 ){
    if( !is.character(from) | !is.character(to)) stop("from and to arguments must both be node labels")
    if( !all(c(from,to) %in% nodes) ) stop("from and to arguments must both be node labels contained in graph")
    if( !isNodesConnected(graph,from,to) ) stop("from and to nodes are not connected in graph")
  } else {
    if(type=="dag"){
      A<-gRbase::as.adjMAT(graph)
      from<-nodes[which(apply(A,2,sum)==0)]
      to<-nodes[which(apply(A,1,sum)==0)]
    }
    if(type=="ug"){
      from<-sample(nodes,1)
      to<-findLeafUg(graph,from)
    }
  }

  P<-suppressWarnings(igraph::get.all.shortest.paths(igraph::igraph.from.graphNEL(graph),from,to)$res)
  P<-as.vector(P[[sample(1:length(P),1)]])
  P<-nodes[P]
  np<-length(P)
  path<-strsplit(paste(left=P[1:(np-1)],rigth=P[2:np],sep="~"),split="~")


  return(path)
}


#' Estimate covariance matrix of a graphical model
#'
#' Fit a Gaussian Graphical Model or a Gaussian Bayesian Network  by maximum likelihood.
#' @param graph A directed or undirected graph represented as a graphNEL object.
#' @param S A sample covariance matrix
#' @details If graph is undirected it uses the Iterative Proprotional Fitting algoritm (\pkg{qpgraph} package). If graph is directed it uses Iterative Conditional Fitting (\pkg{ggm} package).
#' @return A covariance matrix with the independence constraints entailed by the graph.
#' @references Drton, M. & Richardson, T. S. (2003). A new algorithm for maximum likelihood estimation in Gaussian graphical models for marginal independence. Proceedings of the Ninetheen Conference on Uncertainty in Artificial Intelligence, 184-191.
#' @references Whittaker, J. Graphical models in applied multivariate statistics. Wiley, 1990.
#' @seealso \code{\link[ggm]{icfmag}},  \code{\link[qpgraph]{qpIPF}}
#' @export
fitSgraph<-function(graph,S){

  if(!(class(graph)=="graphNEL")) stop("graph argument is not a graphNEL object")
  else {
    if(gRbase::is.DG.graphNEL(graph)) type<-"dag"
    if(gRbase::is.UG.graphNEL(graph)) type<-"ug"
    if(gRbase::is.TUG.graphNEL(graph)) type<-"tug"


    if(is.null(type)) stop("argument graph is not an admitted object")
  }

  if(type=="dag"){
    ICF<-ggm::icfmag(gRbase::as.adjMAT(graph),S)
    S<-ICF$Sigmahat
  }

  if(type=="ug"){
    names<-colnames(S)
    string<-R.utils::captureOutput(S<-qpgraph::qpIPF(vv = S,clqlst = qpgraph::qpGetCliques(graph),verbose=FALSE))
    colnames(S)<-rownames(S)<-names
  }

  if(type=="tug"){
    S<-SMLEdecomposable(S,graph)
  }

  return(S)
}



#' Positive definite matrix
#'
#' Adjust the diagonal of a symmetric square matrix, by the smallest eigenvalue method, in order to make it positive definite.
#' @param M1,M2 A squared numeric matrix, typically a correlation or a covariance matrix. It must be symmetric.
#' @param threshold A correction factor.
#' @details Finds the smallest eigenvalue lambda of \code{M1} (or \code{M1} and \code{M2} if supplied) and adds (threshold-lambda) to the diagonal to make it positive definite.
#' @return A list with the corrected input matrices and the correction \code{threshold}-lambda.
#' @export
makePositiveDefinite<-function(M1,M2=NULL,threshold=0.1){

  if(is.null(M2)) {
    eig<-min( round(eigen(M1)$values,2)  )
  } else {

    if(!(ncol(M1)==nrow(M1)) | !(ncol(M2)==nrow(M2))) stop("the matrices are not square")
    if(!all(rownames(M1)==colnames(M1)) | !all(rownames(M2)==colnames(M2))  ) stop("column and row elements must be equal")
    if( !all(dim(M1)==dim(M2)) | !all(colnames(M1)==colnames(M2))  ) stop("M1 and M2 have different sizes")

    eig<-min(  round(c(eigen(M1)$values,eigen(M2)$values),2)  )
  }

  if( !(eig>0) )
  {
    n<-ncol(M1)

    M1n<-M1+(threshold-eig)*diag(n)
    M2n<-if(!is.null(M2)) M2+(threshold-eig)*diag(n) else NULL

    correction<- TRUE
    value<-threshold-eig

  } else {
    M1n<-M1; M2n<-M2
    correction<-FALSE
    value<-NULL
  }

  return(list(M1=M1n,M2=M2n,correction=correction,value=value))

}

#' Local Maximum Likelihood Estimation
#'
#' Compute a maximum likelihood estimate of a covariance matrix in a decomposable Gaussian graphical model.
#' @param S a covariance matrix.
#' @param graph a decomposable graph represented as a graphNEL object.
#' @references Lauritzen, S. L. (1996). Graphical Models. Clarendon Press, Oxford.
#' @return The MLE of a covarince matrix.
#' @export
SMLEdecomposable<-function(S,graph){

  if( !(class(graph)=="graphNEL") ) stop("graph argument is not a graphNEL object")
  if( !gRbase::is.TUG.graphNEL(graph) ) stop("graph argument must be decomposable")

  nodes<-graph::nodes(graph)


  cond<-all(nodes %in% colnames(S))
  if(!cond) stop("all of the graph variables must be contained in the covariance matrix")
  if( sum(!cond)>0 ) S<-S[cond,cond]


  ripped<-gRbase::rip(graph)
  cliques<-ripped$cliques
  names(cliques)<-paste0("C",1:length(cliques))

  separators<-ripped$separators
  separators<-separators[-which(sapply(separators,length)==0)]
  if(length(separators)>0) names(separators)<-paste0("S",1:length(separators))

  Scc<-lapply(append(cliques,separators),function(cl,ss=S) ss[cl,cl,drop=FALSE] )
  Pcc<-lapply(Scc,function(cl) solve(cl))

  Blocks<-lapply(Pcc,getBlock,nodes=nodes)
  indC<-grep("C",names(Blocks))
  indS<-grep("S",names(Blocks))

  if( length(indS)>0 ) P<-Reduce("+",Blocks[indC])-Reduce("+",Blocks[indS])
  else P<-Reduce("+",Blocks[indC])

  return( solve(P) )

}


#' Simulate data from a graphical model
#'
#' Simulate data in two different conditions with a common structure of dependences. The two different conditions are characterized by different strengths of the links between nodes (dysregulation).
#' @param graph A graphNEL object.
#' @param S The sample covariance matrix.
#' @param n1,n2 Number of observations to generate from the two conditions.
#' @param min,max Vectors of length 1 or of the same length as \code{path} containing the lower and upper limits of a uniform distribution.
# The parameters define the strength of the dysregulation for each edges.
#' The strength of dysregulation is sampled uniformly from the interval [\code{min}, \code{max}]: a value smaller than 1 represents deactivation, a value greater than 1 represents activation. If \code{path=NULL} only the first element is used.
#' @param prob A vector of size 1 or of the same length as path, giving the probability to change the sign of the correlation coefficient for each edge.
#' \code{prob=0} implying that the sign of the dysregulation should be changed, and \code{prob=1} implying that the sign should be left unaltered (default). Values between these two extremes allow for random sign switch: the sign is changed with probability \code{1-prob}.
#' @param path A list of edges in edgesList format (see \pkg{gRbase}).
#' @param digits Integer indicating the number of decimal places to be used.
#' @param mu1,mu2 A vector of size 1 or of the length equal to the number of nodes in the graph. Means of the multivariate normal distributions from which observations are generated. If \code{mu1} (and/or \code{mu2}) is a vector it has to be named in accordance with the names of the nodes of the graph.
#' @param muRandom Logical. If \code{muRandom=TRUE} the means of the variables are randomly generated.
#' @return It returns a list containing:
#' \itemize{
#'  \item{\code{data}}       {random samples generated from multivariate normal distributions with covariance matrices \code{S1} (reference condition) and \code{S2} (dysregulated condition);}
#'  \item{\code{S1,S2}}      {two covariance matrices;}
#'  \item{\code{path}}       {the dysregulated path;}
#'  \item{\code{strength}}   {the dysregulation strength for each edge in the path;}
#'  \item{\code{mu1,mu2}}    {two mean vectors;}
#'  \item{\code{correction}} {correction details.}
#' }
#' @details
#' If the matrix \code{S} does not reflect conditional independence constraints imposed by the graph \code{simPATHy} uses the maximum likelihood estimation of covariance matrices for graphical models via internal function \code{\link[simPATHy]{fitSgraph}}.\cr \cr
#' When the dysregulation of the initial (reference condition) covariance matrix leads to a matrix that is no longer positive definite, the resulting matrix is corrected via internal function \code{\link[simPATHy]{makePositiveDefinite}}.\cr \cr
#' To avoid excessively strong dysregulations, the upper limit for the absolute value of the dysregulated correlation coefficient is set to: \deqn{min( 0.9, 1.25*max(abs(C[upper.tri(C)])) )} where C is the correlation matrix of the reference condition.
#' @seealso \code{\link[simPATHy]{easyLookDys}}, \code{\link[simPATHy]{easyLookShiny}}, \code{\link[simPATHy]{plotCorGraph}},  \code{\link[simPATHy]{plotGraphNELD3}}
#' @examples
#' if(require(gRbase) & require(graph)){
#'
#'   ## Directed graph
#'   ## sub-graph Acute Myel... Leukemia
#'   graph<-gRbase::dag(~867:25+867:613+5295:867+5294:867+
#'                        + 207:5295+207:5294+4193:207+3551:207+
#'                        + 4792:3551+7157:4193+3265:6654+
#'                        + 3845:6654+6654:2885+2885:25+2885:613)
#'   genes<-graph::nodes(graph)
#'
#'   # covariance matrix of the reference condition
#'   data<-t(chimera[genes,colnames(chimera)==1])
#'   S<-cov(data)
#'   S<-fitSgraph(graph,S)
#'
#'   # select path to dysregulate
#'   path<-list(c("613","867"),c("867","5295"),c("5295","207"),
#'              c("207","4193"),c("4193","7157"))
#'   ## ..or select the path in an interactive plot
#'   # path<-getPathShiny(graph)
#'
#'   # select parameters of the dysregulation
#'   min<-c(2,8,2,0.1,0.5)
#'   max<-c(2,10,2,4,0.5)
#'   prob<-c(1,0,0,0.5,1)
#'
#'   # activation, switch, switch, random, deactivation
#'   dys<-cbind(min,max,prob)
#'   rownames(dys)<-sapply(path,paste,collapse = "~")
#'   dys
#'
#'   set.seed(123)
#'   # main function
#'   Result<-simPATHy(graph,path,S,min,max,prob)
#'   class(Result)
#'   names(Result)
#'
#'   # simulated data from two conditions
#'   round(Result$dataset[c(1:3,501:503),1:5],3)
#'
#'   # Summary
#'   easyLookDys(Result)
#'   # ..or interactive summary
#'   # easyLookShiny(resObj=Result,graph=graph)
#'
#'
#'   # Visualization
#'   plotCorGraph(S1=Result$S1,S2 = Result$S2,graph = graph,path = path,colLim = c(-0.3,0.3))
#'   plotGraphNELD3(S1=Result$S1,S2 = Result$S2,graph = graph,colLim = c(-0.3,0.3))
#'   rm(list=ls())
#'
#'
#'   ## Undirected graph
#'   graph <- gRbase::ug(~a:b, ~a:c, ~c:d, ~b:d, ~b:c)
#'   # when reference condition covariance matrix is not supplied simPATHy generate a random one
#'   Result_ug<-simPATHy(graph)
#'   easyLookDys(Result_ug)
#'   plotGraphNELD3(S1=Result_ug$S1,S2 = Result_ug$S2,graph = graph,colLim = c(-0.5,0.5))
#'
#' }
#' @export
simPATHy<-function(graph,path=NULL,S=NULL,min=2,max=3,prob=1,n1=500,n2=n1,digits=5,mu1=0,mu2=mu1,muRandom=FALSE ){

  type<-NULL
  force(mu2);
  #force(min); force(max); force(prob)
  ### step 0: checks ###
  # check1 -> min,max,prob
  if( !all(sapply(c(min,max,prob),class)=="numeric") ) stop("invalid type: min,max and prob must be numeric")
  if( sum(min>max)>0 ) stop("invalid argument: min must be less then max")
  if( !all(c(min,max)>0) ) stop("invalid argument: min and max must be greater than one (activation or deactivation must be set with prob argument)")
  if( sum(prob<0 | prob>1)>0 ) stop("invalid argument: prob argument must be a probability")

  # check2 -> graph
  if(!(class(graph)=="graphNEL")) stop("graph argument is not a graphNEL object")
  else {
    if(gRbase::is.DG.graphNEL(graph)) type<-"dag"
    if(gRbase::is.UG.graphNEL(graph)) type<-"ug"
    if(is.null(type)) stop("graph argument is not a valid object")
  }

  # check3 -> path
  nodes<-graph::nodes(graph)
  indEdges<-which(gRbase::as.adjMAT(graph)==1,arr.ind = TRUE)
  edList<-as.list(data.frame(apply(indEdges,1,function(x,n=nodes) nodes[x] ),stringsAsFactors = FALSE))
  names(edList)<-NULL


  ### step 000: mean ###
  if( !muRandom ){

    ## mu1
    if( length(mu1)==1 ) { mu1<-rep(mu1,length(nodes)); names(mu1)<-nodes }
    else {
      if( !(length(mu1)==length(nodes)) ) stop("all of the graph variables must be contained in mu1")
      indMu1<-match(names(mu1), nodes)
      # repeated names or not in graph
      if( !all( table(indMu1)==1 )  |  all( is.na(indMu1) ) ) stop("insuitable mu1 vector names")
    }

    ## mu2
    if( length(mu2)==1 ) { mu2<-rep(mu2,length(nodes)); names(mu2)<-nodes}
    else {

      if( !(length(mu2)==length(nodes)) ) stop("all of the graph variables must be contained in mu2")
      indMu2<-match(names(mu2), nodes)
      # repeated names or not in graph
      if( !all( table(indMu2)==1 )  |  all( is.na(indMu2) ) ) stop("insuitable mu2 vector names")
    }

  } else { if( !all(mu1==0) | !all(mu2==0) ) warning("Means randomly generated") }


  if( !(class(path) %in% c("list","NULL")) ) stop("path argument must be a list object or set to NULL: see edgeList format")
  else if( !all(unique(unlist(path)) %in% nodes) ) stop("all of the path variables must be contained in the graph")

  # check4 -> path and param compatibility
  param<-list(min=min,max=max,prob=prob)
  if(is.null(path)){
    param<-lapply(param,function(x) x[1])
    path<-generatePath(graph)
  } else param<-lapply(1:length(param),function(i,x=param,p=path){
       if( !(length(x[[i]])==1) && !(length(x[[i]])==length(p)) ) {
         warning(paste0(names(x)[i]," and path must have the same size (only the first element is used)"))
         x[[i]][1]
       } else x[[i]]
     })
  names(param)<-c("min","max","prob")

  # check5 -> invalid path edges
  param<-lapply(param,function(x,p=path) if(length(x)==1) rep(x,length(p)) else x )
  invalid<- which(!(path %in% edList))

  if( length(invalid)>0 ) {
    if( length(invalid)==length(path) ) stop("all of the path variables are not contained in the graph")
    else {
      warning(paste0("n.",length(invalid)," edge(s) is not in graph"))
      param<-lapply(param,function(x,ind=invalid) x[-ind])
      path<-path[-invalid]
    }
  }

  or_path<-path
  # one param for each unique edge
  U<-uniqueParam(nodes,path,param)
  path<-U$path
  param<-U$param
  rm(U)

  # check6 -> S
  if(!is.null(S)){
    cond<-all(nodes %in% colnames(S))
    if(!cond) stop("all of the graph variables must be contained in S matrix")
    if( sum(!cond)>0 ) S<-S[cond,cond]
    mpd<-makePositiveDefinite(S)
    S<-mpd$M1
    if(mpd$correction) warning("S must be positive definitive: forced, see makePositiveDefinitive() ")

    if(muRandom){
      rr<-randomS(nodes,mean = muRandom)
      mu1<-mu2<-rr$mu
      rm(rr)
    }

  } else {
    rr<-randomS(nodes,mean = muRandom)
    S<-rr$S
    if(muRandom) mu1<-mu2<-rr$mu
    rm(rr)
  }




  ### step 1: ###
  # if dag-> icf; if ug->ipf
  S<-fitSgraph(graph,S)

  if(type=="dag") structure<-list(graph=gRbase::triangulate(gRbase::moralize(graph)),Dgraph=graph)
  else structure<-list(graph=gRbase::triangulate(graph),Dgraph=NULL)


  ### step 4: Sregolo....
  D<-Dysregulate(S,path,param$min,param$max,param$prob)

  ### step 5: Rendo Definito positivo
  E<-makePositiveDefinite(S,D$dysS)

  S1n<-round(SMLEdecomposable(E$M1,structure$graph),digits)
  S2n<-round(SMLEdecomposable(E$M2,structure$graph),digits)

  S1n<-fitSgraph(graph,S1n)
  S2n<-fitSgraph(graph,S2n)

  # match name mu1 and mu2
  mu1<-mu1[match(colnames(S1n),names(mu1))]
  mu2<-mu2[match(colnames(S2n),names(mu2))]

  data1<-mvtnorm::rmvnorm(n1,mean=mu1,sigma = S1n)
  data2<-mvtnorm::rmvnorm(n2,mean=mu2,sigma = S2n)
  colnames(data1)<-colnames(S1n)
  colnames(data2)<-colnames(S2n)
  rownames(data1)<-rep("cl1",n1)
  rownames(data2)<-rep("cl2",n2)
  data<-rbind(data1,data2)

  ## return original order path
  ind_order<-sapply(path,function(lp,op=or_path){
    var<-rep(FALSE,length(op))
    for(i in 1:length(op)) var[i]<- all(lp %in% op[[i]])
    which(var)[1]
  })
  path<-or_path[ind_order]

  path<-path[order(ind_order)]
  param<-lapply(param,function(x,ii=ind_order) x[order(ii)])


  result<-list( dataset=data, S1=S1n, S2=S2n, path=path, strength=D$strength[order(ind_order)], param=param, correction=list(isCorrected=E$correction,correction=E$value) , mu1=mu1, mu2=mu2    )
  class(result)<-"simPATHy"
  return(result)
}



#' Dysregulation summary
#'
#' Summary of the result for a quick look of simPATHy function.
#' @return Nice formatted output of simPATHy dysregulation
#' @param resObj The output of simPATHy function (simPATHy class object).
#' @param digits Integer indicating the number of decimal places to be used.
#' @return Nicely formatted output of simPATHy dysregulation.
#' @export
easyLookDys<-function(resObj,digits=4){
  if(!(class(resObj)=="simPATHy") ) stop("invalid class: resObj must be a simPATHy object")

  Type<-function(x) if(x>0) if(x<1) return("deactivation") else if(x==1) return("unchanged") else return("activation") else return("switch")

  E<-sapply(resObj$path,paste,collapse="~")
  ind<-matrix(unlist(resObj$path),byrow=TRUE,ncol=2)

  typ<-sapply(resObj$strength,Type)

  C1<-riscala(resObj$S1)
  C2<-riscala(resObj$S2)

  data.frame(edge=E,type=typ,strength=round(resObj$strength,digits = digits),"cov.S1"=round(resObj$S1[ind],digits=digits),"cov.S2"=round(resObj$S2[ind],digits=digits),"cor..S1"=round(C1[ind],digits = digits),"cor..S2"=round(C2[ind],digits = digits))
}

