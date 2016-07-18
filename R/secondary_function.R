# Scale symmetrix matrix
#
# Scales a covariance matrix into the corresponding correlation matrix efficiently
# @param mat symmetric numeric matrix, usually positive definite such as a covariance matrix
# @return A scaled matrix
riscala<-function(mat) {
  Scala<-sqrt(diag(mat) %*% t(diag(mat)))
  mat/Scala
}

# Connected nodes
#
# Check if two nodes are connected in a graph
# @param graph A graph rapresented as graphNEL (from the graph package)
# @param node1,node2 the nodes to check
isNodesConnected<-function(graph,node1,node2){
  if( !all(c(node1,node2) %in% graph::nodes(graph)) ) stop("node1 or node2 are not in graph")
  ans<-suppressWarnings(igraph::get.shortest.paths(igraph::igraph.from.graphNEL(graph),node1,node2))
  length(paste(ans$vpath[[1]]))>0
}


# Find a leaf nodes in an undirected graph
#
# Given a root node, it finds a leaf node in the assciated decomposable model of the graph
# @param graph an undirected graph represented as a graphNEL object
# @param root a variable
# @return It return a variable.
findLeafUg<-function(graph,root){
  ripped<-gRbase::rip(gRbase::triangulate(graph),root = root)

  if(length(ripped$cliques)>1){
    childs<-ripped$children[1]
    while( !is.na(childs[length(childs)]) ){
      childs<-c(childs,ripped$children[childs[length(childs)]])
    }
    leaf<-childs[length(childs)-1]
    to<-sample(ripped$cliques[[leaf]],1)
  } else to<- sample(ripped$cliques[[1]][-which(ripped$cliques[[1]]==root)],1)

  return(to)
}


# Return 0 matrix with block
#
# to do....
# @param scl a symmetric positive definite submatrix
# @param nodes the row and col elements of the original matrix
# @return A matrix with all elements set to 0 except for the element in scl
getBlock<-function(scl,nodes){
  cl<-colnames(scl)
  if( !all(cl %in% nodes) ) stop("Not all elements are contanined in nodes")

  A<-matrix(0,nrow=length(nodes),ncol=length(nodes),dimnames = list(nodes,nodes))
  A[cl,cl]<-scl
  A
}

# Generate a random covariance matrix
#
# to do...
# @param nodes a vector with the name of the variables
# @param mean if set TRUE return mean nodes
# @details The maximal number of variables supported is ...
# @return A covariance matrix
randomS<-function(nodes,mean=FALSE){

  classes<-colnames(simPATHy::chimera)
  ## Data: dataset di default (come metterli nel pacchetto)
  nr<-nrow(simPATHy::chimera); nc<-ncol(simPATHy::chimera)

  if(length(nodes)>nr) stop("too much nodes to generate a random covariance matrix")

  randNodes<-sample(1:nr,length(nodes),replace=FALSE)

  randData<-simPATHy::chimera[randNodes,classes==1]
  rownames(randData)<-nodes

  if(mean) mu<-apply(randData,1,mean)
  else mu<-NULL

  S<-stats::cov(t(randData))
  S<-makePositiveDefinite(S)$M1

  return(list(S=S,mu=mu))
}

# Check path
#
# It check if there are repeated edges in path and adjust the parameters
# @param nodes A character vector of node labels
# @param path A list of edges in edges format (see graphNEL)
# @param param A list with numeric vector elements named min, max and prob. Min, max and prob must have the same length of path list or length==1
# @return A list with correct parameters
uniqueParam<-function(nodes,path,param){

  Min<-Max<-Prob<-matrix(0,ncol=length(nodes),nrow=length(nodes),dimnames = list(nodes,nodes))
  A<-matrix(0,ncol=length(nodes),nrow=length(nodes),dimnames = list(nodes,nodes))

  ind.path<-t(data.frame(path,stringsAsFactors = FALSE))
  rownames(ind.path)<-NULL

  for(i in 1:length(path)){

    Min[ind.path[i,,drop=FALSE]]<-Min[ind.path[i,,drop=FALSE]]+param$min[i]
    Max[ind.path[i,,drop=FALSE]]<-Max[ind.path[i,,drop=FALSE]]+param$max[i]
    Prob[ind.path[i,,drop=FALSE]]<-Prob[ind.path[i,,drop=FALSE]]+param$prob[i]
    A[ind.path[i,,drop=FALSE]]<-A[ind.path[i,,drop=FALSE]]+1
  }

  A2<- A+t(A)

  nmultiple<-which(A2[upper.tri(A2)]>1)
  if(length(nmultiple)>0) warning("An edge can not be reapeted in the path: param mean computes")

  Min<- (Min+t(Min))/A2
  Max<- (Max+t(Max))/A2
  Prob<- (Prob+t(Prob))/A2

  ind<-which(!(A2==0),arr.ind = TRUE)
  ind<-ind[ind[,1]<ind[,2],,drop=FALSE]

  prob<-Prob[ind]
  min<-Min[ind]
  max<-Max[ind]

  P<-apply(ind,1,function(x,n=nodes) nodes[x] )
  P<-as.list(data.frame(P,stringsAsFactors = FALSE))
  names(P)<-NULL

  list(path=P,param=list(min=min,max=max,prob=prob))
}


# Dysregulate a covariance matrix
#
# It activate or deactivate a path in a graph
# @param S A covariance matrix
# @param path A list of edges in edges format (see graphNEL)
# @param min,max a number or a numeric vector containing the (lower and upper) limits of the uniform distribution. The parameters define the strength of the dysregulation for each edges.
# @param prob prov A vector of probability weights for obtaining an activation (prob) or a deactivation (1-prob)
# @return A list with the dysregulate matrix and the dysregulation strengths.
Dysregulate<-function(S,path,min,max,prob,tolerance=0.25){

  nodes<-colnames(S)
  variance<-diag(S)

  C<-riscala(S)
  bound<-max(abs(stats::quantile(C[upper.tri(C)],c(0,1))))
  bound<-min(0.9,(1+tolerance)*bound)

  index<-t(sapply(path,function(el,nn=nodes) match(el,nodes)))
  strength<-sapply(1:length(path),function(i,m1=min,m2=max,p=prob,q=c(1-prob))  sample(c(-1,1),1,prob = c(q[i],p[i]))*stats::runif(1,m1[i],m2[i])   )

  dys<- C[index]*(strength)

  strength[abs(dys)>bound]<- sign(dys[abs(dys)>bound])* bound / C[index][abs(dys)>bound]
  dys[abs(dys)>bound]<-sign(dys[abs(dys)>bound])*bound

  dysS<-C
  dysS[index[,c(1,2),drop=FALSE]]<-dysS[index[,c(2,1),drop=FALSE]]<-dys

  dysS<-dysS * sqrt(variance %*% t(variance))
  return(list(dysS=dysS,strength=strength))
}












