## ----echo=FALSE,include=FALSE-------------------------------------------------
library(simPATHy)
library(graph)

## ----eval=FALSE---------------------------------------------------------------
#  library(simPATHy)
#  library(graph)

## -----------------------------------------------------------------------------
data("chimera")
dim(chimera)
table(colnames(chimera))

## -----------------------------------------------------------------------------
graph<-gRbase::dag(~867:25+867:613+5295:867+5294:867+
207:5295+207:5294+4193:207+3551:207+
4792:3551+7157:4193+3265:6654+
3845:6654+6654:2885+2885:25+2885:613)

## -----------------------------------------------------------------------------
genes<-graph::nodes(graph)
data<-t(chimera[genes,colnames(chimera)==1])
S<-cov(data) 

## -----------------------------------------------------------------------------
S<-fitSgraph(graph,S)
round(S[1:5,1:5],3)

## ---- eval=FALSE--------------------------------------------------------------
#  plotGraphNELD3(graph,type = "cor",S1 = S)

## ---- fig.height=5, fig.width=7,  fig.align='center'--------------------------
plotCorGraph(S1 = S,type = "cor")

## ---- fig.height=5, fig.width=7,  fig.align='center'--------------------------
lim<-round(max(abs(simPATHy:::riscala(S))[upper.tri(S)]),2)
#plotGraphNELD3(graph,type = "cor",S1 = S,colLim = c(-lim,lim))
plotCorGraph(S1 = S,type ="cor",colLim = c(-lim,lim))

## ---- fig.height=5, fig.width=7,  fig.align='center'--------------------------
plotCorGraph(S1 = S,type = "cor", graph = graph)

## -----------------------------------------------------------------------------
path<-list(c("613","867"),c("867","5295"),c("5295","207"),
             c("207","4193"),c("4193","7157"))

## ---- eval=FALSE--------------------------------------------------------------
#  path <- getPathShiny(graph)

## -----------------------------------------------------------------------------
path

## ---- eval=FALSE--------------------------------------------------------------
#  path <- generatePath(graph,from="613",to="7157")

## ---- eval=FALSE--------------------------------------------------------------
#  path1 <- list(c("613","2885"), c("4193","7157"))

## -----------------------------------------------------------------------------
min<-c(2,8,2,0.1,0.5)
max<-c(2,10,2,4,0.5)

## -----------------------------------------------------------------------------
prob<-c(1,0,0,0.5,1)
dys<-cbind(min,max,prob)
rownames(dys)<-sapply(path,paste,collapse = "~")
dys

## -----------------------------------------------------------------------------
set.seed(123)
Result<-simPATHy(graph,path,S,min,max,prob)

## -----------------------------------------------------------------------------
class(Result)
names(Result)

## -----------------------------------------------------------------------------
round(Result$dataset[c(1:3,501:503),1:5],3)

## -----------------------------------------------------------------------------
Result$param

## -----------------------------------------------------------------------------
Result$strength

## -----------------------------------------------------------------------------
Result$correction

## ----eval = FALSE-------------------------------------------------------------
#  easyLookDys(Result)

## ----echo= FALSE--------------------------------------------------------------
knitr::kable(easyLookDys(Result))

## ---- fig.height=5, fig.width=7,  fig.align='center'--------------------------
#plotGraphNELD3(graph,type = "cor",S1 = Result$S1, S2 = Result$S2, colLim = c(-0.4,0.4))
plotCorGraph(S1 = Result$S1, S2 = Result$S2, type = "cor",
graph = graph, path = Result$path,colLim = c(-0.4,0.4))

## ---- eval=FALSE--------------------------------------------------------------
#  easyLookShiny(Result, graph)

## ---- echo = T, results = 'hide'----------------------------------------------
#library(topologyGSA)
#?pathway.var.test
y1<-Result$dataset[rownames(Result$dataset)=="cl1",]
y2<-Result$dataset[rownames(Result$dataset)=="cl2",]
alpha<-0.05
#pathway.var.test(y1, y2, dag = graph, alpha)

## -----------------------------------------------------------------------------
library(clipper)
#?clipper
expr<-t(Result$dataset)
classes<-as.numeric(gsub("cl","",colnames(expr)))
clipped<-clipper(expr,classes,graph)
clipped

