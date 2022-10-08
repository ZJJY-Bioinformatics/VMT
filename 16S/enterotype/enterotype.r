
library(cluster)
library(clusterSim)
library(ade4)
library(vegan)
library(tidyverse)

data=read.table("all.r.input.txt", header=T, row.names=1, dec=".", sep="	")
data=data[-1,]



dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)

  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) {
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix)
}



data.dist=dist.JSD(data)                                                                                               

sink("all.dist.txt")
data.dist
sink()

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters                                   
  require(cluster)                                                                                                     
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)                                                        
  return(cluster)                                                                                                      
}   


require(clusterSim)                                                                                                                                        
                                                                                                                       
nclusters=NULL                                                                                                         
                                                                                                                       
for (k in 1:20) {                                                                                                      
  if (k==1) {                                                                                                          
    nclusters[k]=NA                                                                                                    
  } else {                                                                                                             
    data.cluster_temp=pam.clustering(data.dist, k)                                                                     
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,                                                   
                          centrotypes = "medoids")                                                                     
  }                                                                                                                    
} 
k = which.max(nclusters)
write.table(file="all.nclusters.txt",nclusters,sep="	",quote=FALSE)

data.cluster=pam.clustering(data.dist, k=k)  
write.table(file=paste0("all.", k, "_clusters.data.cluster.txt"),data.cluster,quote=FALSE)

pdf("all.nclusters.pdf")
plot(nclusters, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters")                        
dev.off()                                                                                                                       
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,k])                                                           
write.table(file=paste0("all.", k, "_clusters.silhouette.txt"),obs.silhouette)                                                                                         
                                                                                                                       
#data=noise.removal(data, percent=0.01)

#colour-set
coul=c("red","blue","green")                                                                                           

## plot 1                                                                                                              
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)                                                                 
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=2)                                                    
write.table(file=paste0("all.", k, "_clusters.tab.txt"),obs.bet$tab,sep="\t",quote=FALSE)

pdf(paste0("all.", k, "_clusters.bca.pdf"))                                                                                                              
s.class(obs.bet$ls, fac=as.factor(data.cluster), col=coul, grid=F,sub="Between-class analysis")                        
dev.off()

obs.pcoa=dudi.pco(data.dist, scannf=F, nf=2)                                                                           
pdf(paste0("all.", k, "_clusters.pcoa.pdf"))                                                                                                              
s.class(obs.pcoa$li, fac=as.factor(data.cluster), col=coul, grid=F,sub="Principal coordiante analysis")        
dev.off()
