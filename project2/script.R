

library(cluster)
library(MASS)
library(corrplot)

#preparazione tabella

pa=read.csv2("tabella.csv", row.names=1)
head(pa)

corrplot(round(cor(pa),2), "number", tl.cex= 0.9)
pa=data.frame(scale(pa))

#pca
pa.pca=princomp(pa)
summary(pa.pca)
# Scree plot
library(ggplot2)
variance = pa.pca$sdev^2 / sum(pa.pca$sdev^2)
qplot(c(1:12),variance, label= round(variance))+geom_col(fill="dodgerblue3")+geom_line()
+geom_text(label = scales::percent(variance),vjust = -0.25,hjust=-0.05)+xlab("Principal Component")+ylab("Variance Explained")+ggtitle("Scree Plot")
plot(cumsum(pa.pca$sdev^2)/sum(pa.pca$sdev^2), type="b",  ylab="Proportion of variance", ylim=c(0,1))
segments(1,0.8,7,0.8,col="red")
pa.ld=loadings(pa.pca)
corrplot(t(pa.ld[,1:3]),"number")
corrplot(t(varimax(pa.ld[,1:3])$loadings),"number")
layout(t(1:2))
biplot(pa.pca,col=c("gray","red"))
biplot(pa.pca,col=c("gray","red"), choices = c(1,3), cex=.8)
layout(1)

#wss
wss=rep(0,10)
for(k in 2:10){
  wss[k]=kmeans(pa,k,nstart=20)$tot.withinss
}
plot(2:10,wss[2:10],type="b",pch=20, xlab= "number of cluster", ylab="wss")
segments(5.5,700,5.5,1600,col="red")
title("within-cluster sum of squares")

#silhouette media kmeans + pam euclidea + pam manhattan

as=rep(0,9)
for(k in 2:10){
  cl=kmeans(pa,k,nstart=20)$cluster
  as[k]=mean(silhouette(cl,dist(pa))[,3])
}
plot(2:10,as[2:10],type="b",pch=20, ylim=c(0.15,0.41), xlab="number of cluster", ylab="average silhouette")
c=rep(0,10)
for(i in 2:10){
  c[i]=pam(pa,i)$silinfo$avg.width
}
points(2:10,c[2:10],type="b",pch=19,col = "blue")
c=rep(0,10)
for(i in 2:10){
  c[i]=pam(pa,i,metric="manhattan")$silinfo$avg.width
}
points(2:10,c[2:10],type="b",pch=19, col= "green3")
legend(8, 0.4, legend=c("kmeans", "pam euclidea", "pam manhattan"), col=c("black", "blue", "green3"), lty=1:2, cex=1)
title("Silhouette media" )
segments(5.78,0.15,5.78,0.40,col="red")

#silhouette media metodi gerarchici

pa.hc=hclust(de) # complete linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,de)[,3])
}
plot(2:10,as[2:10], pch=20, type="b", ylim=c(0.12,0.51), xlab="number of cluster", ylab="average silhouette", col="red2")
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
pa.hc=hclust(de,method="single") # single linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,de)[,3])
}
points(2:10,as[2:10],pch=20, type="b", lty=2, lwd=2, col= "red3")
pa.hc=hclust(de,method="average") # average linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,de)[,3])
}
points(2:10,as[2:10],pch=20, type="b", lty=2, col= "red4")
pa.hc=hclust(dmax) # complete linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dmax)[,3])
}
points(2:10,as[2:10],pch=20, type="b", col= "dodgerblue2")

pa.hc=hclust(dmax,method="single") # single linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dmax)[,3])
}
points(2:10,as[2:10],pch=20, type="b", col= "dodgerblue3", lty=2, lwd=2)

pa.hc=hclust(dmax,method="average") # average linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dmax)[,3])
}
points(2:10,as[2:10],pch=20, type="b", col= "dodgerblue4",lty=2)

pa.hc=hclust(dman) # complete linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dman)[,3])
}
points(2:10,as[2:10],pch=20, type="b", col= "green2")

pa.hc=hclust(dman,method="single") # single linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dman)[,3])
}
points(2:10,as[2:10],pch=20, type="b", col= "green3",lty=2, lwd=2)

pa.hc=hclust(dman,method="average") # average linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dman)[,3])
}
points(2:10,as[2:10],pch=20, type="b", col= "green4",lty=2)

legend("topright", legend=c("euclidean complete", "euclidean simple", "euclidean average", "maximum complete", "maximum simple","maximum average","manhattan complete","manhattan simple","manhattan average"), 
       col=c("red2","red3","red4", "dodgerblue2", "dodgerblue3", "dodgerblue4", "green2","green3","green4"), lty=1:9, inset=c(-0.34,0),
       bty = "n")
title("Silhouette media" )
segments(5.78,0.11,5.78,0.50,col="red")


#metodi gerarchici: qui plotto tutte le silhouette con le varie distanze e inserisco in due matrici la grandezza minima/massima dei cluster
#poi plotto l'stogramma che mi mostra la grandezza minima del cluster nel caso di 3 cluster (come nella relazione)

min.cl=matrix(rep(0,48),12,4)
max.cl=matrix(rep(0,48),12,4)

layout(matrix(1:12,3,4))

pa.hc=hclust(de) # complete linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,de)[,3])
}
plot(2:10,as[2:10],type="b", main="Euclidea-complete")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[1,k-1]=min(l)
  max.cl[1,k-1]=max(l)
}

pa.hc=hclust(de,method="single") # single linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,de)[,3])
}
plot(2:10,as[2:10],type="b", main="Euclidea-single")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[2,k-1]=min(l)
  max.cl[2,k-1]=max(l)
}

pa.hc=hclust(de,method="average") # average linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,de)[,3])
}
plot(2:10,as[2:10],type="b", main="Euclidea-average")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[3,k-1]=min(l)
  max.cl[3,k-1]=max(l)
}


pa.hc=hclust(deq) # complete linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,deq)[,3])
}
plot(2:10,as[2:10],type="b", main="Euclidea^2-complete")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[4,k-1]=min(l)
  max.cl[4,k-1]=max(l)
}

pa.hc=hclust(deq,method="single") # single linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,deq)[,3])
}
plot(2:10,as[2:10],type="b", main="Euclidea^2-single")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[5,k-1]=min(l)
  max.cl[5,k-1]=max(l)
}

pa.hc=hclust(deq,method="average") # average linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,deq)[,3])
}
plot(2:10,as[2:10],type="b", main="Euclidea^2-average")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[6,k-1]=min(l)
  max.cl[6,k-1]=max(l)
}


pa.hc=hclust(dmax) # complete linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dmax)[,3])
}
plot(2:10,as[2:10],type="b", main="Massimo-complete")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[7,k-1]=min(l)
  max.cl[7,k-1]=max(l)
}

pa.hc=hclust(dmax,method="single") # single linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dmax)[,3])
}
plot(2:10,as[2:10],type="b", main="Massimo-single")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[8,k-1]=min(l)
  max.cl[8,k-1]=max(l)
}

pa.hc=hclust(dmax,method="average") # average linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dmax)[,3])
}
plot(2:10,as[2:10],type="b", main="Massimo-average")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[9,k-1]=min(l)
  max.cl[9,k-1]=max(l)
}

pa.hc=hclust(dman) # complete linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dman)[,3])
}
plot(2:10,as[2:10],type="b", main="Manhattan-complete")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[10,k-1]=min(l)
  max.cl[10,k-1]=max(l)
}

pa.hc=hclust(dman,method="single") # single linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dman)[,3])
}
plot(2:10,as[2:10],type="b", main="Manhattan-single")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[11,k-1]=min(l)
  max.cl[11,k-1]=max(l)
}

pa.hc=hclust(dman,method="average") # average linkage
as=rep(0,9)
for(i in 2:10){
  pa.cut=cutree(pa.hc,i)
  as[i]=mean(silhouette(pa.cut,dman)[,3])
}
plot(2:10,as[2:10],type="b", main="Manhattan-average")
# quantità minima/massima in un cluster
for (k in 2:5){
  pa.cut=cutree(pa.hc,k)
  l=rep(0,k)
  for (i in 1:k){
    l[i]=length(pa.cut[which(pa.cut==i)])
  }
  min.cl[12,k-1]=min(l)
  max.cl[12,k-1]=max(l)
}

layout(1)

min.cl
hist(min.cl[-c(4,5,6),2], breaks=15, main="Minimun of cluster", col="dodgerblue3",xlab="minimum number of one cluster")


#3 cluster

#plotto silhouette kmeans e i due pam
pa.km=kmeans(pa,3,nstart=20)
plot(silhouette(kmeans(pa,3,nstart=20)$cluster,dist(pa)),col=2:4, border=NA,cex=0.4)
plot(pam(pa,3))
plot(pam(pa,3,metric="manhattan"))

#plotto cluster sul piano principale kmeans e i due pam
plot(pa.pca$scores,col=1+pa.km$cluster,pch=10)
title("k-means")
#text(pa.pca$scores,labels=as.character(row.names(pa)),col=1+pa.km$cluster,pos=3) 
pa.pam=pam(pa,3)
plot(pa.pca$scores,pch=10,col=pa.pam$cluster)
title("pam")

#plotto silhouette metodi gerarchici buoni
pa.hc=hclust(de) # complete linkage
pa.cut=cutree(pa.hc,3)
plot(silhouette(pa.cut,de),col=heat.colors(3), border=par("fg"))

pa.hc=hclust(dman) # complete linkage
pa.cut=cutree(pa.hc,3)
plot(silhouette(pa.cut,dman),col=heat.colors(3), border=par("fg"))

#plotto cluster sul piano principale metodi gerarchici buoni
pa.hc.de=hclust(de) # complete linkage
pa.hc.dman=hclust(dman) # complete linkage
pa.cut.de=cutree(pa.hc.de,3)
pa.cut.dman=cutree(pa.hc.dman,3)
plot(pa.pca$scores,col=1+pa.cut.de, pch=10, main="distanza euclidea")
plot(pa.pca$scores,col=1+pa.cut.dman,pch=10, main= "distanza manhattan")

#parcoord kmeans
parcoord(pa.pca$scores[,1:3],col=pa.km$cluster)


#4o5cluster

#plotto silhouette e cluster sul piano principale kmeans
layout(matrix(1:4,2,2))
for (k in 4:5){
pa.km=kmeans(pa,k,nstart=20)
plot(silhouette(kmeans(pa,k,nstart=20)$cluster,dist(pa)),col=2:(k+1), border=NA,cex=0.4)
plot(pa.pca$scores,col=1+pa.km$cluster,pch=10)
}
#plotto silhouette e cluster sul piano principale pam euclidea
for (k in 4:5){
plot(pam(pa,k))
}
#plotto silhouette e cluster sul piano principale pam manhattan
for (k in 4:5){
plot(pam(pa,k,metric="manhattan"))
}

#plotto silhouette e cluster sul piano principale complete linkage distanza euclidea
pa.hc=hclust(de) # complete linkage de
for (k in 4:5){
pa.cut=cutree(pa.hc,k)
plot(silhouette(pa.cut,de),col=heat.colors(k), border=par("fg"))
plot(pa.pca$scores,col=1+pa.cut.de, pch=10, main="distanza euclidea")
}

#plotto silhouette e cluster sul piano principale complete linkage distanza manhattan
pa.hc=hclust(dman) # complete linkage dman
for (k in 4:5){
pa.cut=cutree(pa.hc,k)
plot(silhouette(pa.cut,dman),col=heat.colors(k), border=par("fg"))
plot(pa.pca$scores,col=1+pa.cut,pch=10, main= "distanza manhattan")
}

layout(1)

#grafici definitivi e inseriti nella relazione
layout(matrix(1:6,3,2))
k=4
pa.km=kmeans(pa,k,nstart=20)
plot(silhouette(pa.km$cluster,dist(pa)),col=2:5, border=NA,cex=0.4)
plot(pa.pca$scores,col=1+pa.km$cluster,pch=10)
parcoord(pa.pca$scores[,1:3],col=1+pa.km$cluster)
k=5
pa.km=kmeans(pa,k,nstart=20)
plot(silhouette(pa.km$cluster,dist(pa)),col=6:10, border=NA,cex=0.4)
plot(pa.pca$scores,col=5+pa.km$cluster,pch=10)
parcoord(pa.pca$scores[,1:3],col=5+pa.km$cluster)

layout(1)

