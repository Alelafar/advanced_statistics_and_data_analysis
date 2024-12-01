
#preparazione tabella

pa=read.csv2("tabella.csv",row.names=1)
head(pa)
library(corrplot)
pa[,1:16]=pa[,1:16]/(pa[,17]/100)
pa=pa[,-c(17)]
corrplot(round(cor(pa),2), "number", tl.cex= 0.6, number.cex=0.55)

#analisi pca

pa.pca=princomp(scale(pa))
summary(pa.pca)
plot(cumsum(pa.pca$sdev^2)/sum(pa.pca$sdev^2),type="b",ylim=c(0,1))
segments(1,0.8,7,0.8,col="red")

#piani principali e loadings

biplot(pa.pca,col=c("gray49","red"), cex=.8)
pa.ld=loadings(pa.pca)
corrplot(pa.ld,"number",tl.cex= 0.7, number.cex=0.87)
corrplot(varimax(pa.ld[,1:5])$loadings[c(2,3,7,12,13,14,15),],"number", number.cex=1.5)
layout(t(1:2))
biplot(pa.pca,col=c("gray49","red"), choices = c(2,3), cex=.8)
biplot(pa.pca,col=c("gray49","red"), choices = c(1,5), cex=.8)
layout(1)
biplot(pa.pca,col=c("gray49","red"), choices = c(3,4), cex=.8)


#analisi di stabilità

err=rep(0,35)
for(j in 1:35){
  pa_s=data.frame(scale(pa))
  pa_sm=pa_s[-j,]
  pa_sm.pca=princomp(pa_sm)
  err[j]=sqrt(mean(predict(pa_sm.pca,newdata=pa_s[j,])[c(1,2,3,4,5)]-predict(pa.pca)[j,c(1,2,3,4,5)])^2)
}
mean(err)
err