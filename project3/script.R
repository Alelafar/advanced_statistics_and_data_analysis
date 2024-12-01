

# Preparazione tabella
france.data=read.csv("tabella.csv",header=F,sep=";",dec=",")
head(france.data)
tp.data=france.data[,2]
head(tp.data)
tp=ts(tp.data,frequency=12,start=c(1994))
start(tp)
end(tp)
frequency(tp)

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
#Analisi autocorrelazione
plot(tp, main="Time Series")
acf(tp,30, main="Series ACF")
# vedere le differenze potrebbe aiutare a cercare la stagionalità
acf(diff(tp),30, main="Series ACF-DIFF")

# rappresentiamo i grafici dei vari periodi sovrapposti
m_tp=matrix(tp.data,12,27)
k=length(m_tp[1,])
k
par(bg="black")
layout(matrix(1:2,2,1))
ts.plot(m_tp,col=heat.colors(k))  # (non sovrapposti, ogni periodo ha la sua media)
ts.plot(scale(m_tp,scale=F),col=heat.colors(k))
par(bg="white")
layout(1)

layout(matrix(1:2,1,2))
#decompose additiva
tp.d=decompose(tp)
plot(tp.d)
#decomposizione moltiplicativa e confronto con la additiva
tp.dm=decompose(tp,type="multiplicative")
plot(tp.dm)
layout(1)
#decomposizione con stagionalità non uniforme
layout(matrix(1:4,2,2))
for (k in 1:4){
  plot(stl(tp,4+(2*k+1)))
}
layout(1)
k=9
plot(stl(tp,9),main="Stl - Decomposition")

#confronto tra trend, stagionalità, residui dei tre modelli
layout((matrix(1:3,3,1)))
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
plot(decompose(tp)$trend,col="blue",main="trend",ylab="trend",ylim=c(min(stl(tp,7)$time.series[,2])-10,max(stl(tp,7)$time.series[,2])+10))
lines(stl(tp,7)$time.series[,2],col="red")
legend("topright", legend=c("decompose", "stl"), 
       col=c("blue","red"), lty=1:9, inset=c(-0.09,0),
       bty = "n")
plot(decompose(tp)$seasonal,col="green3",main="seasonal",ylab="seasonal")
lines(mean(decompose(tp,type="multiplicative")$trend,na.rm=T)*(decompose(tp,type="multiplicative")$seasonal-1),col="blue")
lines(stl(tp,9)$time.series[,1],col="red")
legend("topright", legend=c("decompose ad", "decompose mult", "stl"), 
       col=c("green3","blue","red"), lty=1:9, inset=c(-0.12,0),
       bty = "n")
plot(decompose(tp)$random,col="green3",ylab="resid",main="residuals")
lines(mean(decompose(tp,type="multiplicative")$trend,na.rm=TRUE)*(decompose(tp,type="multiplicative")$random-1),col="blue")
lines(stl(tp,9)$time.series[,3],col="red")
legend("topright", legend=c("decompose ad", "decompose mult", "stl"), 
       col=c("green3","blue","red"), lty=1:9, inset=c(-0.12,0),
       bty = "n")
layout(1)


layout(matrix(1:2,2,1))
#struttura temporale sui residui
#additivo
tp.d=decompose(tp)
tp.dr=as.vector(window(tp.d$random,c(1994,7),c(2020,6)))
#plot(tp.dr,pch=20)
#var(tp.dr)/var(window(tp.data,c(1949,7),c(1960,6)))
acf(as.ts(tp.dr),main="Decompose ad")
#stl
tp.stl=stl(tp,9)
tp.stlr=as.vector(window(tp.stl$time.series[,3],c(1994,7),c(2020,6)))
acf(tp.stlr,main="Stl")

library(gridExtra)
library(grid)
tb=data.frame(matrix(c(round(sd(acf(tp.dr,plot=F)$acf),3),round(sd(acf(tp.stlr,plot=F)$acf),3)),1,2))
colnames(tb)<-c("sd acf decompose resid","sd acf stl resid")
tb
grid.newpage()
grid.table(tb)

layout(t(1:2))
#vera analisi dei residui (additivo + stl)
#additivo
hist(tp.dr,20,freq=F, main="Densità Decompose")
lines(density(tp.dr),col="blue")
lines(sort(tp.dr),dnorm(sort(tp.dr),mean(tp.dr),sd(tp.dr)),col="red")
legend("topright", legend=c("Empirica","Teorica"), col=c("blue","red"), lty=1:2)
qqnorm(tp.dr,main="Q-Q Plot Decompose")
qqline(tp.dr)
shapiro.test(tp.dr)
#stl
hist(tp.stlr,20,freq=F, main="Densità Stl")
lines(density(tp.stlr),col="blue")
lines(sort(tp.stlr),dnorm(sort(tp.stlr),mean(tp.stlr),sd(tp.stlr)),col="red")
legend("topright", legend=c("Empirica","Teorica"), col=c("blue","red"), lty=1:2)
qqnorm(tp.stlr,main="Q-Q Plot Stl")
qqline(tp.stlr)
shapiro.test(tp.stlr)
#indice numerico che può confermarci quale sia meglio (residui con sd più bassa)
sd(acf(tp.dr,plot=F)$acf)
sd(acf(tp.stlr,plot=F)$acf)

# Metodi di Holt-Winters

tp.hw=HoltWinters(tp)
plot(tp.hw)
c(tp.hw$alpha,tp.hw$beta,tp.hw$gamma)
plot(tp.hw$fitted)

layout(t(1:2))
tp.stl=stl(tp,9)
# confronto grafico (delle intercette) con il trend di stl
ts.plot(tp.stl$time.series[,2],tp.hw$fitted[,2],col=c("black","red"))
# confronto grafico con la stagionalità di stl
ts.plot(tp.stl$time.series[,1],tp.hw$fitted[,4],col=c("black","red"))
layout(1)
#mini regressione lineare per calcolare i coefficienti iniziali
k=round(length(tp.data)/7)
k
x=1:k
coefficients(lm(tp.data[1:k]~x))
plot(HoltWinters(tp,l.start=coefficients(lm(tp.data[1:k]~x))[1],
                 b.start=coefficients(lm(tp.data[1:k]~x))[2]))

#variare i parametri alpha.beta,gamma per avere dei modelli diversi che possono essere migliore da poi confrontare
layout(matrix(1:3,3,1))
for (a in 1:3){
  for (b in 1:3){
    for (c in 1:3){
      plot(HoltWinters(tp,alpha=(3*a-0.5)/10,beta=(3*b-0.5)/10,gamma=(3*c-0.5)/10), ylim=c(min(tp.data)-1,max(tp.data)+1),
           xlab=paste("alpha=",(3*a-0.5)/10," - beta=", (3*b-0.5)/10, " - gamma=", (3*c-0.5)/10))
    }
  }
}

layout(1)
#metodi a confronto
a1=0.341     #alpha soft
b1=0.013   #beta soft
c1=0.466     #gamma soft
a2=0.5     #alpha migliore
b=0.1
c2=0.55     #gamma migliore

#1
plot(tp.hw,xlab=paste("alpha=",round(a1,2)," - beta=", round(b1,2), " - gamma=", round(c1,2)))
#2
plot(tp.hw1,xlab=paste("alpha=",a2," - beta=", b, " - gamma=", c2))

layout(matrix(1:2,2,1))
# predizione (di un intero periodo)
#1
plot(tp.hw,predict(tp.hw,12),xlab=paste("alpha=",round(a1,2)," - beta=", round(b1,2), " - gamma=", round(c1,2)), xlim=c(2017,2022))
#2
plot(tp.hw1,predict(tp.hw1,12),xlab=paste("alpha=",a2," - beta=", b, " - gamma=", c2), xlim=c(2017,2022))

#modello moltiplicativo
tp.hw.m=HoltWinters(tp,l.start=coefficients(lm(tp.data[1:k]~x))[1],
                    b.start=coefficients(lm(tp.data[1:k]~x))[2],seasonal="multiplicative")
tp.hw.m1=HoltWinters(tp,alpha=0.5,beta=0.1,gamma=0.55,l.start=coefficients(lm(tp.data[1:k]~x))[1],
                     b.start=coefficients(lm(tp.data[1:k]~x))[2],seasonal="multiplicative")
plot(tp.hw.m,predict(tp.hw.m,12),xlab=paste("alpha=",round(a1,2)," - beta=", round(b1,2), " - gamma=", round(c1,2)))
plot(tp.hw.m1,predict(tp.hw.m1,12),xlab=paste("alpha=",a2," - beta=", b, " - gamma=", c2))
#vediamolo solo sugli ultimi anni
plot(tp.hw.m,predict(tp.hw.m,12),xlab=paste("alpha=",round(a1,2)," - beta=", round(b1,2), " - gamma=", round(c1,2)), xlim=c(2017,2022))
plot(tp.hw.m1,predict(tp.hw.m1,12),xlab=paste("alpha=",a2," - beta=", b, " - gamma=", c2), xlim=c(2017,2022))

#metodi di autoregressione

acf(tp)
pacf(tp)
L = length(tp.data)
l = 13 # numero di lag in ingresso
mtp = matrix(nrow = L - l, ncol = l + 1)
for (i in 1:(l + 1)) {
  mtp[, i] = tp.data[i:(L - l - 1 + i)]
}
mtp <- data.frame(mtp)
head(mtp)
tp.lm <- lm(X14 ~ ., data = mtp)  # X14 perché 13 lag in ingresso
summary(tp.lm) #possiamo decidere se escludere dei fattori con la riduzione del modello

#previsione modello completo (sul futuro: 1 periodo)
anni = 1
pt= rep(0, L + 12 * anni)
pt[1:L] = tp.data 
for (i in 1:(12 * anni)) {
  pt[L + i] = coef(tp.lm) %*% c(1,rev(pt[L + i - 1:l]))
}
tp.lm.pt = ts(pt, frequency = 12, start = c(1994))
length(window(tp,c(1995,2)))-length(resid(tp.lm))
tp.lm.a = window(tp,c(1995,2)) - resid(tp.lm)


layout(matrix(1:3,3,1))
ts.plot(tp, tp.lm.a, window(tp.lm.pt, start=c(2020, 12)), col = c("black","blue", "red"))
#Holtwinters
tp.hw=HoltWinters(tp)
c(tp.hw$alpha,tp.hw$beta,tp.hw$gamma)
plot(tp.hw,predict(tp.hw,12))
tp.lm.pt = window(tp.lm.pt, c(2020, 12))
tp.hw.pt = predict(tp.hw, 24)
ts.plot(tp, tp.lm.pt, tp.hw.pt, col = c("black", "red", "blue"))
layout(1)

ts.plot(window(tp,c(2017,1)), tp.lm.pt, tp.hw.pt, col = c("black", "red", "blue"))

#riduzione del modello
r=matrix(ncol=2,nrow=11)
tp.lm1= lm(X14 ~ ., data = mtp) 
summary(tp.lm1)
r[1,]=c(summary(tp.lm1)$r.squared,summary(tp.lm1)$adj.r.squared)
tp.lm2= lm(X14 ~ .-X9, data = mtp) 
summary(tp.lm2)
r[2,]=c(summary(tp.lm2)$r.squared,summary(tp.lm2)$adj.r.squared)
tp.lm3= lm(X14 ~ .-X9-X7, data = mtp) 
summary(tp.lm3)
r[3,]=c(summary(tp.lm3)$r.squared,summary(tp.lm3)$adj.r.squared)
tp.lm4= lm(X14 ~ .-X9-X7-X11, data = mtp) 
summary(tp.lm4)
r[4,]=c(summary(tp.lm4)$r.squared,summary(tp.lm4)$adj.r.squared)
tp.lm5= lm(X14 ~ .-X9-X7-X11-X8, data = mtp) 
summary(tp.lm5)
r[5,]=c(summary(tp.lm5)$r.squared,summary(tp.lm5)$adj.r.squared)
tp.lm6= lm(X14 ~ .-X9-X7-X11-X8-X4, data = mtp) 
summary(tp.lm6)
r[6,]=c(summary(tp.lm6)$r.squared,summary(tp.lm6)$adj.r.squared)
tp.lm7= lm(X14 ~ .-X9-X7-X11-X8-X4-X5, data = mtp) 
summary(tp.lm7)
r[7,]=c(summary(tp.lm7)$r.squared,summary(tp.lm7)$adj.r.squared)
tp.lm8= lm(X14 ~ .-X9-X7-X11-X8-X4-X5-X3, data = mtp) 
summary(tp.lm8)
r[8,]=c(summary(tp.lm8)$r.squared,summary(tp.lm8)$adj.r.squared)
tp.lm9= lm(X14 ~ .-X9-X7-X11-X8-X4-X5-X3-X6, data = mtp) 
summary(tp.lm9)
r[9,]=c(summary(tp.lm9)$r.squared,summary(tp.lm9)$adj.r.squared)
tp.lm10= lm(X14 ~ .-X9-X7-X11-X8-X4-X5-X3-X6-X10, data = mtp) 
summary(tp.lm10)
r[10,]=c(summary(tp.lm10)$r.squared,summary(tp.lm10)$adj.r.squared)
tp.lm11= lm(X14 ~ .-X9-X7-X11-X8-X4-X5-X3-X6-X10-X12, data = mtp) 
summary(tp.lm11)
r[11,]=c(summary(tp.lm11)$r.squared,summary(tp.lm11)$adj.r.squared)
ymin=min(r)
ymax=max(r)
plot(r[,1],pch=19,type="b",col="red",ylim=c(0.845,ymax))
lines(r[,2],pch=19,type="b",col="blue")


#Metodo minimi quadrati (per serie stazionarie, magari con trend pronunciati)

tp.ls = ar(tp, method = "ols")
tp.ls$order
ts.plot(tp, tp - tp.ls$resid, col = c("black", "red"))

#confronto tra modelli

#modelli
tp.hw=HoltWinters(tp,l.start=coefficients(lm(tp.data[1:k]~x))[1],
                  b.start=coefficients(lm(tp.data[1:k]~x))[2])
tp.hw1=HoltWinters(tp,alpha=0.5,beta=0.1,gamma=0.55,l.start=coefficients(lm(tp.data[1:k]~x))[1],
                   b.start=coefficients(lm(tp.data[1:k]~x))[2])
tp.lm <- lm(X14 ~ ., data = mtp)
tp.lm1=lm(X14 ~ .-X9-X7-X11-X8-X4-X5-X3-X6-X10-X12, data = mtp)
tp.ls = ar(tp, method = "ols")

#analisi dei residui
layout(matrix(1:4, 2, 2))
# estrazione dei residui
tp.lm.r=resid(tp.lm)
tp.lm1.r=resid(tp.lm1)
tp.ls.r = as.double(na.omit(tp.ls$resid))
tp.ls.fitted = as.double(na.omit(tp - tp.ls$resid))
tp.hw1.r=resid(tp.hw1)
# proporzione di varianza non spiegata
varianze= c(round(var(tp.hw1.r)/var(window(tp,1995)),2),
            round(var(tp.lm.r)/var(window(tp,c(1995,2))),2),
            round(var(tp.lm1.r)/var(window(tp,c(1995,2))),2),
            round(var(tp.ls.r)/var(tp.ls.r + tp.ls.fitted),2))
# autocorrelazione
acf(tp.hw1.r,main="HW-best")
acf(tp.lm.r,main="lm-compl")
acf(tp.lm1.r,main="lm-rid")
acf(tp.ls.r,main="min-quad")
sd.acf.res= c(round(sd(acf(tp.hw1.r,plot=F)$acf),3),round(sd(acf(tp.lm.r,plot=F)$acf),3),
              round(sd(acf(tp.lm1.r,plot=F)$acf),3),round(sd(acf(tp.ls.r,plot=F)$acf),3))
# densità empiriche
hist(tp.hw1.r, 20, freq = F,main="HW-best")
lines(density(tp.hw1.r),col="blue")
lines(sort(tp.hw1.r), dnorm(sort(tp.hw1.r), mean(tp.hw1.r), sd(tp.hw1.r)), col = "red")
legend("topright", legend=c("Empirica","Teorica"), col=c("blue","red"), lty=1:2)
hist(tp.lm.r, 20, freq = F,main="lm-compl")
lines(density(tp.lm.r),col="blue")
lines(sort(tp.lm.r), dnorm(sort(tp.lm.r), mean(tp.lm.r), sd(tp.lm.r)), col = "red")
legend("topright", legend=c("Empirica","Teorica"), col=c("blue","red"), lty=1:2)
hist(tp.lm1.r, 20, freq = F,main="lm-rid")
lines(density(tp.lm1.r),col="blue")
lines(sort(tp.lm1.r), dnorm(sort(tp.lm1.r), mean(tp.lm1.r), sd(tp.lm1.r)), col = "red")
legend("topright", legend=c("Empirica","Teorica"), col=c("blue","red"), lty=1:2)
hist(tp.ls.r, 20, freq = F,main="min-quad")
lines(density(tp.ls.r),col="blue")
lines(sort(tp.ls.r), dnorm(sort(tp.ls.r), mean(tp.ls.r), sd(tp.ls.r)), col = "red")
legend("topright", legend=c("Empirica","Teorica"), col=c("blue","red"), lty=1:2)
# grafico quantile-quantile
qqnorm(tp.hw1.r, pch = 20, main="HW-best")
qqline(tp.hw1.r)
qqnorm(tp.lm.r, pch = 20, main="lm-compl")
qqline(tp.lm.r)
qqnorm(tp.lm1.r, pch = 20, main="lm-rid")
qqline(tp.lm1.r)
qqnorm(tp.ls.r, pch = 20, main="min-quad")
qqline(tp.ls.r)
# test
shapiro.test(tp.hw1.r)
shapiro.test(tp.lm.r)
shapiro.test(tp.lm1.r)
shapiro.test(tp.ls.r)
shapiro= c(0.52,0.04,0.04,0.01)
layout(1)
library(gridExtra)
library(grid)
tb=data.frame(matrix(c(varianze,sd.acf.res,shapiro),4,3))
row.names(tb)<-c("HW-best", "lm-compl", "lm-rid", "min-quad")
colnames(tb)<-c("Varianza non spiegata", "sd acf.residuals","p-value Shapiro.test")
tb
grid.newpage()
grid.table(tb)

#autovalidazione1
train = window(tp, end = c(2018, 12))
test = window(tp, 2019)
#Hw e hw1
tp.hw=HoltWinters(train,l.start=coefficients(lm(tp.data[1:k]~x))[1],
                  b.start=coefficients(lm(tp.data[1:k]~x))[2])
tp.hw1=HoltWinters(train,alpha=0.5,beta=0.1,gamma=0.55,l.start=coefficients(lm(tp.data[1:k]~x))[1],
                   b.start=coefficients(lm(tp.data[1:k]~x))[2])
tp.hw.p=predict(tp.hw,12)
tp.hw1.p=predict(tp.hw1,12)
#lm totale
L=length(train)
pt= rep(0, L + 12)
pt[1:L] = tp.data[1:L] 
for (i in 1:(12)) {
  pt[L + i] = coef(tp.lm) %*% c(1,rev(pt[L + i - 1:l]))
}
tp.lm.pt = ts(pt, frequency = 12, start = c(1994))
#lm ristretto
ptr= rep(0, L + 12)
ptr[1:L] = tp.data[1:L]
for (i in 1:(12)) {
  ptr[L + i] = coef(tp.lm1) %*% c(1,  ptr[L + i - 13], ptr[L + i - 12], ptr[L + i - 1])
}
tpr.lm.pt = ts(ptr, frequency = 12, start = c(1994,1))
#minimiquadrati
tp.ls = ar(train, method = "ols")
tp.ls.pt = predict(tp.ls, n.ahead = 12, se.fit = FALSE)

layout(1)
# disegno previsione (2019)
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
#Hw1
ts.plot(window(tp,2019),tp.hw1.p, col=c("black","red"),lwd=c(2,2))
#lm totale
lines(window(tp.lm.pt,2019),col="blue",lwd=2)
#lm ristretto
lines(window(tpr.lm.pt,2019), col= "green3",lwd=2)
#minimiquadrati

err.19=c(round(sqrt(mean((tp.hw1.p - test)^2)),2),
         round(sqrt(mean((window(tp.lm.pt,2019) - test)^2)),2),
         round(sqrt(mean((window(tpr.lm.pt,2019) - test)^2)),2),
         round(sqrt(mean((tp.ls.pt - test)^2)),2))

#autovalidazione2
train = window(tp, end = c(2019, 12))
test = window(tp, 2020)
#Hw e hw1
tp.hw=HoltWinters(train,l.start=coefficients(lm(tp.data[1:k]~x))[1],
                  b.start=coefficients(lm(tp.data[1:k]~x))[2])
tp.hw1=HoltWinters(train,alpha=0.5,beta=0.1,gamma=0.55,l.start=coefficients(lm(tp.data[1:k]~x))[1],
                   b.start=coefficients(lm(tp.data[1:k]~x))[2])
tp.hw.p=predict(tp.hw,12)
tp.hw1.p=predict(tp.hw1,12)
#lm totale
L=length(train)
pt= rep(0, L + 12)
pt[1:L] = tp.data[1:L] 
for (i in 1:(12)) {
  pt[L + i] = coef(tp.lm) %*% c(1,rev(pt[L + i - 1:l]))
}
tp.lm.pt = ts(pt, frequency = 12, start = c(1994))
#lm ristretto
ptr= rep(0, L + 12)
ptr[1:L] = tp.data[1:L]
for (i in 1:(12)) {
  ptr[L + i] = coef(tp.lm1) %*% c(1,  ptr[L + i - 13], ptr[L + i - 12], ptr[L + i - 1])
}
tpr.lm.pt = ts(ptr, frequency = 12, start = c(1994,1))
#minimiquadrati
tp.ls = ar(train, method = "ols")
tp.ls.pt = predict(tp.ls, n.ahead = 12, se.fit = FALSE)

layout(1)
# disegno previsione (2020)
#HW1
lines(tp.hw1.p,col="red",lwd=2)
#lm totale
lines(window(tp.lm.pt,2020),col="blue",lwd=2)
#lm ristretto
lines(window(tpr.lm.pt,2020), col= "green3",lwd=2)
#minimiquadrati
lines(tp.ls.pt, col = "dodgerblue3",lwd=2)
legend("topright", legend=c("HW-best", "lm-compl", "lm-rid", "min-quad"), 
       col=c("red","blue", "green3", "dodgerblue3"), lty=1:5, inset=c(-0.11,0),
       bty = "n")

#sqm errori
library(gridExtra)
library(grid)
err.20=c(round(sqrt(mean((tp.hw1.p - test)^2)),2),
         round(sqrt(mean((window(tp.lm.pt,2020) - test)^2)),2),
         round(sqrt(mean((window(tpr.lm.pt,2020) - test)^2)),2),
         round(sqrt(mean((tp.ls.pt - test)^2)),2))
tb=data.frame(matrix(c(err.19,err.20),4,2))
row.names(tb)<-c("HW-best", "lm-compl", "lm-rid", "min-quad")
colnames(tb)<-c("2019","2020")
tb
grid.newpage()
grid.table(tb)

#2018
train = window(tp, end = c(2017, 12))
test = window(tp, 2018)
#hw1
tp.hw1=HoltWinters(train,alpha=0.5,beta=0.1,gamma=0.55,l.start=coefficients(lm(tp.data[1:k]~x))[1],
                   b.start=coefficients(lm(tp.data[1:k]~x))[2])
tp.hw1.p=predict(tp.hw1,12)
ts.plot(window(tp,start=c(2018,1),end=c(2018,12)),tp.hw1.p, col=c("black","red"),lwd=c(2,2),main="Predizione 2018")
sqrt(mean(tp.hw1.p - test)^2)

#autovalutazione più robusta 
l=length(tp)
k=24
res.hw1=rep(0,k)
res.lm=rep(0,k)
j=1
for(i in (l-k):(l-1)){
  tp_cv=ts(tp[1:i],frequency=12,start=c(1994,1))
  tp.hw1=HoltWinters(tp_cv,alpha=0.5,beta=0.1,gamma=0.55)
  L=length(tp_cv)
  mtp = matrix(nrow = L - 13, ncol = 13 + 1)
  for (h in 1:(13 + 1)) {
    mtp[, h] = tp_cv[h:(L - 13 - 1 + h)]
  }
  mtp <- data.frame(mtp)
  tp.lm <- lm(X14 ~ ., data = mtp)
  tp.hw1.p=predict(tp.hw1,1)
  tp.lm.p = coef(tp.lm) %*% c(1,rev(tp.data[L + 1 - 1:13]))
  res.hw1[j]=tp.hw1.p - tp[i+1]
  res.lm[j]=tp.lm.p - tp[i+1]
  j=j+1
}
sqrt(mean(res.hw1[1:12]^2))
sqrt(mean(res.hw1[13:24]^2))
sqrt(mean(res.lm[1:12]^2))
sqrt(mean(res.lm[12:24]^2))
plot(res.hw1,type="b",pch=20,col="blue")
lines(res.lm,type="b",pch=20,col="green3")

layout(1)
#2021

data=c(58395,54000,52300,61500,60400,61200,60900,66000,65800,65200,67000)
tp2021=ts(data,frequency=12,start=c(2020,12))
start(tp2021)
end(tp2021)

tp.hw1=HoltWinters(tp,alpha=0.5,beta=0.1,gamma=0.55)
tp.hw1.r=resid(tp.hw1)  
ts.plot(window(tp,2020), predict(tp.hw1,12), col=c("black","red"), lwd=c(2,2), main="Previsione 2021")
lines(predict(tp.hw1,12)+qnorm(0.05,mean(tp.hw1.r),sd(tp.hw1.r)),col="blue")
lines(predict(tp.hw1,12)+qnorm(0.95,mean(tp.hw1.r),sd(tp.hw1.r)),col="blue")
lines(predict(tp.hw1,12)+quantile(tp.hw1.r,0.05),col="green3")
lines(predict(tp.hw1,12)+quantile(tp.hw1.r,0.95),col="green3")
legend("topright", legend=c("non parametriche", "parametriche"),
       col=c("blue", "green3"), lty=1:2, cex=1.1)

ts.plot(window(tp,2020),tp2021, col=c("black","black"), lwd=c(2,2), main="Previsione 2021")
lines(predict(tp.hw1,10),col="red",lwd=2)
train = window(tp, end = c(2019, 12))
#Hw e hw1
tpr.hw1=HoltWinters(train,alpha=0.5,beta=0.1,gamma=0.55)
tpr.hw1.p=predict(tpr.hw1,24)
lines(window(tpr.hw1.p,c(2021,1),c(2021,10)),col="blue",lwd=2)
legend("bottomleft", legend=c("HW-best-2020", "HW-best-2019"),
       col=c("red", "blue"), lty=1:2, cex=1.3)














