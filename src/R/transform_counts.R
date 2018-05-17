library(statmod)
library(gdata)
library(genefilter)
library(EBImage)
library(rhdf5)
library(DESeq)
library(statmod)
library("hom.Hs.inp.db")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(preprocessCore)
source('~/Code/hipsci-singlecell/R/utils.R')


#1. REad in data
#iPS
dataAll_0 = h5dump('../data_files/iPS/data_raw.hdf5')
dataAll_2 = h5dump('../data_files/day2/data_raw.hdf5')
dataAll_3 = h5dump('../data_files/day3/data_raw.hdf5')


times = c('0','2','3')
#cell_phases = c('G2M')
nCountsHum=c()
nCountsERCC=c()
countsHumAll=c()
countsERCCAll=c()
labels = c()
labKey = c()

for(i in 1:length(times)){
  
eval(parse(text=paste("data = dataAll_",times[i],sep="")))
dataHum = load_data(data)

geneTypes <- factor( c( ENSG="ENSG", ERCC="ERCC" )[
  substr( rownames(dataHum), 1, 4 ) ] )

#2. calculate different normalisations for counts
countsHum <- dataHum[ which( geneTypes=="ENSG" ), ]
countsERCC <- dataHum[ which( geneTypes=="ERCC" ), ]

sfHum <- estimateSizeFactorsForMatrix( countsHum )
sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
#we do not want to normalize for differences in biological starting material... 
#as this can be due to different stages in cell cycle
#sfHum <- sfERCC #!ojo!aufgemerkt!
rbind( sfHum, sfERCC )

nCountsERCC_i <- t( t(countsERCC) / sfERCC )
nCountsHum_i <- t( t(countsHum) / sfHum )

nCountsERCC_i <- normalize.quantiles(countsERCC)
nCountsHum_i <- normalize.quantiles(countsHum)
                    
nCountsERCC <- cbind(nCountsERCC,nCountsERCC_i)
nCountsHum <- cbind(nCountsHum,nCountsHum_i)

countsERCCAll <- cbind(countsERCCAll,countsERCC)
countsHumAll <- cbind(countsHumAll,countsHum)


labels = c(labels, rep(i, dim(countsHum)[2]))
labKey = c(labKey, colnames(countsHum))
}
labKey = as.factor(labKey)



sfHum <- estimateSizeFactorsForMatrix( countsHumAll )
sfERCC <- estimateSizeFactorsForMatrix( countsERCCAll )

ratioERCC <- colSums(countsHumAll)/(colSums(countsERCCAll)+colSums(countsHumAll))

#sfERCC <- rowSums(t(countsERCCAll))
#we do not want to normalize for differences in biological starting material... 
#as this can be due to different stages in cell cycle

sfHum <- sfERCC #!ojo!aufgemerkt!
rbind( sfHum, sfERCC )
nCountsERCC <- t( t(countsERCCAll) / sfERCC )
nCountsHum <- t( t(countsHumAll) / sfHum )

rownames(nCountsHum)=rownames(countsHumAll)

barplot(c(mean(log10(nCountsHum[,labels==1]+1)),mean(log10(nCountsHum[,labels==2]+1)),mean(log10(nCountsHum[,labels==3]+1))))
#pdf('./fraction_ERCC.pdf', height=6, width=6)
histogram(~c(ratioERCC[labels==1],ratioERCC[labels==2],ratioERCC[labels==3],ratioERCC[labels==4])| labKey,xlab='Fraction of ERCCs of total mapped reads')
#dev.off()

x11();histogram(~c(apply(log10(nCountsERCC[,labels==1]+1),2,mean),apply(log10(nCountsERCC[,labels==2]+1),2,mean),apply(log10(nCountsERCC[,labels==3]+1),2,mean),apply(log10(nCountsERCC[,labels==4]+1),2,mean))| labels,xlab='log10(ERCC)')

x11();histogram(~c(apply(log10(nCountsHum[,labels==1]+1),2,mean),apply(log10(nCountsHum[,labels==2]+1),2,mean),apply(log10(nCountsHum[,labels==3]+1),2,mean),apply(log10(nCountsHum[,labels==4]+1),2,mean))| labels,xlab='log10(Hum)')



meansERCC <- rowMeans( nCountsERCC )
varsERCC <- rowVars( nCountsERCC )
cv2ERCC <- varsERCC / meansERCC^2



#normalized counts (brennecke)
meansHum <- rowMeans( nCountsHum )
varsHum <- rowVars( nCountsHum )
cv2Hum <- varsHum / meansHum^2

#absolute mRNA
totCERCC <- rowSums(t(nCountsERCC)) #use counts corrected for sequencing depth
totCHum <- rowSums(t(nCountsHum))
fac <- 1
totRHum <- fac*(totCHum/totCERCC)
sfabs=totCHum/totRHum #(=fac/totCERCC)
anCountsHum <- t( t(countsHumAll) / sfabs )
anCountsERCC <- t( t(countsERCCAll) / sfabs )

#log counts
LCountsHum <- log10(nCountsHum+1)
LCountsERCC <- log10(nCountsERCC+1)
 LmeansERCC <- rowMeans( LCountsERCC )
 LvarsERCC <- rowVars( LCountsERCC )
 Lcv2ERCC <- LvarsERCC / LmeansERCC^2
 LmeansHum <- rowMeans( LCountsHum )
 LvarsHum <- rowVars( LCountsHum )
 Lcv2Hum <- LvarsHum / LmeansHum^2

anmeansERCC <- rowMeans( anCountsERCC )
anvarsERCC <- rowVars( anCountsERCC )
ancv2ERCC <- anvarsERCC / anmeansERCC^2
anmeansHum <- rowMeans( anCountsHum )
anvarsHum <- rowVars( anCountsHum )
ancv2Hum <- anvarsHum / anmeansHum^2

#3. fit technical noise

#normalised counts (with size factor)
minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .8 ) )
useForFitA <- meansERCC >= minMeanForFitA
fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
                    cv2ERCC[useForFitA] )
#variance explained
residualA <- var( log( fitted.values(fitA) ) - log( cv2ERCC[useForFitA] ) )
totalA <- var( log( cv2ERCC[useForFitA] ) )
#plot fit
pdf('./technical_noiseESCKedar.pdf',width=5,height=5)
plot( meansERCC, cv2ERCC, log="xy", col=1+useForFitA)
xg <- 10^seq( -3, 5, length.out=100 )
lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg )
segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
          meansERCC[useForFitA], fitA$fitted.values, col="gray" )
title(paste('coeffs',as.character(round(fitA$coefficients,3)),sep=""))
dev.off()

minBiolDisp <- .5^2
xi <- mean( 1 / sfERCC )
m <- ncol(countsHum)
psia1thetaA <- mean( 1 / sfERCC ) +
  ( coefficients(fitA)["a1tilde"] - xi ) * mean( sfERCC / sfHum )
cv2thA <- coefficients(fitA)["a0"] + minBiolDisp + coefficients(fitA)["a0"] * minBiolDisp
testDenomA <- ( meansHum * psia1thetaA + meansHum^2 * cv2thA ) / ( 1 + cv2thA/m )
pA <- 1 - pchisq( varsHum * (m-1) / testDenomA, m-1 )
padjA <- p.adjust( pA, "BH" )
table( padjA < .1 )

#4. Transform to log-space and propagate error
eps=1
LogNcountsHum=log10(nCountsHum+eps)
dLogNcountsHum=1/(meansHum+eps)
var_techHum=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansHum)*meansHum^2
LogVar_techHum=(dLogNcountsHum*sqrt(var_techHum))^2 #error propagation 

#check in plot...
xg <- 10^seq( -3, 5, length.out=100 )
var_tech_fit=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg)*xg^2 

#pdf('./technical_noise.pdf',width=8,height=4)
x11()
plot( LmeansHum, Lcv2Hum, log="y", col=1+(padjA<0.1),xlab='meansLogHum',ylab='cv2LogHum',ylim=c(1e-3,1e2))  
#plot( LmeansERCC, Lcv2ERCC, log="y", col=1+(useForFitA))  
#if std(xg)=s_a then: std(log(xg))=s_a/(xg*log(10))
lines( log10(xg+1), ((sqrt(var_tech_fit)/(xg+1))^2)/log10(xg+1)^2 ,lwd=2,col='green' )
points(LmeansERCC, Lcv2ERCC,col='blue',pch=15,cex=1.1)
legend('topright',c('Hum (padj >= 0.1)','Hum (padj<0.1)','ERCC'),pch=c(1,1,15),col=c('black','red','blue'),cex=0.8)

  #plot in natural space  
x11()
#pdf('./Hum_mean_CV.pdf',width=4.5,height=4.5)
plot( meansHum, cv2Hum, log="xy", col=1+(padjA<0.1),ylim=c(0.1,95), xlab='Mean Counts', ylab='CV2 Counts')
#plot( meansHum, cv2Hum*meansHum*2, log="xy", col=1+(padjA<0.1), xlab='Mean Counts', ylab='CV2 Counts')
xg <- 10^seq( -3, 5, length.out=100 )
lines( xg, (coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg),lwd=2,col='blue' )
points(meansERCC, cv2ERCC,col='blue',pch=15,cex=1.1)
points(meansHum[cc_gene_indices], cv2Hum[cc_gene_indices],col=rgb(0,255,0,50,maxColorValue=255),pch=2,cex=0.75)
points(meansHum[ccCBall_gene_indices], cv2Hum[ccCBall_gene_indices],col=rgb(0,255,0,50,maxColorValue=255),pch=2,cex=0.8)
#legend('bottomleft',c('T-cells (padj >= 0.1)','T-cells (padj<0.1)','ERCC','Cell Cycle genes'),pch=c(1,1,15),col=c('black','red','blue','green'),cex=0.7)
dev.off()  

#alternatively: fit mean-CV relation of ERCC in log-space
#LogNcountsERCC=log10(nCountsERCC+eps)
LogNcountsList=list()
useForFitL=LmeansERCC>2
LogNcountsList$mean=LmeansERCC[useForFitL]
LogNcountsList$cv2=Lcv2ERCC[useForFitL]
fit_loglin=nls(cv2 ~ a* 10^(-k*mean), LogNcountsList,start=c(a=20,k=1))
LogVar_techHum_logfit <- coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*LmeansHum)
x11()
plot( LmeansERCC, Lcv2ERCC, log="y", col=1+useForFitL)

#pdf('./technical_noise_logfit.pdf',width=4,height=4)

plot( LmeansHum, Lcv2Hum, log="y", col=1+(padjA<0.1),ylim=c(1e-3,1e2),xlab='meansLogHum',ylab='cv2LogHum')
xg <- seq( 0, 5.5, length.out=100 )
lines( xg, coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*xg ),lwd=2,col='green' )
points(LmeansERCC, Lcv2ERCC,col='blue',pch=15,cex=1.1)
legend('topright',c('Hum (padj >= 0.1)','Hum (padj<0.1)','ERCC'),pch=c(1,1,15),col=c('black','red','blue'))
#dev.off()

#mean-variance
varERCC=LvarsERCC
x11()
#pdf('./mean_varTcellsNocc.pdf',width=6,height=6)
plot( LmeansHum, Lcv2Hum*LmeansHum^2, col=1+(padjA<0.1),ylim=c(1e-3,150.5),log="y",xlab='meansLogHum',ylab='VarLogHum')
xg <- seq( 0, 5.5, length.out=100 )
xg2 <- 10^seq( -3, 5.5, length.out=100 )
lines( xg, coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*xg )*xg^2,lwd=2,col='green' )
lines( log10(xg2+1),((((sqrt(var_tech_fit)/(xg2+1))^2))) ,lwd=2,col='green',lty=2)
datamv<-list()
datamv$mean=LmeansERCC[LmeansERCC>0.5]
datamv$var=varERCC[LmeansERCC>0.5]
fit_var=nls(var ~ mean^2*a* 10^(-k*mean), datamv,start=c(a=20,k=1))
fit_var2 = loess(varERCC[useForFitA] ~ LmeansERCC[useForFitA], span=0.8, control=loess.control(surface="direct"))
fit_varSmall = loess(varERCC ~ LmeansERCC, span=0.8)
Var_techHum_logfit <- xg^2*coefficients(fit_var)["a"] *10^(-coefficients(fit_var)["k"]*xg)
Var_techHum_logfit_loess <-  predict(fit_var2, xg)
Var_techHum_logfit_loess[xg>3.73] = 4e-3

points(LmeansERCC, Lcv2ERCC*LmeansERCC^2,col='blue',pch=15,cex=1.1)
lines(xg, Var_techHum_logfit_loess,lwd=3,col='blue',lty=3)
legend('topright',c('Hum (padj >= 0.1)','Hum (padj<0.1)','ERCC'),pch=c(1,1,15),col=c('black','red','blue'))
legend('bottomright',c('Fit in log-space','Fit in count space','Fit in count space'),lty=c(1,2,3),col=c('green','green','blue'))
dev.off()
xge=LmeansHum
Var_techHum_logfit <- xge^2*coefficients(fit_var)["a"] *10^(-coefficients(fit_var)["k"]*xge)
Var_techHum_logfit_loess <-  predict(fit_var2, xge)
Var_techHum_logfit_loess[xge>3.73] = 4e-3


#now save everything to hdf5 file.
#gene names 
gene_names=rownames(nCountsHum)
gene_names_het=gene_names[which(padjA<0.1)]



#all Cycle base genes homologs
dataCB=read.table(file='~/Downloads/genes_cycle_base2.txt', header=T)
hu2musAll=as.character(dataCB[1:600,3])
ccCBall_gene_indices=match(unlist(hu2musAll),rownames(nCountsHum))

  #get cell cycle genes from GO 
#GO:0045064#Th2 diff
#GO:0030217#T cell diff
#GO:0042110# T cell activ
#GO:0045630 pos reg Th2 cell 
#GO:0034976 ER stress
#GO:0006979 ox. Stress

  xxGO <- as.list(org.Hs.egGO2EG)
  cell_cycleEG <-unlist(xxGO['GO:0007049'])

#get ENSEMBLE ids
  x <- org.Hs.egENSEMBL
  mapped_genes <- mappedkeys(x)
  xxE <- as.list(x[mapped_genes])
  ens_ids_cc<-unlist(xxE[cell_cycleEG])


gene_ids = row.names(countsHum)
gene_ids_het=gene_ids[which(padjA<0.1)] 

#get indices of cell cycle genes (GO)
  cc_gene_indices <- na.omit(match(ens_ids_cc, rownames(dataHum)))
  #rename a few variables...
  cellcyclegenes <- ens_ids_cc
  cellcyclegenes_filter <- cc_gene_indices
  cell_names <- colnames(nCountsHum)
  Y <- nCountsHum
  genes_heterogen <- (padjA<0.1)*1
  countsERCC_mat=as.matrix(countsERCC * 1)
  countsERCCAll=as.matrix(countsERCCAll * 1)
  countsHum_mat = as.matrix(countsHum * 1)
ct = c('iPS', 'Day2', 'Day3')
cell_type = ct[labels]
geno_id = cell_names

  h5save(cell_type,geno_id,labels,nCountsERCC,countsERCCAll,countsHumAll,labels,ccCBall_gene_indices,gene_ids,gene_ids_het,cellcyclegenes_filter,cellcyclegenes,nCountsHum,genes_heterogen,LogVar_techHum,LogNcountsHum,Var_techHum_logfit,file='normCountsHum_all.h5f')
