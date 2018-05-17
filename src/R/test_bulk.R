library(statmod)
library(gdata)
library(genefilter)
library(EBImage)
library(rhdf5)
#library(DESeq)
library(statmod)
library("hom.Hs.inp.db")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(preprocessCore)
source('~/Code/hipsci-singlecell/R/utils.R')
library(DESeq2)
library(gtools)

####load bulk data
#1. REad in bulk data

times = c('0','3')

cell_types = c('ffdm_2','xavk', 'ffdj_1', 'coxy')

labels = c()
labT = c()
labC = c()
labKey = c()

indx = 1
countsHum = c()
for(i in 1:length(times)){
  for(j in cell_types){
    fn = paste('"../data_files/bulk/sc_',as.character(indx), '.bam.quant.gene.txt"', sep="")
    eval(parse(text=paste("data = read.table(",fn,', row.names = 1)',sep="")))
    rn  = row.names(data)
    rn = strsplit(rn, '.', fixed = T)
    rn = unlist(lapply(rn, function(x)unlist(x)[1]))
    row.names(data) <- rn
    
    geneTypes <- factor( c( ENSG="ENSG", ERCC="ERCC" )[
      substr( rownames(data), 1, 4 ) ] )
    
    labels = c(labels, paste(j,times[i], sep='.'))
    labT = c(labT, times[i])
    labC = c(labC, j)
    #2. calculate different normalisations for counts
    countsHum <- cbind(countsHum,data[ which( geneTypes=="ENSG" ), ])
    indx= indx+1
  }
}

sfMmus = estimateSizeFactorsForMatrix(countsHum)
RawCountsHum = countsHum
countsHum = RawCountsHum/sfMmus

labKey = as.factor(labels)
facT = as.factor(labT)
facC = as.factor(labC)
rownames(countsHum)<-rn[which( geneTypes=="ENSG" )]
rownames(RawCountsHum)<-rn[which( geneTypes=="ENSG" )]

gene_names = rownames(countsHum)


###load DE genes form SC analysis 
dataSC = h5dump('./normCountsHum_all2.h5f')
dataCorr= h5dump('../analysis/results/CcGenoPluriInt/varsHipsciCLVMfac_CcGenoPluriIntGPl.hdf5')
is_het = dataSC$genes_heterogen
#is_03 = dataSC$cell_type=='Day3'
is_03 = dataSC$cell_type=='iPS' | dataSC$cell_type=='Day3'
cell_type = dataSC$cell_type
ctF = as.factor(cell_type)
geno03 = dataSC$geno_id[is_03]
cell_type03 = dataSC$cell_type[is_03]
gtF = as.factor( dataSC$geno_id)
data_het = t(dataSC$LogNcountsHum[which(is_het==1),])

data_hetCorr = dataCorr$Ycorr
dataC_hetCorr = data_hetCorr[is_03,]

countsHumSC = dataSC$countsHumAll
rownames(countsHumSC) = dataSC$gene_ids

dataC_het = t(dataSC$nCountsHum[which(is_het==1),])
dataC_het = dataC_het[is_03,]




####now do test and check for concordance
#concordance between genotypes

gt1 = 'ffdj_1'
gt2 = 'ffdm_2'

ct12 = 'iPS'

do_rlog = 0
do_DESeq = 1
threshold = 0.1

pvGiPS = c()
for(i in 1:dim(data_het)[2]){
  res_ = wilcox.test(dataC_het[geno03==gt1 & cell_type03==ct12,i], dataC_het[geno03==gt2 & cell_type03==ct12,i])
  pvGiPS = c(pvGiPS, res_$p.value)
}


L2FC = foldchange2logratio(foldchange(apply(dataC_het[geno03==gt1 & cell_type03==ct12,],2,mean), apply(dataC_het[geno03==gt2 & cell_type03==ct12,],2,mean)))

pvadjGiPS = p.adjust(pvGiPS, 'fdr')
DEsceHetGiPS = dataSC$gene_ids_het[which(pvadjGiPS<threshold & is.finite(L2FC))]

genes_plot = DEsceHetGiPS

indsSC = match(genes_plot,dataSC$gene_ids_het)



if(ct12=='Day3'){
  idx_gt = which(cell_types %in% c(gt1,gt2))+4
}else{
  idx_gt = which(cell_types %in% c(gt1,gt2))
}

#DE bulk 
if(do_DESeq==1){
  colData_iPS <- DataFrame(cell_type = facC[idx_gt ])
  dds_iPS <- DESeqDataSetFromMatrix(RawCountsHum[,idx_gt], colData_iPS, formula(~ cell_type))
  
  dds_iPS <- DESeq(dds_iPS)
  res_iPS <- results(dds_iPS)
  
  idx_bulkhet = (rownames(res_iPS) %in% dataSC$gene_ids_het)
  indsBulk = match(genes_plot,rownames(res_iPS)[idx_bulkhet])
  if(do_rlog==1){
    L2FCbulk = foldchange2logratio(foldchange(rlogMat[idx_bulkhet,1], rlogMat[idx_bulkhet,2]))  
  }else{L2FCbulk = res_iPS[idx_bulkhet,2]}
  
}else{    
  idx_bulkhet = (rownames(countsHum) %in% dataSC$gene_ids_het)
  data_hetBulk = countsHum[idx_bulkhet,idx_gt]
  L2FCbulk = foldchange2logratio(foldchange(data_hetBulk[,1], data_hetBulk[,2]))  
}

#plot
indsBulk = match(genes_plot,rownames(countsHum)[idx_bulkhet])
#pdf(paste('./bulk/scatterFC',gt1,gt2,ct12,as.character(do_rlog),as.character(do_DESeq),'2.pdf',sep='_'),width=4.5,height=4.5)
plot(L2FC[indsSC],-L2FCbulk[indsBulk], type = "p",xlab='log2FC SC',ylab='log2FC Bulk', xlim=c(-8,8), ylim = c(-4,4))
#dev.off()

r = cor(x=L2FC[indsSC], y=L2FCbulk[indsBulk], method = "s", use = "complete.obs")
r


plot_fdr = c(0.25,0.1, 0.01,0.001, 1e-6)
fdr_all0 = pvadjGiPS[indsSC]
fdr_all = fdr_all0

for(pval in plot_fdr){
  fdr_all[fdr_all0<pval] = pval  
}

#pdf('./F2FCGTalpha.pdf',width = 5, height = 5)
dfr = data.frame(L2FC[indsSC],L2FCbulk[indsBulk])
d <- ggplot(dfr,aes(L2FC.indsSC.,L2FCbulk.indsBulk.))
d <- d+geom_point(aes(alpha = as.factor(1-fdr_all))) +scale_alpha_discrete(labels=as.character(plot_fdr))+labs(alpha="FDR")
d <- d+xlab("Log2 Fold Change Single Cells")
d <- d+ylab("Log2 Fold Change Bulk")
d + ggtitle('iPS vs Day3')


#concordance between days
ct1 = 'iPS'
ct2 = 'Day3'

do_rlog = 0
do_DESeq = 1
threshold = 0.25

pvGiPS = c()
for(i in 1:dim(data_het)[2]){
  res_ = wilcox.test(dataC_het[cell_type03==ct1,i], dataC_het[cell_type03==ct2,i])
  pvGiPS = c(pvGiPS, res_$p.value)
}


L2FC = foldchange2logratio(foldchange(apply(dataC_het[cell_type03==ct1,],2,mean), apply(dataC_het[cell_type03==ct2,],2,mean)))

pvadjGiPS = p.adjust(pvGiPS, 'fdr')
DEsceHetGiPS = dataSC$gene_ids_het[which(pvadjGiPS<threshold & is.finite(L2FC))]

genes_plot = DEsceHetGiPS

indsSC = match(genes_plot,dataSC$gene_ids_het)

#DE bulk 
if(do_DESeq==1){
  colData_iPS <- DataFrame(cell_type = facT)
  dds_iPS <- DESeqDataSetFromMatrix(RawCountsHum, colData_iPS, formula(~ cell_type))
  
  dds_iPS <- DESeq(dds_iPS)
  res_iPS <- results(dds_iPS)    
  
  idx_bulkhet = (rownames(res_iPS) %in% dataSC$gene_ids_het)
  indsBulk = match(genes_plot,rownames(res_iPS)[idx_bulkhet])
  if(do_rlog==1){
    L2FCbulk = foldchange2logratio(foldchange(rlogMat[idx_bulkhet,1], rlogMat[idx_bulkhet,2]))  
  }else{L2FCbulk = res_iPS[idx_bulkhet,2]}
  
}else{    
  idx_bulkhet = (rownames(countsHum) %in% dataSC$gene_ids_het)
  data_hetBulk = countsHum[idx_bulkhet,]
  L2FCbulk = foldchange2logratio(foldchange(data_hetBulk[,facT==0], data_hetBulk[,facT==3]))  
}

#plot
indsBulk = match(genes_plot,rownames(countsHum)[idx_bulkhet])
pdf(paste('./bulk/scatterFC',ct1,ct2,as.character(do_rlog),as.character(do_DESeq),'.pdf',sep='_'),width=4.5,height=4.5)
plot(L2FC[indsSC],-L2FCbulk[indsBulk], type = "p",xlab='log2FC SC',ylab='log2FC Bulk', xlim=c(-8,8), ylim = c(-4,4))
dev.off()

r = cor(x=L2FC[indsSC], y=L2FCbulk[indsBulk], method = "p", use = "complete.obs")
r


#FPs 
fdr_vec = sort(c(2^seq(-20,-3,1),0.001,0.01,0.1))

r_vec = c()

indsSCList = list()
indsBulkList = list()
j=1
for(i in fdr_vec){
  DEsc_ = dataSC$gene_ids_het[which(pvadjGiPS<i & is.finite(L2FC))]
  genes_plot = intersect(DEsc_,row.names(res_iPS))
  
  indsSCList[[j]] = match(genes_plot,dataSC$gene_ids_het)
  indsBulkList[[j]] = match(genes_plot,rownames(res_iPS))
  r_ = cor(x=L2FC[indsSCList[[j]]], y=-res_iPS[indsBulkList[[j]],2], method = "p", use = "complete.obs")
  r_vec = c(r_vec, r_)
  j = j+1
}

pdf('./bulk/corrFC.pdf',width=4.5,height=4.5)
plot(fdr_vec, r_vec, type = "p", log="x",xlab='FDR',ylab='correlation')
dev.off()


plot_fdr = c(0.25,0.1, 0.01,0.001, 1e-6)
fdr_all0 = pvadjGiPS[indsSC]
fdr_all = fdr_all0

for(pval in plot_fdr){
  fdr_all[fdr_all0<pval] = pval  
}

#pdf('./L2FCiPSDayr_NoCorr.pdf',width = 5, height = 5)
dfr = data.frame(L2FC[indsSC],-L2FCbulk[indsBulk])
d <- ggplot(dfr,aes(L2FC.indsSC.,X.L2FCbulk.indsBulk.))
d <- d+geom_point(aes(alpha = as.factor(1-fdr_all))) +scale_alpha_discrete(labels=as.character(plot_fdr))+labs(alpha="FDR")
d <- d+xlab("Log2 Fold Change Single Cells")
d <- d+ylab("Log2 Fold Change Bulk")
d + ggtitle('iPS vs Day3')
#dev.off()

