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


#1. REad in data
#iPS


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

labKey = as.factor(labels)
facT = as.factor(labT)
facC = as.factor(labC)
rownames(countsHum)<-rn[which( geneTypes=="ENSG" )]

gene_names = rownames(countsHum)
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
xxenseg <- as.list(org.Hs.egENSEMBL2EG)
gene_syms=unlist(xx[unlist(xxenseg[gene_names])])
gene_names_list<-(lapply(xxenseg[gene_names],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
sym_names=unlist(lapply(xx[unlist(gene_names_list)],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
sym_names[is.na(sym_names)]=gene_names[is.na(sym_names)]


colData <- DataFrame(time = facT, cell_type =facC)
dds <- DESeqDataSetFromMatrix(countsHum, colData, formula(~ time+cell_type+time:cell_type))

dds <- DESeq(dds)
res <- results(dds)
resType <- results(dds, contrast = list(
                          c("cell_typeffdm_2", "time0.cell_typeffdm_2"),
                          c("cell_typeffdj_1", "time0.cell_typeffdj_1"))) 


resType <- results(dds, contrast = c("cell_type", "ffdm_2", "ffdj_1"))
resType$syms = as.character(sym_names)


res$syms = as.character(sym_names)

resOrdered <- res[order(res$padj),]

head(resOrdered)


resTypeOrdered <- resType[order(resType$padj),]
head(resTypeOrdered)

pdf('./bulk/MAbulk.pdf',width=4.5,height=4.5)
plotMA(dds, main="MA plot", ylim=c(-2,2))
dev.off()

rld <- rlogTransformation(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

distsRL <- dist(t(assay(rld)))

mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- labKey
hc <- hclust(distsRL)
pdf('./bulk/Heatmapbulk.pdf',width=4.5,height=4.5)
heatmap.2(mat, symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13),Rowv = TRUE,
          Colv=T)
dev.off()



pdf('./bulk/PCAbulk.pdf',width=4.5,height=4.5)
plotPCA(rld, intgroup=c("time", "cell_type"))
dev.off()




###load DE genes form SC analysis 
dataSC = h5dump('./normCountsHum_all2.h5f')
dataCorr= h5dump('./varsHipsciCLVM_CcGenoPluriIntGCC.hdf5')
is_het = dataSC$genes_heterogen
is_03 = dataSC$cell_type=='iPS' | dataSC$cell_type=='Day3'
cell_type03 = dataSC$cell_type[is_03]
ctF = as.factor(cell_type03)
geno03 = dataSC$geno_id[is_03]
gtF = as.factor(geno03)
data_het = t(dataSC$LogNcountsHum[which(is_het==1),])
data_het = data_het[is_03,]

data_hetCorr = dataCorr$Ycorr
data_hetCorr = data_hetCorr[is_03,]

countsHumSC = dataSC$countsHumAll[,is_03]
rownames(countsHumSC) = dataSC$gene_ids


colData <- DataFrame(time = ctF, cell_type =gtF)
dds <- DESeqDataSetFromMatrix(countsHumSC, colData, formula(~ cell_type+time+cell_type:condition))
ddsSC <- DESeq(dds)
resSC <- results(ddsSC)

resSC$syms = as.character(dataSC$gene_ids)

resScOrdered <- resSC[order(resSC$padj),]

head(resScOrdered)

ensDESc<-resScOrdered$syms[which(resScOrdered$padj<0.1)]

genes_plot = intersect(resScOrdered$syms,DEbulkHet)

indsSC = match(genes_plot,resScOrdered$syms)
indsBulk = match(genes_plot,rownames(resOrdered))

pdf('./bulk/scatterFC.pdf',width=4.5,height=4.5)
plot(-resOrdered$log2FoldChange[indsBulk],resScOrdered$log2FoldChange[indsSC], type = "p",xlab='log2FC Bulk',ylab='log2FC SC')
dev.off()

r = cor(x=resScOrdered$log2FoldChange[indsSC], y=-resOrdered$log2FoldChange[indsBulk], method = "s", use = "complete.obs")


fdr_vec = c(2^seq(-15,-2,1))
r_vec = c()

for(i in fdr_vec){
  DEbulk_ = row.names(resOrdered)[which(resOrdered$padj<i)]
  DEbulkHet_ = intersect(DEbulk_, dataSC$gene_ids_het)
  genes_plot = intersect(resScOrdered$syms,DEbulkHet_)
  
  indsSC_ = match(genes_plot,resScOrdered$syms)
  indsBulk_ = match(genes_plot,rownames(resOrdered))
  r_ = cor(x=resScOrdered$log2FoldChange[indsSC_], y=-resOrdered$log2FoldChange[indsBulk_], method = "s", use = "complete.obs")
  r_vec = c(r_vec, r_)
}

pdf('./bulk/corrFC.pdf',width=4.5,height=4.5)
plot(fdr_vec, r_vec, type = "p", log="x",xlab='FDR',ylab='correlation')
dev.off()



#get het DE genes from bulk
DEbulk = row.names(resOrdered)[which(resOrdered$padj<.1)]
DEbulkHet = intersect(DEbulk, dataSC$gene_ids_het)

DEbulkType = row.names(resTypeOrdered)[which(resTypeOrdered$padj<.1)]
DEbulkTypeHet = intersect(DEbulkType, dataSC$gene_ids_het)




pv = c()
#lfc = c()
for(i in 1:dim(data_het)[2]){
  res = wilcox.test(data_het[cell_type03=='iPS',i], data_het[cell_type03!='iPS',i])
  pv = c(pv, res$p.value)
  #lfc = c(lfc,log2(data_het[mean(cell_type03!='iPS',i])/mean(data_het[cell_type03!='iPS',i]))
}

pvadj = p.adjust(pv, 'fdr')
  
DEsceHet = dataSC$gene_ids_het[which(pvadj<.1)]

#DE SC 
pvG = c()
for(i in 1:dim(data_het)[2]){
  res_ = wilcox.test(data_het[geno03=='coxy',i], data_het[geno03=='ffdj_1',i])
  pvG = c(pvG, res_$p.value)
}

pvadjG = p.adjust(pvG, 'fdr')



DEsceHetG = dataSC$gene_ids_het[which(pvadjG<.1)]




#ffdm2 ffdj_1
pvGiPS = c()
for(i in 1:dim(data_het)[2]){
  res_ = wilcox.test(data_het[geno03=='ffdm_2' & cell_type03=='iPS',i], data_het[geno03=='ffdj_1' & cell_type03=='iPS',i])
  pvGiPS = c(pvGiPS, res_$p.value)
}


L2FC = foldchange2logratio(foldchange(apply(data_het[geno03=='ffdm_2' & cell_type03=='iPS',],2,mean), apply(data_het[geno03=='ffdj_1' & cell_type03=='iPS',],2,mean)))

pvadjGiPS = p.adjust(pvGiPS, 'fdr')
DEsceHetGiPS = dataSC$gene_ids_het[which(pvadjGiPS<.1 & is.finite(L2FC))]

genes_plot = DEsceHetGiPS

indsSC = match(genes_plot,dataSC$gene_ids_het)

idx_bulkhet = (rownames(countsHum) %in% dataSC$gene_ids_het)
indsBulk = match(genes_plot,rownames(countsHum)[idx_bulkhet])

data_hetBulk = log10(countsHum+1)[idx_bulkhet,c(1,3)]
L2FCbulk = foldchange2logratio(foldchange(data_hetBulk[,1], data_hetBulk[,2]))

pdf('./bulk/scatterFC_ffdm2vsffdj1.pdf',width=4.5,height=4.5)
plot(L2FC[indsSC],L2FCbulk[indsBulk], type = "p",xlab='log2FC Bulk',ylab='log2FC SC', xlim=c(-8,8), ylim = c(-4,4))
dev.off()

r = cor(x=L2FC[indsSC], y=L2FCbulk[indsBulk], method = "s", use = "complete.obs")


#corrected data
pvGiPS = c()
for(i in 1:dim(data_het)[2]){
  res_ = wilcox.test(data_hetCorr[geno03=='ffdm_2' & cell_type03=='iPS',i], data_hetCorr[geno03=='ffdj_1' & cell_type03=='iPS',i])
  FC_ = foldchange(mean(data_hetCorr[geno03=='ffdm_2' & cell_type03=='iPS',i]), mean(data_hetCorr[geno03=='ffdj_1' & cell_type03=='iPS',i]))
  
  pvGiPS = c(pvGiPS, res_$p.value)
}

L2FC = foldchange2logratio(foldchange(apply(data_hetCorr[geno03=='ffdm_2' & cell_type03=='iPS',],2,mean), apply(data_hetCorr[geno03=='ffdj_1' & cell_type03=='iPS',],2,mean)))

pvadjGiPS = p.adjust(pvGiPS, 'fdr')
DEsceHetGiPS = dataSC$gene_ids_het[which(pvadjGiPS<.1)]


genes_plot = DEsceHetGiPS

indsSC = match(genes_plot,dataSC$gene_ids_het)

idx_bulkhet = (rownames(countsHum) %in% dataSC$gene_ids_het)
indsBulk = match(genes_plot,rownames(countsHum)[idx_bulkhet])

data_hetBulk = log10(countsHum+1)[idx_bulkhet,c(1,3)]

pdf('./bulk/scatterFC_ffdm2vsffdj1Corr.pdf',width=4.5,height=4.5)
plot(L2FCbulk[indsBulk],L2FC[indsSC], type = "p",xlab='log2FC Bulk',ylab='log2FC SC', xlim=c(-4,4), ylim = c(-8,8))
dev.off()





#DE bulk 
colData_iPS <- DataFrame(cell_type = facC[facT==0 & (facC=='ffdj_1' | facC=='ffdm_2') ])
dds_iPS <- DESeqDataSetFromMatrix(countsHum[,facT==0  & (facC=='ffdj_1' | facC=='ffdm_2')], colData_iPS, formula(~ cell_type))
rld <- rlogTransformation(dds_iPS)
rlogMat <- assay(rld)


dds_iPS <- DESeq(dds_iPS)
#res_iPS <- results(dds_iPS)
#DEbulk_iPS = row.names(res_iPS)[which(res_iPS$padj<.1)]
res_iPSType <- results(dds_iPS, contrast = c("cell_type", "ffdm_2", "ffdj_1"))

DEbulk_iPSType = row.names(res_iPSType)[which(res_iPSType$padj<.1)]


###
#DE in iPS good vs bad
idx_bad = (geno03=='ffdm_2' | geno03=='ffdj_1') & cell_type03=='iPS'
colData <- DataFrame(cell_type =gtF[idx_bad])
dds_SC <- DESeqDataSetFromMatrix(countsHumSC[,idx_bad], colData, formula(~ cell_type))
rldSC <- rlogTransformation(dds_SC)
rlogMatSC <- assay(rldSC)










