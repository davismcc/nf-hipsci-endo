library(diffusionmap)
library(monocle)
library(rhdf5)
setwd('~/projects/hipsci-singlecell//R')

###load DE genes form SC analysis 
dataCorr = h5dump('../analysis/varsHipsciCLVMfac_CcGenoPluriIntGPl.hdf5')
dataSC = h5dump('./normCountsHum_all2.h5f')
is_het = dataSC$genes_heterogen
is_03 = dataSC$cell_type=='Day3'
cell_type = dataSC$cell_type
ctF = as.factor(cell_type)
geno03 = dataSC$geno_id[is_03]
gtF = as.factor( dataSC$geno_id)
data_het = t(dataSC$LogNcountsHum[which(is_het==1),])
gene_short_name = dataSC$sym_names_het
data_hetCorr = dataCorr$Ycorr

countsHumSC = dataSC$countsHumAll
rownames(countsHumSC) = dataSC$gene_ids
#Psi <- diffuse(data_het, 10^1.8)
PsiCorr <- diffuse(data_hetCorr, 10^1.3)

#X11()
#plot(Psi$vector[,1], Psi$vector[,2])

#X11()
psi = PsiCorr$vector
DC1=psi[,1]
DC2=psi[,2]

pdf('./diffmap.pdf',width = 5, height = 5)
dfr = data.frame(DC1,DC2)
d <- qplot(DC1,DC2,data=dfr,colour=ctF)
d + ggtitle('Diffusion map, corrected data')
dev.off()



###uncorrected data
Psi <- diffuse(data_het, 10^1.8)
psi = Psi$vector
DC1=psi[,1]
DC2=psi[,2]

pdf('./diffmap_nocorr.pdf',width = 5, height = 5)
dfr = data.frame(DC1,DC2)
d <- qplot(DC1,DC2,data=dfr,colour=ctF)
d + ggtitle('Diffusion map, uncorrected data')
dev.off()

###ICA
ICs <- fastICA(data_het, 2,method = 'C')



###run monocle###

times0 = cell_type

Y = data_hetCorr

genSym = as.data.frame(unlist(dataSC$gene_ids_het))


rownames(genSym) = as.character(genSym[,])
colnames(genSym) = 'GenSym'
times = as.data.frame(times0)
colnames(times) = 'Time'

fd <- new("AnnotatedDataFrame", data = genSym)
pd <- new("AnnotatedDataFrame", data = times)


cnts = as.matrix(t(Y))
colnames(cnts)<-as.character(seq(1,dim(Y)[1]))
rownames(cnts)<-as.character(seq(1,dim(genSym)[1]))

cds = new("CellDataSet",exprs=matrix(10^t(Y), nrow=dim(Y)[2], ncol=dim(Y)[1]), featureData=fd, phenoData = pd,
          expressionFamily = VGAM::tobit(Lower =0.1))

#cds <- reduceDimension(cds, use_irlba=F)
#cds <- orderCells(cds, num_paths=1, reverse=T)

reducedDimS(cds)<-t(psi)
cds <- orderCells(cds, num_paths=1, reverse=T)

pdf('./spanning_treeNoCorr.pdf', height = 5, width = 5)
plot_spanning_tree(cds)
dev.off()

genes_inds = which(as.character(genSym[,]) %in% c('ENSG00000175387'))
pdf('./genes_pseudotime_diffmapCorr_Smad2.pdf', height=4, width=4)
plot_genes_in_pseudotime(cds[genes_inds,ctF != 'iPS'], color_by="Time")
dev.off()



#test a few genes

genes_DE <- read.table('../analysis/results/CcGenoPluriInt/lists/genes_CcGenoPluriInt_DE_Day2.txt')
genes_DE = as.character(genes_DE[,1])
to_be_tested = which(dataSC$sym_names_het %in% genes_DE)#[1:100])
tested_genes = dataSC$sym_names_het[to_be_tested]


cds_subset <- cds[to_be_tested, ctF != 'iPS']

diff_test_res <-differentialGeneTest(cds_subset, fullModelFormulaStr="expression~sm.ns(Pseudotime)")
diff_test_res <- merge(fData(cds), diff_test_res, by="row.names")
row.names(diff_test_res) = gene_short_name[to_be_tested]
diff_test_res[diff_test_res[,"qval"]<.001,]

pdf('./pseudotime.pdf', height=3, width=5)
plot_genes_in_pseudotime(cds_subset, color_by="Time")
dev.off()

#cluster genes
full_model_fits <- fitModel(cds_subset, modelFormulaStr="expression~sm.ns(Pseudotime, df=3)")
expression_curve_matrix <- responseMatrix(full_model_fits)

clusters <- clusterGenes(expression_curve_matrix, k=3)

pdf('./pseudotime_clusters3AllGenes.pdf', height=4, width=6)
plot_clusters(cds_subset, clusters)
dev.off()

pdf('./pseudotime_clusters_medoids3AllGenes.pdf', height=4, width=4)
plot_genes_in_pseudotime(cds_subset[clusters$id.med,], color_by="Time")
dev.off()

tested_genes[clusters$clustering==2]

clusteringAll = cbind(tested_genes,clusters$clustering)
write.table(clusteringAll, file = './clustering_timeAll.txt', quote = F)
