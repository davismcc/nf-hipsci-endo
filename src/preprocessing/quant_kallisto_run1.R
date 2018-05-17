# Quantify expression with kallisto for run1 data
# Davis McCarthy
# November 2015

## setup the directory and index
setwd('../data/run1/quantKallisto/')
kallisto_idx <- "/nfs/research/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.rel75.cdna.all.ERCC92.kallisto_idx"
library(scater)

## returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

## Define targets file
fastq_left <- dir("../raw/")[grep("_1.fastq$", dir("../raw/"))]
sample_names <- gsub("#", "_", gsub("_[12].fastq", "", fastq_left))
fastq_left <-  paste0("../raw/", fastq_left)
fastq_left <- trim(fastq_left)
fastq_right <- gsub("_1.fastq", "_2.fastq", fastq_left)
identical(gsub("_[12].fastq", "", fastq_left),
          gsub("_[12].fastq", "", fastq_right))
targets_df <- data.frame(Sample = sample_names, FASTQ_1 = fastq_left,
                         FASTQ_2 = fastq_right)
head(targets_df)
write.table(targets_df, row.names = FALSE, file = "targets_run1_kallisto.txt",
            quote = TRUE, sep = "\t")


## Run kallisto on the fastq files
kallisto_run1 <- runKallisto("targets_run1_kallisto.txt",
                             kallisto_idx, output_prefix = "out",
                             n_cores = 18, single_end = FALSE,
                             correct_bias = TRUE, verbose = TRUE,
                             dry_run = TRUE)

## Submit kallisto jobs to cluster
cluster_prefix <- 'bsub -n 4 -q research-rh6 -R "rusage[mem=16000,tmp=50000]" -M 16000 -o ./cluster_out'
for( i in seq_len(length(kallisto_run1)) ) {
    cmd <- paste0(cluster_prefix, "_", i, " ",  kallisto_run1[[i]]$kallisto_call)
    cmd <- gsub(" kallisto ",
                " /homes/davis/davis/src/kallisto_linux-v0.42.4/kallisto ",
                cmd)
    print(cmd)
    system(cmd)
    Sys.sleep(1)
}

## Read kallisto results in R
sce_run1_kallisto_tx <- readKallistoResults(kallisto_run1, read_h5 = TRUE)

## Look at results
which(apply(exprs(sce_run1_kallisto_tx), 2, function(x) {any(is.na(x))}))
sce_run1_kallisto_tx

# load('scesets_run1_kallisto_preqc.RData')

## Get feature annotations from biomaRt
sce_run1_kallisto_tx<- getBMFeatureAnnos(
    sce_run1_kallisto_tx, filters = "ensembl_transcript_id",
       attributes = c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol",
                      "chromosome_name", "transcript_biotype", "transcript_start",
                      "transcript_end", "transcript_count"),
    feature_symbol = "hgnc_symbol", feature_id = "ensembl_gene_id",
    dataset = "hsapiens_gene_ensembl", host = "feb2014.archive.ensembl.org")

head(fData(sce_run1_kallisto_tx))
sum(is.na(fData(sce_run1_kallisto_tx)$hgnc_symbol))

## Summarise expression at gene level
sce_run1_kallisto_gene <- summariseExprsAcrossFeatures(sce_run1_kallisto_tx)

sce_run1_kallisto_gene <- getBMFeatureAnnos(
    sce_run1_kallisto_gene, filters = "ensembl_gene_id",
       attributes = c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol",
                      "chromosome_name", "gene_biotype"),
    feature_symbol = "hgnc_symbol", feature_id = "ensembl_gene_id",
    dataset = "hsapiens_gene_ensembl", host = "feb2014.archive.ensembl.org")

head(fData(sce_run1_kallisto_gene))
sum(is.na(fData(sce_run1_kallisto_gene)$hgnc_symbol))

## Save SCESet objects
save(sce_run1_kallisto_tx, sce_run1_kallisto_gene,
     file = 'scesets_run1_kallisto_preqc.RData', compress = 'gzip',
     compression_level = 9)
