'
Make an RDS file with QTL results data from aggregated output

Usage:
make_rds_qtl_results.R --input_prefix <IN> --output_prefix <OUT> [-hv]

Options:
-i --input_prefix <IN>        input filename prefix
-o --output_prefix <OUT>      output filename prefix
-h --help                        show this
-v --version                     print version and stop


This program requires the R packages "docopt", "magrittr", "dplyr", "feather"
(CRAN) and "biomaRt" (Bioconductor).

Davis McCarthy
March 2017
' -> doc

library(magrittr)
library(dplyr)
library(feather)
library(biomaRt)


## Define main function
main <- function(input_prefix, output_prefix) {
    ## load results
    eqtl <- read_feather(paste0(input_prefix, '.aggregated.feather'))
    gene_info <- read_feather(paste0(input_prefix, '.agg_gene_info.feather'))
    eqtl <- eqtl %>% group_by(assoc_gene) %>%
        dplyr::mutate(fdr = p.adjust(pv, method = "BH"),
                      holm = p.adjust(pv, method = "holm"),
                      bonferroni = p.adjust(pv, method = "bonferroni"))
    eqtl_lead <- eqtl %>% group_by(assoc_gene) %>% dplyr::slice(which.min(pv))
    eqtl <- eqtl %>% group_by(assoc_gene) %>%
        dplyr::mutate(lead = (gdid %in% eqtl_lead[["gdid"]][eqtl_lead[["gdid"]]
                                                            == assoc_gene[1]]))
    ## annotate genes
    mm <- match(eqtl_lead[["assoc_gene"]], gene_info[["hgnc_symbol"]])
    eqtl_lead2 <- eqtl_lead %>% ungroup %>%
        mutate(gene_start = gene_info[["start_position"]][mm],
               gene_end = gene_info[["end_position"]][mm],
               gene_strand = gene_info[["strand"]][mm],
               gene_biotype = gene_info[["gene_biotype"]][mm],
               ensembl_gene_id = gene_info[["ensembl_gene_id"]][mm],
               hgnc_symbol = gene_info[["hgnc_symbol"]][mm],
               distance_to_tss = pos - gene_start)
    ens <- useMart("ENSEMBL_MART_ENSEMBL", host = "feb2014.archive.ensembl.org",
                   dataset="hsapiens_gene_ensembl")
    ens87 <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
    tss87 <- getBM(attributes = c("ensembl_gene_id", "transcription_start_site",
                                  "transcript_start", "start_position"),
                   filters = "ensembl_gene_id",
                   values = eqtl_lead2[["ensembl_gene_id"]], mart = ens87)
    tss87 <- as_data_frame(tss87)
    tss87 <- mutate(
        tss87, gene_start_to_tss = start_position - transcription_start_site)
    df_gs_to_tss <- tss87 %>% group_by(ensembl_gene_id) %>%
        summarise(mean_gene_start_to_tss = mean(gene_start_to_tss))
    eqtl_lead2[["gene_start_to_tss"]] <- 0
    mm <- match(df_gs_to_tss[["ensembl_gene_id"]],
                eqtl_lead2[["ensembl_gene_id"]])
    eqtl_lead2[["gene_start_to_tss"]][mm] <-
        df_gs_to_tss[["mean_gene_start_to_tss"]]
    eqtl_lead2 <- mutate(eqtl_lead2, tss = gene_start - gene_start_to_tss)
    eqtl_lead2 <- mutate(eqtl_lead2, corrected_dist_to_tss = pos - tss)
    ## save results to RDS
    saveRDS(eqtl, paste0(output_prefix, '.df_all.rds'))
    saveRDS(eqtl_lead2, paste0(output_prefix, '.df_leads.rds'))
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

## Run main function
main(opt$input_prefix, opt$output_prefix)
