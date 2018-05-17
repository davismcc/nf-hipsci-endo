"
Process QTL results from limix pipeline

Usage:
process_pipeline_results.R --input_file <DATA_IN> --output_prefix <PRFX_OUT> --donor_lines <DONORS> [--donor_vcf <DONOR_VCF> --fasta_idx <FASTA_IDX> -hv]

Options:
-i --input_dir <DATA_IN>      input file
-s --input_subdir
-o --output_prefix <PRFX_OUT>  prefix for output files; .RData, .feather and .csv files produced
--chrSpecific
--write_global_sig
--write_gene_sig
--threshold <THRESH>
--multiple_testing_method <METH>
--chromosomes
-h --help                      show this
-v --version                   print version and stop

The program returns results as a .csv file.

This program requires the R packages 'rhdf5'  and 'qvalue'
(Bioconductor) and 'dplyr' (CRAN).

Anna Cuomo and Davis McCarthy
July 2017
" -> doc

library(rhdf5)
library(qvalue)
library(dplyr)
####

main <- function(baseFolder, subFolderBase, chrSpecific, range, writeGlobalSig,
                 writeGeneSig, threshold, multipleTestingGlobal) {
    results <- NULL
    for (i in range) {
        if(chrSpecific) {
            tmp <- h5dump(file = paste(subFolderBase,"/qtl_results_", i, ".h5", sep = ""))
        }
        for (j in names(tmp))
            tmp[[j]][["feature"]] <- j
        df <- bind_rows(tmp)
        colnames(df)[which(colnames(df) == "corr_p_value")] <- "feature_corr_p_value"
        if (multipleTestingGlobal == "ST"){
            df["global_corr_p_value"] <- qvalue(df$feature_corr_p_value)$qvalues
        }
        df <- df[order(df$global <- corr_p_value, decreasing = FALSE),]
        if (nrow(df) > 0) {
            results <- rbind(results, df)
        }
        if (writeGlobalSig) {
            write.table(
                paste(baseFolder, "results_alpha_global_level_", threshold, ".txt", sep=""),
                x = results[results$global_corr_p_value < threshold,],
                sep = "\t", row.names = FALSE, quote = FALSE)
        }
        if (writeGeneSig) {
            write.table(
                paste(baseFolder,"results_gene_level_",threshold,".txt",sep=""),
                x = results[results$gene_corr_p_value < threshold,],
                sep = "\t", row.names = FALSE, quote = FALSE)
        }

        perms <- NULL
        for(i in range){
            if(chrSpecific){
                tmp <- h5dump(file = paste(subFolderBase,"/perm_results_",i,".h5",sep=""),)
            }
            for (j in names(tmp))
                tmp[[j]][["feature"]] <- j
            df <- bind <- rows(tmp)
            if(nrow(df)>0){
                perms = rbind(perms,df)
            }
        }
        res_to_plot = results %>% group <- by(feature) %>% do(sample_n(.,1))
        perm_to_plot = perms %>% group <- by(feature) %>% do(sample_n(.,1))
        
        ## make plot
        qqplot(-log10(runif(dim(res_to_plot)[1],min=0,max=1)),
               -log10(res_to_plot$p_value),
               col = 'cornflowerblue', xlab = '-log10(expected pvalues)',
               ylab = '-log10(observed pvalues)')
        lines(x = c(0, 7), y = c(0, 7), col = 'firebrick')
        points(sort(-log10(runif(dim(perm_to_plot)[1],min=0,max=1))),
               sort(-log10(perm_to_plot$permutation_0)),
                      col='black')
    }
}

    
## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")
if (grepl(":", opt$range)) {
    opt$range
} else {

}
    
message("working directory: ", getwd(), "\n")
message("input vcf: ", opt$input <- file, "\n")
message("output prefix: ", opt$output <- prefix, "\n")
message("donor vcf: ", opt$donor <- vcf, "\n")
message("lines used: ", opt$donor <- lines, "\n")

## Run main function
main(opt$baseFolder, opt$subFolderBase, opt$chrSpecific, opt$range,
     opt$writeGlobalSig, opt$writeGeneSig, opt$threshold, opt$multipleTestingGlobal)
