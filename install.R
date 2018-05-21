source("https://bioconductor.org/biocLite.R")
library(BiocInstaller) # shouldn't be necessary


pkgs <- c(
  "apcluster",
  "batchtools",
  "bayesplot",
  "bblme",
  "clustermq",
  "coda",
  "cowplot",
  "d3heatmap",
  "data.table",
  "devtools",
  "docopt",
  "DT",
  "dynamicTreeCut",
  "e1071",
  "future",
  "future.batchtools",
  "flexmix",
  "formatR",
  "fpc",
  "GGally",
  "ggbeeswarm",
  "ggdendro",
  "ggmcmc",
  "ggpubr",
  "ggrepel",
  "ggridges",
  "ggsci",
  "ggthemes",
  "ggtree",
  "glmnet",
  "gdata",
  "gplots",
  "gtools",
  "greta",
  "keras",
  "lattice",
  "latticeExtra",
  "lintr",
  "lme4",
  "Matrix",
  "MatrixModels",
  "matrixStats",
  "microbenchmark",
  "mvoutlier",
  "NMF",
  "packrat",
  "pheatmap",
  "pryr",
  "RColorBrewer",
  "reshape2",
  "roxygen2",
  "rprojroot",
  "scales",
  "scater",
  "snpStats",
  "superheat",
  "tensorflow",
  "testthat",
  "tidyverse",
  "tufte",
  "UpSetR",
  "vcfR",
  "VariantAnnotation",
  "VGAM",
  "viridis",
  "wesanderson",
  "xtable"
)

ap.db <- available.packages(contrib.url(biocinstallRepos()))
ap <- rownames(ap.db)

pkgs_to_install <- pkgs[pkgs %in% ap]

biocLite(pkgs_to_install)
# install.packages("rmote", repos = c(getOption("repos"), "http://cloudyr.github.io/drat"))
# greta::install_tensorflow()

# Single-cell relevant packages
sc_pkgs <- c(
  "beachmat",
  "Canopy",
  "clusterExperiment",
  "DESeq2",
  "destiny",
  "edgeR",
  "GO.db",
  "goseq",
  "limma",
  "MAST",
  "MultiAssayExperiment",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "pcaMethods",
  "pheatmap",
  "preprocessCore",
  "rhdf5",
  "Rsamtools",
  "Rsubread",
  "Rtsne",
  "scater",
  "scran",
  "slalom",
  "snpStats",
  "tximport",
  "variancePartition",
  "vcfR",
  "zinbwave"
)

pkgs_to_install <- sc_pkgs[sc_pkgs %in% ap]

biocLite(pkgs_to_install)

devtools::install_github("davismcc/cardelino")

## just in case there were warnings, we want to see them
## without having to scroll up:
warnings()

if (!is.null(warnings()))
{
  w <- capture.output(warnings())
  if (length(grep("is not available|had non-zero exit status", w)))
    quit("no", 1L)
}

suppressWarnings(BiocInstaller::biocValid(fix=TRUE, ask=FALSE))

