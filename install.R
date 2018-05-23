source("https://bioconductor.org/biocLite.R")
library(BiocInstaller) # shouldn't be necessary


pkgs <- c(
  "apcluster",
  "bblme",
  "cowplot",
  "data.table",
  "devtools",
  "docopt",
  "DT",
  "dynamicTreeCut",
  "ggridges",
  "ggthemes",
  "glmnet",
  "gdata",
  "gplots",
  "Matrix",
  "MatrixModels",
  "matrixStats",
  "microbenchmark",
  "mvoutlier",
  "packrat",
  "pheatmap",
  "pryr",
  "qvalue",
  "RColorBrewer",
  "reshape2",
  "roxygen2",
  "rprojroot",
  "scales",
  "superheat",
  "VGAM",
  "viridis"
)

# check that desired packages are available
ap.db <- available.packages(contrib.url(biocinstallRepos()))
ap <- rownames(ap.db)
pkgs_to_install <- pkgs[pkgs %in% ap]

# do not reinstall packages that are already installed in the image
ip.db <- installed.packages()
ip <- rownames(ip.db)
pkgs_to_install <- pkgs_to_install[!(pkgs_to_install %in% ip)]

biocLite(pkgs_to_install)
# install.packages("rmote", repos = c(getOption("repos"), "http://cloudyr.github.io/drat"))
# greta::install_tensorflow()

# Single-cell relevant packages
sc_pkgs <- c(
  "beachmat",
  "Canopy",
  "destiny",
  "edgeR",
  "limma",
  "preprocessCore",
  "Rtsne",
  "scater",
  "scran",
  "slalom",
  "snpStats",
  "tximport",
  "VariantAnnotation",
  "vcfR"
)

pkgs_to_install <- sc_pkgs[sc_pkgs %in% ap]
# do not reinstall packages that are already installed in the image
ip.db <- installed.packages()
ip <- rownames(ip.db)
pkgs_to_install <- pkgs_to_install[!(pkgs_to_install %in% ip)]

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

