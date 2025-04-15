# ---- 00. Install & Load libraries ----
# cran
cran_pkgs <- c("here", "dplyr", "tidyr", "pander", "devtools","ggplot2")

installed <- rownames(installed.packages())
to_install_cran <- setdiff(cran_pkgs, installed)
if (length(to_install_cran) > 0) install.packages(to_install_cran)

# load cran pkg
lapply(cran_pkgs, library, character.only = TRUE)

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("BiocParallel", "mzR", "CompoundDb", "MsBackendMgf", 
               "AnnotationHub", "Spectra","CluMSID")

for (pkg in bioc_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, force = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Load AnnotationHub and set global object
ah <- AnnotationHub()

# GitHub packages (if not already installed)
if (!require("MetaboAnnotation")) {
  devtools::install_github("rformassspectrometry/MetaboAnnotation")
  library(MetaboAnnotation)
}

if (!require("MsBackendMassbank")) {
  devtools::install_github("rformassspectrometry/MsBackendMassbank")
  library(MsBackendMassbank)
}

# ---- Set working directory ----
setwd(here())
message("Working directory set to: ", getwd())

# ---- Create folders ----
dir.create(here::here("data", "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("output"), recursive = TRUE, showWarnings = FALSE)
#dir.create(here::here("data", "xcms_SERRF"), recursive = TRUE, showWarnings = FALSE)
#dir.create(here::here("data", "Metaboanalyst"), showWarnings = FALSE)
# dir.create(here::here("output", "plots"), recursive = TRUE, showWarnings = FALSE)
# dir.create(here::here("output", "xcms_SERRF"), showWarnings = FALSE)
# dir.create(here::here("output", "Metaboanalyst_input"), showWarnings = FALSE)

# ---- Save session info ----
writeLines(capture.output(sessionInfo()), here::here("session_info.txt"))