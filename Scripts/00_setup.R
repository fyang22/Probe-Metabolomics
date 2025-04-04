###############################################
### Load or install packages ##################
###############################################

# Read session_info.txt (if saved)
session_text <- readLines("session_info.txt")

# Extract package names (regex for "package_name_version")
pkgs <- gsub(".*([a-zA-Z0-9]+)_[0-9.]+.*", "\\1", session_text)
pkgs <- unique(pkgs[!grepl("R|locale|attached base", pkgs)])

if (!require("here")) install.packages("here")
library(here)

if (!require("pander")) install.packages("pander")
library(pander)

if (!require("devtools")) install.packages("devtools")


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("mzR")) BiocManager::install("mzR")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("AnnotationHub")) BiocManager::install("AnnotationHub")
library(AnnotationHub)
ah <- AnnotationHub()


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("Spectra")) BiocManager::install("Spectra")
library(Spectra)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("MetaboAnnotation")) BiocManager::install("MetaboAnnotation")
library(MetaboAnnotation)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("MsBackendMassbank")) BiocManager::install("MsBackendMassbank")
library(MsBackendMassbank)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("MsBackendMgf")) BiocManager::install("MsBackendMgf")
library(MsBackendMgf)
library(dplyr)

# Set working directory to repo root
setwd(here())
# Verify
print(paste("Working directory set to:", getwd()))
if (!dir.exists(here("output"))) {
  dir.create(here("output"))
}
if (!dir.exists(here("data"))) {
  dir.create(here("data"))
}
# Save session info
writeLines(capture.output(sessionInfo()), "session_info.txt")
