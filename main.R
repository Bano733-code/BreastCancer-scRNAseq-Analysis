# ==================================
# Single-cell RNA-seq Analysis
# Breast Cancer GSE176078
# ==================================


# Load libraries
library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
# Check working directory

getwd()
# Check dataset files

list.files("data/GSE176078/GSE176078")
counts <- ReadMtx(
  mtx = "data/GSE176078/GSE176078/matrix.mtx",
  cells = "data/GSE176078/GSE176078/barcodes.tsv",
  features = "data/GSE176078/GSE176078/genes.tsv",
  feature.column = 1
)
