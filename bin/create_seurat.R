#!/usr/bin/env Rscript

args <-commandArgs(TRUE)
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))


library(AnnotationDbi)
library(AnnotationHub)
library(rtracklayer)




if(length(args) == 0){
  ## Default Annotation Hub is for Homo_sapiens
  genome <- "AH75393"
} else {
   ah_index <- args[1]
}
  

ah <- AnnotationHub()


test <- query(ah, c("Macaca mulatta", "release-98"))
test2 <- query(ah, c("Homo sapiens", "release-98"))


gtf_chr <- ah[["AH75429"]]

gtf <- ah[["AH75430"]]