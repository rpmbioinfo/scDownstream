#!/usr/bin/env Rscript

args <-commandArgs(TRUE)
library(yaml)



ymlfile <- read_yaml("_quarto_full.yml")




if(length(args) == 0){
	valid_scripts <- list.dirs(path = "_freeze", recursive = F, full.names = F)
	ymlfile$book$chapters <- ymlfile$book$chapters[grepl(paste(c("index.qmd", valid_scripts), collapse = "|"),ymlfile$book$chapters)]
} else {
	ymlfile$book$chapters <- c("index.qmd", args[1])
}




write_yaml(ymlfile, file = "_quarto.yml", handlers = list(logical=verbatim_logical))




