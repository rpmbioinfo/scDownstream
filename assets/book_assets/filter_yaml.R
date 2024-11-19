#!/usr/bin/env Rscript

args <-commandArgs(TRUE)
library(yaml)

test <- TRUE

if(test){
	ymlfile <- read_yaml("_quarto_full.yml")
} else {
	ymlfile <- 	read_yaml("../../../sle_cgamp/nanostring2/_quarto.yml")
}



if(length(args) == 0){
	if(test){
		valid_scripts <- list.dirs(path = "_freeze", recursive = F, full.names = F)
	} else {
		valid_scripts <- list.dirs(path = "../../../sle_cgamp/nanostring2/_freeze", recursive = F, full.names = F)
	}
	ymlfile$book$chapters <- ymlfile$book$chapters[grepl(paste(c("index.qmd", valid_scripts), collapse = "|"),ymlfile$book$chapters)]
} else {
	ymlfile$book$chapters <- c("index.qmd", args[1])
}



if(test){
	write_yaml(ymlfile, file = "_quarto.yml", handlers = list(logical=verbatim_logical))
}



