#!/usr/bin/env Rscript
# get rid of everything there
rm(list=ls())

source("stress_file_functions.r")

# obtain a list of all files
all.files <- list.files(path="."
        ,recursive=F
        ,pattern="^stressL.*K\\d.*\\.txt")


summary.data <- NULL

for (filename_i in all.files)
{
    params <- read.parameters(filename_i)

    attack.file.name <- gsub(
            pattern="stress"
            ,replacement="simAttacks"
            ,x=strategy.file.name)

    data.attack <- read.attack.file(attack.file.name)

    max.h <- max(data.attack$hormone)

    summary.data <- rbind(summary.data,params)
}
