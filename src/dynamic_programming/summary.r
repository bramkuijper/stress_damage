#!/usr/bin/env Rscript
# get rid of everything there
rm(list=ls())

library("here")

script.dir <- here()

source(file.path(script.dir,"src/dynamic_programming/stress_file_functions.r"))
# obtain a list of all files
all.files <- list.files(path="."
        ,recursive=T
        ,pattern="^stressL.*K\\d.*\\.txt")


summary.data <- NULL

for (filename_i in all.files)
{
    print("filename:")
    print(filename_i)
    params <- read.parameters(filename_i)

    attack.file.name <- gsub(
            pattern="stress"
            ,replacement="simAttacks"
            ,x=filename_i)

    data.attack <- read.attack.file(attack.file.name)

    # maximum hormone after stressor: chronic
    max.h <- data.attack[
                    data.attack$time == 10
                    & data.attack$type == "chronic"
                    ,"hormone"]

    # baseline hormone after stressor 
    h.base <- data.attack[
                    data.attack$time == 45
                    & data.attack$type == "chronic"
                    ,"hormone"]
    
     # maximum hormone after stressor: acute
    max.h.a <- data.attack[
                    data.attack$time == 10
                    & data.attack$type == "acute"
                    ,"hormone"]

    # baseline hormone after stressor: acute
    h.base.a <- data.attack[
                    data.attack$time == 45
                    & data.attack$type == "acute"
                    ,"hormone"]

    params["h_max_c"] <- max.h
    params["h_base_c"] <- h.base
    params["h_max_a"] <- max.h.a
    params["h_base_a"] <- h.base.a
    params["file"] <- filename_i

    summary.data <- rbind(summary.data,params)
}

write_csv(x=summary.data
        ,file="summary.csv")
