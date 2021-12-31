#!/usr/bin/env Rscript
# get rid of everything there
rm(list=ls())

library("here")

script.dir <- here()

source(file.path(script.dir,"src/dynamic_programming/stress_file_functions.r"))
# obtain a list of all files
all.files <- list.files(path="."
        ,recursive=F
        ,pattern="^stressL.*K\\d.*\\.txt")


summary.data <- NULL

for (filename_idx in 1:length(all.files))
{
    filename_i <- all.files[[filename_idx]]

    #    if (filename_idx > 10)
    #    {
    #        break;
    #    }

    print("filename:")
    print(filename_i)
    params <- read.parameters(filename_i)

    attack.file.name <- gsub(
            pattern="stress"
            ,replacement="simAttacks"
            ,x=filename_i)

    data.attack <- read.attack.file(attack.file.name)
   
    # get the fwd calc data 
    fwdcalc.file.name <- gsub(
            pattern="stress"
            ,replacement="fwdCalc"
            ,x=filename_i)
    
    fwdcalc.data <- read.fwdcalc(fwdcalc.file.name)

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

    # maximum damage level
    max.damage <- max(data.attack$damage)

    params["h_max_c"] <- max.h
    params["h_base_c"] <- h.base
    params["h_max_a"] <- max.h.a
    params["h_base_a"] <- h.base.a
    params["max_damage"] <- max.damage
    params["file"] <- filename_i

    # add fwdcalc data
    params <- cbind(params,fwdcalc.data)

    summary.data <- rbind(summary.data,params)
}

write_csv(x=summary.data
        ,file="summary.csv")
