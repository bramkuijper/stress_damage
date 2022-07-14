#!/usr/bin/env Rscript
# get rid of everything there
rm(list=ls())

library("here")

script.dir <- here()

source(file.path(script.dir,"src/dynamic_programming/stress_file_functions.r"))
# obtain a list of all files
all.files <- list.files(path="."
        ,recursive=T
        ,pattern="^stress_strategy.*\\.txt")
        #        ,pattern="^stressL.*Kfec\\d.*\\.txt")


summary.data <- NULL

fwdcalc.data <- T

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
            pattern="stress_strategy"
            ,replacement="simAttacks"
            ,x=filename_i)

    data.attack <- read.attack.file(filename=attack.file.name)
   
    # get the fwd calc data 
    fwdcalc.file.name <- gsub(
            pattern="stress_strategy"
            ,replacement="fwdCalc"
            ,x=filename_i)
    
    if (fwdcalc.data && file.exists(fwdcalc.file.name))
    {
        fwdcalc.data <- read.fwdcalc(fwdcalc.file.name)
    }
    else
    {
        fwdcalc.data <- F
    }
    

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
    max.damage.a <- max(data.attack[
                    data.attack$type == "acute"
                    ,"damage"])
    
    max.damage.c <- max(data.attack[
                    data.attack$time <= 50 &
                    data.attack$type == "chronic"
                    ,"damage"])

    params["h_max_c"] <- max.h
    params["h_base_c"] <- h.base
    params["h_max_a"] <- max.h.a
    params["h_base_a"] <- h.base.a
    params["max_damage_a"] <- max.damage.a
    params["max_damage_c"] <- max.damage.c
    params["file"] <- filename_i

    # add fwdcalc data
    if (fwdcalc.data)
    {
        params <- cbind(params,fwdcalc.data)
    }

    summary.data <- rbind(summary.data,params)
}

write_csv(x=summary.data
        ,file="summary.csv")
