library("ggplot2")
library("tidyverse")
library("patchwork")
library("here")

# get general stress functionality script
proj.dir <- here()
stress_functions_file <- file.path(
              proj.dir
              ,"src/dynamic_programming/stress_file_functions.r")
source(stress_functions_file)


#' Function that obtains a summarized time 
#' profile by averaging over several simulations
summarize_profile <- function(filenames)
{
    all.data <- NULL
    for (file.idx in 1:length(filenames))
    {
        filename.i <- filenames[[file.idx]]

        basename.i <- basename(filename.i)
        dirname.i <- dirname(filename.i)
        
        # get the file name of the associated simAttacks file
        # focus on basename only and replace stress by simAttacks
        basename.i <- gsub(pattern = "^stress"
                            ,replacement = "simAttacks"
                            ,x=basename.i)

        # put dirname back in front
        attack.file <- file.path(dirname.i,basename.i)

        # get the hormone profile data
        hormone.data.i <- read.attack.file(attack.file)
        
        hormone.data.i[,"idx"] <- file.idx
        
        # then rbind it 
        all.data <- rbind(all.data, hormone.data.i)
        
        print(head(all.data))
    }
    return(all.data)
} #end function summarize_profile

# get summary data from the replicate simulations
stress.summary <- read_csv(file="summary_hormone_profile.csv")

# filter data on 
stress.summary$risk <- with(stress.summary,
        pArrive / (pLeave + pArrive))

stress.summary$autocorrelation <- with(stress.summary,
        1.0 - pArrive - pLeave)

autocorr.list <- c(0.3,0.7)

for (autocorr.i in autocorr.list)
{
    stress.summary.sub <- filter(stress.summary, 
            risk == 0.05 & abs(autocorrelation - autocorr.i) < 0.01) 

    stopifnot(nrow(stress.summary.sub) > 0)
    readline(prompt=paste0("There are ",nrow(stress.summary.sub)," simulations here, OK? [Enter]"))

    # the two repair values we want to plot
    repair.vals <- c(0,1)

    # make a list of plots to print later using patchwork
    plots <- list()

    hormone_time_data <- NULL

    for (idx in 1:2)
    {
        filenames.i <- stress.summary.sub %>% 
            filter(repair == repair.vals[[idx]]) %>%
            select(file) %>%
            pull(file)
        
        profile.data <- summarize_profile(filenames.i)
        
        hormone_by_time <- profile.data %>%
            filter(type == "acute") %>%
            group_by(time) %>%
            summarise(
                mean_damage = mean(damage)
                ,mean_hormone = mean(hormone)
            )
        
        hormone_by_time[,"repair"] <- repair.vals[[idx]]
        
        hormone_time_data <-bind_rows(hormone_time_data, hormone_by_time)
    } # end for idx

    # add a column for the panel label
    hormone_time_data[,"repair_text"] <- with(
        hormone_time_data
        ,paste("Repair: ",repair)
    )

    # scale damage between 0 and 1
    hormone_time_data <- hormone_time_data %>% 
        mutate(
            mean_damage_scaled = mean_damage / max(hormone_time_data$mean_damage)
        )

    p1 <- hormone_time_data %>% 
        ggplot(mapping=aes(x=time)) +
        geom_line(mapping=aes(y=mean_hormone)) +
        facet_grid(.~ repair_text) +
        theme_classic() +
        theme(
            strip.background = element_rect(
                color="transparent"
            )
        ) +
        ylab("Mean hormone") +
        xlab("") 

    p2 <- hormone_time_data %>% 
        ggplot(mapping=aes(x=time)) +
        geom_line(mapping=aes(y=mean_damage)) +
        facet_grid(.~ repair_text) +
        theme_classic() +
        theme(
            strip.background = element_rect(
                color="transparent"
            ),
            strip.text.x = element_text(color="transparent")
        ) +
        ylab("Mean damage") +
        xlab("Time")

    (p1/p2)

    ggsave(filename=paste0("hormone_profile_autocorr",autocorr.i,".pdf")
           ,device=cairo_pdf
           ,height=6
           ,width=5)
} # end autocorr.i
