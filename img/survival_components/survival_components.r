library("tidyverse")

type = "pdf"

data.death.comps <- read_csv(file="../vary_repair/summary.csv")


# calculate risk
data.death.comps[,"risk"] <- with(data.death.comps
                             ,round(pArrive/(pLeave+pArrive)
                                    ,digits=2))
                             
# calculate autocorrelation
data.death.comps[,"autocorrelation"] <- with(data.death.comps
                                     ,round(1.0 - pArrive - pLeave
                                       ,digits=2))


# autocorrelation values used in this plot
autocorr.u <- c(0,0.3,0.9)
risk.u <- c(0.02,0.05,0.1)

data.death.comps.f <- filter(data.death.comps
                        ,autocorrelation %in% autocorr.u & 
                            risk %in% risk.u
                        )

# make the risk as factor so that we can plot differently
# colored lines in each plot highlighting the level of risk
data.death.comps.f <- mutate(data.death.comps.f, 
                      risk_cat=as_factor(data.death.comps.f$risk))

data.death.comps.f$autocorr.text <- with(
  data.death.comps.f
  ,paste("Autocorrelation: ",autocorrelation)
)

# labels for the first row
labels.row.1 <- data.frame(
        label=LETTERS[1:length(autocorr.u)]
        ,autocorr.text=factor(sort(unique(data.death.comps.f$autocorr.text)))
        ,x=rep(0.1,times=length(autocorr.u))
        ,y=rep(0.2,times=length(autocorr.u))
        )

p_death <- ggplot(mapping=aes(x=repair, y=damageDeaths)
             ,data=data.death.comps.f) +
  geom_line(mapping=aes(colour=risk_cat)) +
  facet_grid(. ~autocorr.text) +
  theme_classic() +
  geom_text(
          data = labels.row.1
          ,mapping = aes(x=x, y=y, label=label)
          ) +
  xlab("Repair") +
    theme(
        strip.background = element_rect(
            color="transparent"
        )
    ,panel.spacing = unit(1,"lines")
    ) +
  ylab("Death due to stressor") +
  labs(colour="Risk")

plot.name = "mortality_components"

if (type == "png")
{
    ggsave(filename=paste0(plot.name,".png"))
} else
{
    ggsave(filename=paste0(plot.name,".pdf")
            ,device=cairo_pdf)
}


