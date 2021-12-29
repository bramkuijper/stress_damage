library("ggplot2")
library("patchwork")
library("tidyverse")

type="png"

#data.repair <- read_csv(file="summary_repair2.csv")
data.repair <- read_csv(file="summary.csv")

data.repair <- filter(data.repair, repair <= 10)

# calculate risk
data.repair[,"risk"] <- with(data.repair
                             ,round(pArrive/(pLeave+pArrive)
                                    ,digits=2))
                             
# calculate autocorrelation
data.repair[,"autocorrelation"] <- with(data.repair
                                     ,round(1.0 - pArrive - pLeave
                                       ,digits=2))

# autocorrelation values used in this plot
autocorr.u <- c(0,0.3,0.9)
risk.u <- c(0.02,0.05, 0.1)

data.repair.f <- filter(data.repair
                        ,autocorrelation %in% autocorr.u & 
                            risk %in% risk.u
                        )

# make the risk as factor so that we can plot differently
# colored lines in each plot highlighting the level of risk
data.repair.f <- mutate(data.repair.f, 
                      risk_cat=as_factor(data.repair.f$risk))

data.repair.f$autocorr.text <- with(
  data.repair.f
  ,paste("Autocorrelation: ",autocorrelation)
)

data.repair.f$max_damage <- data.repair.f$max_damage / 1000

# labels for the first row
labels.row.1 <- data.frame(
        label=LETTERS[1:length(autocorr.u)]
        ,autocorr.text=factor(sort(unique(data.repair.f$autocorr.text)))
        ,x=rep(0.1,times=length(autocorr.u))
        ,y=rep(0.2,times=length(autocorr.u))
        )

labels.row.2 <- data.frame(
        label=LETTERS[(length(autocorr.u)+1):(2*length(autocorr.u))]
        ,autocorr.text=factor(sort(unique(data.repair.f$autocorr.text)))
        ,x=rep(0.1,times=length(autocorr.u))
        ,y=rep(0.2,times=length(autocorr.u))
        )

labels.row.3 <- data.frame(
        label=LETTERS[(2*length(autocorr.u)+1):(3*length(autocorr.u))]
        ,autocorr.text=factor(sort(unique(data.repair.f$autocorr.text)))
        ,x=rep(1,times=length(autocorr.u))
        ,y=rep(0.75,times=length(autocorr.u))
        )

p_baseline <- ggplot(mapping=aes(x=repair, y=h_base_a)
             ,data=data.repair.f) +
  geom_line(mapping=aes(colour=risk_cat)) +
  facet_grid(. ~autocorr.text) +
  theme_classic() +
  geom_text(
          data = labels.row.1
          ,mapping = aes(x=x, y=y, label=label)
          ) +
  xlab("") +
        theme(
            strip.background = element_rect(
                color="transparent"
            )
        ,panel.spacing = unit(1,"lines")
        ) +
  ylab("Baseline hormone level") +
  labs(colour="Risk")

# plot the maximum hormone level for 
# various amounts of repair
p_peak <- ggplot(mapping=aes(x=repair, y=h_max_a)
       ,data=data.repair.f) +
  geom_line(mapping=aes(colour=risk_cat)) +
  facet_grid(. ~autocorr.text) +
  theme_classic() +
  geom_text(
          data = labels.row.2
          ,mapping = aes(x=x, y=y, label = label)
          ) +
  xlab("") +
    theme(
        strip.text.x = element_blank()
        ,legend.position="none"
        ,panel.spacing = unit(1,"lines")
    ) +
  ylab("Peak hormone level") +
  labs(colour="Risk")

# plot the maximum hormone level for 
# various amounts of repair
p_damage <- ggplot(mapping=aes(x=repair, y=max_damage)
       ,data=data.repair.f) +
  geom_line(mapping=aes(colour=risk_cat)) +
  facet_grid(. ~autocorr.text) +
  theme_classic() +
  geom_text(
          data = labels.row.3
          ,mapping = aes(x=x, y=y, label = label)
          ) +
  xlab("Repair") +
    theme(
        strip.text.x = element_blank()
        ,legend.position="none"
        ,panel.spacing = unit(1,"lines")
    ) +
  ylab("Maximum damage") +
  labs(colour="Risk")
(p_baseline/p_peak/p_damage)

if (type == "png")
{
    ggsave(filename="vary_repair.png")
} else
{
    ggsave(filename="vary_repair.pdf",device=cairo_pdf)
}

