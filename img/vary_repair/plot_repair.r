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

# risk values used
risk.u <- c(0.02,0.05, 0.1)

# decrement in both baseline and peak hormone level
# to discern lines in the plot
risk.decrement.hormone <- c(0,-0.005,0.005)

# slightly change hormone levels to make it graphically clearer
for (risk.idx in 1:length(risk.u))
{
    risk.u.i <- risk.u[[risk.idx]]
    data.repair[data.repair$risk == risk.u.i,"h_max_a"] <-
        data.repair[data.repair$risk == risk.u.i,"h_max_a"] + risk.decrement.hormone[[risk.idx]]

    data.repair[data.repair$risk == risk.u.i,"h_base_a"] <-
        data.repair[data.repair$risk == risk.u.i,"h_base_a"] + risk.decrement.hormone[[risk.idx]]
}

data.repair.f <- filter(data.repair
                        ,autocorrelation %in% autocorr.u & 
                            risk %in% risk.u)

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
        ,x=rep(0.5,times=length(autocorr.u))
        ,y=rep(0.01,times=length(autocorr.u))
        )

labels.row.2 <- data.frame(
        label=LETTERS[(length(autocorr.u)+1):(2*length(autocorr.u))]
        ,autocorr.text=factor(sort(unique(data.repair.f$autocorr.text)))
        ,x=rep(1,times=length(autocorr.u))
        ,y=rep(0.75,times=length(autocorr.u))
        )

labels.row.3 <- data.frame(
        label=LETTERS[(2*length(autocorr.u)+1):(3*length(autocorr.u))]
        ,autocorr.text=factor(sort(unique(data.repair.f$autocorr.text)))
        ,x=rep(1,times=length(autocorr.u))
        ,y=rep(0.75,times=length(autocorr.u))
        )


# gather the data so that baseline and peak hormone levels
# are summarized in the same type column
data.repair.recast <- data.repair.f %>%
       select(repair,autocorr.text, risk_cat, h_max_a, h_base_a) %>%
       gather(key="hormone_type", value="hormone_level", h_max_a, h_base_a)

# make text values of the different hormone types for a readable legend
data.repair.recast <- data.repair.recast %>%
        mutate(
        hormone_type_txt=with(data.repair.recast, 
                ifelse(hormone_type == "h_max_a","Peak","Baseline")
                )
        )

# make the plot of the hormone level, baseline + peak
p_baseline <- ggplot(mapping=aes(x=repair, y=hormone_level)
             ,data=data.repair.recast) +
  geom_line(mapping=aes(colour=risk_cat,linetype=hormone_type_txt),size=0.4) +
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
  ylab("Hormone level") +
  ylim(0,1.05) +
  labs(colour="Risk",linetype="Hormone")

  ## plot the maximum hormone level for 
  ## various amounts of repair
  #p_peak <- ggplot(mapping=aes(x=repair, y=h_max_a)
  #       ,data=data.repair.f) +
  #  geom_line(mapping=aes(colour=risk_cat)) +
  #  facet_grid(. ~autocorr.text) +
  #  theme_classic() +
  #  geom_text(
  #          data = labels.row.2
  #          ,mapping = aes(x=x, y=y, label = label)
  #          ) +
  #  xlab("") +
  #    theme(
  #        strip.text.x = element_blank()
  #        ,legend.position="none"
  #        ,panel.spacing = unit(1,"lines")
  #    ) +
  #  ylab("Peak hormone level") +
  #  labs(colour="Risk")

# plot the maximum hormone level for 
# various amounts of repair
p_damage <- ggplot(mapping=aes(x=repair, y=max_damage)
       ,data=data.repair.f) +
  geom_line(mapping=aes(colour=risk_cat)) +
  facet_grid(. ~autocorr.text) +
  theme_classic() +
  geom_text(
          data = labels.row.2
          ,mapping = aes(x=x, y=y, label = label)
          ) +
  xlab("Repair") +
    theme(
        strip.text.x = element_blank()
        #      ,legend.position="none"
        ,panel.spacing = unit(1,"lines")
    ) +
  ylab("Maximum damage") +
  labs(colour="Risk")

(p_baseline/p_damage)

if (type == "png")
{
    ggsave(filename="vary_repair.png")
} else
{
    ggsave(filename="vary_repair.pdf",device=cairo_pdf)
}

