library("ggplot2")
library("patchwork")
library("tidyverse")

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

# make the risk as factor so that we can plot differently
# colored lines in each plot highlighting the level of risk
data.repair <- mutate(data.repair, 
                      risk_cat=as_factor(data.repair$risk))


data.repair.f <- filter(data.repair
                        ,autocorrelation %in% c(0,0.3,0.7) & 
                            risk %in% c(0.05, 0.1, 0.2, 0.3)
                        )

data.repair.f$autocorr.text <- with(
  data.repair.f
  ,paste("Autocorrelation: ",autocorrelation)
)

# plot the maximum hormone level for 
# various amounts of repair
p1 <- ggplot(mapping=aes(x=repair, y=h_max_a)
       ,data=data.repair.f) +
  geom_line(mapping=aes(colour=risk_cat)) +
  facet_grid(. ~autocorr.text) +
  theme_classic() +
  xlab("Repair") +
  ylab("Hormone level at \u03c4 = 1") +
  labs(colour="Risk")

ggsave(filename="vary_repair_hormone_max.pdf", device=cairo_pdf, height=3, width=7)


p2 <- ggplot(mapping=aes(x=repair, y=h_base_a)
             ,data=data.repair.f) +
  geom_line(mapping=aes(colour=risk_cat)) +
  facet_grid(. ~autocorr.text) +
  theme_classic() +
  xlab("Repair") +
  ylab("Baseline hormone level at \u03c4 = 40") +
  labs(colour="Risk")

ggsave(filename="vary_repair_hormone_baseline.pdf", device=cairo_pdf)
