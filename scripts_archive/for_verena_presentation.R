
# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiscovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA")

# domestic species to label
domestic = read.csv("./data/clover/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)

# associations data and no humans
clover = read.csv("./data/clover/Clover_v1.0_NBCIreconciled_20201218.csv", stringsAsFactors = FALSE) %>%
  #dplyr::filter(Host != "Homo sapiens") %>%
  dplyr::mutate(Domestic = ifelse(Host %in% domestic$Host, TRUE, FALSE)) %>%
  dplyr::filter(DetectionMethod != "Not specified") %>%
  # dplyr::filter(!DetectionMethod %in% c("Not specified", "Antibodies")) %>%
  #dplyr::filter(!DetectionMethod == "Isolation/Observation") %>%
  dplyr::filter(!is.na(Year)) %>%
  dplyr::filter(Year <= 2010)

# viral discovery for specified species/order
dd = clover[ clover$HostOrder == "Cetartiodactyla", ] %>%
  dplyr::group_by(Virus) %>%
  dplyr::summarize(Year = min(Year)) %>%
  dplyr::arrange(Year) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(Discovered = n_distinct(Virus)) %>%
  dplyr::mutate(Cumulative = cumsum(Discovered),
                Cumulative_Nminus1 = Cumulative - Discovered)

p0 = ggplot(dd) + 
  geom_line(aes(Year, Cumulative), col="darkred", size=1) + 
  theme_classic() +
  xlab("Year") + ylab("Cumulative num. virus species") +
  scale_x_continuous(breaks=seq(1930, 2010, by=20), labels=seq(1930, 2010, by=20)) +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12))

ggsave(p0, file="./output/figures/Cetartiodactyla_trends.png", device="png", units="in", width=4.5, height=4, dpi=600, scale=0.8)


# curves
curve_preds = read.csv("./output/order_models/curves_allspecies_byorder_inlabrupois_20200106.csv")
curve_preds = curve_preds[ curve_preds$HostOrder == "Cetartiodactyla", ]
curve_preds$model2 = paste(curve_preds$model, "-2010", sep="")
p1 = ggplot(curve_preds[ curve_preds$model %in% c(1930, 1960, 1990),  ]) +
  geom_point(aes(Year, Discovered), size=1.5, col="grey60", alpha=0.8) +
  geom_line(aes(Year, median, group=factor(model2))) +
  geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=factor(model2)), alpha=0.3) +
  theme_classic() +
  scale_fill_viridis_d( name="Time epoch", begin=0.2, end=0.7) +
  #ggtitle(Hmisc::capitalize(daty$Host[1])) +
  #ylab(expression(lambda)) +
  ylab("Viral discovery rate (viruses/year)") +
  xlab("Year") +
  ylim(0, 12) +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=11),
        legend.title = element_text(size=12), 
        strip.text = element_text(size=13), 
        legend.text = element_text(size=11), 
        axis.text = element_text(size=11))

ggsave(p1, file="./output/figures/Cetartiodactyla_discovery_points_3_leg.png", device="png", units="in", width=6, height=4, dpi=600, scale=0.8)


px = gridExtra::grid.arrange()


cl = clover %>%
  group_by(HostOrder) %>%
  dplyr::summarise(vrich = n_distinct(Virus))

fixed = read.csv("./output/order_models/fixedeffects_allspecies_byorder_inlabrupois_20200106.csv")
fixed = left_join(fixed, cl)
fixed = fixed[ order(fixed$vrich, decreasing = TRUE), ]
fixed = fixed[! fixed$HostOrder %in% c("Didelphimorphia", "Diprotodontia", "Eulipotyphla"), ]
fixed$HostOrder = factor(fixed$HostOrder, levels=unique(fixed$HostOrder), ordered=TRUE)
fixed$model2 = paste(fixed$model, "-2010", sep="")
#ggplot(fixed[ fixed$param %in% c("yearx") & fixed$model != "2000", ]) + 
p2 = ggplot(fixed[ fixed$param %in% c("yearx"), ]) + 
  geom_point(aes(HostOrder, median, col=model2, group=model2), position=position_dodge(width=0.5)) +
  geom_linerange(aes(HostOrder, ymin=q0.025, ymax=q0.975, col=model2, group=model2), position=position_dodge(width=0.5)) + 
  geom_hline(yintercept=0, lty=2) +
  theme_minimal() +
  ylab(expression(beta[year])) +
  scale_color_viridis_d( name="Time epoch", begin=0.2, end=0.7) +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=11), 
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11, angle=90))
ggsave(p2, file="./output/figures/Order_Betaestimates_20201203_forpresentation.png", device="png", dpi=300, width=8, height=5, units="in", scale=0.85)



cl1 = clover[ clover$HostOrder == "Cetartiodactyla" & clover$Domestic == FALSE, ] %>%
  group_by(Host) %>%
  dplyr::summarise(vrich = n_distinct(Virus)) %>%
  dplyr::arrange(desc(vrich))

cl2 = clover[ clover$HostOrder == "Cetartiodactyla" &  clover$Domestic == FALSE & clover$Year <= 1995, ] %>%
  group_by(Host) %>%
  dplyr::summarise(vrich = n_distinct(Virus)) %>%
  dplyr::arrange(desc(vrich))

cl3 = clover[clover$HostOrder == "Cetartiodactyla" &  clover$Domestic == FALSE & clover$Year <= 1980, ] %>%
  group_by(Host) %>%
  dplyr::summarise(vrich = n_distinct(Virus)) %>%
  dplyr::arrange(desc(vrich))


table = data.frame(
  `Vdiv 1980` = cl3$Host[1:20], `Vdiv 1995`=cl2$Host[1:20], `Vdiv 2010` = cl1$Host[1:20]
)
library(ggpubr)
tt = ggtexttable(table)
ggsave(tt, file="./output/figures/Cetartio_table_rank.png", device="png", dpi=300, width=8, height=9, units="in")
