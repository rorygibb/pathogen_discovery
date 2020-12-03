


# ====================== Generate species viral discovery curves for mammals using Shaw, GMPD2, HP3 =====================

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiscovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA")

# associations data
assoc = read.csv("./data/clover_v1/Clover_reconciledassociations_v1_20201120.csv", stringsAsFactors = FALSE) %>%
  #dplyr::filter(Database != "EID2") %>%
  dplyr::filter(YearType != "Nucleotide")

# domestic
load("./data/domestic/domestic_species.R")
domestic = c(domestic, "canis familiaris", "bos frontalis")
assoc$Domestic = ifelse(assoc$Host_Harmonised %in% domestic, "Domestic", "Wild")

# total viral richness by species
tr = assoc %>%
  group_by(Host_Harmonised) %>%
  dplyr::summarise(VRichness = n_distinct(Pathogen_Harmonised))

# unique associations by year
dd = assoc %>%
  dplyr::group_by(Host_Harmonised, Pathogen_Harmonised) %>%
  dplyr::summarise(HostOrder = head(HostOrder, 1), 
                   HostFamily = head(HostFamily, 1), 
                   Database = paste(unique(Database), collapse=", "),
                   NumRecords = length(Year),
                   Year = min(Year, na.rm=TRUE),
                   Domestic = head(Domestic, 1)) %>%
  left_join(tr) %>%
  dplyr::arrange(desc(VRichness), Host_Harmonised, Year) %>%
  #dplyr::filter(Host_Harmonised != "homo sapiens") %>%
  dplyr::filter(Year != Inf) # temporary fix for pathogens with no year

# create cumulative curves of all years
curves = expand.grid(unique(dd$Host_Harmonised), 1930:2016) %>%
  dplyr::rename("Host_Harmonised" = 1, "Year" = 2) %>%
  left_join(dd[ , c("Host_Harmonised", "Year", "Pathogen_Harmonised")]) %>%
  dplyr::group_by(Host_Harmonised, Year) %>%
  dplyr::summarise(Discovered = sum(!is.na(Pathogen_Harmonised)),
                   Virus = paste(unique(Pathogen_Harmonised), collapse=", ")) %>%
  left_join(dd[ !duplicated(dd$Host_Harmonised), c("Host_Harmonised", "HostOrder", "HostFamily", "VRichness", "Domestic") ]) %>%
  dplyr::rename("HostSpecies" = Host_Harmonised) %>%
  dplyr::arrange(desc(VRichness), HostSpecies, Year) %>%
  dplyr::group_by(HostSpecies) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered))
  

# # some initial viz
# ggplot(curves[ curves$Host %in% unique(curves$Host)[1:40], ]) + 
#   geom_line(aes(Year, VirusCumulative, group=Host, col=Host), size=1)
# ggplot(curves[ curves$Host %in% unique(curves$Host)[20:40], ]) + 
#   geom_line(aes(Year, VirusCumulative, group=Host, col=Host), size=1)
# ggplot(curves[ curves$Host %in% unique(curves$Host)[20:40], ]) + 
#   geom_bar(aes(Year, Discovered), stat="identity") + 
#   facet_wrap(~Host)
# 
# # mean viruses discovered per year per species over time
# # interesting taper down from the mid-2000s: reflects shift towards GenBank and delay in publishing time?
# ggplot(curves) + 
#   geom_smooth(aes(Year, Discovered), method="gam", method.args=list(family="poisson"))
# ggplot(curves) + 
#   geom_smooth(aes(Year, Discovered, group=Domestic, col=Domestic), method="gam", method.args=list(family="poisson"))
# ggplot(curves[ curves$Domestic == "Wild", ]) + 
#   #geom_point(aes(Year, Discovered)) +
#   geom_smooth(aes(Year, Discovered), method="gam", method.args=list(family="poisson"))
# 
# # mean cumulative viruses
# ggplot(curves) + 
#   geom_smooth(aes(Year, VirusCumulative), method="gam", method.args=list(family="poisson"))
# ggplot(curves) + 
#   geom_smooth(aes(Year, VirusCumulative, group=Domestic, col=Domestic), method="gam", method.args=list(family="poisson"))
# ggplot(curves[ curves$Domestic == "Wild", ]) + 
#   geom_smooth(aes(Year, VirusCumulative), method="gam", method.args=list(family="poisson"))





# ================ for each species, fit poisson model to discovery rates ===================

# fits Poisson model of discovery rates (novel virus counts per year) for a specified order with specified data
# equivalent to fitting inhomogenous 1D Poisson process but without smudging event times
# for 3 time epochs: 1930 to present, 1960 to present, and 1990 to present
# currently effect of year is linear (i.e. exponential with log-link); could explore SPDE but may be more intractable if we're interested in broad trends

fitDiscoveryRateCurve = function(species, data){
  
  # data with 2012 cutoff
  print(species)
  dx = data[ data$HostSpecies == species & data$Year <= 2012, ]
  
  # 1. fit model from 1930 to present
  dx$yearx = 1:nrow(dx)
  form = Discovered ~ yearx + Intercept
  bru_mod = bru(form, dx, family = "poisson")

  # extract fixed effects
  fx = bru_mod$summary.fixed
  fx$param = row.names(fx)
  names(fx)[ 3:5 ] = c("q0.025", "median", "q0.975")
  fx$model = "1930"
  row.names(fx) = c()
  
  # predict curve
  x4pred = data.frame(yearx = 1:nrow(dx))
  predx_bru1 = predict(bru_mod, x4pred, ~ exp(yearx + Intercept), n.samples=2000)
  resx = left_join(dx, predx_bru1, by="yearx")
  resx$model = "1930"
  
  # 2. fit model from 1960 to present
  dx = dx[ dx$Year >= 1960, ]
  dx$yearx = 1:nrow(dx)
  form = Discovered ~ yearx + Intercept
  bru_mod = bru(form, dx, family = "poisson")
  
  # extract fixed effects
  fy = bru_mod$summary.fixed
  fy$param = row.names(fy)
  names(fy)[ 3:5 ] = c("q0.025", "median", "q0.975")
  fy$model = "1960"
  row.names(fy) = c()
  
  # predict curve
  x4pred = data.frame(yearx = 1:nrow(dx))
  predx_bru1 = predict(bru_mod, x4pred, ~ exp(yearx + Intercept), n.samples=2000)
  resy = left_join(dx, predx_bru1, by="yearx")
  resy$model = "1960"
  
  # 2. fit model from 1990 to present
  dx = dx[ dx$Year >= 1990, ]
  dx$yearx = 1:nrow(dx)
  form = Discovered ~ yearx + Intercept
  bru_mod = bru(form, dx, family = "poisson")
  
  # extract fixed effects
  fz = bru_mod$summary.fixed
  fz$param = row.names(fz)
  names(fz)[ 3:5 ] = c("q0.025", "median", "q0.975")
  fz$model = "1990"
  row.names(fz) = c()
  
  # predict curve
  x4pred = data.frame(yearx = 1:nrow(dx))
  predx_bru1 = predict(bru_mod, x4pred, ~ exp(yearx + Intercept), n.samples=1995)
  resz = left_join(dx, predx_bru1, by="yearx")
  resz$model = "1990"
  
  # create results
  ff = rbind(fx, fy, fz)
  ff$HostSpecies = species
  ff$HostOrder = dx$HostOrder[1]
  ff$Domestic = dx$Domestic[1]
  res = rbind(resx, resy, resz)
  
  return(list(
    fixed = ff,
    pred_curve = res
  ))
}

# 1. run models for all species with viral richness of >8 (n=165)
species_for_model = curves[ !duplicated(curves$HostSpecies) & curves$VRichness >= 8, c("HostSpecies", "HostOrder", "HostFamily", "Domestic", "VRichness")]
#table(species_for_model$HostOrder)
result = lapply(species_for_model$HostSpecies, fitDiscoveryRateCurve, data=curves)

# extract estimates and include total viral richness
fixed = do.call(rbind.data.frame, lapply(result, "[[", 1)) %>%
  left_join(curves[ !duplicated(curves$HostSpecies), c("HostSpecies", "VRichness")])
curve_preds = do.call(rbind.data.frame, lapply(result, "[[", 2)) %>%
  left_join(curves[ !duplicated(curves$HostSpecies), c("HostSpecies", "VRichness")])

write.csv(fixed, "./output/order_models/fixedeffects_byspecies_inlabrupois.csv", row.names=FALSE)
write.csv(curve_preds, "./output/order_models/curves_byspecies_inlabrupois.csv", row.names=FALSE)




# ==================== some viz =========================

curve_preds = curve_preds[ curve_preds$VRichness > ]
curve_preds$HostSpecies = factor(curve_preds$HostSpecies, levels=unique(curve_preds$HostSpecies), ordered=TRUE)

ggplot(curve_preds[ curve_preds$model == 1930 & curve_preds$HostSpecies %in% unique(curve_preds$HostSpecies)[1:60], ]) +
  #geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
  geom_line(aes(Year, median)) +
  geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=HostOrder), alpha=0.25) +
  theme_minimal() +
  facet_wrap(~HostSpecies, scales="free_y") +
  #ggtitle(Hmisc::capitalize(daty$Host[1])) +
  ylab(expression(lambda)) +
  xlab("Year") +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=12),
        axis.text = element_text(size=11))

ggplot(curve_preds[ curve_preds$model == 1960 & curve_preds$HostSpecies %in% unique(curve_preds$HostSpecies)[1:60], ]) +
  #geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
  geom_line(aes(Year, median)) +
  geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=HostOrder), alpha=0.25) +
  theme_minimal() +
  facet_wrap(~HostSpecies, scales="free_y") +
  #ggtitle(Hmisc::capitalize(daty$Host[1])) +
  ylab(expression(lambda)) +
  xlab("Year") +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=12),
        axis.text = element_text(size=11))

# 
# ggplot(curve_preds[ curve_preds$model == 1960, ]) +
#   geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
#   geom_line(aes(Year, median)) +
#   geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=HostOrder), alpha=0.25) +
#   theme_minimal() +
#   facet_wrap(~HostOrder) +
#   #ggtitle(Hmisc::capitalize(daty$Host[1])) +
#   ylab(expression(lambda)) +
#   xlab("Year") +
#   theme(plot.title=element_text(size=14, hjust=0.5),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11))
# 
ggplot(curve_preds[ curve_preds$model == 1990 & curve_preds$HostSpecies %in% unique(curve_preds$HostSpecies)[1:40], ]) +
  #geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
  geom_line(aes(Year, median)) +
  geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=HostOrder), alpha=0.25) +
  theme_minimal() +
  facet_wrap(~HostSpecies, scales="free_y") +
  #ggtitle(Hmisc::capitalize(daty$Host[1])) +
  ylab("Viral discovery rate (viruses/year)") +
  xlab("Year") +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=12),
        axis.text = element_text(size=11))

curve_preds$model2 = paste(curve_preds$model, "-2012", sep="")
p1 = ggplot(curve_preds[ curve_preds$model %in% c(1930, 1960) & curve_preds$HostSpecies %in% unique(curve_preds$HostSpecies)[1:60], ]) +
  geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
  geom_line(aes(Year, median, group=factor(model2))) +
  geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=factor(model2)), alpha=0.3) +
  theme_minimal() +
  facet_wrap(~HostSpecies, scales="free_y") +
  scale_fill_viridis_d( name="Time epoch", begin=0.2, end=0.7) +
  #ggtitle(Hmisc::capitalize(daty$Host[1])) +
  #ylab(expression(lambda)) +
  ylab("Viral discovery rate (viruses/year)") +
  xlab("Year") +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=11),
        legend.title = element_text(size=12), 
        strip.text = element_text(size=13), 
        legend.text = element_text(size=11), 
        axis.text = element_text(size=11))


fixed$model2 = paste(fixed$model, "-2012", sep="")
fixed = fixed[ fixed$median > -10, ]
fixed = fixed[ fixed$q0.025>-10, ]
#ggplot(fixed[ fixed$param %in% c("yearx") & fixed$model != "2000", ]) + 
ggplot(fixed[ fixed$param %in% c("yearx"), ]) + 
  #geom_point(aes(HostOrder, median, group=HostSpecies), alpha=0.1, position=position_dodge(width=0.2)) +
  geom_linerange(aes(HostOrder, ymin=q0.025, ymax=q0.975, group=HostSpecies, col=factor(HostOrder)), position=position_dodge(width=0.8)) + 
  geom_hline(yintercept=0, lty=2) + 
  theme_minimal() +
  facet_wrap(~model2, scales="free_y") +
  scale_color_viridis_d( name="Order", begin=0.2, end=0.7) +
  theme(legend.position = "none",
        plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11, angle=90))

fixed$sig_pos = fixed$q0.025 > 0
fixed$sig_neg = fixed$q0.975 < 0
fixed$non_sig = fixed$sig_pos == FALSE & fixed$sig_neg == FALSE
resx = fixed[ fixed$param == "yearx", ] %>%
  dplyr::group_by(model2, HostOrder) %>%
  dplyr::summarise(Increasing = sum(sig_pos) / length(sig_pos),
                   Flat = sum(non_sig) / length(non_sig),
                   Decreasing = sum(sig_neg) / length(sig_neg))




# is there a negative correlation between trend in recent years and overall richness? (i.e. are species with more viruses decreasing?)
fcor = fixed %>%
  dplyr::filter(param == "yearx" & model == "1990") %>%
  dplyr::filter(median > -10) %>%
  left_join(curves[ !duplicated(curves$HostSpecies), c("HostSpecies", "VRichness")])

# pretty much no relationship between slope and total richness
ggplot(fcor[ fcor$Domestic == "Wild" & fcor$HostSpecies != "homo sapiens", ]) + 
  geom_point(aes(VRichness, median, group=HostSpecies), position=position_dodge(width=0.1), alpha=0.25) + 
  geom_linerange(aes(VRichness, ymin=q0.025, ymax=q0.975, group=HostSpecies), position=position_dodge(width=0.1), alpha=0.25) + 
  geom_hline(yintercept=0, lty=2) + 
  theme_minimal()

cor(fcor$median, fcor$VRichness)







# ======================= How stable are ranks of species and order level viral richness over time? ============================

# 1. Order level comparison of rank of mean viral richness

getOrderVRich = function(year){
  resx = curves %>%
    dplyr::filter(Year <= year) %>%
    dplyr::filter(HostSpecies != "homo sapiens") %>%
    dplyr::group_by(HostSpecies) %>%
    dplyr::summarise(VRich = sum(Discovered),
                     HostOrder = head(HostOrder, 1)) %>%
    dplyr::group_by(HostOrder) %>%
    dplyr::summarise(VRich_Mean = mean(VRich),
                     VRich_SEM = plotrix::std.error(VRich)) %>%
    dplyr::arrange(desc(VRich_Mean)) %>%
    dplyr::mutate(Year = year)
  return(resx)
}
ord_over_time = do.call(rbind.data.frame, lapply(seq(1960, 2015, by=5), getOrderVRich))
oot = reshape2::dcast(ord_over_time[ , c("HostOrder", "Year", "VRich_Mean")], HostOrder ~ Year, value.var="VRich_Mean")

# calculate spearman coefficient over time
corx = function(x){ cor(oot$`2015`, x, method="spearman") }
spearman_o = as.data.frame(apply(oot[ , 2:ncol(oot)], 2, corx))
spearman_o$Year = row.names(spearman_o)
names(spearman_o)[1] = "SpearmanRho"

ggplot(spearman_o) + 
  geom_point(aes(Year, SpearmanRho)) +
  geom_line(aes(Year, SpearmanRho))

# ggplot(disc_beta_time[ disc_beta_time$HostOrder %in% t1$HostOrder[ 1:15], ]) +
#   geom_point(aes(Year, VRich_Mean, group=HostOrder, col=HostOrder), size=1.5, position=position_dodge(width=5)) +
#   geom_line(aes(Year, VRich_Mean, group=HostOrder, col=HostOrder), position=position_dodge(width=5)) +
#   geom_linerange(aes(Year, ymin=VRich_Mean - VRich_SEM, ymax = VRich_Mean + VRich_SEM, group=HostOrder, col=HostOrder), position=position_dodge(width=5))
#   

# =========== at the species level ===============

# diversity estimates in 5 year increments
getSpeciesVRich = function(year){
  resx = curves %>%
    dplyr::filter(Year <= year) %>%
    dplyr::filter(HostSpecies != "homo sapiens") %>%
    dplyr::group_by(HostSpecies) %>%
    dplyr::summarise(VRich = sum(Discovered),
                     HostOrder = head(HostOrder, 1),
                     Domestic = head(Domestic, 1)) %>%
    dplyr::mutate(Year = year)
  return(resx)
}
rich_over_time = do.call(rbind.data.frame, lapply(seq(1960, 2015, by=5), getSpeciesVRich))
rot = rich_over_time[ rich_over_time$Domestic == "Wild", ]

# ggplot(rot[ rot$HostOrder == "Rodentia", ]) + 
#   geom_point(aes(Year, VRich, group=HostSpecies, col=HostSpecies)) + 
#   geom_line(aes(Year, VRich, group=HostSpecies, col=HostSpecies)) + 
#   theme(legend.position = "none")
# 
# oo = "Carnivora"
# px = ggplot(rot[ rot$HostOrder == oo, ]) + 
#   geom_point(aes(Year, VRich, group=HostSpecies, col=HostSpecies)) + 
#   geom_line(aes(Year, VRich, group=HostSpecies, col=HostSpecies)) + 
#   theme_minimal() +
#   theme(legend.position = "none") + 
#   ylab("Viral richness") + 
#   ggtitle(paste(oo, " (", n_distinct(rot$HostSpecies[ rot$HostOrder == oo ]), " species)", sep=""))
# ggsave(px, "./output/figures/Carnivora_spaghetti.png", device="png", )

# over time
rot2 = reshape2::dcast(rich_over_time[ , c("HostSpecies", "Year", "VRich")], HostSpecies ~ Year, value.var="VRich") %>%
  left_join(rich_over_time[ !duplicated(rich_over_time$HostSpecies), c("HostSpecies", "HostOrder") ])

# run for all orders
result = data.frame()
orders = c("Carnivora", "Chiroptera", "Cetartiodactyla", "Lagomorpha", "Perissodactyla", "Primates", "Rodentia", "Diprotodontia", "Didelphimorphia")
for(i in 1:(length(orders)+1)){

  if(i == 1){ 
    datx = rot2
    orderx = "Mammalia"
  } else{
    orderx = orders[i-1]
    datx = rot2[ rot2$HostOrder == orderx, ]
  }
  
  # run 
  corx = function(x){ cor(datx$`2015`, x, method="spearman") }
  resx = as.data.frame(apply(datx[ , 2:(ncol(datx)-1)], 2, corx)) %>%
    rename("rho" = 1) %>%
    dplyr::mutate(Order = orderx)
  resx$Year = seq(1960, 2015, by=5)
  result = rbind(result, resx)
}

#
result$Order = factor(result$Order, levels=c("Mammalia", orders), ordered=TRUE)
p0 = ggplot(result) + 
  geom_point(aes(Year, rho), size=0.75) +
  geom_line(aes(Year, rho), col="darkred") + 
  facet_wrap(~Order, nrow=2) +
  theme_minimal() + 
  geom_vline(xintercept=2015, lty=2) +
  scale_x_reverse() +
  ylab("Spearman rho (rank correlation)")
ggsave(p0, file="./output/figures/Order_viralrichnessranks_temporaldecay.png", device="png", dpi=300, height=4.5, width=10, units="in")


# ggplot(rich_over_time) + geom_abline() + geom_point(aes(`2010`, `2015`), size=2, alpha=0.35, col="skyblue4") + theme_minimal()
# ggplot(rich_over_time) + geom_abline() + geom_point(aes(`2000`, `2015`), size=2, alpha=0.35, col="skyblue4") + theme_minimal()
# ggplot(rich_over_time) + geom_abline() + geom_point(aes(`1990`, `2015`), size=2, alpha=0.35, col="skyblue4") + theme_minimal()
# ggplot(rich_over_time) + geom_abline() + geom_point(aes(`1980`, `2015`), size=2, alpha=0.35, col="skyblue4") + theme_minimal()
# ggplot(rich_over_time) + geom_abline() + geom_point(aes(`1970`, `2015`), size=2, alpha=0.35, col="skyblue4") + theme_minimal()
# ggplot(rich_over_time) + geom_abline() + geom_point(aes(`1960`, `2015`), size=2, alpha=0.35, col="skyblue4") + theme_minimal()
