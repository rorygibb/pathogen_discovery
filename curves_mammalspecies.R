


# ====================== Generate species viral discovery curves for mammals using Shaw, GMPD2, HP3, EID2 =====================

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

# publication effort by year
# this needs some thought and work if it's worth including
pubs = read.csv("./output/host_effort/PubMed_Hosts_PubsByYear_Final.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Note == "") %>%
  dplyr::select(1:3) %>%
  dplyr::rename("HostSpecies" = Host)

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
  
# publications
curves = left_join(curves, pubs) %>%
  dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0))


# ======================== initial =============================

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



# ================ 1. Model of viral discovery rate at species level ===================

# fits Poisson model of discovery rates (novel virus counts per year) for a specified species with specified data
# equivalent to fitting nonhomogenous 1D Poisson process but without smudging event times
# for 3 time epochs: 1930 to present, 1960 to present, and 1990 to present
# currently effect of year is linear (i.e. exponential with log-link); could explore SPDE for gam-like fits but may be more intractable if we're interested in broad trends

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

# write.csv(fixed, "./output/order_models/fixedeffects_byspecies_inlabrupois.csv", row.names=FALSE)
# write.csv(curve_preds, "./output/order_models/curves_byspecies_inlabrupois.csv", row.names=FALSE)




# ======================= 2. Summary statistics of model outputs ============================

# model outputs 
fixed = read.csv("./output/order_models/fixedeffects_byspecies_inlabrupois.csv", stringsAsFactors = FALSE)
curve_preds = read.csv("./output/order_models/curves_byspecies_inlabrupois.csv", stringsAsFactors=FALSE)

# 1. proportion of species with strong evidence of +ve or -ve trend
fixed$pos = fixed$q0.025 > 0
fixed$neg = fixed$q0.975 < 0
fixed$flat = fixed$pos == FALSE & fixed$neg == FALSE

# across all wild species
resx = fixed[ fixed$param == "yearx" & fixed$Domestic == "Wild", ] %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(Increasing = sum(pos) / length(pos),
                   Flat = sum(flat) / length(flat),
                   Decreasing = sum(neg) / length(neg))

# split by order
resy = fixed[ fixed$param == "yearx" & fixed$Domestic == "Wild", ] %>%
  dplyr::group_by(model, HostOrder) %>%
  dplyr::summarise(Increasing = sum(pos) / length(pos),
                   Flat = sum(flat) / length(flat),
                   Decreasing = sum(neg) / length(neg))

# visualisation of fixed effects
fixed$model2 = paste(fixed$model, "-2012", sep="")
fixed = fixed[ fixed$median > -10, ]
fixed = fixed[ fixed$q0.025>-10, ]
fixed = fixed[ fixed$q0.975< 2, ]
#ggplot(fixed[ fixed$param %in% c("yearx") & fixed$model != "2000", ]) +
ggplot(fixed[ fixed$param %in% c("yearx"), ]) +
  geom_linerange(aes(HostOrder, ymin=q0.025, ymax=q0.975, group=HostSpecies, col=factor(HostOrder)), position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0, lty=2) +
  theme_minimal() +
  facet_wrap(~model2, scales="free_y", nrow=3) +
  scale_color_viridis_d( name="Order", begin=0.2, end=0.7) +
  theme(legend.position = "none",
        plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11, angle=90))

# is there evidence that species with more known viruses are more likely to be decreasing trends?
fcor = fixed %>%
  dplyr::filter(param == "yearx" & model == "1990") %>%
  dplyr::filter(median > -10) %>%
  left_join(curves[ !duplicated(curves$HostSpecies), c("HostSpecies", "VRichness")])

# pretty much no relationship between slope and total richness
ggplot(fcor[ fcor$Domestic == "Wild" & fcor$HostSpecies != "homo sapiens" & fcor$flat == FALSE, ]) +
  geom_point(aes(VRichness, median, group=HostSpecies), position=position_dodge(width=0.1), alpha=0.25) +
  geom_linerange(aes(VRichness, ymin=q0.025, ymax=q0.975, group=HostSpecies), position=position_dodge(width=0.1), alpha=0.25) +
  geom_hline(yintercept=0, lty=2) +
  theme_minimal()





# ====================== 3. Analysis of stability of relative viral richness estimates over time ===========================

# 1. Order-level comparison of mean viral richness across species

# mean viral richness across species within each order
getOrderVRich = function(year){
  resx = curves %>%
    dplyr::filter(Year <= year) %>%
    dplyr::filter(HostSpecies != "homo sapiens") %>%
    dplyr::filter(Domestic == "Wild") %>%
    dplyr::group_by(HostSpecies) %>%
    dplyr::summarise(VRich = sum(Discovered),
                     HostOrder = head(HostOrder, 1)) %>%
    dplyr::group_by(HostOrder) %>%
    dplyr::summarise(VRich_Mean = mean(VRich),
                     VRich_Median = median(VRich),
                     VRich_SEM = plotrix::std.error(VRich)) %>%
    dplyr::arrange(desc(VRich_Mean)) %>%
    dplyr::mutate(Year = year)
  return(resx)
}
ord_over_time = do.call(rbind.data.frame, lapply(seq(1960, 2015, by=5), getOrderVRich))
oot = reshape2::dcast(ord_over_time[ , c("HostOrder", "Year", "VRich_Mean")], HostOrder ~ Year, value.var="VRich_Mean")

# calculate spearman rank correlation between 2015 and all other years (examine temporal decay)
cor_rho = function(x){ return( as.vector(cor.test(oot$`2015`, x, method="spearman")$estimate) ) }
stab1 = as.data.frame(apply(oot[ , 2:ncol(oot)], 2, cor_rho)) %>%
  rename("rho" = 1) %>%
  dplyr::mutate(year = seq(1960, 2015, by=5),
                model = "Order")


# 2. Family-level comparison of mean viral richness across species

# mean viral richness across species within each order
getFamilyVRich = function(year){
  resx = curves %>%
    dplyr::filter(Year <= year) %>%
    dplyr::filter(HostSpecies != "homo sapiens") %>%
    dplyr::filter(Domestic == "Wild") %>%
    dplyr::group_by(HostSpecies) %>%
    dplyr::summarise(VRich = sum(Discovered),
                     HostFamily = head(HostFamily, 1),
                     HostOrder = head(HostOrder, 1)) %>%
    dplyr::group_by(HostFamily) %>%
    dplyr::summarise(VRich_Mean = mean(VRich),
                     VRich_Median = median(VRich),
                     VRich_SEM = plotrix::std.error(VRich),
                     HostOrder = head(HostOrder, 1)) %>%
    dplyr::arrange(desc(VRich_Mean)) %>%
    dplyr::mutate(Year = year)
  return(resx)
}
fam_over_time = do.call(rbind.data.frame, lapply(seq(1960, 2015, by=5), getFamilyVRich))
fot = reshape2::dcast(fam_over_time[ , c("HostFamily", "Year", "VRich_Mean")], HostFamily ~ Year, value.var="VRich_Mean")

# calculate spearman rank correlation between 2015 and all other years (examine temporal decay)
cor_rho = function(x){ return( as.vector(cor.test(fot$`2015`, x, method="spearman")$estimate) ) }
stab2 = as.data.frame(apply(fot[ , 2:ncol(fot)], 2, cor_rho)) %>%
  rename("rho" = 1) %>%
  dplyr::mutate(year = seq(1960, 2015, by=5),
                model = "Family")

ggplot(rbind(stab1, stab2)) +
  geom_line(aes(year, rho, col=model)) +
  geom_point(aes(year, rho, col=model))


# 3. Species-level comparison of ranking of viral richness across years

# diversity estimates in 5 year increments
getSpeciesVRich = function(year){
  resx = curves %>%
    dplyr::filter(Year <= year) %>%
    dplyr::filter(HostSpecies != "homo sapiens") %>%
    dplyr::filter(Domestic == "Wild") %>%
    dplyr::group_by(HostSpecies) %>%
    dplyr::summarise(VRich = sum(Discovered),
                     HostOrder = head(HostOrder, 1),
                     Domestic = head(Domestic, 1),
                     Pubs = sum(NumPubs, na.rm=TRUE),
                     VRich_2015 = head(VRichness, 1)) %>%
    dplyr::mutate(Year = year)
  return(resx)
}
rich_over_time = do.call(rbind.data.frame, lapply(seq(1960, 2015, by=5), getSpeciesVRich))

# over time
rot2 = reshape2::dcast(rich_over_time[ , c("HostSpecies", "Year", "VRich")], HostSpecies ~ Year, value.var="VRich") %>%
  left_join(rich_over_time[ !duplicated(rich_over_time$HostSpecies), c("HostSpecies", "HostOrder", "VRich_2015") ])

# run for all orders
result = data.frame()
orders = c("Carnivora", "Chiroptera", "Cetartiodactyla", "Lagomorpha", "Perissodactyla", "Primates", "Rodentia", "Diprotodontia", "Didelphimorphia")
for(i in 1:(length(orders)+2)){

  if(i == 1){ 
    datx = rot2
    orderx = "Species"
  } else if(i == 2){
      datx = rot2[ rot2$VRich_2015 >= 5, ]
      orderx = "Species (>5)"
    } else{
      orderx = orders[i-2]
      datx = rot2[ rot2$HostOrder == orderx, ]
    }
  
  # run 
  cor_rho = function(x){ return( as.vector(cor.test(datx$`2015`, x, method="spearman")$estimate) ) }
  resx = as.data.frame(apply(datx[ , 2:(ncol(datx)-2)], 2, cor_rho)) %>%
    rename("rho" = 1) %>%
    dplyr::mutate(model = orderx)
  resx$year = seq(1960, 2015, by=5)
  result = rbind(result, resx)
}

# stab3 
stab3 = result


# ==================== combine all results ======================

stab_agg = do.call(rbind.data.frame, list(stab1, stab2, stab3[ grep("Species", stab3$model), ]))
stab_orders = stab3[ !grepl("Species", stab3$model), ]

# plots
ss = stab_agg[ stab_agg$model != "Species (>5)", ]
ss$model = factor(ss$model, levels=c("Order", "Family", "Species"), ordered=TRUE)
p1 = ggplot(ss) + 
  geom_line(aes(year, rho, col=model), size=1) + 
  ylim(0, 1) + 
  scale_color_viridis_d(begin=0, end=0.8) +
  theme_minimal() + 
  geom_vline(xintercept=2015, lty=2) + 
  ylab("Spearman rho (rank correlation)") + 
  xlab("Year") + 
  theme(legend.title = element_blank(), 
        axis.title = element_text(size=13),
        axis.text = element_text(size=11),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12),
        legend.position = c(0.2, 0.9))
p2 = ggplot(stab_orders) + 
  geom_line(aes(year, rho, col=model), size=1) + 
  ylim(0, 1) + 
  scale_color_viridis_d(begin=0, end=0.8) +
  facet_wrap(~model) + 
  theme_minimal() + 
  theme(legend.position="none",
        strip.text = element_text(size=12),
        axis.title.x = element_text(size=13),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=11),
        axis.title.y = element_blank()) + 
  geom_vline(xintercept=2015, lty=2) + 
  ylab("Spearman rho (rank correlation)") + 
  xlab("Year")

p_comb = gridExtra::grid.arrange(p1, p2, nrow=1)
ggsave(p_comb, file="./output/figures/TemporalDecayFigure.png", device="png", dpi=600, units="in", width=12, height=6)
  