
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon", "mgcv", "vroom")

library(sf)
library(dplyr)

mam1 = sf::st_read("C:/Users/roryj/Documents/PhD/202008_discovery/data/iucn_range/mammals_freshwater/MAMMALS_FRESHWATER.shp")
mam2 = sf::st_read("C:/Users/roryj/Documents/PhD/202008_discovery/data/iucn_range/mammals_terrestrial/MAMMALS_TERRESTRIAL_ONLY.shp")
mam3 = sf::st_read("C:/Users/roryj/Documents/PhD/202008_discovery/data/iucn_range/mammals_marine/MAMMALS_MARINE_ONLY.shp")
mam4 = sf::st_read("C:/Users/roryj/Documents/PhD/202008_discovery/data/iucn_range/mammals_terrmar/MAMMALS_MARINE_AND_TERRESTRIAL/MAMMALS_MARINE_AND_TERRESTRIAL.shp")

mam1 = mam1 %>% st_drop_geometry() %>% dplyr::select(binomial, order_, family, genus) %>% distinct()
mam2 = mam2 %>% st_drop_geometry() %>% dplyr::select(binomial, order_, family, genus) %>% distinct()
mam3 = mam3 %>% st_drop_geometry() %>% dplyr::select(binomial, order_, family, genus) %>% distinct()
mam4 = mam4 %>% st_drop_geometry() %>% dplyr::select(binomial, order_, family, genus) %>% distinct()

# recode to cetacea
mam1$order_[ mam1$family %in% c("BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE", 
                                "ESCHRICHTIIDAE", "KOGIIDAE", "PHOCOENIDAE", "PONTOPORIIDAE", "PHYSETERIDAE", "INIIDAE", "LIPOTIDAE", "PLATANISTIDAE") ] = "CETACEA"
mam3$order_[ mam3$family %in% c("BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE", 
                                "ESCHRICHTIIDAE", "KOGIIDAE", "PHOCOENIDAE", "PONTOPORIIDAE", "PHYSETERIDAE", "INIIDAE", "LIPOTIDAE", "PLATANISTIDAE") ] = "CETACEA"

# order level species richness
mam_o = rbind(mam1, mam2) %>%
  rbind(mam3) %>%
  rbind(mam4) %>%
  dplyr::group_by(order_) %>%
  dplyr::summarise(SR = n_distinct(binomial)) %>%
  dplyr::rename("Order"=order_) %>%
  dplyr::mutate(Order = Hmisc::capitalize(tolower(Order)),
                Order = replace(Order, Order == "Cetartiodactyla", "Artiodactyla"))

# family level species richness
mam_f = rbind(mam1, mam2) %>%
  rbind(mam3) %>%
  rbind(mam4) %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(SR = n_distinct(binomial),
                   Order = head(order_, 1)) %>%
  dplyr::rename("Family"=family) %>%
  dplyr::mutate(Family = Hmisc::capitalize(tolower(Family)),
                Order = Hmisc::capitalize(tolower(Order)),
                Order = replace(Order, Order == "Cetartiodactyla", "Artiodactyla"))




# ---------------- build VIRION database ----------------

ictv_flag = "ictvpredict"

# domestic species to label
domestic = read.csv("./data/clovert/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)

# additional laboratory species to exlude
lab = c("macaca mulatta", "macaca fasicularis")

# publiations
pubs = read.csv("./output/host_effort/PubMed_HostsEffort_PerYear_VirusRelated_VIRION.csv") %>%
  dplyr::filter(Note != "Lookup error") %>%
  dplyr::filter(!Host %in% c(domestic$Host, lab)) %>% 
  dplyr::mutate(Host = Hmisc::capitalize(Host))

# VIRION
#virion = vroom::vroom("./data/virion/pre_push/Virion.csv.gz")
virion = vroom::vroom("./data/virion/Virion.csv.gz")

# build virion and add flag for whether resolved to PREDICT internal tax
vir = virion %>%
  dplyr::filter(HostClass == "mammalia" & Host != "homo sapiens") %>%
  dplyr::filter(!(Host %in% lab)) %>%
  dplyr::filter(Hmisc::capitalize(Host) %in% pubs$Host) %>%
  dplyr::mutate(Domestic = ifelse(Host %in% tolower(domestic$Host), TRUE, FALSE)) %>%
  dplyr::filter(!DetectionMethod %in% c("kmer", "Not specified", "Antibodies")) %>%
  dplyr::filter(!is.na(ReleaseYear) | !is.na(PublicationYear) | !is.na(CollectionYear)) %>%
  dplyr::mutate(PredictFlag = str_detect(Virus, "predict\\_"),
                CLOVERflag = ifelse(Database %in% c("EID2", "Shaw", "GMPD2", "HP3"), TRUE , FALSE)) %>%
  dplyr::filter(ICTVRatified | PredictFlag)

# set years for separate subsets of the data
vp = vir %>%
  dplyr::filter(Database == "PREDICT") %>%
  dplyr::mutate(Year = CollectionYear)
vg = vir %>%
  dplyr::filter(Database == "GenBank") %>%
  dplyr::mutate(Year = ReleaseYear)
vc = vir %>%
  dplyr::filter(CLOVERflag) %>%
  dplyr::mutate(Year = PublicationYear)
vc$Year[ is.na(vc$Year) ] = vc$ReleaseYear[ is.na(vc$Year) ]

# combine
vir = do.call(rbind.data.frame, list(vp, vg, vc))

# subset to keep only years pre endyear
endyear = 2020
vir = vir[ vir$Year <= endyear, ]

# virus or family
#vir$Virus = vir$VirusGenus




# -------------- 1. Order level -----------------

# total viral richness by order
tr = vir %>%
  group_by(HostOrder) %>%
  dplyr::summarise(VRichness = n_distinct(Virus))

# 1. all species including both wild and domestic

# unique associations by order and year
v1 = vir %>%
  dplyr::filter(Domestic == FALSE) %>%
  dplyr::group_by(HostOrder, Virus) %>%
  dplyr::summarise(Database = paste(unique(Database), collapse=", "),
                   HostSpecies = paste(unique(Host), collapse=", "),
                   NumRecords = length(Year),
                   YearEarliest = min(Year, na.rm=TRUE),
                   YearLatest = max(Year, na.rm=TRUE)) %>%
  left_join(tr) %>%
  dplyr::arrange(desc(VRichness), HostOrder, YearEarliest)

# create cumulative curves of all years
curveso = expand.grid(unique(v1$HostOrder), 1930:endyear) %>%
  dplyr::rename("HostOrder" = 1, "YearEarliest" = 2) %>%
  left_join(v1[ , c("HostOrder", "YearEarliest", "Virus")]) %>%
  dplyr::group_by(HostOrder, YearEarliest) %>%
  dplyr::summarise(Discovered = sum(!is.na(Virus)),
                   Virus = paste(unique(Virus), collapse=", ")) %>%
  left_join(v1[ !duplicated(v1$HostOrder), c("HostOrder", "VRichness") ]) %>%
  dplyr::arrange(desc(VRichness), HostOrder, YearEarliest) %>%
  dplyr::group_by(HostOrder) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered)) %>%
  dplyr::rename("Year" = YearEarliest)

# calculate cum pubs at each step
pc = expand.grid(unique(pubs$Host), 1930:endyear) %>%
  dplyr::rename("Host" = 1, "Year" = 2) %>%
  left_join(pubs[ , c("Host", "Year", "NumPubs")]) %>%
  dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0)) %>%
  dplyr::left_join(vir %>% dplyr::select(Host, HostOrder) %>% dplyr::mutate(Host = Hmisc::capitalize(Host)) %>% distinct()) %>%
  dplyr::mutate(HostOrder = Hmisc::capitalize(HostOrder)) %>%
  dplyr::group_by(HostOrder, Year) %>%
  dplyr::summarise(NumPubs = sum(NumPubs)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(HostOrder, Year) %>%
  dplyr::group_by(HostOrder) %>%
  dplyr::mutate(CumPubs = cumsum(NumPubs))

# combine with virus curves
vir_pubs = left_join(pc, curveso %>% dplyr::select(HostOrder, Year, VirusCumulative) %>% dplyr::mutate(HostOrder = Hmisc::capitalize(HostOrder)))


# ------------------ subset to 3 time epochs (1990, 2005, 2018) ----------------------

# data from three time epochs
vir_pubs = vir_pubs %>%
  dplyr::filter(Year %in% c(1980, 1990, 2000, 2010, 2020)) %>%
  dplyr::mutate(logPubs = log(CumPubs+1))

# combine with order richness data
vir_pubs = vir_pubs %>%
  dplyr::left_join(mam_o %>% dplyr::mutate(Order = replace(Order, Order=="Cetartiodactyla", "Artiodactyla")), by=c("HostOrder"="Order"))




# ------------ glm for each year --------------

# negbinom
m1.1 = MASS::glm.nb(VirusCumulative ~ logPubs + log(SR), data=vir_pubs %>% dplyr::filter(Year == 1990))
m1.2 = MASS::glm.nb(VirusCumulative ~ logPubs + log(SR), data=vir_pubs %>% dplyr::filter(Year == 2000))
m1.3 = MASS::glm.nb(VirusCumulative ~ logPubs + log(SR), data=vir_pubs %>% dplyr::filter(Year == 2010))
m1.4 = MASS::glm.nb(VirusCumulative ~ logPubs  + log(SR), data=vir_pubs %>% dplyr::filter(Year == 2020))

# pois
# m1.0 = glm(VirusCumulative ~ logPubs + log(SR), family="poisson", data=vir_pubs %>% dplyr::filter(Year == 1980))
# m1.1 = glm(VirusCumulative ~ logPubs + log(SR), family="poisson", data=vir_pubs %>% dplyr::filter(Year == 1990))
# m1.2 = glm(VirusCumulative ~ logPubs + log(SR), family="poisson", data=vir_pubs %>% dplyr::filter(Year == 2000))
# m1.3 = glm(VirusCumulative ~ logPubs + log(SR), family="poisson", data=vir_pubs %>% dplyr::filter(Year == 2010))
# m1.4 = glm(VirusCumulative ~ logPubs + log(SR), family="poisson", data=vir_pubs %>% dplyr::filter(Year == 2020))

getCoef = function(m, year){
  x = as.data.frame(coef(summary(m)))
  x$param = row.names(x)
  names(x) = c("estimate", "se", "z", "p", "param")
  row.names(x) = c()
  x$upper = x$estimate + (1.96 * x$se)
  x$lower = x$estimate - (1.96 * x$se)
  x$year = year
  return(x)
}

coefs_ord = do.call(
  rbind.data.frame,
    list(
      #getCoef(m1.0, year=1980),
      getCoef(m1.1, year=1990),
      getCoef(m1.2, year=2000), 
      getCoef(m1.3, year=2010),
      getCoef(m1.4, year=2020)
    )
  )




# ------------------------ family level ----------------------

# total viral richness by family
tr = vir %>%
  group_by(HostFamily) %>%
  dplyr::summarise(VRichness = n_distinct(Virus))

# unique associations by order and year
v1 = vir %>%
  dplyr::filter(Domestic == FALSE) %>%
  dplyr::group_by(HostFamily, Virus) %>%
  dplyr::summarise(Database = paste(unique(Database), collapse=", "),
                   HostSpecies = paste(unique(Host), collapse=", "),
                   NumRecords = length(Year),
                   YearEarliest = min(Year, na.rm=TRUE),
                   YearLatest = max(Year, na.rm=TRUE)) %>%
  left_join(tr) %>%
  dplyr::arrange(desc(VRichness), HostFamily, YearEarliest)

# create cumulative curves of all years
curvesf = expand.grid(unique(v1$HostFamily), 1930:endyear) %>%
  dplyr::rename("HostFamily" = 1, "YearEarliest" = 2) %>%
  left_join(v1[ , c("HostFamily", "YearEarliest", "Virus")]) %>%
  dplyr::group_by(HostFamily, YearEarliest) %>%
  dplyr::summarise(Discovered = sum(!is.na(Virus)),
                   Virus = paste(unique(Virus), collapse=", ")) %>%
  left_join(v1[ !duplicated(v1$HostFamily), c("HostFamily", "VRichness") ]) %>%
  dplyr::arrange(desc(VRichness), HostFamily, YearEarliest) %>%
  dplyr::group_by(HostFamily) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered)) %>%
  dplyr::rename("Year" = YearEarliest)

# calculate cum pubs at each step
pc = expand.grid(unique(pubs$Host), 1930:endyear) %>%
  dplyr::rename("Host" = 1, "Year" = 2) %>%
  left_join(pubs[ , c("Host", "Year", "NumPubs")]) %>%
  dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0)) %>%
  dplyr::left_join(vir %>% dplyr::select(Host, HostFamily) %>% dplyr::mutate(Host = Hmisc::capitalize(Host)) %>% distinct()) %>%
  dplyr::mutate(HostFamily = Hmisc::capitalize(HostFamily)) %>%
  dplyr::group_by(HostFamily, Year) %>%
  dplyr::summarise(NumPubs = sum(NumPubs)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(HostFamily, Year) %>%
  dplyr::group_by(HostFamily) %>%
  dplyr::mutate(CumPubs = cumsum(NumPubs))

# combine with virus curves
vir_pubs = left_join(pc, curvesf %>% dplyr::select(HostFamily, Year, VirusCumulative) %>% dplyr::mutate(HostFamily = Hmisc::capitalize(HostFamily)))



# ------------------ subset to 3 time epochs (1990, 2005, 2018) ----------------------

# data from three time epochs
vir_pubs = vir_pubs %>%
  dplyr::filter(Year %in% c(1980, 1990, 2000, 2010, 2020)) %>%
  dplyr::mutate(logPubs = log(CumPubs+1))

# combine with order richness data
vir_pubs = vir_pubs %>%
  dplyr::left_join(mam_f, by=c("HostFamily"="Family"))

# view distribtuion over time and how it relates to expected pub counts given SR
# vir_pubs = vir_pubs %>%
#   dplyr::filter(!is.na(Order)) %>%
#   group_by(Year) %>%
#   dplyr::mutate(expectedpubs = sum(CumPubs) * (SR / sum(SR)) )

# vir_pubs %>% 
#   dplyr::filter(Order %in% c("Cetartiodactyla", "Lagomorpha", "Carnivora", "Primates", "Chiroptera", "Rodentia")) %>%
#   ggplot() + 
#   geom_point(aes(log(SR), logPubs, col=Order), size=3) + 
#   geom_line(aes(log(SR), log(expectedpubs+1))) + 
#   facet_wrap(~Year) + 
#   theme_bw()
# 
# vir_pubs %>% 
#   dplyr::filter(Order %in% c("Cetartiodactyla", "Lagomorpha", "Carnivora", "Primates", "Chiroptera", "Rodentia")) %>%
#   dplyr::filter(Year == 2020) %>%
#   ggplot() + 
#   geom_point(aes(log(SR), logPubs, col=Order), size=3) + 
#   geom_line(aes(log(SR), log(expectedpubs+1))) + 
#   facet_wrap(~Order) + 
#   theme_bw()
# 
# vir_pubs %>%
#   dplyr::filter(Order %in% c("Cetartiodactyla", "Lagomorpha", "Carnivora", "Primates", "Chiroptera", "Rodentia")) %>%
#   dplyr::mutate(bias = logPubs - log(expectedpubs+1)) %>%
#   dplyr::filter(Year %in% c(1990, 2020)) %>%
#   dplyr::mutate(Year = factor(Year, levels=c(2020, 1990), ordered=TRUE),
#                 HostFamily = factor(HostFamily, levels=rev(unique(HostFamily)), ordered=TRUE)) %>%
#   ggplot() + 
#   geom_bar(aes(HostFamily, bias, fill=Year, group=Year), stat="identity", position=position_dodge(width=0.5), width=0.5) +
#   geom_hline(yintercept=0, lty=2) + 
#   theme_minimal() +
#   coord_flip() +
#   facet_wrap(~Order, scales="free_y") +
#   ylab("Log virus-related publication bias (Observed - [Expected | Species Richness])") + xlab("Host Family") + 
#   scale_fill_viridis_d(option="F", begin=0.2, end=0.8,guide = guide_legend(reverse = TRUE)) +
#   theme(strip.text = element_text(size=14),
#         axis.text.x = element_text(size=11.5),
#         axis.title = element_text(size=12))
# 
# vir_pubs %>%
#   dplyr::filter(Order %in% c("Cetartiodactyla", "Lagomorpha", "Carnivora", "Primates", "Chiroptera", "Rodentia")) %>%
#   dplyr::filter(Year %in% c(1990, 2020)) %>%
#   dplyr::group_by(Year) %>%
#   dplyr::mutate(bias = (CumPubs - expectedpubs)/sum(CumPubs)) %>%
#   ungroup() %>%
#   dplyr::mutate(Year = factor(Year, levels=c(2020, 1990), ordered=TRUE),
#                 HostFamily = factor(HostFamily, levels=rev(unique(HostFamily)), ordered=TRUE)) %>%
#   ggplot() + 
#   geom_bar(aes(HostFamily, bias, fill=Year, group=Year), stat="identity", position=position_dodge(width=0.5), width=0.5) +
#   geom_hline(yintercept=0, lty=2) + 
#   theme_minimal() +
#   coord_flip() +
#   facet_wrap(~Order, scales="free_y") +
#   ylab("Virus-related publication bias (Observed - [Expected | Species Richness])") + xlab("Host Family") + 
#   scale_fill_viridis_d(option="F", begin=0.2, end=0.8,guide = guide_legend(reverse = TRUE)) +
#   theme(strip.text = element_text(size=14),
#         axis.text.x = element_text(size=11.5),
#         axis.title = element_text(size=12))



# ------------ glm for each year --------------

vir_pubs %>%
  ggplot() + 
  geom_point(aes(logPubs, VirusCumulative)) + 
  facet_wrap(~Year, nrow=1)
vir_pubs %>%
  ggplot() + 
  geom_point(aes(log(SR), log(VirusCumulative))) + 
  facet_wrap(~Year, nrow=1)
vir_pubs %>%
  ggplot() + 
  geom_point(aes(logPubs, log(SR))) + 
  facet_wrap(~Year, nrow=1)

# ng (for virus species)
#m1.0 = MASS::glm.nb(VirusCumulative ~ logPubs + log(SR) + Order, data=vir_pubs %>% dplyr::filter(Year == 1980))
m1.1 = MASS::glm.nb(VirusCumulative ~ logPubs + log(SR) + Order, data=vir_pubs %>% dplyr::filter(Year == 1990))
m1.2 = MASS::glm.nb(VirusCumulative ~ logPubs + log(SR) + Order, data=vir_pubs %>% dplyr::filter(Year == 2000))
m1.3 = MASS::glm.nb(VirusCumulative ~ logPubs + log(SR) + Order, data=vir_pubs %>% dplyr::filter(Year == 2010))
m1.4 = MASS::glm.nb(VirusCumulative ~ logPubs + log(SR) + Order, data=vir_pubs %>% dplyr::filter(Year == 2020))

# # pois (for virus families)
# m1.0 = glm(VirusCumulative ~ logPubs + log(SR), family="poisson", data=vir_pubs %>% dplyr::filter(Year == 1980))
# m1.1 = glm(VirusCumulative ~ logPubs + log(SR), family="poisson", data=vir_pubs %>% dplyr::filter(Year == 1990))
# m1.2 = glm(VirusCumulative ~ logPubs + log(SR), family="poisson", data=vir_pubs %>% dplyr::filter(Year == 2000))
# m1.3 = glm(VirusCumulative ~ logPubs + log(SR), family="poisson", data=vir_pubs %>% dplyr::filter(Year == 2010))
# m1.4 = glm(VirusCumulative ~ logPubs + log(SR), family="poisson", data=vir_pubs %>% dplyr::filter(Year == 2020))

getCoef = function(m, year){
  x = as.data.frame(coef(summary(m)))
  x$param = row.names(x)
  names(x) = c("estimate", "se", "z", "p", "param")
  row.names(x) = c()
  x$upper = x$estimate + (1.96 * x$se)
  x$lower = x$estimate - (1.96 * x$se)
  x$year = year
  return(x)
}

coefs_fam = do.call(
  rbind.data.frame,
  list(
    #getCoef(m1.0, year=1980),
    getCoef(m1.1, year=1990),
    getCoef(m1.2, year=2000), 
    getCoef(m1.3, year=2010),
    getCoef(m1.4, year=2020)
  )
)



# ============ viz ==============

plot = rbind(
  coefs_ord %>% dplyr::mutate(level="Order-level"),
  coefs_fam %>% dplyr::mutate(level="Family-level")
) %>%
  dplyr::filter(param %in% c("logPubs", "log(SR)")) %>%
  dplyr::mutate(param = replace(param, param == "logPubs", "Virus-related\ncitations(log)"),
                param = replace(param, param == "log(SR)", "Species richness\n(log)")) %>%
  dplyr::mutate(sig_05 = p < 0.05,
                sig_01 = p < 0.01,
                level = factor(level, levels=c("Order-level", "Family-level"), ordered=TRUE)) %>%
  ggplot() + 
  geom_point(aes(param, estimate, group=year, col=factor(year)), size=3, position=position_dodge(width=0.5)) + 
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=year, col=factor(year)), position=position_dodge(width=0.5)) +
  theme_minimal() + 
  geom_hline(yintercept=0, lty=2) + 
  facet_wrap(~level, nrow=2) + 
  ylab("Estimate (mean ± 95% CI)") + 
  xlab("") + 
  ggtitle("Predictors of mammalian viral richness") +
  theme(legend.title=element_blank(),
        plot.title=element_text(hjust=0.5, size=13),
        strip.text = element_text(size=11),
        axis.text.x = element_text(size=10.5),
        legend.text = element_text(size=10.5))

ggsave(plot, 
       file="./output/timewarp_srplot.png",
       device="png", units="in", width=4.5, height=6.6, dpi=600)

