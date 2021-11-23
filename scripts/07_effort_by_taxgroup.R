

# ========================================= ======================================================

# Gibb et al., "Mammal virus diversity estimates are unstable due to accelerating discovery effort"
# Script 7: Visualise the relationship of publication counts by taxonomic group and relation to 
# species richness

# ================================================================================================

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon", "mgcv", "vroom", "sf")




# ----------------- Species richness estimates for mammal groups from IUCN ----------------

# mam1 = sf::st_read("C:/Users/roryj/Documents/PhD/202008_discovery/data/iucn_range/mammals_freshwater/MAMMALS_FRESHWATER.shp")
# mam2 = sf::st_read("C:/Users/roryj/Documents/PhD/202008_discovery/data/iucn_range/mammals_terrestrial/MAMMALS_TERRESTRIAL_ONLY.shp")
# mam3 = sf::st_read("C:/Users/roryj/Documents/PhD/202008_discovery/data/iucn_range/mammals_marine/MAMMALS_MARINE_ONLY.shp")
# mam4 = sf::st_read("C:/Users/roryj/Documents/PhD/202008_discovery/data/iucn_range/mammals_terrmar/MAMMALS_MARINE_AND_TERRESTRIAL/MAMMALS_MARINE_AND_TERRESTRIAL.shp")
# 
# mam1 = mam1 %>% st_drop_geometry() %>% dplyr::select(binomial, order_, family, genus) %>% distinct()
# mam2 = mam2 %>% st_drop_geometry() %>% dplyr::select(binomial, order_, family, genus) %>% distinct()
# mam3 = mam3 %>% st_drop_geometry() %>% dplyr::select(binomial, order_, family, genus) %>% distinct()
# mam4 = mam4 %>% st_drop_geometry() %>% dplyr::select(binomial, order_, family, genus) %>% distinct()
# 
# # recode to cetacea
# mam1$order_[ mam1$family %in% c("BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE", 
#                                 "ESCHRICHTIIDAE", "KOGIIDAE", "PHOCOENIDAE", "PONTOPORIIDAE", "PHYSETERIDAE", "INIIDAE", "LIPOTIDAE", "PLATANISTIDAE") ] = "CETACEA"
# mam3$order_[ mam3$family %in% c("BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE", 
#                                 "ESCHRICHTIIDAE", "KOGIIDAE", "PHOCOENIDAE", "PONTOPORIIDAE", "PHYSETERIDAE", "INIIDAE", "LIPOTIDAE", "PLATANISTIDAE") ] = "CETACEA"
# 
# # order level species richness
# mam_o = rbind(mam1, mam2) %>%
#   rbind(mam3) %>%
#   rbind(mam4) %>%
#   dplyr::group_by(order_) %>%
#   dplyr::summarise(SR = n_distinct(binomial)) %>%
#   dplyr::rename("Order"=order_) %>%
#   dplyr::mutate(Order = Hmisc::capitalize(tolower(Order)),
#                 Order = replace(Order, Order == "Cetartiodactyla", "Artiodactyla"))
# 
# # family level species richness
# mam_f = rbind(mam1, mam2) %>%
#   rbind(mam3) %>%
#   rbind(mam4) %>%
#   dplyr::group_by(family) %>%
#   dplyr::summarise(SR = n_distinct(binomial),
#                    Order = head(order_, 1)) %>%
#   dplyr::rename("Family"=family) %>%
#   dplyr::mutate(Family = Hmisc::capitalize(tolower(Family)),
#                 Order = Hmisc::capitalize(tolower(Order)),
#                 Order = replace(Order, Order == "Cetartiodactyla", "Artiodactyla"))

# write.csv(mam_o, "./data/speciesrichness/order_speciesrichness_iucn.csv", row.names=FALSE)
# write.csv(mam_f, "./data/speciesrichness/family_speciesrichness_iucn.csv", row.names=FALSE)

mam_o = read.csv("./data/speciesrichness/order_speciesrichness_iucn.csv")
mam_f = read.csv("./data/speciesrichness/family_speciesrichness_iucn.csv")



# ================== build VIRION database ====================

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

# add host 
pubs = dplyr::left_join(
  pubs %>% dplyr::mutate(Host = tolower(Host)), 
  vir %>% dplyr::select(Host, HostOrder, HostFamily) %>% distinct()
)




# ========================== distribution of observed and expected effort at Order and Family level =================================

# order level summary
viro = vir %>%
  dplyr::filter(Domestic == FALSE) %>%
  dplyr::group_by(HostOrder) %>%
  dplyr::summarise(VRich = n_distinct(Virus),
                   VGGenusRich = n_distinct(VirusGenus)) %>%
  dplyr::right_join(vir %>% dplyr::select(HostOrder) %>% distinct) %>%
  dplyr::left_join(mam_o %>% dplyr::select(Order, SR) %>% dplyr::rename("HostOrder"=Order) %>% dplyr::mutate(HostOrder = tolower(HostOrder))) %>%
  dplyr::left_join(
    pubs %>% dplyr::group_by(HostOrder) %>% dplyr::summarise(NumPubs = sum(NumPubs))
  ) %>%
  dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0),
                VRich  = replace(VRich , is.na(VRich ), 0),
                VGGenusRich = replace(VGGenusRich, is.na(VGGenusRich), 0)) %>%
  dplyr::mutate(NumPubs_Expected = sum(NumPubs) * (SR/sum(SR)))

vv = viro %>% 
  dplyr::mutate(logPubs = log(NumPubs + 1),
                logExpPubs = log(NumPubs_Expected + 1),
                HostOrder = Hmisc::capitalize(HostOrder))
  
p1 = vv %>%
  ggplot() + 
  geom_point(aes(log(SR), logPubs, col=HostOrder), size=4) +
  geom_text(aes(log(SR), logPubs, label=HostOrder),  check_overlap=TRUE, size=4, nudge_y=-0.2) +
  geom_text(data=vv[ vv$HostOrder == "Primates", ], aes(log(SR), logPubs, label=HostOrder), size=4, nudge_y=0.3, nudge_x = 0.2) +
  geom_text(data=vv[ vv$HostOrder == "Carnivora", ], aes(log(SR), logPubs, label=HostOrder), size=4, nudge_y=0.2, nudge_x = -0.4) +
  geom_line(aes(log(SR), logExpPubs)) +
  theme_bw() + 
  xlab("Order species richness (log)") + ylab("Virus-related publications (log)") + 
  theme(legend.position="none", legend.title = element_blank()) + 
  xlim(0, 8.2) + 
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8), labels=c(0, 2, 4, 6, 8)) +
  theme(axis.title = element_text(size=13))


# family level summary
virf = vir %>%
  dplyr::filter(Domestic == FALSE) %>%
  dplyr::group_by(HostFamily) %>%
  dplyr::summarise(VRich = n_distinct(Virus),
                   VGGenusRich = n_distinct(VirusGenus)) %>%
  dplyr::right_join(vir %>% dplyr::select(HostFamily) %>% distinct) %>%
  dplyr::left_join(mam_f %>% dplyr::select(Family, Order, SR) %>% dplyr::rename("HostFamily"=Family) %>% dplyr::mutate(HostFamily = tolower(HostFamily))) %>%
  dplyr::rename("HostOrder" = Order) %>%
  dplyr::left_join(
    pubs %>% dplyr::group_by(HostFamily) %>% dplyr::summarise(NumPubs = sum(NumPubs))
  ) %>%
  dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0),
                VRich  = replace(VRich , is.na(VRich ), 0),
                VGGenusRich = replace(VGGenusRich, is.na(VGGenusRich), 0)) %>%
  dplyr::filter(!is.na(SR)) %>%
  dplyr::mutate(NumPubs_Expected = sum(NumPubs) * (SR/sum(SR)))

ff = virf %>% 
  dplyr::filter(HostOrder %in% c("Carnivora", "Artiodactyla", "Chiroptera", "Eulipotyphla", "Lagomorpha", "Perissodactyla",
                                 "Primates", "Rodentia", "Diprotodontia")) %>%
  dplyr::mutate(logPubs = log(NumPubs + 1),
                logExpPubs = log(NumPubs_Expected + 1)) 
p2 = ggplot() + 
  geom_point(data=ff, aes(log(SR), logPubs), size=3, alpha=0.6) +
  geom_line(aes(log(SR), logExpPubs), data=ff %>% dplyr::select(-HostOrder), size=0.6, alpha=0.7) +
  theme_bw() +
  facet_wrap(~HostOrder) + 
  xlab("Family species richness (log)") + ylab("Virus-related publications (log)") + 
  theme(legend.position="none", strip.background = element_blank(),
        strip.text = element_text(size=13), axis.title = element_text(size=13))


# combine all of these
p_comb = gridExtra::grid.arrange(p1, p2, nrow=1)
ggsave(p_comb, file="./output/figures_2021/SI_EffortDistribution_byTaxonomy.png", device="png", units="in", width=11, height=5.5, dpi=600)


