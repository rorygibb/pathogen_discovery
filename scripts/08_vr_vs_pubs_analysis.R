
# ================================================================================================

# Gibb et al., "Mammal virus diversity estimates are unstable due to accelerating discovery effort"
# Script 8: Visualise the correlation of viral richness and cumulative publication counts over time
# at several taxonomic levels

# ================================================================================================


# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon", "mgcv", "vroom")



# ---------------- build VIRION database ----------------

ictv_flag = "ictvpredict"

# domestic species to label
domestic = read.csv("./data/clovert/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)

# additional laboratory species to exlude
lab = c("macaca mulatta", "macaca fasicularis")

# publiations
pubs = read.csv("./output/host_effort/PubMed_HostsEffort_PerYear_VirusRelated_VIRION.csv") %>%
  dplyr::filter(!Host %in% c(domestic$Host, lab))

# VIRION
virion = vroom::vroom("./data/virion/Virion.csv.gz")

# build virion and add flag for whether resolved to PREDICT internal tax
vir = virion %>%
  dplyr::filter(HostClass == "mammalia" & Host != "homo sapiens") %>%
  dplyr::filter(!(Host %in% lab)) %>%
  dplyr::filter(Host %in% pubs$Host) %>%
  dplyr::mutate(Domestic = ifelse(Host %in% tolower(domestic$Host), TRUE, FALSE)) %>%
  dplyr::filter(!DetectionMethod %in% c("kmer", "Not specified")) %>%
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




# -------------- 1. Order level: Compare correlations in cumulative viral richness and cumulative publications over time -----------------

# total viral richness by order
tr = vir %>%
  group_by(HostOrder) %>%
  dplyr::summarise(VRichness = n_distinct(Virus))

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
  dplyr::mutate(Host = Hmisc::capitalize(Host), 
                NumPubs = replace(NumPubs, is.na(NumPubs), 0)) %>%
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

reso = data.frame()
for(i in 1970:endyear){
  rr = vir_pubs[ vir_pubs$Year == i, ] %>%
    dplyr::filter(VirusCumulative > 0)
  resx = data.frame(Year = i, rho = cor.test(rr$CumPubs, rr$VirusCumulative, method = "spearman")$estimate)
  resx$prho = cor.test(log(rr$CumPubs+1), log(rr$VirusCumulative+1), method = "pearson")$estimate
  calcTies = function(x){
    tx = as.data.frame(table(rank(x)))
    return(sum(tx$Freq[ tx$Freq > 1]))
  }
  resx$n_obs = nrow(rr)
  resx$ties_vr = calcTies(rr$VirusCumulative)
  resx$ties_pubs = calcTies(rr$CumPubs)
  reso = rbind(reso, resx)
}
reso$model = "Order"



# ---------------------- 2. Family level comparison of pubs versus vr ------------------------

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
  dplyr::mutate(Host = Hmisc::capitalize(Host),
                NumPubs = replace(NumPubs, is.na(NumPubs), 0)) %>%
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

resf = data.frame()
for(i in 1970:endyear){
  rr = vir_pubs[ vir_pubs$Year == i, ] %>%
    dplyr::filter(VirusCumulative > 0)
  resx = data.frame(Year = i, rho = cor.test(rr$CumPubs, rr$VirusCumulative, method = "spearman")$estimate)
  resx$prho = cor.test(log(rr$CumPubs+1), log(rr$VirusCumulative+1), method = "pearson")$estimate
  calcTies = function(x){
    tx = as.data.frame(table(rank(x)))
    return(sum(tx$Freq[ tx$Freq > 1]))
  }
  resx$n_obs = nrow(rr)
  resx$ties_vr = calcTies(rr$VirusCumulative)
  resx$ties_pubs = calcTies(rr$CumPubs)
  resf = rbind(resf, resx)
}
resf$model = "Family"




# ---------------------- 3. Species level comparison of pubs versus cumulative vr ---------------------

# Calculate discovery information for species

# total richness
tr = vir %>%
  group_by(HostOrder) %>%
  dplyr::summarise(VRichness = n_distinct(Virus))

# unique associations by year
dd = vir %>%
  dplyr::filter(Domestic == FALSE) %>%
  dplyr::group_by(Host, Virus) %>%
  dplyr::summarise(HostOrder = head(HostOrder, 1), 
                   HostFamily = head(HostFamily, 1), 
                   Database = paste(unique(Database), collapse=", "),
                   NumRecords = length(Year),
                   Year = min(Year, na.rm=TRUE),
                   YearLatest = max(Year, na.rm=TRUE),
                   Domestic = head(Domestic, 1)) %>%
  left_join(tr) %>%
  dplyr::arrange(desc(VRichness), Host, Year) 

# only keep species who had any viruses by end year of study 
dd = dd[ dd$Year <= endyear, ]

# create cumulative curves of all years
curves = expand.grid(unique(dd$Host), 1930:endyear) %>%
  dplyr::rename("Host" = 1, "Year" = 2) %>%
  left_join(dd[ , c("Host", "Year", "Virus")]) %>%
  dplyr::group_by(Host, Year) %>%
  dplyr::summarise(Discovered = sum(!is.na(Virus)),
                   Virus = paste(unique(Virus), collapse=", ")) %>%
  left_join(dd[ !duplicated(dd$Host), c("Host", "HostOrder", "HostFamily", "VRichness", "Domestic") ]) %>%
  dplyr::arrange(desc(VRichness), Host, Year) %>%
  dplyr::group_by(Host) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered))

# calculate cum pubs at each step
pc = expand.grid(unique(pubs$Host), 1930:endyear) %>%
  dplyr::rename("Host" = 1, "Year" = 2) %>%
  left_join(pubs[ , c("Host", "Year", "NumPubs")]) %>%
  dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0)) %>%
  dplyr::arrange(Host, Year) %>%
  dplyr::group_by(Host) %>%
  dplyr::mutate(CumPubs = cumsum(NumPubs))

# combine with virus curves
vir_pubs = left_join(pc, curves %>% dplyr::select(Host, Year, VirusCumulative))

ress = data.frame()
for(i in 1970:endyear){
  rr = vir_pubs[ vir_pubs$Year == i, ] %>%
    dplyr::filter(VirusCumulative > 0)
  resx = data.frame(Year = i, rho = cor.test(rr$CumPubs, rr$VirusCumulative, method = "spearman")$estimate)
  resx$prho = cor.test(log(rr$CumPubs+1), log(rr$VirusCumulative+1), method = "pearson")$estimate
  calcTies = function(x){
    tx = as.data.frame(table(rank(x)))
    return(sum(tx$Freq[ tx$Freq > 1]))
  }
  resx$n_obs = nrow(rr)
  resx$ties_vr = calcTies(rr$VirusCumulative)
  resx$ties_pubs = calcTies(rr$CumPubs)
  ress = rbind(ress, resx)
}
ress$model = "Species"


# ========== visualise these results ==============

lims = range(c(ress$rho, resf$rho, reso$rho))
lims = c(lims[1], 1)

p1 = rbind(reso, resf, ress) %>%
  dplyr::mutate(model = factor(model, levels=c("Order", "Family", "Species"), ordered=TRUE)) %>%
  ggplot() + 
  geom_line(aes(Year, prho, col=model), size=0.8) +
  theme_minimal() +
  xlab("Year") + ylab("Correlation between\ncitation count and viral richness") +
  theme(legend.title=element_blank(),
        legend.position=c(0.85, 0.15),
        legend.text = element_text(size=12),
        axis.title = element_text(size=14),
        axis.text = element_text(size=11)) + 
  scale_y_continuous(limits=lims) +
  #ggtitle("Order and species-level (mammals)") +
  theme(plot.title = element_text(size=14, hjust=0.5)) + 
  scale_color_viridis_d(option="F", end=0.8)

ggsave(p1, file="./output/figures_2021/SI_PublicationViralRichnessCorrelation.png", device="png", units="in", width=6, height=5)
