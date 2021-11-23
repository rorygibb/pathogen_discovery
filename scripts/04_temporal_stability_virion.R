


# ========================================= ======================================================

# Gibb et al., "Mammal virus diversity estimates are unstable due to accelerating discovery effort"
# Script 4: Examine temporal stability of relative virus diversity estimates across species/groups

# ================================================================================================



# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon", "mgcv")

# domestic and lab species to label
domestic = read.csv("./data/clovert/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)
lab = c("macaca mulatta", "macaca fasicularis")

# VIRION
virion = vroom::vroom("./data/virion/Virion.csv.gz")

# build virion and add flag for whether resolved to PREDICT internal tax
vir = virion %>%
  dplyr::filter(HostClass == "mammalia" & Host != "homo sapiens") %>%
  dplyr::filter(!(Host %in% lab)) %>%
  dplyr::mutate(Domestic = ifelse(Host %in% tolower(domestic$Host), TRUE, FALSE)) %>%
  dplyr::filter(!DetectionMethod %in% c("kmer")) %>%
  dplyr::filter(!is.na(ReleaseYear) | !is.na(PublicationYear) | !is.na(CollectionYear)) %>%
  dplyr::mutate(PredictFlag = str_detect(Virus, "predict\\_"),
                CLOVERflag = ifelse(Database %in% c("EID2", "Shaw", "GMPD2", "HP3"), TRUE , FALSE)) %>%
  dplyr::filter(ICTVRatified | PredictFlag)

# temporary fix
# vir$HostOrder[ vir$HostFamily %in% c("balaenidae", "balaenopteridae", "delphinidae",
#                                      "eschrichtiidae", "monodontidae", "phocoenidae",
#                                      "physeteridae", "ziphiidae")] = "cetacea"

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
endyear = 2018
vir = vir[ vir$Year <= endyear, ]




# ------------------ Calculate discovery information for species --------------------

# total richness
tr = vir %>%
group_by(HostOrder) %>%
  dplyr::summarise(VRichness = n_distinct(Virus))

# unique associations by year
dd = vir %>%
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






# ===================== Stability of viral richness estimates over time ===========================

# species-level: 1271 hosts
# function that gets viral richness of all species in a specified year
getYearVRichness = function(year, data){
  
  yrich = data %>%
    dplyr::filter(Year <= year) %>%
    dplyr::filter(Host != "Homo sapiens") %>%
    dplyr::filter(Domestic == FALSE) %>%
    dplyr::group_by(Host) %>%
    dplyr::summarise(ViralRichness = sum(Discovered),
                     #HostOrder = head(HostOrder, 1),
                     Domestic = head(Domestic, 1)) %>%
    dplyr::mutate(Year = year)
  
}

# function that compares a given year to the future year (input is output of getYearVRichness function)
spearmanStability = function(year, data){
  
  # dataframe of current and future end step viral richness
  spx = data[ data$Year == year, ]
  spx = left_join(spx, 
                  data[ data$Year == endyear, c("Host", "ViralRichness")] %>% dplyr::rename("ViralRichness_End" = ViralRichness))
  
  # calculate coef
  rho = as.vector(cor.test(spx$ViralRichness, spx$ViralRichness_End, method="spearman")$estimate)
  return(rho)
}

# start and end years
startyear = 1960
endyear = endyear



# --------------- total viral richness ----------------

# run across all wild species (n=935)
sps = do.call(rbind.data.frame, lapply(startyear:endyear, getYearVRichness, data=curves))
stab = data.frame(Year = startyear:endyear, SpearmanCoef = sapply(startyear:endyear, spearmanStability, data=sps))
stab$Group = "Species"

# run at the level of orders and families
curvesx = curves %>%
  dplyr::mutate(Host = HostOrder)
spso = do.call(rbind.data.frame, lapply(startyear:endyear, getYearVRichness, data=curvesx))
stabo = data.frame(Year = startyear:endyear, SpearmanCoef = sapply(startyear:endyear, spearmanStability, data=spso))
stabo$Group = "Order"

curvesx = curves %>%
  dplyr::mutate(Host = HostFamily)
spsf = do.call(rbind.data.frame, lapply(startyear:endyear, getYearVRichness, data=curvesx))
stabf = data.frame(Year = startyear:endyear, SpearmanCoef = sapply(startyear:endyear, spearmanStability, data=spsf))
stabf$Group = "Family"

# run at the species-level across Orders
orders = c("artiodactyla", "rodentia", "lagomorpha", "chiroptera", "carnivora", "primates")
stabos = data.frame()
for(i in orders){
  
  vrx = do.call(rbind.data.frame, lapply(startyear:endyear, getYearVRichness, data=curves[ curves$HostOrder == i, ]))
  stabx = data.frame(Year = startyear:endyear, SpearmanCoef = sapply(startyear:endyear, spearmanStability, data=vrx))
  stabx$Group = i
  stabos = rbind(stabos, stabx)
}



# ---------------- mean species level richness at family and order level -------------------

# species-level: 935 species that had any viruses discovered by 2010
# function that gets viral richness of all species in a specified year
getYearVRichnessAggregate = function(year, data, aggregate = "order"){
  
  yrich = data %>%
    dplyr::filter(Year <= year) %>%
    dplyr::filter(Host != "Homo sapiens") %>%
    dplyr::filter(Domestic == FALSE) %>%
    dplyr::group_by(Host) %>%
    dplyr::summarise(ViralRichness = sum(Discovered),
                     #HostOrder = head(HostOrder, 1),
                     Domestic = head(Domestic, 1),
                     HostOrderx = head(HostOrder, 1),
                     HostFamilyx = head(HostFamily, 1)) %>%
    dplyr::mutate(Year = year)
  
  # aggregate
  if(aggregate == "order"){
    yrich = yrich %>%
      dplyr::group_by(HostOrderx) %>%
      dplyr::summarise(Host = head(HostOrderx, 1),
                       ViralRichness = mean(ViralRichness),
                       Year = head(Year, 1))
  }
  if(aggregate == "family"){
    yrich = yrich %>%
      dplyr::group_by(HostFamilyx) %>%
      dplyr::summarise(Host = head(HostFamilyx, 1),
                       ViralRichness = mean(ViralRichness),
                       Year = head(Year, 1))
  }
  
  return(yrich)
}


# run for orders
sps_oa = do.call(rbind.data.frame, lapply(startyear:endyear, getYearVRichnessAggregate, data=curves, aggregate = "order"))
stab_oa = data.frame(Year = startyear:endyear, SpearmanCoef = sapply(startyear:endyear, spearmanStability, data=sps_oa))
stab_oa$Group = "Order"

# run for families
sps_fa = do.call(rbind.data.frame, lapply(startyear:endyear, getYearVRichnessAggregate, data=curves, aggregate = "family"))
stab_fa = data.frame(Year = startyear:endyear, SpearmanCoef = sapply(startyear:endyear, spearmanStability, data=sps_fa))
stab_fa$Group = "Family"






# ================= Combine and plot all with total viral richness ================

# combine and plot all
stab_all = rbind(stab, stabo, stabf, stabos)

# combine mean species-level in tax groups
sp_mean = rbind(stab_oa, stab_fa)
sp_mean$Group = factor(sp_mean$Group, levels = c("Order", "Family"), ordered=TRUE)

# Across mammals
d1 = stab_all[ stab_all$Group %in% c("Species", "Order", "Family"), ] %>%
  dplyr::mutate(Group = factor(Group, levels = c("Order", "Family", "Species"), ordered=TRUE))
p1 = ggplot() + 
  geom_vline(xintercept=2018, lty=1) +
  #geom_hline(yintercept=0.5, lty=2) + 
  geom_line(data = d1, aes(Year, SpearmanCoef, col=Group), size=0.8) +
  geom_line(data = sp_mean, aes(Year, SpearmanCoef, col=Group), size=0.7, lty=5, alpha=0.9) +
  geom_line(data = sp_mean, aes(Year, SpearmanCoef, col=Group), size=0.25, alpha=0.2) +
  theme_minimal() +
  ylim(0.1, 1.0) +
  ylab("") +
  scale_color_viridis_d(option="magma", begin=0, end=0.85) +
  scale_x_continuous(breaks=c(1960, 1980, 2000, 2018), labels=c(1960, 1980, 2000, 2018)) +
  theme(legend.position=c(0.8, 0.2),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size=11),
        axis.text.x = element_blank(),
        legend.text = element_text(size=10),
        axis.title.y = element_text(hjust = -8, size=13))


# At the Order level with colourblind friendly palette
pal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
         "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(pal)
scale = pal[ c(1, 4:8) ]
stab_ord = stab_all[ !stab_all$Group %in% c("Species", "Order", "Family"), ]
stab_ord$Group = Hmisc::capitalize(stab_ord$Group)

p2 = ggplot(stab_ord) +
  geom_vline(xintercept=2018, lty=1) +
  #geom_hline(yintercept=0.5, lty=2) + 
  geom_line(aes(Year, SpearmanCoef, col=Group), size=0.8, alpha=0.9) +
  theme_minimal() +
  ylim(0.1, 1.0) +
  ylab("Correlation with viral richness rank in 2018") +
  scale_color_manual(values=scale) +
  scale_x_continuous(breaks=c(1960, 1980, 2000, 2018), labels=c(1960, 1980, 2000, 2018)) +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=13),
        panel.grid.minor.y = element_blank(),
        legend.text = element_text(size=10.5),
        axis.title.y = element_text(hjust = -4.6, size=13, color="black"))

# combine and annotate
p_comb = gridExtra::grid.arrange(p1, p2, ncol=1, heights=c(0.8, 1)) 
pp = ggpubr::as_ggplot(p_comb)  +
  cowplot::draw_plot_label(label = c("a", "b"), size = 20, 
                           x = c(0.15, 0.15), y = c(0.99, 0.55))
ggsave(pp, file="./output/figures_2021/MS_Figure2_Stability_types.png", dpi=600, height=8, width=4.4, units="in")



