


# ============================= Use "time slice" correlation coefficients to examine temporal stability of virus diversity ==========================

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon", "mgcv")

# domestic species to label
domestic = read.csv("./data/clovert/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)

# associations data and no humans
clover = read.csv("./data/clovert/CLOVERMammalViruses1.0_AssociationsFlatFile.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Host != "Homo sapiens") %>%
  dplyr::mutate(Domestic = ifelse(Host %in% domestic$Host, TRUE, FALSE)) %>%
  #dplyr::filter(DetectionMethod != "Antibodies") %>%
  dplyr::filter(DetectionMethod != "Not specified")

# create combined year column
clover$Year = clover$PublicationYear
clover$Year[ is.na(clover$Year) ] = clover$ReleaseYear[ is.na(clover$Year) ]
clover = clover[ !is.na(clover$Year), ]

# total viral richness by species
tr = clover %>%
  group_by(Host) %>%
  dplyr::summarise(VRichness = n_distinct(Virus))



# Calculate discovery information for species

# unique associations by year
dd = clover %>%
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
endyear = 2010
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

# species-level: 935 species that had any viruses discovered by 2010
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
orders = c("Cetartiodactyla", "Rodentia", "Lagomorpha", "Chiroptera", "Carnivora", "Primates")
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
  geom_line(data = d1, aes(Year, SpearmanCoef, col=Group), size=0.8) +
  geom_line(data = sp_mean, aes(Year, SpearmanCoef, col=Group), size=0.7, lty=5, alpha=0.9) +
  geom_line(data = sp_mean, aes(Year, SpearmanCoef, col=Group), size=0.25, alpha=0.2) +
  geom_vline(xintercept=2010, lty=2) + 
  theme_minimal() +
  ylim(0, 1.0) +
  ylab("Correlation with viral richness in 2010") +
  scale_color_viridis_d(option="magma", begin=0, end=0.85) +
  theme(legend.position=c(0.8, 0.2),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=11),
        axis.text.x = element_blank(),
        legend.text = element_text(size=10),
        axis.title.y = element_text(hjust = -9, size=13))


# At the Order level with colourblind friendly palette
pal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
         "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(pal)
scale = pal[ c(1, 4:8) ]
p2 = ggplot(stab_all[ !stab_all$Group %in% c("Species", "Order", "Family"), ]) +
  geom_line(aes(Year, SpearmanCoef, col=Group), size=0.8, alpha=0.9) +
  geom_vline(xintercept=2010, lty=2) +
  theme_minimal() +
  ylim(0.1, 1.0) +
  ylab("") +
  scale_color_manual(values=scale) +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=13),
        legend.text = element_text(size=10.5),
        axis.title.y = element_text(size=13, color="white"))

# p2 = ggplot(stab_all[ !stab_all$Group %in% c("Species", "Order", "Family"), ]) + 
#   geom_line(aes(Year, SpearmanCoef, col=Group), size=0.8) +
#   geom_vline(xintercept=2010, lty=2) + 
#   theme_minimal() +
#   ylim(0, 1.0) +
#   ylab("") +
#   scale_color_manual(values=scale) +
#   theme(legend.position=c(0.5, 0.1),
#         legend.title = element_blank(),
#         axis.text = element_text(size=11),
#         axis.title.x = element_text(size=13),
#         legend.text = element_text(size=10),
#         axis.title.y = element_text(size=13, color="white")) + 
#   guides(color=guide_legend(nrow=2))

# combine and annotate
p_comb = gridExtra::grid.arrange(p1, p2, ncol=1, heights=c(0.8, 1)) 
pp = ggpubr::as_ggplot(p_comb)  +
  cowplot::draw_plot_label(label = c("a", "b"), size = 20, 
                           x = c(0.15, 0.15), y = c(0.99, 0.55))
ggsave(pp, file="./output/figures/MS_Figure2_Stability_types.png", dpi=600, height=8, width=4.4, units="in")




# ===================== Combine and plot with mean species level viral richness =========================

# combine and plot all
stab_all = rbind(stab, stab_oa, stab_fa, stabos)

# Across mammals
d1 = stab_all[ stab_all$Group %in% c("Species", "Order", "Family"), ] %>%
  dplyr::mutate(Group = factor(Group, levels = c("Order", "Family", "Species"), ordered=TRUE))
p1 = ggplot(d1) + 
  geom_line(aes(Year, SpearmanCoef, col=Group), size=0.8) +
  geom_vline(xintercept=2010, lty=2) + 
  theme_minimal() +
  ylim(0, 1.0) +
  ylab("Correlation with viral richness in 2010") +
  scale_color_viridis_d(option="magma", begin=0, end=0.8) +
  theme(legend.position=c(0.8, 0.2),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=11),
        legend.text = element_text(size=10.5),
        axis.title.y = element_text(hjust = -9, size=13))


# At the Order level with colourblind friendly palette
pal <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
         "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(pal)
scale = pal[ c(1, 4:8) ]
p2 = ggplot(stab_all[ !stab_all$Group %in% c("Species", "Order", "Family"), ]) + 
  geom_line(aes(Year, SpearmanCoef, col=Group), size=0.8) +
  geom_vline(xintercept=2010, lty=2) + 
  theme_minimal() +
  ylim(0.2, 1.0) +
  ylab("") +
  scale_color_manual(values=scale) +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.text = element_text(size=11),
        axis.title.x = element_text(size=13),
        legend.text = element_text(size=10.5),
        axis.title.y = element_text(size=13, color="white"))

p2 = ggplot(stab_all[ !stab_all$Group %in% c("Species", "Order", "Family"), ]) + 
  geom_line(aes(Year, SpearmanCoef, col=Group), size=0.8) +
  geom_vline(xintercept=2010, lty=2) + 
  theme_minimal() +
  ylim(0, 1.0) +
  ylab("") +
  scale_color_manual(values=scale) +
  theme(legend.position=c(0.5, 0.1),
        legend.title = element_blank(),
        axis.text = element_text(size=11),
        axis.title.x = element_text(size=13),
        legend.text = element_text(size=10.5),
        axis.title.y = element_text(size=13, color="white")) + 
  guides(color=guide_legend(nrow=2))

# combine and annotate
p_comb = gridExtra::grid.arrange(p1, p2, ncol=1, heights=c(0.8, 1)) 
pp = ggpubr::as_ggplot(p_comb)  +
  cowplot::draw_plot_label(label = c("a", "b"), size = 20, 
                           x = c(0.15, 0.15), y = c(0.99, 0.55))
ggsave(pp, file="./output/figures/MS_Figure2_Stability_MeanVR.png", dpi=600, height=8, width=4.3, units="in")




# =================== Bump plot of Order ranks =====================

# get ranks in 5 year intervals
yearrank = function(x){
  bx = spso[ spso$Year == x, ] %>%
    dplyr::arrange(desc(ViralRichness)) 
  bx$Rank = rank(bx$ViralRichness)
  bx
}
bpd = do.call(rbind.data.frame, lapply(seq(1960, 2010, by=5), yearrank))

bp = ggplot(bpd) + 
  geom_point(aes(Year, Rank, col=Host, group=Host), size=3) + 
  geom_line(aes(Year, Rank, col=Host, group=Host), alpha=0.5, size=1) + 
  theme_minimal() +
  scale_y_continuous(breaks=1:20, labels=1:20) + 
  geom_text(data = bpd[ bpd$Year == 2010, ], aes(x = 2012, y=Rank, label=Host), size = 4.5, hjust = 0) +
  theme(legend.position="none") + 
  scale_x_continuous(limits=c(1960, 2022), breaks=seq(1960, 2010, by=10), labels=seq(1960, 2010, by=10)) +
  theme(axis.text = element_text(size=11.5),
        axis.title = element_text(size=13)) +
  ylab("Viral diversity rank")
ggsave(bp, file="./output/figures/SI_OrderRanks_BumpChart.png", dpi=600, height=8, width=9.5, units="in", scale=0.9)




# get ranks in 5 year intervals
yearrank = function(x){
  bx = sps_oa[ sps_oa$Year == x, ] %>%
    dplyr::arrange(desc(ViralRichness)) 
  bx$Rank = rank(bx$ViralRichness)
  bx
}
bpd = do.call(rbind.data.frame, lapply(seq(1960, 2010, by=5), yearrank))

bp = ggplot(bpd) + 
  geom_point(aes(Year, Rank, col=Host, group=Host), size=3) + 
  geom_line(aes(Year, Rank, col=Host, group=Host), alpha=0.5, size=1) + 
  theme_minimal() +
  scale_y_continuous(breaks=1:20, labels=1:20) + 
  geom_text(data = bpd[ bpd$Year == 2010, ], aes(x = 2012, y=Rank, label=Host), size = 4.5, hjust = 0) +
  theme(legend.position="none") + 
  scale_x_continuous(limits=c(1960, 2022), breaks=seq(1960, 2010, by=10), labels=seq(1960, 2010, by=10)) +
  theme(axis.text = element_text(size=11.5),
        axis.title = element_text(size=13)) +
  ylab("Viral diversity rank")
ggsave(bp, file="./output/figures/SI_OrderRanks_BumpChart_MeanVR.png", dpi=600, height=8, width=9.5, units="in", scale=0.9)


