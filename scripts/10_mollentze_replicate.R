

# =================== Heuristic examination of sensitity of ecological inferences to changes in knowledge =====================

# We replicate the Order-level viral richness analyses of Mollentze et al which show that zoonotic viral richness is mainly predicted 
# by species richness of each order (i.e. a neutral explanation for patterns of viral diversity in nature)
# Fit model to order-level total viral richness data (notably, from our analyses in Figure 2, the most stable metric over time in terms of correlation)
# Fit a negative binomial GLM considering either logCitations or logCitations + logSR at Order-level (n=17 orders)
# Examine change in parameter over time

# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon", "mgcv", "vroom")

library(sf)
library(dplyr)




# ------------------ Species richness estimates at Order and Family level from IUCN -------------------

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




# --------------- VIRION viral discovery data -------------------

ictv_flag = "ictvpredict"

# domestic species and lab species
domestic = read.csv("./data/clovert/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)
lab = c("macaca mulatta", "macaca fasicularis")

# VIRION
virion = vroom::vroom("./data/virion/Virion.csv.gz")

# publiations
pubs = read.csv("./output/host_effort/PubMed_HostsEffort_PerYear_VirusRelated_VIRION.csv") %>%
  dplyr::filter(Note != "Lookup error") %>%
  dplyr::filter(!Host %in% c(domestic$Host, lab)) %>% 
  dplyr::mutate(Host = Hmisc::capitalize(Host)) %>%
  dplyr::left_join(virion %>% dplyr::select(HostOrder, Host) %>% distinct() %>% dplyr::mutate(Host = Hmisc::capitalize(Host)))

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




# -------------------- Order-level viral richness and publication effort estimates at each time epoch ------------------------

epochs = c(1990, 1995, 2000, 2005, 2010, 2015, 2020)
dd = data.frame()

for(i in epochs){
  
  # summarise
  vir_i = vir %>%
    dplyr::filter(Year <= i) %>%
    dplyr::filter(Domestic == FALSE) %>%
    dplyr::group_by(HostOrder) %>%
    dplyr::summarise(VirRich = n_distinct(Virus),
                     NHosts = n_distinct(Host))
  
  # combine with other orders with no records yet
  toadd = vir %>%
    dplyr::select(HostOrder) %>%
    dplyr::filter(!HostOrder %in% vir_i$HostOrder) %>%
    distinct() %>%
    dplyr::mutate(VirRich = 0, NHosts = 0)
  vir_i = rbind(vir_i, toadd)
  
  # add publications
  pubs_i = pubs %>%
    dplyr::filter(Year <= i) %>%
    dplyr::group_by(HostOrder) %>%
    dplyr::summarise(Citations = sum(NumPubs, na.rm=TRUE))
  vir_i = vir_i %>%
    dplyr::left_join(pubs_i) %>%
    dplyr::mutate(Citations = replace(Citations, is.na(Citations), 0),
                  Year = i)
  
  dd = rbind(dd, vir_i)
}

# ggplot(dd) + 
#   geom_line(aes(Year, VirRich, col=HostOrder, group=HostOrder))

# combine with species richness data
dd = left_join(
  dd,
  mam_o %>% dplyr::mutate(Order = tolower(Order)),
  by=c("HostOrder"="Order")
)

# log transform
dd = dd %>% dplyr::mutate(logSR = log(SR), logCitations = log(Citations+1))


host_levels =  dd %>% 
  dplyr::filter(Year == 1990) %>%
  dplyr::arrange(logSR)

p1 = dd %>% 
  dplyr::filter(Year == 1990) %>%
  dplyr::arrange(desc(logSR)) %>%
  dplyr::mutate(HostOrder =Hmisc::capitalize(HostOrder), 
                HostOrder = factor(HostOrder, levels=Hmisc::capitalize(host_levels$HostOrder), ordered=TRUE)) %>%
  ggplot() + 
  geom_bar(aes(HostOrder, logSR, fill=HostOrder), stat="identity", width=0.65) + 
  coord_flip() +
  theme_classic() +
  theme(legend.position="none") +
  #scale_fill_viridis_d(direction=-1) + 
  scale_fill_discrete(viridisLite::turbo(17)) +
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=10.5)) +
  xlab("Host Order") + ylab("Species richness (log)")

p3 = dd %>% 
  dplyr::mutate(HostOrder = factor(HostOrder, levels=host_levels$HostOrder, ordered=TRUE)) %>%
  ggplot() + 
  geom_line(aes(Year, logCitations, group=HostOrder, col=HostOrder), size=0.7) + 
  theme_classic() +
  theme(legend.position="none") + 
  #scale_color_viridis_d(direction=-1) + 
  scale_color_discrete(viridisLite::turbo(17)) +
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=10.5)) +
  xlab("Year") + ylab("Virus-related citations (log)")
  
p2 = dd %>% 
  dplyr::mutate(HostOrder = factor(HostOrder, levels=host_levels$HostOrder, ordered=TRUE)) %>%
  ggplot() + 
  geom_line(aes(Year, log(VirRich+1), group=HostOrder, col=HostOrder), size=0.7) + 
  theme_classic() + 
  theme(legend.position="none") +
  #scale_color_viridis_d(direction=-1) + 
  scale_color_discrete(viridisLite::turbo(17)) + 
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=10.5)) +
  xlab("Year") + ylab("Viral richness (log)")

pc1 = gridExtra::grid.arrange(p1, p2, p3, nrow=1, widths=c(1, 0.8, 0.8))



# -------------------- Fit models comparing AIC and R2 ------------------

fitNBModel = function(year, data){
  
  data = data %>% dplyr::filter(Year == year)
  m1 = MASS::glm.nb(VirRich ~ logCitations, data=data)
  m2 = MASS::glm.nb(VirRich ~ logCitations + logSR, data=data)
  
  # calculate GoF (R2)
  r2.1 = with(summary(m1), 1 - deviance/null.deviance)
  r2.2 = with(summary(m2), 1 - deviance/null.deviance)
  
  # extract coefs
  getCoef = function(m){
    x = as.data.frame(coef(summary(m)))
    x$param = row.names(x)
    names(x) = c("estimate", "se", "z", "p", "param")
    row.names(x) = c()
    x$upper = x$estimate + (1.96 * x$se)
    x$lower = x$estimate - (1.96 * x$se)
    return(x)
  }
  cc = getCoef(m2)
  
  # add stats
  cc$R2_withoutSR = r2.1
  cc$R2_withSR = r2.2
  cc$year = year
  
  return(cc)
}

# run for each year combination
results = do.call(
  rbind.data.frame,
  lapply(c(1990, 1995, 2000, 2005, 2010, 2015, 2020), fitNBModel, data=dd)
)

# visualise
p4 = results %>%
  dplyr::filter(param != "(Intercept)") %>%
  dplyr::mutate(param = replace(param, param == "logCitations", "Virus-related citations (log)"),
                param = replace(param, param == "logSR", "Species richness (log)")) %>%
  dplyr::mutate(sig_05 = p < 0.05,
                sig_01 = p < 0.01) %>%
  ggplot() + 
  geom_point(aes(param, estimate, group=year, col=factor(year)), size=3, position=position_dodge(width=0.5)) + 
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=year, col=factor(year)), position=position_dodge(width=0.5)) +
  theme_classic() + 
  geom_hline(yintercept=0, lty=2) + 
  ylab("Slope estimate (mean + 95% CI)") + 
  xlab("") + 
  ggtitle("Predictors of Order-level viral richness") +
  theme(legend.title=element_blank(),
        legend.position="left",
        #plot.title=element_text(hjust=0.5, size=13),
        plot.title = element_blank(),
        strip.text = element_text(size=11),
        axis.text = element_text(size=11),
        axis.text.x = element_text(size=11.5),
        axis.title = element_text(size=11),
        legend.text = element_text(size=10.5)) +
  scale_color_viridis_d(option="F", direction=-1, end=0.8)

# p5 = results %>%
#   dplyr::select(year, R2_withoutSR, R2_withSR) %>%
#   distinct() %>%
#   reshape2::melt(id.vars=1) %>%
#   dplyr::mutate(variable = as.vector(variable), 
#                 variable = replace(variable, variable == "R2_withoutSR", "logCitations"),
#                 variable = replace(variable, variable == "R2_withSR", "logCitations + logSR"),
#                 variable = factor(variable, levels=c("logCitations + logSR", "logCitations"), ordered=TRUE),
#                 year = factor(year, levels=seq(from=2020, to=1990, by=-5), ordered=TRUE)) %>%
#   dplyr::rename("Model"=variable) %>%
#   ggplot() + 
#   #geom_point(aes(factor(year), value, group=Model, col=Model), position=position_dodge(width=0.5), size=3.5) + 
#   geom_point(aes(factor(year), value, group=Model, pch=Model), position=position_dodge(width=0.6), size=4, col="skyblue4", fill="skyblue4") + 
#   ylab("R-squared") + xlab("Year")  + 
#   theme_classic() +
#   theme(legend.position=c(0.33, 0.86)) + 
#   #scale_color_viridis_d(begin=0.2, end=0.7, guide=guide_legend(reverse=TRUE)) + 
#   scale_shape_discrete(guide=guide_legend(reverse=TRUE)) + 
#   theme(legend.title=element_text(size=10),
#         panel.grid.major = element_line(color="grey92"),
#         plot.title=element_text(hjust=0.5, size=13),
#         strip.text = element_text(size=11),
#         axis.text = element_text(size=10.5),
#         axis.title = element_text(size=11),
#         legend.text = element_text(size=9.5),
#         legend.background = element_blank()) + 
#   coord_flip()
  
p5 = results %>%
  dplyr::select(year, R2_withoutSR, R2_withSR) %>%
  distinct() %>%
  reshape2::melt(id.vars=1) %>%
  dplyr::mutate(variable = as.vector(variable),
                variable = replace(variable, variable == "R2_withoutSR", "logCites"),
                variable = replace(variable, variable == "R2_withSR", "logCites + logSR"),
                variable = factor(variable, levels=c("logCites", "logCites + logSR"), ordered=TRUE),
                year = factor(year, levels=seq(from=1990, to=2020, by=5), ordered=TRUE)) %>%
  dplyr::rename("Model"=variable) %>%
  ggplot() +
  #geom_point(aes(factor(year), value, group=Model, col=Model), position=position_dodge(width=0.5), size=3.5) +
  geom_point(aes(factor(year), value, group=Model, pch=Model), position=position_dodge(width=0.6), size=3.2, col="skyblue4", fill="skyblue4") +
  ylab("R-squared") + xlab("Year")  +
  theme_classic() +
  theme(legend.position=c(0.3, 0.25)) +
  #scale_color_viridis_d(begin=0.2, end=0.7, guide=guide_legend(reverse=TRUE)) +
  scale_shape_discrete(guide=guide_legend(reverse=TRUE)) +
  theme(legend.title=element_text(size=10),
        panel.grid.major = element_line(color="grey92"),
        plot.title=element_text(hjust=0.5, size=13),
        strip.text = element_text(size=11),
        axis.text = element_text(size=10.5),
        axis.title = element_text(size=11),
        legend.text = element_text(size=9.5),
        legend.background = element_blank()) 

pc2 = gridExtra::grid.arrange(p4, p5, ncol=2, widths=c(0.9, 0.5))
#ggsave(pc2, file="./output/figures_2021/SI_VirRichbySR.png", device="png", units="in", width=10.5, height=4, dpi=600)


# combine
pc_full = gridExtra::grid.arrange(pc1, pc2, nrow=2, heights=c(0.7, 0.9))
pc_full = ggpubr::as_ggplot(pc_full)  +
  cowplot::draw_plot_label(label = c("a", "b", "c", "d", "e"), 
                           fontface = "bold", size = 23, 
                           x = c(0.01, 0.425, 0.745, 0.08, 0.67), y = c(1, 1, 1, 0.58, 0.58))
ggsave(pc_full, file="./output/figures_2021/SI_VirRichbySR.png", device="png", units="in", width=10.5, height=6.8, dpi=600)

