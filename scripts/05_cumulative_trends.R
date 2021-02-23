
# ====================== Overall cumulative trends in species-level viral discovery and effort =====================

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

# create cumulative curves of all years
curves = expand.grid(unique(dd$Host), 1930:2010) %>%
  dplyr::rename("Host" = 1, "Year" = 2) %>%
  left_join(dd[ , c("Host", "Year", "Virus")]) %>%
  dplyr::group_by(Host, Year) %>%
  dplyr::summarise(Discovered = sum(!is.na(Virus)),
                   Virus = paste(unique(Virus), collapse=", ")) %>%
  left_join(dd[ !duplicated(dd$Host), c("Host", "HostOrder", "HostFamily", "VRichness", "Domestic") ]) %>%
  dplyr::arrange(desc(VRichness), Host, Year) %>%
  dplyr::group_by(Host) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered))




# =================== set-up for fitting GAMs ============================

# GAMs with smoothing effect of year fitted in mgcv
# additional functions for calculating spline derivatives accessed from Gavin Simpson's articles
# https://fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/

# deriv function
dest1 = "C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/scripts/mgcv_simpson/derivFun.R"
#download.file("https://github.com/gavinsimpson/random_code/raw/master/derivFun.R", destfile=dest1, method = "wget")
source(dest1)

dest2 = "C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/scripts/mgcv_simpson/tsDiagGamm.R"
#download.file("https://github.com/gavinsimpson/random_code/raw/master/tsDiagGamm.R", dest2, method = "wget")
source(dest2)



# ================== fit overall cumulative curve of virus richness at the wild species-level ===================

# GAM across all species fitted using REML
dd = curves[ curves$Domestic == FALSE & curves$Year < 2011, ]
m_sp = mgcv::gamm(VirusCumulative ~ s(Year), family=nb(link="log"), data=dd, method="REML")

# fitted spline and pointwise confidence intervals, transform to natural scale
# n.b. simultaneous intervals can be obtained via predictive sim https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
#preds = data.frame(Year = dd$Year)
preds = data.frame(Year = seq(min(dd$Year), max(dd$Year), length.out=n_distinct(dd$Year)))
preds = cbind(preds, predict(m_sp$gam, preds, type = "link", se.fit = TRUE))
preds$upper = exp(preds$fit + (1.96*preds$se.fit))
preds$lower = exp(preds$fit - (1.96*preds$se.fit))
preds$fitted = exp(preds$fit)

#
plot_cumul = ggplot() +
  #geom_ribbon(data=r2[ -which(r2$Year < 1960 & r2$Order %in% c("Perissodactyla")), ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.3, fill="grey50", col=NA, size=0.05) +
  geom_ribbon(data=preds, aes(x=Year, ymin=lower, ymax=upper), alpha=0.3, fill="skyblue4", col=NA, size=0.05) +
  geom_line(data=preds, aes(x=Year, y=fitted), size=0.75, col="skyblue4") +
  #geom_line(data=preds[ preds$sig_decr == TRUE, ], aes(x=Year, y=fitted), col="red", size=1.2) +
  theme_classic() +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  #facet_wrap(~Order, ncol=1, scales="free_y", strip.position = "right") +
  theme(strip.background = element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust=0.5),
        strip.text = element_text(size=13),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title=element_text(size=13)) +
  ylab("Mean viral richness (all species)") +
  ggtitle("Viral richness") +
  scale_x_continuous(breaks=seq(1930, 2010, by=20), seq(1940, 2000, by=20), name="Year")




# =================== overall cumulative publications at the species level ===================

# read effort
hx = read.csv("./output/host_effort/PubMed_HostsEffort_PerYear_19302020.csv") %>%
  dplyr::select(-Note) %>%
  dplyr::mutate(Type = "All publications")
hy = read.csv("./output/host_effort/PubMed_HostsEffort_PerYear_VirusRelated_19302020.csv") %>%
  dplyr::select(-Note, -VirusRelated) %>%
  dplyr::mutate(Type = "Virus-related")
effort = rbind(hx, hy) %>%
  dplyr::left_join(curves[ !duplicated(curves$Host), c("Host", "HostOrder", "HostFamily")])

# remove domestics
dom = read.csv("./data/clover/domestic_status/HostLookup_Domestic.csv")
effort$Domestic = ifelse(effort$Host %in% dom$Host, TRUE, FALSE)

# create publication curves (virus related)
pcurves = expand.grid(unique(effort$Host), 1930:2010) %>%
  dplyr::rename("Host" = 1, "Year" = 2) %>%
  left_join(effort[ effort$Type == "Virus-related", c("Host", "Year", "NumPubs")])
pcurves$NumPubs[ is.na(pcurves$NumPubs) ] = 0
pcurves = pcurves %>%
  dplyr::arrange(Host, Year) %>%
  dplyr::group_by(Host) %>%
  dplyr::mutate(CumPubs = cumsum(NumPubs)) 
pcurves = left_join(pcurves, effort[ !duplicated(effort$Host), c("Host", "HostOrder", "HostFamily", "Domestic")])
  
# GAM across all species fitted using REML
dd = pcurves[ pcurves$Domestic == FALSE, ]
m_ef = mgcv::gamm(CumPubs ~ s(Year), family=nb(link="log"), data=dd, method="REML")

# fitted spline and pointwise confidence intervals, transform to natural scale
# n.b. simultaneous intervals can be obtained via predictive sim https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
#preds = data.frame(Year = dd$Year)
preds = data.frame(Year = seq(min(dd$Year), max(dd$Year), length.out=n_distinct(dd$Year)))
preds = cbind(preds, predict(m_ef$gam, preds, type = "link", se.fit = TRUE))
preds$upper = exp(preds$fit + (1.96*preds$se.fit))
preds$lower = exp(preds$fit - (1.96*preds$se.fit))
preds$fitted = exp(preds$fit)

#
plot_eff = ggplot() +
  geom_ribbon(data=preds, aes(x=Year, ymin=lower, ymax=upper), alpha=0.3, fill="skyblue4", col=NA, size=0.05) +
  geom_line(data=preds, aes(x=Year, y=fitted), size=0.75, col="skyblue4") +
  #geom_line(data=preds[ preds$sig_decr == TRUE, ], aes(x=Year, y=fitted), col="red", size=1.2) +
  theme_classic() +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  #facet_wrap(~Order, ncol=1, scales="free_y", strip.position = "right") +
  theme(strip.background = element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust=0.5),
        strip.text = element_text(size=13),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title=element_text(size=13)) +
  ylab("Mean virus-related publications (all species)") +
  ggtitle("Virus-related publications") +
  scale_x_continuous(breaks=seq(1930, 2010, by=20), seq(1940, 2000, by=20), name="Year")

# creat combined plot
p_comb = gridExtra::grid.arrange(plot_cumul, plot_eff, nrow=1)
ggsave(p_comb, file="./output/figures/SI_CumulativeCurves.png", dpi=600, width=8, height=3.8, units="in", device="png")





# ================== Publication trends at the Order level =====================

# # read effort
# hx = read.csv("./output/host_effort/PubMed_HostsEffort_PerYear_19302020.csv") %>%
#   dplyr::select(-Note) %>%
#   dplyr::mutate(Type = "All publications")
# hy = read.csv("./output/host_effort/PubMed_HostsEffort_PerYear_VirusRelated_19302020.csv") %>%
#   dplyr::select(-Note, -VirusRelated) %>%
#   dplyr::mutate(Type = "Virus-related")
# 
# effort = rbind(hx, hy) %>%
#   dplyr::left_join(curves[ !duplicated(curves$Host), c("Host", "HostOrder", "HostFamily")])
# 
# # remove domestics
# dom = read.csv("./data/clover/domestic_status/HostLookup_Domestic.csv")
# effort$Domestic = ifelse(effort$Host %in% dom$Host, TRUE, FALSE)
# 
# # plot effort over time
# eff_all = effort[ effort$Year < 2021 & effort$Domestic == FALSE, ] %>%
#   dplyr::group_by(Type, Year) %>%
#   dplyr::summarise(Publications = sum(NumPubs))
# ggplot(eff_all) +
#   geom_point(aes(Year, Publications)) +
#   facet_wrap(~Type, nrow=1, scales="free_y")
# 
# # ef2 = rbind(eff_all, data.frame(Year = 2020.01, Publications=0, Type = unique(eff_all$Type)))
# # ggplot(ef2) +
# #   geom_polygon(aes(Year, Publications, fill=Type)) +
# #   theme_minimal()
# 
# 
# # effort by Order
# eff_all = effort[ effort$Year < 2016 & effort$Domestic == FALSE, ] %>%
#   dplyr::filter(NumPubs > 0) %>%
#   dplyr::group_by(Type, HostOrder, Year) %>%
#   dplyr::summarise(Publications = sum(NumPubs),
#                    `Number of species` = n_distinct(Host),
#                    `Number of families` = n_distinct(HostFamily))
# fac_order = eff_all[ eff_all$Type == "Virus-related", ] %>%
#   dplyr::group_by(HostOrder) %>%
#   dplyr::summarise(TotalPubs = sum(Publications)) %>%
#   dplyr::arrange(desc(TotalPubs))
# fac_order$facet = paste(fac_order$HostOrder, " (", fac_order$TotalPubs, ")", sep="")
# eff_all = left_join(eff_all, fac_order[ , c("HostOrder", "facet") ])
# eff_all$facet = factor(eff_all$facet, levels=fac_order$facet, ordered=TRUE)
# eff_all = eff_all[ eff_all$HostOrder %in% fac_order$HostOrder[ fac_order$TotalPubs > 16 ], ]
# 
# 
# # plot virus related effort across all orders
# effort_by_order = ggplot(eff_all[ eff_all$Type == "Virus-related", ]) +
#   geom_point(aes(Year, Publications), pch=21, fill="coral2") +
#   geom_line(stat="smooth", aes(Year, Publications, group=facet), method="gam", se=FALSE, alpha=0.6, size=0.5, col="blue") +
#   lemon::facet_rep_wrap(~facet, scales="free_y") +
#   theme_classic() +
#   theme(strip.background = element_blank()) +
#   ylab("Virus-related publications") +
#   theme(strip.background = element_blank(),
#         legend.position="none",
#         strip.text = element_text(size=12),
#         axis.text.y = element_text(size=11),
#         axis.text.x = element_text(size=12),
#         axis.title=element_text(size=13.5))
# 
# ggsave(effort_by_order, file="./output/figures/MS_SIFigure_PublicationEffortByOrder.png", device="png", units="in", width=9, height=6, dpi=600, scale=0.95)
# 
# 
# # size by number of species
# effort_by_order2 = ggplot(eff_all[ eff_all$Type == "Virus-related", ]) +
#   geom_point(aes(Year, Publications, size=`Number of species`), pch=21, fill="coral2", alpha=0.5) +
#   geom_line(stat="smooth", aes(Year, Publications, group=facet), method="gam", se=FALSE, alpha=0.6, size=0.6, col="blue") +
#   lemon::facet_rep_wrap(~facet, scales="free_y") +
#   theme_classic() +
#   theme(strip.background = element_blank()) +
#   scale_size_continuous(range=c(0.2, 4)) +
#   ylab("Virus-related publications") +
#   theme(strip.background = element_blank(),
#         legend.position="bottom",
#         strip.text = element_text(size=12),
#         axis.text.y = element_text(size=11),
#         axis.text.x = element_text(size=12),
#         axis.title=element_text(size=13.5))
# 
# ggsave(effort_by_order2, file="./output/figures/MS_SIFigure_PublicationEffortByOrder_NSpecies.png", device="png", units="in", width=9, height=6.8, dpi=600, scale=0.95)
# 
# 
# effort_by_order3 = ggplot(eff_all[ eff_all$Type == "Virus-related", ]) +
#   geom_point(aes(Year, Publications, size=`Number of families`), pch=21, fill="coral2", alpha=0.5) +
#   geom_line(stat="smooth", aes(Year, Publications, group=facet), method="gam", se=FALSE, alpha=0.6, size=0.6, col="blue") +
#   lemon::facet_rep_wrap(~facet, scales="free_y") +
#   theme_classic() +
#   theme(strip.background = element_blank()) +
#   scale_size_continuous(range=c(0.2, 4)) +
#   ylab("Virus-related publications") +
#   theme(strip.background = element_blank(),
#         legend.position="bottom",
#         strip.text = element_text(size=12),
#         axis.text.y = element_text(size=11),
#         axis.text.x = element_text(size=12),
#         axis.title=element_text(size=13.5))
# 
# 
# ##method.args=list(family="poisson"))



# ================= Taxonomic breadth of discovery ====================

# 
tb1 = curves %>%
  dplyr::filter(Discovered > 0) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(HostRange_Species = n_distinct(Host),
                   HostRange_Family = n_distinct(HostFamily))
ggplot(tb1) + 
  geom_point(aes(Year, HostRange_Family), pch=21, fill=NA) + 
  theme_classic() +
  geom_line(stat="smooth", aes(Year, HostRange_Family), method="gam", method.args=list(family="poisson"))

ggplot(tb1) + 
  geom_point(aes(Year, HostRange_Family), pch=21, fill=NA) + 
  theme_classic() +
  geom_smooth(aes(Year, HostRange_Family), method="gam", method.args=list(family="poisson"))
ggplot(tb1) + 
  geom_point(aes(Year, HostRange_Species), pch=21, fill=NA) + 
  theme_classic() +
  geom_smooth(aes(Year, HostRange_Species), method="glm", method.args=list(family="poisson"))


tb2 = curves[ curves$HostOrder %in% c("Primates", "Rodentia", "Carnivora", "Chiroptera", "Cetartiodactyla", "Lagomorpha"), ] %>%
  dplyr::filter(Discovered > 0) %>%
  dplyr::group_by(HostOrder, Year) %>%
  dplyr::summarise(HostRange_Species = n_distinct(Host),
                   HostRange_Family = n_distinct(HostFamily))

ggplot(tb2) + 
  geom_point(aes(Year, HostRange_Species), pch=21, fill=NA) + 
  theme_classic() +
  geom_smooth(aes(Year, HostRange_Species), method="gam", method.args=list(family="poisson")) + 
  facet_wrap(~HostOrder, scales="free_y")

ggplot(tb2) + 
  geom_point(aes(Year, HostRange_Family), pch=21, fill=NA) + 
  theme_classic() +
  geom_smooth(aes(Year, HostRange_Family), method="gam", method.args=list(family="poisson")) + 
  facet_wrap(~HostOrder, scales="free_y")

