

# ========================================= ======================================================

# Gibb et al., "Mammal virus diversity estimates are unstable due to accelerating discovery effort"
# Script 3: Fit species-level viral discovery GAMs for the top 50 most virus-rich mammal species

# ================================================================================================




# fit GAM curves across time epochs to examine changes in trends 

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon", "mgcv", "vroom")



# ---------------- build VIRION database ----------------

ictv_flag = "ictvpredict"

# domestic species to label
domestic = read.csv("./data/clovert/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)

# additional laboratory species to exlude
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

# total viral richness by species
tr = vir %>%
  group_by(Host) %>%
  dplyr::summarise(VRichness = n_distinct(Virus))



# ---------------------- Calculate discovery information for species -------------------------

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

# create cumulative curves of all years
curves = expand.grid(unique(dd$Host), 1930:2018) %>%
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



# ================== fit overall cumulative curve at the wild species-level ===================

# GAM across all species fitted using REML
dd = curves[ curves$Domestic == FALSE & curves$Year <= 2018, ]
m_sp = mgcv::gamm(VirusCumulative ~ s(Year), family=nb(), data=dd, method="REML")

# fitted spline and pointwise confidence intervals, transform to natural scale
# n.b. simultaneous intervals can be obtained via predictive sim https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
#preds = data.frame(Year = dd$Year)
preds = data.frame(Year = seq(min(dd$Year), max(dd$Year), length.out=n_distinct(dd$Year)))
preds = cbind(preds, predict(m_sp$gam, preds, type = "link", se.fit = TRUE))
preds$upper = exp(preds$fit + (1.96*preds$se.fit))
preds$lower = exp(preds$fit - (1.96*preds$se.fit))
preds$fitted = exp(preds$fit)

#
ggplot() +
  #geom_ribbon(data=r2[ -which(r2$Year < 1960 & r2$Order %in% c("Perissodactyla")), ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.3, fill="grey50", col=NA, size=0.05) +
  geom_ribbon(data=preds, aes(x=Year, ymin=lower, ymax=upper), alpha=0.25, fill="skyblue4", col=NA, size=0.05) +
  geom_line(data=preds, aes(x=Year, y=fitted), size=0.8, col="skyblue4") +
  #geom_line(data=preds[ preds$sig_decr == TRUE, ], aes(x=Year, y=fitted), col="red", size=1.2) +
  theme_classic() +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  #facet_wrap(~Order, ncol=1, scales="free_y", strip.position = "right") +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text = element_text(size=13),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title=element_text(size=16)) +
  ylab("Cumulative viral richness") +
  scale_x_continuous(breaks=seq(1930, 2010, by=20), seq(1940, 2000, by=20), name="Year")






# =================== fit GAMs and estimate curve derivative for top n species ============================

# data frame for results and specified orders
result = data.frame()

# get first n species
n = 50
cc = curves[! curves$Host %in% c("triaenops persicus", "hipposideros ruber"), ] # model for h ruber does not converge; too recent uptick
spp = unique(cc$Host)[1:n]
datax = curves

spp2 = spp

# for each order 
for(i in 1:length(spp2)){
  
  # data (either wild only, or all) with cutoff of 2018
  spx = spp2[i]
  dd = datax[ datax$Host == spx & datax$Year <= 2018, ]
  
  # specify time threshold for inclusion of zeros prior to the first virus discovered in that taxa
  # either start from the first year of virus discovery, or n years beforehand, or can exclude this and run from 0
  time_thresh = 0
  first_year = min(dd$Year[ dd$Discovered > 0 ])
  dd = dd[ dd$Year >= first_year - time_thresh, ]
  
  # GAM fit with spline on Year and using ML, Poisson likelihood
  if(dd$HostOrder == "chiroptera"){
    m1 = mgcv::gamm(Discovered ~ s(Year), family=nb(), data=dd, method="REML")
  } else{
    m1 = mgcv::gamm(Discovered ~ s(Year), family="poisson", data=dd, method="REML")
  }
  #plot(m1$gam)
  
  # acf(resid(m1$lme, type = "normalized"))
  # pacf(resid(m1$lme, type = "normalized"))
  # hist(resid(m1$lme, type = "normalized"), 15)
  
  # fitted spline and pointwise confidence intervals, transform to natural scale
  # n.b. simultaneous intervals can be obtained via predictive sim https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
  #preds = data.frame(Year = dd$Year)
  preds = data.frame(Year = seq(min(dd$Year), max(dd$Year), length.out=100))
  preds = cbind(preds, predict(m1$gam, preds, type = "link", se.fit = TRUE))
  preds$upper = exp(preds$fit + (1.96*preds$se.fit))
  preds$lower = exp(preds$fit - (1.96*preds$se.fit))
  preds$fitted = exp(preds$fit)
  
  # use Deriv function to estimate the first derivative of split at each year
  m1.d = Deriv(m1$gam, n=length(preds$Year))
  plot(m1.d, sizer=TRUE)
  
  # estimate CIs on first derivative and identify areas where CIs indicate signifiance of trend
  CI = confint(m1.d, alpha = 0.01)
  S = signifD(preds$fit, m1.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
  preds$sig_incr = !is.na(S$incr)
  preds$sig_decr = !is.na(S$decr)
  preds$sig = !is.na(S$incr) | !is.na(S$decr)
  
  preds$Host = spx
  preds$HostOrder = dd$HostOrder[1]
  preds$VRichness = dd$VRichness[1]
  preds$Domestic = dd$Domestic[1]
  result = rbind(result, preds)
}





# ======================= visualise fitted curves ===========================

# plot
r2 = result
r2$Host2 = Hmisc::capitalize(r2$Host)
r2$Host2 = lapply(strsplit(r2$Host2, " "), function(x) paste(x, collapse="\n"))
r2$Host2 = paste(r2$Host2, " (", r2$VRichness, ")", sep="")
r2$Host2 = factor(r2$Host2, levels=unique(r2$Host2), ordered=TRUE)

r2$signif_col = NA
r2$signif_col[ r2$sig_incr == TRUE ] = "Increase";
r2$signif_col[ r2$sig_incr == FALSE ] = "Decrease"

#raw_data = curvesw[ curvesw$HostOrder %in% r2$Order & curvesw$Domestic == FALSE & curvesw$Year <= 2010, ]
raw_data = datax[ datax$Host %in% r2$Host & datax$Year <= 2018, ]
raw_data = left_join(raw_data, r2[ !duplicated(r2$Host), c("Host", "Host2")])
raw_data$Host2 = factor(raw_data$Host2, levels=unique(r2$Host2), ordered=TRUE)

# first 25 species
curve_plot = ggplot() + 
  geom_ribbon(data=r2[ r2$Host2 %in% unique(r2$Host2)[1:25], ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.25, fill="skyblue4", col=NA, size=0.05) +
  geom_point(data=raw_data[ raw_data$Host2 %in% unique(r2$Host2)[1:25], ], aes(x=Year, y=Discovered), col="grey55", size=0.3) +
  geom_line(data=r2[ r2$Host2 %in% unique(r2$Host2)[1:25], ], aes(x=Year, y=fitted, group=Host2, col=sig_incr), size=1.2) +
  geom_line(data=r2[ r2$Host2 %in% unique(r2$Host2)[1:25] & r2$sig_decr == TRUE, ], aes(x=Year, y=fitted, group=Host2), col="yellow", size=1.2) +
  theme_classic() +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  lemon::facet_rep_wrap(~Host2, ncol=5, nrow=5, scales="free_y") +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text = element_text(size=13),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title=element_text(size=16)) +
  ylab("Virus discovery rate (viruses per year)") +
  scale_x_continuous(breaks=seq(1940, 2000, by=20), seq(1940, 2000, by=20), name="Year")
ggsave(curve_plot, file="./output/figures_2021/SI_Figure_SpeciesGAMs_A.png", device="png", units="in", width=11, height=10, dpi=300)

# next 25 species
curve_plot2 = ggplot() + 
  geom_ribbon(data=r2[ r2$Host2 %in% unique(r2$Host2)[26:50], ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.25, fill="skyblue4", col=NA, size=0.05) +
  geom_point(data=raw_data[ raw_data$Host2 %in% unique(r2$Host2)[26:50], ], aes(x=Year, y=Discovered), col="grey55", size=0.3) +
  geom_line(data=r2[ r2$Host2 %in% unique(r2$Host2)[26:50], ], aes(x=Year, y=fitted, group=Host2, col=sig_incr), size=1.2) +
  geom_line(data=r2[ r2$Host2 %in% unique(r2$Host2)[26:50] & r2$sig_decr == TRUE, ], aes(x=Year, y=fitted, group=Host2), col="yellow", size=1.2) +
  theme_classic() +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  lemon::facet_rep_wrap(~Host2, ncol=5, nrow=5, scales="free_y") +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text = element_text(size=13),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title=element_text(size=16)) +
  ylab("Virus discovery rate (viruses per year)") +
  scale_x_continuous(breaks=seq(1940, 2000, by=20), seq(1940, 2000, by=20), name="Year")
ggsave(curve_plot2, file="./output/figures_2021/SI_Figure_SpeciesGAMs_B.png", device="png", units="in", width=11, height=10, dpi=300)


