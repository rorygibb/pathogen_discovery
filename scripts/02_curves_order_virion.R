

# ===================== Run discovery curves at Order-level ====================

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
#virion = vroom::vroom("./data/virion/pre_push/Virion.csv.gz")
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
# vir$HostOrder[ vir$HostOrder == "cetartiodactyla" ] = "artiodactyla"
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



# --------------- specify data subset: include or exclude ICTV ratified viruses --------------

# ictv_filter = FALSE
# 
# if(ictv_filter){ 
#   vir = vir[ vir$ICTVRatified == TRUE, ]
#   ictv_flag = "ictvratified"
# } else{
#     ictv_flag = "allviruses_genus"
#     vir$Virus = vir$VirusGenus
#   }



# -------------- build discovery curves -----------------

# total viral richness by order
tr = vir %>%
  group_by(HostOrder) %>%
  dplyr::summarise(VRichness = n_distinct(Virus))


# 1. all species including both wild and domestic

# unique associations by order and year
v1 = vir %>%
  dplyr::group_by(HostOrder, Virus) %>%
  dplyr::summarise(Database = paste(unique(Database), collapse=", "),
                   HostSpecies = paste(unique(Host), collapse=", "),
                   NumRecords = length(Year),
                   YearEarliest = min(Year, na.rm=TRUE),
                   YearLatest = max(Year, na.rm=TRUE)) %>%
  left_join(tr) %>%
  dplyr::arrange(desc(VRichness), HostOrder, YearEarliest)

# create cumulative curves of all years
curves_all = expand.grid(unique(v1$HostOrder), 1930:endyear) %>%
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


# 2. split out by wild/domestic

# unique associations by order, year
v2 = vir %>%
  dplyr::filter(!is.na(Year)) %>%
  dplyr::group_by(HostOrder, Domestic, Virus) %>%
  dplyr::summarise(Database = paste(unique(Database), collapse=", "),
                   HostSpecies = paste(unique(Host), collapse=", "),
                   NumRecords = length(Year),
                   YearEarliest = min(Year, na.rm=TRUE),
                   YearLatest = max(Year, na.rm=TRUE)) %>%
  left_join(tr) %>%
  dplyr::arrange(desc(VRichness), HostOrder, YearEarliest) 

# create cumulative curves of all years
curvesw = expand.grid(unique(v2$HostOrder), unique(v2$Domestic), 1930:endyear) %>%
  dplyr::rename("HostOrder" = 1, "Domestic" = 2, "YearEarliest" = 3) %>%
  left_join(v2[ , c("HostOrder", "Domestic", "YearEarliest", "Virus")]) %>%
  dplyr::group_by(HostOrder, Domestic, YearEarliest) %>%
  dplyr::summarise(Discovered = sum(!is.na(Virus)),
                   Virus = paste(unique(Virus), collapse=", ")) %>%
  left_join(v2[ !duplicated(v2$HostOrder), c("HostOrder", "VRichness") ]) %>%
  dplyr::arrange(desc(VRichness), HostOrder, Domestic, YearEarliest) %>%
  dplyr::group_by(HostOrder, Domestic) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered)) %>%
  dplyr::rename("Year" = YearEarliest)

# separate out into domestic and wild
curves_wild = curvesw[ curvesw$Domestic == FALSE, ]
curves_dom = curvesw[ curvesw$Domestic == FALSE, ]




# ==================== Plot cumulative virus discovery curves for all Orders, split out by domestic/wild ==================================

# set up data for plotting 
plot_data = curvesw
plot_data$facet = Hmisc::capitalize(paste(plot_data$HostOrder, " (", plot_data$VRichness, ")", sep=""))
fac_order = plot_data[ !duplicated(plot_data$HostOrder), c("HostOrder", "facet", "VRichness")] %>% dplyr::arrange(desc(VRichness))
plot_data$facet = factor(plot_data$facet, levels=fac_order$facet, ordered=TRUE)
plot_data$DW = ifelse(plot_data$Domestic == TRUE, "Domestic", "Wild")

# remove domestic curves from groups with no domesitcs
any_dom_orders = unique(vir$HostOrder[ vir$Domestic == TRUE])
plot_data = plot_data[ -which(plot_data$Domestic == TRUE & !plot_data$HostOrder %in% any_dom_orders), ]

# plot for MS / SI 
px = ggplot(plot_data) + 
  geom_line(aes(Year, VirusCumulative, group=DW, col=DW), size=1) + 
  #facet_wrap(~facet, scales="free_y") +
  lemon::facet_rep_wrap(~facet, scales="free_y") +
  theme_minimal() +
  scale_color_viridis_d(option="viridis", begin=0.25, end=0.7, direction=-1) +
  scale_x_continuous(breaks=seq(1950, 2010, by=30), labels=seq(1950, 2010, by=30)) +
  theme(strip.background = element_blank(), 
        panel.grid = element_blank(),
        strip.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.title = element_blank(), 
        legend.position = "bottom", 
        legend.text = element_text(size=13),
        axis.text = element_text(size=10.5),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.ticks = element_line()) + 
  ylab("Cumulative viral richness")
ggsave(px, file=paste("./output/figures_2021/MS_MammalOrders_CumulativeCurves_VIRION_", ictv_flag, ".png", sep=""), device="png", units="in", width=10, height=6.5, dpi=300)




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





# =================== fit GAMs and estimate curve derivative for specified Orders ============================

# data frame for results and specified orders
result = data.frame()

# specify data for inclusion
# model = "all"
# model = "domestic"
model = "wild"

if(model=="all"){
  datax = curves_all
  modname = "allspecies"
  orders = c("primates", "lagomorpha", "perissodactyla", "rodentia", "carnivora", "artiodactyla", "chiroptera", "eulipotyphla")
  fac_order = c("artiodactyla", "rodentia", "chiroptera", "primates", "carnivora", "perissodactyla", "eulipotyphla", "lagomorpha")
}
if(model=="domestic"){
  datax = curves_dom
  modname = "domestic"
  orders = c("Lagomorpha", "Perissodactyla", "Rodentia", "Carnivora", "Cetartiodactyla")
  fac_order = c("Cetartiodactyla", "Rodentia", "Carnivora", "Perissodactyla", "Lagomorpha")
}
if(model == "wild"){
  datax = curves_wild
  modname = "wild"
  # orders = c("primates", "lagomorpha", "rodentia", "carnivora", "artiodactyla", "chiroptera")
  # fac_order = c("Cetartiodactyla", "Rodentia", "Primates", "Carnivora", "Chiroptera", "Lagomorpha")
  orders = c("primates", "lagomorpha", "perissodactyla", "rodentia", "carnivora", "artiodactyla", "chiroptera", "eulipotyphla")
  fac_order = c("artiodactyla", "rodentia", "chiroptera", "primates", "carnivora", "perissodactyla", "eulipotyphla", "lagomorpha")
}

# for each order 
for(i in 1:length(orders)){
  
  # data (either wild only, or all)
  ord = orders[i]
  dd = datax %>%
    dplyr::filter(HostOrder == ord)
  
  # specify time threshold for inclusion of zeros prior to the first virus discovered in that taxa
  # either start from the first year of virus discovery, or n years beforehand, or can exclude this and run from 0
  if(model == "domestic"){
    time_thresh = 10
    first_year = min(dd$Year[ dd$Discovered > 0 ])
    dd = dd[ dd$Year >= first_year - time_thresh, ]
  }
  
  # GAM fit with spline on Year and using ML, Poisson likelihood
  # bats highly overdispersed so use nb() likelihood
  m1 = mgcv::gamm(Discovered ~ s(Year), family="poisson", data=dd, method="REML")
  m1 = m1$gam
  if(dd$HostOrder[1]%in%c("chiroptera", "rodentia", "primates")){  m1 = mgcv::gam(Discovered ~ s(Year), family=nb(), data=dd, method="REML") }

  # acf(resid(m1$lme, type = "normalized"))
  # pacf(resid(m1$lme, type = "normalized"))
  # hist(resid(m1$lme, type = "normalized"), 15)
  
  # fitted spline and pointwise confidence intervals, transform to natural scale
  # n.b. simultaneous intervals can be obtained via predictive sim https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
  #preds = data.frame(Year = dd$Year)
  preds = data.frame(Year = seq(min(dd$Year), max(dd$Year), length.out=100))
  preds = cbind(preds, predict(m1, preds, type = "link", se.fit = TRUE))
  preds$upper = exp(preds$fit + (1.96*preds$se.fit))
  preds$lower = exp(preds$fit - (1.96*preds$se.fit))
  preds$fitted = exp(preds$fit)
  
  # use Deriv function to estimate the first derivative of split at each year
  m1.d = Deriv(m1, n=length(preds$Year))
  plot(m1.d, sizer=TRUE)
  
  # estimate CIs on first derivative and identify areas where CIs indicate signifiance of trend
  CI = confint(m1.d, alpha = 0.01)
  S = signifD(preds$fit, m1.d$Year$deriv, CI$Year$upper, CI$Year$lower, eval = 0)
  preds$sig_incr = !is.na(S$incr)
  preds$sig_decr = !is.na(S$decr)
  preds$sig = !is.na(S$incr) | !is.na(S$decr)
  
  preds$Order = ord
  result = rbind(result, preds)
}

# ======================= visualise fitted curves ===========================

# plot
r2 = result
r2$Order = factor(Hmisc::capitalize(r2$Order), levels=Hmisc::capitalize(fac_order), ordered=TRUE)
r2$signif_col = NA
r2$signif_col[ r2$sig_incr == TRUE ] = "Increase"
r2$signif_col[ r2$sig_incr == FALSE ] = "Decrease"

#raw_data = curvesw[ curvesw$HostOrder %in% r2$Order & curvesw$Domestic == FALSE & curvesw$Year <= 2010, ]
raw_data = datax[ Hmisc::capitalize(datax$HostOrder) %in% r2$Order, ]
raw_data$Order = Hmisc::capitalize(raw_data$HostOrder)
raw_data$Order = factor(raw_data$Order, levels=Hmisc::capitalize(fac_order), ordered=TRUE)

curve_plot = ggplot() + 
  #geom_ribbon(data=r2[ -which(r2$Year < 1960 & r2$Order %in% c("Perissodactyla")), ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.3, fill="grey50", col=NA, size=0.05) +
  geom_ribbon(data=r2, aes(x=Year, ymin=lower, ymax=upper), alpha=0.2, fill="skyblue4", col=NA, size=0.05) +
  geom_point(data=raw_data, aes(x=Year, y=Discovered), col="grey55", size=0.3) +
  geom_line(data=r2, aes(x=Year, y=fitted, group=Order, col=sig_incr), size=1.2) +
  geom_line(data=r2[ r2$sig_decr == TRUE, ], aes(x=Year, y=fitted, group=Order), color="yellow", size=1.2) +
  theme_classic() +
  #scale_color_manual(values=c("FALSE"="grey70", "TRUE"=viridis::viridis(200)[150])) +
  scale_color_viridis_d(begin=0.3, end=0.75) +
  lemon::facet_rep_wrap(~Order, ncol=2, scales="free_y", strip.position = "top") +
  #facet_wrap(~Order, ncol=1, scales="free_y", strip.position = "right") +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text = element_text(size=12),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title=element_text(size=13.5)) +
  ylab("Virus discovery rate (viruses per year)") +
  scale_x_continuous(breaks=seq(1930, 2020, by=20), seq(1930, 2020, by=20), name="Year")
curve_plot

ggsave(curve_plot, file=paste("./output/figures_2021/MS_Figure1_OrderGAMs_", ictv_flag,".png", sep=""), device="png", units="in", width=8, height=8, dpi=600, scale=0.9)

