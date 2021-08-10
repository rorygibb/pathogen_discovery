


# ====================== Generate species viral discovery rate curves for mammals at the Order level =====================

# fit GAM curves across time epochs to examine changes in trends 

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon", "mgcv")

# domestic species to label
domestic = read.csv("./data/clovert/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)

# associations data and no humans
clover = read.csv("./data/virion/VirionPartialDiscovery.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Host != "homo sapiens") %>%
  dplyr::mutate(Domestic = ifelse(Host %in% tolower(domestic$Host), TRUE, FALSE)) %>%
  #dplyr::filter(DetectionMethod != "Antibodies") %>%
  dplyr::filter(DetectionMethod != "Not specified")

# create combined year column
clover$Year = clover$PublicationYear
clover$Year[ is.na(clover$Year) ] = clover$CollectionYear[ is.na(clover$Year) ]
clover$Year[ is.na(clover$Year) ] = clover$ReleaseYear[ is.na(clover$Year) ]
clover = clover[ !is.na(clover$Year), ]

# mammals only
clover = clover[ clover$HostClass == "mammalia", ]
clover$HostOrder[ clover$HostOrder == "artiodactyla" ] = "cetartiodactyla"

tr = clover %>%
  group_by(HostOrder) %>%
  dplyr::summarise(VGenusRichness = n_distinct(VirusGenus))

# # total viral richness by order
# tr1 = clover[ clover$Year < 2011, ] %>%
#   group_by(HostOrder) %>%
#   dplyr::summarise(VRichness = n_distinct(Virus),
#                    VGenusRichness = n_distinct(VirusGenus),
#                    HostRichness = n_distinct(Host)) %>%
#   reshape2::melt(id.vars=1) %>%
#   dplyr::mutate(Year = 2010)
# # 
# tr2 = clover %>%
#   group_by(HostOrder) %>%
#   dplyr::summarise(VRichness = n_distinct(Virus),
#                    VGenusRichness = n_distinct(VirusGenus),
#                    HostRichness = n_distinct(Host)) %>%
#   reshape2::melt(id.vars=1) %>%
#   dplyr::mutate(Year = 2021)
#  
# tr = rbind(tr1, tr2)
# ggplot(tr[ tr$HostOrder %in% tr$HostOrder[ tr$variable == "VRichness" & tr$value > 100 ], ]) + 
#   geom_point(aes(HostOrder, value, col=factor(Year)), position=position_dodge(width=0.5)) + 
#   facet_wrap(~variable, scales="free_y", ncol=1) + 
#   theme_minimal() + 
#   theme(axis.text.x = element_text(angle=90))

# tr$Percent_Change = (tr$VRichness_2021 / tr$VRichness_2010) * 100
#write.csv(tr, ".yikes.csv", row.names=FALSE)



# 1. all species including both wild and domestic

# unique associations by order and year
cl1 = clover %>%
  dplyr::group_by(HostOrder, VirusGenus) %>%
  dplyr::summarise(Database = paste(unique(Database), collapse=", "),
                   HostSpecies = paste(unique(Host), collapse=", "),
                   NumRecords = length(Year),
                   YearEarliest = min(Year, na.rm=TRUE),
                   YearLatest = max(Year, na.rm=TRUE)) %>%
  left_join(tr) %>%
  dplyr::arrange(desc(VGenusRichness), HostOrder, YearEarliest)

# create cumulative curves of all years
curves_all = expand.grid(unique(cl1$HostOrder), 1930:2016) %>%
  dplyr::rename("HostOrder" = 1, "YearEarliest" = 2) %>%
  left_join(cl1[ , c("HostOrder", "YearEarliest", "VirusGenus")]) %>%
  dplyr::group_by(HostOrder, YearEarliest) %>%
  dplyr::summarise(Discovered = sum(!is.na(VirusGenus)),
                   VirusGenus = paste(unique(VirusGenus), collapse=", ")) %>%
  left_join(cl1[ !duplicated(cl1$HostOrder), c("HostOrder", "VGenusRichness") ]) %>%
  dplyr::arrange(desc(VGenusRichness), HostOrder, YearEarliest) %>%
  dplyr::group_by(HostOrder) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered)) %>%
  dplyr::rename("Year" = YearEarliest)


# 2. split out by wild/domestic

# unique associations by order, year
cl2 = clover %>%
  dplyr::filter(!is.na(Year)) %>%
  dplyr::group_by(HostOrder, Domestic, VirusGenus) %>%
  dplyr::summarise(Database = paste(unique(Database), collapse=", "),
                   HostSpecies = paste(unique(Host), collapse=", "),
                   NumRecords = length(Year),
                   YearEarliest = min(Year, na.rm=TRUE),
                   YearLatest = max(Year, na.rm=TRUE)) %>%
  left_join(tr) %>%
  dplyr::arrange(desc(VGenusRichness), HostOrder, YearEarliest) 

# create cumulative curves of all years
curvesw = expand.grid(unique(cl2$HostOrder), unique(cl2$Domestic), 1930:2016) %>%
  dplyr::rename("HostOrder" = 1, "Domestic" = 2, "YearEarliest" = 3) %>%
  left_join(cl2[ , c("HostOrder", "Domestic", "YearEarliest", "VirusGenus")]) %>%
  dplyr::group_by(HostOrder, Domestic, YearEarliest) %>%
  dplyr::summarise(Discovered = sum(!is.na(VirusGenus)),
                   VirusGenus = paste(unique(VirusGenus), collapse=", ")) %>%
  left_join(cl2[ !duplicated(cl2$HostOrder), c("HostOrder", "VGenusRichness") ]) %>%
  dplyr::arrange(desc(VGenusRichness), HostOrder, Domestic, YearEarliest) %>%
  dplyr::group_by(HostOrder, Domestic) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered)) %>%
  dplyr::rename("Year" = YearEarliest)

# separate out into domestic and wild
curves_wild = curvesw[ curvesw$Domestic == FALSE, ]
curves_dom = curvesw[ curvesw$Domestic == FALSE, ]




# ==================== Plot cumulative virus discovery curves for all Orders, split out by domestic/wild ==================================

# set up data for plotting 
plot_data = curvesw
plot_data$facet = paste(plot_data$HostOrder, " (", plot_data$VGenusRichness, ")", sep="")
fac_order = plot_data[ !duplicated(plot_data$HostOrder), c("HostOrder", "facet", "VGenusRichness")] %>% dplyr::arrange(desc(VGenusRichness))
plot_data$facet = factor(plot_data$facet, levels=fac_order$facet, ordered=TRUE)
plot_data$DW = ifelse(plot_data$Domestic == TRUE, "Domestic", "Wild")

# remove domestic curves from groups with no domesitcs
any_dom_orders = unique(clover$HostOrder[ clover$Domestic == TRUE])
plot_data = plot_data[ -which(plot_data$Domestic == TRUE & !plot_data$HostOrder %in% any_dom_orders), ]

# plot for MS / SI up to 2010
px = ggplot(plot_data[ plot_data$Year <= 2020, ]) + 
  geom_line(aes(Year, VirusCumulative, group=DW, col=DW), size=1) + 
  #facet_wrap(~facet, scales="free_y") +
  lemon::facet_rep_wrap(~facet, scales="free_y") +
  theme_minimal() +
  scale_color_viridis_d(option="viridis", begin=0.25, end=0.7, direction=-1) +
  scale_x_continuous(breaks=seq(1960, 2020, by=30), labels=seq(1960, 2020, by=30)) +
  geom_vline(xintercept=2010, lty=2) +
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
  ylab("Cumulative virus genus richness")
ggsave(px, file="./output/figures/MS_MammalOrders_CumulativeCurves_VIRION_genus.png", device="png", units="in", width=10, height=8, dpi=300)






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
  orders = c("Primates", "Lagomorpha", "Perissodactyla", "Rodentia", "Carnivora", "Cetartiodactyla", "Chiroptera", "Didelphimorphia")
  fac_order = c("Cetartiodactyla", "Rodentia", "Primates", "Carnivora", "Chiroptera", "Perissodactyla", "Lagomorpha", "Didelphimorphia")
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
  orders = c("Primates", "Lagomorpha", "Rodentia", "Carnivora", "Cetartiodactyla", "Chiroptera")
  fac_order = c("Cetartiodactyla", "Rodentia", "Primates", "Carnivora", "Chiroptera", "Lagomorpha")
}

# for each order 
for(i in 1:length(orders)){
  
  # data (either wild only, or all) with cutoff of 2010
  ord = orders[i]
  dd = datax[ datax$HostOrder == ord & datax$Year <= 2010, ]
  
  # specify time threshold for inclusion of zeros prior to the first virus discovered in that taxa
  # either start from the first year of virus discovery, or n years beforehand, or can exclude this and run from 0
  if(model == "domestic"){
    time_thresh = 10
    first_year = min(dd$Year[ dd$Discovered > 0 ])
    dd = dd[ dd$Year >= first_year - time_thresh, ]
  }
  
  # GAM fit with spline on Year and using ML, Poisson likelihood
  m1 = mgcv::gamm(Discovered ~ s(Year), family="poisson", data=dd, method="REML")
  
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
  
  preds$Order = ord
  result = rbind(result, preds)
}

# ======================= visualise fitted curves ===========================

# plot
r2 = result
r2$Order = factor(r2$Order, levels=fac_order, ordered=TRUE)
r2$signif_col = NA
r2$signif_col[ r2$sig_incr == TRUE ] = "Increase"
r2$signif_col[ r2$sig_incr == FALSE ] = "Decrease"

#raw_data = curvesw[ curvesw$HostOrder %in% r2$Order & curvesw$Domestic == FALSE & curvesw$Year <= 2010, ]
raw_data = datax[ datax$HostOrder %in% r2$Order & datax$Year <= 2010, ]
raw_data$Order = raw_data$HostOrder
raw_data$Order = factor(raw_data$Order, levels=fac_order, ordered=TRUE)

curve_plot = ggplot() + 
  #geom_ribbon(data=r2[ -which(r2$Year < 1960 & r2$Order %in% c("Perissodactyla")), ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.3, fill="grey50", col=NA, size=0.05) +
  geom_ribbon(data=r2, aes(x=Year, ymin=lower, ymax=upper), alpha=0.25, fill="skyblue4", col=NA, size=0.05) +
  geom_point(data=raw_data, aes(x=Year, y=Discovered), col="grey55", size=0.3) +
  geom_line(data=r2, aes(x=Year, y=fitted, group=Order, col=sig_incr), size=1.2) +
  geom_line(data=r2[ r2$sig_decr == TRUE, ], aes(x=Year, y=fitted, group=Order), color="red", size=1.2) +
  theme_classic() +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  lemon::facet_rep_wrap(~Order, ncol=2, scales="free_y", strip.position = "top") +
  #facet_wrap(~Order, ncol=1, scales="free_y", strip.position = "right") +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text = element_text(size=12),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=12),
        axis.title=element_text(size=13.5)) +
  ylab("Virus discovery rate (viruses per year)") +
  scale_x_continuous(breaks=seq(1930, 2010, by=20), seq(1930, 2010, by=20), name="Year")
curve_plot

ggsave(curve_plot, file="./output/figures/MS_Figure1_OrderGAMs.png", device="png", units="in", width=8, height=6, dpi=600, scale=0.9)










# ### example
# library(mgcv)
# 
# dd = curvesw[ curvesw$HostOrder == "Chiroptera" & curvesw$Domestic == FALSE & curvesw$Year <= 2010, ]
# m1 = gam(Discovered ~ s(Year), family="poisson", data=dd, method="ML")
# #m2 = gam(Discovered ~ s(Year), data = dd, correlation = corARMA(form = ~ Year, p = 1))
# 
# #hist((dd$Discovered - m1$fitted.values) / sqrt(m1$fitted.values), 20)
# 
# # output = confint(m1, parm="s(Year)", level=0.95, transform=TRUE, type="confidence")
# # ggplot(output) + 
# #   geom_line(aes(Year, est)) +
# #   geom_ribbon(aes(Year, ymin=lower, ymax=upper), alpha=0.2, fill="skyblue4") +
# #   geom_point(data = dd, aes(Year, Discovered))
#   
# 
# plot.gam(m1, residuals = TRUE, pch = 19, cex = 0.75)
# summary(m1)
# hist(resid(m1), 20)
# 
# # calculates first derivative and identifies periods of statistically significant change in slope
# m1.d = Deriv(m1, n=length(dd$Year))
# plot(m1.d, sizer = TRUE, alpha = 0.01)
# 
# # plot
# preds = data.frame(Year = dd$Year)
# preds = cbind(preds, predict(m1, preds, type = "link", se.fit = TRUE))
# preds$upper = exp(preds$fit + (1.96*preds$se.fit))
# preds$lower = exp(preds$fit - (1.96*preds$se.fit))
# preds$fitted = exp(preds$fit)
# 
# # denote significant change
# CI <- confint(m1.d, alpha = 0.01)
# S <- signifD(preds$fit, m1.d$Year$deriv, CI$Year$upper, CI$Year$lower,
#              eval = 0)
# preds$sig = !is.na(S$incr) | !is.na(S$decr)
# preds$all = 1
# ggplot() +
#   geom_line(data=preds, aes(Year, fitted, color=sig, group=all), size=1) + 
#   scale_colour_viridis_d(option="viridis", begin=0.2, end=0.8) +
#   geom_ribbon(data=preds, aes(Year, ymin=lower, ymax=upper), alpha=0.2, fill="skyblue4") +
#   geom_point(data = dd, aes(Year, Discovered)) +
#   theme_classic() +
#   ylab("Viral discovery rate (viruses/year)") 





# # ============== facetted plot ===============]
# 
# r2 = result
# r2$Order = factor(r2$Order, levels=c("Cetartiodactyla", "Rodentia", "Primates", "Carnivora", "Chiroptera", "Perissodactyla", "Lagomorpha"), ordered=TRUE)
# r2$signif_col = NA
# r2$signif_col[ r2$sig_incr == TRUE ] = "Increase";
# r2$signif_col[ r2$sig_incr == FALSE ] = "Decrease"
# 
# #raw_data = curvesw[ curvesw$HostOrder %in% r2$Order & curvesw$Domestic == FALSE & curvesw$Year <= 2010, ]
# raw_data = curves[ curves$HostOrder %in% r2$Order & curves$Year <= 2010, ]
# raw_data$Order = raw_data$HostOrder
# raw_data$Order = factor(raw_data$Order, levels=c("Cetartiodactyla", "Rodentia", "Primates", "Carnivora", "Chiroptera", "Perissodactyla", "Lagomorpha"), ordered=TRUE)
# 
# pt2 = ggplot() + 
#   #geom_ribbon(data=r2[ -which(r2$Year < 1960 & r2$Order %in% c("Perissodactyla")), ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.3, fill="grey50", col=NA, size=0.05) +
#   geom_ribbon(data=r2, aes(x=Year, ymin=lower, ymax=upper), alpha=0.25, fill="skyblue4", col=NA, size=0.05) +
#   geom_point(data=raw_data, aes(x=Year, y=Discovered), col="grey60", size=0.15) +
#   geom_line(data=r2, aes(x=Year, y=fitted, group=Order, col=sig), size=1.2) +
#   theme_classic() +
#   scale_color_viridis_d(begin=0.25, end=0.75) +
#   lemon::facet_rep_wrap(~Order, ncol=1, scales="free_y", strip.position = "right") +
#   #facet_wrap(~Order, ncol=1, scales="free_y", strip.position = "right") +
#   theme(strip.background = element_blank(),
#         legend.position="none",
#         axis.text.y = element_text(size=11),
#         axis.text.x = element_text(size=13),
#         axis.title=element_text(size=14)) +
#   ylab("Virus discovery rate (viruses/year") +
#   scale_x_continuous(breaks=seq(1930, 2010, by=20), seq(1930, 2010, by=20), name="Year")
# pt2
# 
# ggsave(pt2, file="./test3_gam_figure1.png", device="png", units="in", width=5, height=7.5, dpi=300, scale=0.9)
# 
# install.packages("ggjoy")
# ggplot() + 
#   ggjoy::geom_joy(data=r2, aes(x=Year, y=fitted, group=Order, col=sig), size=1.2)
# 
# ggplot(iris, aes(x = Sepal.Length, y = Species)) + 
#   geom_joy() + 
#   theme_joy()



















# # ================= function to calculate windows of significance for heatmap plotting ====================
# 
# get_sig_boxes = function(x){
#   dat = result[ result$Order == unique(result$Order)[x], ]
#   dat = dat[ dat$sig == TRUE, ]
#   
#   if(nrow(dat) == 0){
#     to_return = data.frame(Order =  unique(result$Order)[x],
#                            YearStart = NA,
#                            YearEnd = NA, 
#                            type = NA)
#     
#   } else{
#     dat$type = ifelse(dat$sig_incr == TRUE, "Increase", "Decrease")
#     dat$gap = c(10, dat$Year[2:nrow(dat)] - dat$Year[1:(nrow(dat)-1)])
#     startyear = dat$Year[ dat$gap > 1 ]
#     endyear = c(dat$Year[ which(dat$gap > 1)-1 ], dat$Year[ nrow(dat)])
#     type = dat$type[ dat$gap > 1 ]
#     to_return = data.frame(Order =  unique(result$Order)[x],
#                            YearStart = startyear,
#                            YearEnd = endyear, 
#                            type = type)
#   }
#   
#   return(to_return)
# }
# boxx = do.call(rbind.data.frame, lapply(1:n_distinct(result$Order), get_sig_boxes))
# boxx = boxx[ !is.na(boxx$YearStart), ]


# =========== heat map =============

# colRamp = colorRampPalette(RColorBrewer::brewer.pal("YlGnBu", n=9))(400)
# 
# result$Order = factor(result$Order, levels=rev(c("Cetartiodactyla", "Rodentia", "Primates", "Carnivora", "Chiroptera", "Perissodactyla", "Lagomorpha")), ordered=TRUE)
# p_test = ggplot() +
#   geom_tile(data = result, aes(x=Year, y=Order, fill=fitted), width=1) +
#   #scale_fill_viridis_c(option="magma", name="Virus\ndiscovery\nrate") + theme_classic() + 
#   scale_fill_gradientn(colors = rev(colRamp), na.value="grey90", name="Virus\ndiscovery\nrate") +
#   geom_errorbar(data=boxx, aes(y=Order, xmin=YearStart, xmax = YearEnd, col=type), width=0.15, size=0.9) +
#   scale_color_manual(values=c("Decrease" = "red", "Increase" = "green"), name="Trend\ndirection") +
#   theme_classic() +
#   scale_x_continuous(breaks=seq(1930, 2010, by=20), seq(1930, 2010, by=20), name="Year") +
#   ylab("") + 
#   theme(axis.text.x = element_text(size=11),
#         axis.text.y = element_text(size=12),
#         axis.title.x = element_text(size=12),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank())
# ggsave(p_test, file="./test_gam_figure1.png", device="png", units="in", width=5, height=3.5, dpi=300)


