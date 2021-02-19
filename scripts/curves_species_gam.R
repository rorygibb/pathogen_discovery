


# ====================== Generate species viral discovery rate curves for mammals at the Order level =====================

# fit GAM curves across time epochs to examine changes in trends 

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_discovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon", "mgcv")

# domestic species to label
domestic = read.csv("./data/clover/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)

# associations data and no humans
clover = read.csv("./data/clover/Clover_v1.0_NBCIreconciled_20201218.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Host != "Homo sapiens") %>%
  dplyr::mutate(Domestic = ifelse(Host %in% domestic$Host, TRUE, FALSE)) %>%
  #dplyr::filter(DetectionMethod != "Antibodies") %>%
  dplyr::filter(DetectionMethod != "Not specified")

# total viral richness by order
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




# =================== fit GAMs and estimate curve derivative for specified Orders ============================

# data frame for results and specified orders
result = data.frame()

# get first n species
n = 50
spp = unique(curves$Host)[1:n]
datax = curves

# for each order 
for(i in 1:length(spp)){
  
  # data (either wild only, or all) with cutoff of 2010
  spx = spp[i]
  dd = datax[ datax$Host == spx & datax$Year <= 2010, ]
  
  # specify time threshold for inclusion of zeros prior to the first virus discovered in that taxa
  # either start from the first year of virus discovery, or n years beforehand, or can exclude this and run from 0
  time_thresh = 10
  first_year = min(dd$Year[ dd$Discovered > 0 ])
  dd = dd[ dd$Year >= first_year - time_thresh, ]
  
  # GAM fit with spline on Year and using ML, Poisson likelihood
  m1 = mgcv::gamm(Discovered ~ s(Year), family="poisson", data=dd, method="ML")
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
r2$Host2 = r2$Host
r2$Host2[ r2$Host == "Canis lupus familiaris" ] = "Canis familiaris"
r2$Host2 = lapply(strsplit(r2$Host2, " "), function(x) paste(x, collapse="\n"))
r2$Host2 = paste(r2$Host2, " (", r2$VRichness, ")", sep="")
r2$Host2 = factor(r2$Host2, levels=unique(r2$Host2), ordered=TRUE)

r2$signif_col = NA
r2$signif_col[ r2$sig_incr == TRUE ] = "Increase";
r2$signif_col[ r2$sig_incr == FALSE ] = "Decrease"

#raw_data = curvesw[ curvesw$HostOrder %in% r2$Order & curvesw$Domestic == FALSE & curvesw$Year <= 2010, ]
raw_data = datax[ datax$Host %in% r2$Host & datax$Year <= 2010, ]
raw_data = left_join(raw_data, r2[ !duplicated(r2$Host), c("Host", "Host2")])
raw_data$Host2 = factor(raw_data$Host2, levels=unique(r2$Host2), ordered=TRUE)

curve_plot = ggplot() + 
  #geom_ribbon(data=r2[ -which(r2$Year < 1960 & r2$Order %in% c("Perissodactyla")), ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.3, fill="grey50", col=NA, size=0.05) +
  geom_ribbon(data=r2[ r2$Host2 %in% unique(r2$Host2)[1:25], ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.25, fill="skyblue4", col=NA, size=0.05) +
  geom_point(data=raw_data[ raw_data$Host2 %in% unique(r2$Host2)[1:25], ], aes(x=Year, y=Discovered), col="grey55", size=0.3) +
  geom_line(data=r2[ r2$Host2 %in% unique(r2$Host2)[1:25], ], aes(x=Year, y=fitted, group=Host2, col=sig_incr), size=1.2) +
  geom_line(data=r2[ r2$Host2 %in% unique(r2$Host2)[1:25] & r2$sig_decr == TRUE, ], aes(x=Year, y=fitted, group=Host2), col="red", size=1.2) +
  theme_classic() +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  lemon::facet_rep_wrap(~Host2, ncol=5, nrow=5, scales="free_y") +
  #facet_wrap(~Order, ncol=1, scales="free_y", strip.position = "right") +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text = element_text(size=13),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title=element_text(size=16)) +
  ylab("Virus discovery rate (viruses per year)") +
  scale_x_continuous(breaks=seq(1940, 2000, by=20), seq(1940, 2000, by=20), name="Year")
ggsave(curve_plot, file="./output/figures/SI_Figure_SpeciesGAMs_A.png", device="png", units="in", width=10, height=10, dpi=300)

curve_plot2 = ggplot() + 
  #geom_ribbon(data=r2[ -which(r2$Year < 1960 & r2$Order %in% c("Perissodactyla")), ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.3, fill="grey50", col=NA, size=0.05) +
  geom_ribbon(data=r2[ r2$Host2 %in% unique(r2$Host2)[26:50], ], aes(x=Year, ymin=lower, ymax=upper), alpha=0.25, fill="skyblue4", col=NA, size=0.05) +
  geom_point(data=raw_data[ raw_data$Host2 %in% unique(r2$Host2)[26:50], ], aes(x=Year, y=Discovered), col="grey55", size=0.3) +
  geom_line(data=r2[ r2$Host2 %in% unique(r2$Host2)[26:50], ], aes(x=Year, y=fitted, group=Host2, col=sig_incr), size=1.2) +
  geom_line(data=r2[ r2$Host2 %in% unique(r2$Host2)[26:50] & r2$sig_decr == TRUE, ], aes(x=Year, y=fitted, group=Host2), col="red", size=1.2) +
  theme_classic() +
  scale_color_viridis_d(begin=0.25, end=0.75) +
  lemon::facet_rep_wrap(~Host2, ncol=5, nrow=5, scales="free_y") +
  #facet_wrap(~Order, ncol=1, scales="free_y", strip.position = "right") +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text = element_text(size=13),
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title=element_text(size=16)) +
  ylab("Virus discovery rate (viruses per year)") +
  scale_x_continuous(breaks=seq(1940, 2000, by=20), seq(1940, 2000, by=20), name="Year")
ggsave(curve_plot2, file="./output/figures/SI_Figure_SpeciesGAMs_B.png", device="png", units="in", width=10, height=10, dpi=300)



ggplot(preds) +
  geom_line(aes(Year, fitted, group=Host, col=sig)) +
  geom_ribbon(aes(Year, ymin=lower, ymax=upper), alpha=0.2, fill="skyblue4")
  



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



# ================= function to calculate windows of significance for heatmap plotting ====================

# get_sig_boxes = function(x){
#   
#   dat = result[ result$Host == unique(result$Host)[x], ]
#   dat = dat[ dat$sig == TRUE, ]
# 
#   if(nrow(dat) == 0){
#     to_return = data.frame(Host =  unique(result$Host)[x],
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
#     to_return = data.frame(Host =  unique(result$Host)[x],
#                            YearStart = startyear,
#                            YearEnd = endyear,
#                            type = type)
#   }
# 
#   return(to_return)
# }
# boxx = do.call(rbind.data.frame, lapply(1:n_distinct(result$Host), get_sig_boxes))
# boxx = boxx[ !is.na(boxx$YearStart), ]



# =========== heat map =============

# colRamp = colorRampPalette(RColorBrewer::brewer.pal("YlGnBu", n=9))(400)
# 
# result$Host = factor(result$Host, levels=rev(spp), ordered=TRUE)
# result$Host2 = paste(result$Host, " (", result$VRichness, ")", sep="")
# result$Host2 = factor(result$Host2, levels=rev(unique(result$Host2)), ordered=TRUE)
# 
# boxx = left_join(boxx, result[ !duplicated(result$Host), c("Host", "Host2") ])
# 
# p_test = ggplot() +
#   geom_tile(data = result, aes(x=Year, y=Host2, fill=fitted), width=1) +
#   #scale_fill_viridis_c(option="magma", name="Virus\ndiscovery\nrate") + theme_classic() +
#   scale_fill_gradientn(colors = rev(colRamp), na.value="grey90", name="Virus\ndiscovery\nrate") +
#   geom_errorbar(data=boxx, aes(y=Host2, xmin=YearStart, xmax = YearEnd, col=type), width=0.15, size=0.9) +
#   scale_color_manual(values=c("Decrease" = "red", "Increase" = "green"), name="Trend\ndirection") +
#   theme_classic() +
#   scale_x_continuous(breaks=seq(1930, 2010, by=20), seq(1930, 2010, by=20), name="Year") +
#   ylab("") +
#   theme(axis.text.x = element_text(size=14),
#         axis.text.y = element_text(size=12),
#         axis.title.x = element_text(size=14),
#         axis.line.y = element_blank(),
#         legend.title = element_text(size=13),
#         legend.text = element_text(size=13),
#         axis.ticks.y = element_blank())
# p_test
# ggsave(p_test, file="./test_gam_heatmap_spp.png", device="png", units="in", width=8, height=8.5, dpi=300, scale=0.95)
# 

