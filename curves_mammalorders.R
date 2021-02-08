


# ====================== Generate species viral discovery rate curves for mammals at the Order level =====================

# fit curves across time epochs to examine changes in trends 

# root dir and dependencies
# dependencies and basedir
setwd(here::here())
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon")

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
  group_by(HostOrder) %>%
  dplyr::summarise(VRichness = n_distinct(Virus))



# 1. all species

# unique associations by order and year
dd = clover %>%
  dplyr::group_by(HostOrder, Virus) %>%
  dplyr::summarise(Database = paste(unique(Database), collapse=", "),
                   HostSpecies = paste(unique(Host), collapse=", "),
                   NumRecords = length(Year),
                   YearEarliest = min(Year, na.rm=TRUE),
                   YearLatest = max(Year, na.rm=TRUE)) %>%
  left_join(tr) %>%
  dplyr::arrange(desc(VRichness), HostOrder, YearEarliest)

# create cumulative curves of all years
curves = expand.grid(unique(dd$HostOrder), 1930:2016) %>%
  dplyr::rename("HostOrder" = 1, "YearEarliest" = 2) %>%
  left_join(dd[ , c("HostOrder", "YearEarliest", "Virus")]) %>%
  dplyr::group_by(HostOrder, YearEarliest) %>%
  dplyr::summarise(Discovered = sum(!is.na(Virus)),
                   Virus = paste(unique(Virus), collapse=", ")) %>%
  left_join(dd[ !duplicated(dd$HostOrder), c("HostOrder", "VRichness") ]) %>%
  dplyr::arrange(desc(VRichness), HostOrder, YearEarliest) %>%
  dplyr::group_by(HostOrder) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered)) %>%
  dplyr::rename("Year" = YearEarliest)

# ggplot(curves[ curves$Year <= 2010, ]) +
#   geom_line(aes(Year, VirusCumulative)) +
#   facet_wrap(~HostOrder, scales="free_y")

# combine with publication effort
# curves = left_join(curves, pubs_order_all) %>%
#   dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0))

# 2. split out by wild/domestic
  
# unique associations by order, year
ddw = clover %>%
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
curvesw = expand.grid(unique(ddw$HostOrder), unique(ddw$Domestic), 1930:2016) %>%
  dplyr::rename("HostOrder" = 1, "Domestic" = 2, "YearEarliest" = 3) %>%
  left_join(ddw[ , c("HostOrder", "Domestic", "YearEarliest", "Virus")]) %>%
  dplyr::group_by(HostOrder, Domestic, YearEarliest) %>%
  dplyr::summarise(Discovered = sum(!is.na(Virus)),
                   Virus = paste(unique(Virus), collapse=", ")) %>%
  left_join(ddw[ !duplicated(ddw$HostOrder), c("HostOrder", "VRichness") ]) %>%
  dplyr::arrange(desc(VRichness), HostOrder, Domestic, YearEarliest) %>%
  dplyr::group_by(HostOrder, Domestic) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered)) %>%
  dplyr::rename("Year" = YearEarliest)

# combine with publication effort
# # need to debug this
# curvesw = left_join(curves, pubs_order_dom) %>%
#   dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0))




# ==================== Plot cumulative virus discovery curves for all Orders, split out by domestic/wild ==================================

# set up data for plotting 
plot_data = curvesw
plot_data$facet = paste(plot_data$HostOrder, " (", plot_data$VRichness, ")", sep="")
fac_order = plot_data[ !duplicated(plot_data$HostOrder), c("HostOrder", "facet", "VRichness")] %>% dplyr::arrange(desc(VRichness))
plot_data$facet = factor(plot_data$facet, levels=fac_order$facet, ordered=TRUE)
plot_data$DW = ifelse(plot_data$Domestic == TRUE, "Domestic", "Wild")

# plot for MS / SI up to 2010
px = ggplot(plot_data[ plot_data$Year <= 2010, ]) + 
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
ggsave(px, file="./output/figures/MS_MammalOrders_CumulativeCurves.png", device="png", units="in", width=10, height=6.5, dpi=300)




# ================ for each Order, fit poisson model to discovery rates ===================

# fits Poisson model of discovery rates (novel virus counts per year) for a specified order with specified data
# similar to fitting inhomogenous 1D Poisson process but without smudging event times
# for 3 time epochs: 1930 to present, 1960 to present, and 2000 to present
# currently effect of year is linear (i.e. exponential with log-link); could explore SPDE but may be more intractable if we're interested in broad trends

fitDiscoveryRateCurve = function(order, data){
  
  # data with hard cut-off at 2010 (from Reconciliation paper, this is when reports taper off)
  print(order)
  dx = data[ data$HostOrder == order & data$Year <= 2010, ]
  
  # 1. fit model from 1930 to present
  dx$yearx = 1:nrow(dx)
  form = Discovered ~ yearx + Intercept
  bru_mod = bru(form, dx, family = "poisson")
  #bru_mod = bru(form, dx, family = "nbinomial")
  summary(bru_mod)
  
  # extract fixed effects
  fx = bru_mod$summary.fixed
  fx$param = row.names(fx)
  names(fx)[ 3:5 ] = c("q0.025", "median", "q0.975")
  fx$model = "1930"
  fx$HostOrder = order
  row.names(fx) = c()
  
  # predict curve
  x4pred = data.frame(yearx = 1:nrow(dx))
  predx_bru1 = predict(bru_mod, x4pred, ~ exp(yearx + Intercept), n.samples=2000)
  resx = left_join(dx, predx_bru1, by="yearx")
  resx$model = "1930"
  
  # 2. fit model from 1960 to present
  dx = dx[ dx$Year >= 1960, ]
  dx$yearx = 1:nrow(dx)
  form = Discovered ~ yearx + Intercept
  bru_mod = bru(form, dx, family = "poisson")
  
  # extract fixed effects
  fy = bru_mod$summary.fixed
  fy$param = row.names(fy)
  names(fy)[ 3:5 ] = c("q0.025", "median", "q0.975")
  fy$model = "1960"
  fy$HostOrder = order
  row.names(fy) = c()
  
  # predict curve
  x4pred = data.frame(yearx = 1:nrow(dx))
  predx_bru1 = predict(bru_mod, x4pred, ~ exp(yearx + Intercept), n.samples=2000)
  resy = left_join(dx, predx_bru1, by="yearx")
  resy$model = "1960"
  
  # 2. fit model from 1990 to present
  dx = dx[ dx$Year >= 1990, ]
  dx$yearx = 1:nrow(dx)
  form = Discovered ~ yearx + Intercept
  bru_mod = bru(form, dx, family = "poisson")
  
  # extract fixed effects
  fz = bru_mod$summary.fixed
  fz$param = row.names(fz)
  names(fz)[ 3:5 ] = c("q0.025", "median", "q0.975")
  fz$model = "1990"
  fz$HostOrder = order
  row.names(fz) = c()
  
  # predict curve
  x4pred = data.frame(yearx = 1:nrow(dx))
  predx_bru1 = predict(bru_mod, x4pred, ~ exp(yearx + Intercept), n.samples=1995)
  resz = left_join(dx, predx_bru1, by="yearx")
  resz$model = "1990"
  
  # create results
  ff = rbind(fx, fy, fz)
  res = rbind(resx, resy, resz)
  return(list(
    fixed = ff,
    pred_curve = res
  ))
}

# 1. run models for all species
# result = lapply(unique(curves$HostOrder)[ 1:10 ], fitDiscoveryRateCurve, data=curves)
# 
# # extract estimates
# fixed = do.call(rbind.data.frame, lapply(result, "[[", 1))
# curve_preds = do.call(rbind.data.frame, lapply(result, "[[", 2))
# write.csv(fixed, "./output/order_models/fixedeffects_allspecies_byorder_inlabrupois_20200106.csv", row.names=FALSE)
# write.csv(curve_preds, "./output/order_models/curves_allspecies_byorder_inlabrupois_20200106.csv", row.names=FALSE)

# # 2. run models for wild species only
result = lapply(unique(curvesw$HostOrder)[ 1:7 ], fitDiscoveryRateCurve, data=curvesw[ curvesw$Domestic == FALSE, ])

# extract estimates
fixed = do.call(rbind.data.frame, lapply(result, "[[", 1))
curve_preds = do.call(rbind.data.frame, lapply(result, "[[", 2))
# write.csv(fixed, "./output/order_models/fixedeffects_wild_byorder_inlabrupois_Feb2021.csv", row.names=FALSE)
# write.csv(curve_preds, "./output/order_models/curves_wild_byorder_inlabrupois_Feb2021.csv", row.names=FALSE)

fixed = read.csv("./output/order_models/fixedeffects_wild_byorder_inlabrupois_Feb2021.csv")
curve_preds = read.csv("./output/order_models/curves_wild_byorder_inlabrupois_Feb2021.csv")




# ======================== visualise =====================

# plot curves
curve_preds$HO2 = paste(curve_preds$HostOrder, "\n(", curve_preds$VRichness, ")", sep="")
fac_order = curve_preds[ !duplicated(curve_preds$HostOrder), c("HostOrder", "VRichness", "HO2")] %>% 
  dplyr::arrange(desc(VRichness))
curve_preds$HostOrder = factor(curve_preds$HostOrder, levels=fac_order$HostOrder, ordered=TRUE)
curve_preds$model2 = paste(curve_preds$model, "-2010", sep="")
p1 = ggplot(curve_preds[ curve_preds$model %in% c(1930, 1960, 1990), ]) +
  geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.6) +
  geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=factor(model2)), alpha=0.25) +
  geom_line(aes(Year, median, group=factor(model2), col=factor(model2))) +
  theme_minimal() +
  lemon::facet_rep_wrap(~HostOrder, scales="free_y", nrow=2) +
  scale_colour_viridis_d( name="Time epoch", begin=0.2, end=0.5) +
  scale_fill_viridis_d( name="Time epoch", begin=0.2, end=0.7) +
  #ggtitle(Hmisc::capitalize(daty$Host[1])) +
  #ylab(expression(lambda)) +
  ylab("Viral discovery rate (viruses yer year)") +
  xlab("Year") +
  scale_x_continuous(breaks=seq(1950, 2010, by=30), labels=seq(1950, 2010, by=30)) +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        legend.title = element_text(size=12), 
        panel.grid=element_blank(),
        strip.text = element_text(size=13), 
        legend.position=c(0.88, 0.2),
        legend.text = element_text(size=11), 
        axis.text = element_text(size=11),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.ticks = element_line())
ggsave(p1, file="./output/figures/MS_Order_ViralDiscoveryRates.png", device="png", dpi=300, width=11, height=5.5, units="in", scale=0.95)

fixed$model2 = paste(fixed$model, "-2010", sep="")
fixed = left_join(fixed, fac_order)
fixed$HO2 = factor(fixed$HO2, levels=rev(fac_order$HO2), ordered=TRUE)
fixed$model2 = factor(fixed$model2, levels=rev(unique(fixed$model2)), ordered=TRUE)
range_lim = max(abs(range(c(fixed$q0.025[fixed$param=="yearx"], fixed$q0.975[fixed$param=="yearx"]))))
#ggplot(fixed[ fixed$param %in% c("yearx") & fixed$model != "2000", ]) + 
p2 = ggplot(fixed[ fixed$param %in% c("yearx"), ]) + 
  geom_point(aes(HO2, median, col=model2, group=model2), position=position_dodge(width=0.5), size=2) +
  geom_linerange(aes(HO2, ymin=q0.025, ymax=q0.975, col=model2, group=model2), position=position_dodge(width=0.5)) + 
  geom_hline(yintercept=0, lty=2) +
  theme_classic() +
  #ylab(expression(beta[year])) +
  ylab("Trend in virus discovery rate") +
  scale_color_viridis_d( name="Time epoch", begin=0.2, end=0.7, direction=-1, guide=guide_legend(reverse=TRUE)) +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.x = element_text(size=12),
        panel.border = element_rect(fill=NA),
        axis.title.y = element_blank(),
        legend.position = c(0.18, 0.88),
        legend.title = element_blank(), 
        legend.text = element_text(size=10), 
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11)) + 
  coord_flip() + 
  ylim(-range_lim, range_lim)
ggsave(p2, file="./output/figures/MS_Order_BetaYear_Slopes.png", device="png", dpi=300, width=5.5, height=5.1, units="in", scale=0.9)

# exponentiate to natural scale so 
ff_exp = fixed[ fixed$param == "yearx", ]
ff_exp$median = (exp(ff_exp$median)-1)*100
ff_exp$q0.025 = (exp(ff_exp$q0.025)-1)*100
ff_exp$q0.975 = (exp(ff_exp$q0.975)-1)*100
range_lim = max(abs(range(c(ff_exp$q0.025, ff_exp$q0.975)))) + 3
p2_exp = ggplot(ff_exp) + 
  geom_point(aes(HO2, median, col=model2, group=model2), position=position_dodge(width=0.5), size=2) +
  geom_linerange(aes(HO2, ymin=q0.025, ymax=q0.975, col=model2, group=model2), position=position_dodge(width=0.5)) + 
  geom_hline(yintercept=0, lty=2) +
  theme_classic() +
  #ylab(expression(beta[year])) +
  ylab("Trend in virus discovery rate (% change per year)") +
  scale_color_viridis_d( name="Time epoch", begin=0.2, end=0.7, direction=-1, guide=guide_legend(reverse=TRUE)) +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.x = element_text(size=12),
        panel.border = element_rect(fill=NA),
        axis.title.y = element_blank(),
        legend.position = c(0.18, 0.88),
        legend.title = element_blank(), 
        legend.text = element_text(size=10), 
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11)) + 
  coord_flip() + 
  ylim(-range_lim, range_lim)
ggsave(p2_exp, file="./output/figures/MS_Order_BetaYear_Slopes_Exp.png", device="png", dpi=300, width=5.5, height=5.1, units="in", scale=0.9)




# plot cumulative curves
curves2 = curves
curves2$HostOrder = factor(curves2$HostOrder, levels=unique(curves2$HostOrder), ordered=TRUE)
p3 =ggplot(curves2[ curves2$Year <= 2010, ]) + 
  geom_line(aes(Year, VirusCumulative, group=HostOrder), size=0.8, col="skyblue4") +
  facet_wrap(~HostOrder, scales="free_y") + 
  theme_minimal() + 
  ylab("Viral richness") +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=11),
        legend.title = element_text(size=12), 
        strip.text = element_text(size=12.5), 
        legend.text = element_text(size=11), 
        axis.text = element_text(size=11))
ggsave(p3, file="./output/figures/Order_CumulativeCurves_20201203.png", device="png", dpi=300, width=12, height=8, units="in", scale=0.95)

  
  
# 
# 
# ggplot(curve_preds[ curve_preds$model == 1930, ]) +
#   geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
#   geom_line(aes(Year, median)) +
#   geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=HostOrder), alpha=0.25) +
#   theme_minimal() +
#   facet_wrap(~HostOrder, scales="free_y") +
#   #ggtitle(Hmisc::capitalize(daty$Host[1])) +
#   ylab(expression(lambda)) +
#   xlab("Year") +
#   theme(plot.title=element_text(size=14, hjust=0.5),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11))
# 
# 
# 
# ggplot() +
#   geom_point(data = resy, aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
#   geom_line(data = resy, aes(Year, median)) +
#   geom_ribbon(data = resy, aes(Year, ymin=q0.025, ymax=q0.975), fill="skyblue4", alpha=0.25) +
#   theme_minimal() +
#   #ggtitle(Hmisc::capitalize(daty$Host[1])) +
#   ylab(expression(lambda)) +
#   xlab("Year") +
#   theme(plot.title=element_text(size=14, hjust=0.5),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11))
# 
# 
# # species
# spp = "pan troglodytes"
# datx = curves[ curves$Host == spp & curves$Year <= 2015, ]
# #datx = curves[ curves$Year <= 2015 & curves$Host %in% unique(curves$Host)[1:250] & curves$Host != "homo sapiens", ]
# 
# # formula: linear effect of year + intercept
# datx$x = 1:nrow(datx)
# form = Discovered ~ x + Intercept
# bru_mod = bru(form, datx, family = "poisson")
# summary(bru_mod)
# 
# # predict field for each year
# x4pred = data.frame(x = 1:nrow(datx))
# predx_bru1 = predict(bru_mod, x4pred, ~ exp(x + Intercept), n.samples=2000)
# daty = left_join(datx, predx_bru1)
# 
# ggplot() +
#   geom_point(data = daty, aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
#   geom_line(data = daty, aes(Year, median)) + 
#   geom_ribbon(data = daty, aes(Year, ymin=q0.025, ymax=q0.975), fill="skyblue4", alpha=0.25) +
#   theme_minimal() + 
#   #ggtitle(Hmisc::capitalize(daty$Host[1])) + 
#   ylab(expression(lambda)) + 
#   xlab("Year") +
#   theme(plot.title=element_text(size=14, hjust=0.5),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11))
# 
# # formula: linear effect of year + intercept > 1995
# datx = datx[ datx$Year >= 1995 ,]
# datx$x = 1:nrow(datx)
# form = Discovered ~ x + Intercept
# bru_mod = bru(form, datx, family = "poisson")
# summary(bru_mod)
# 
# # predict field for each year
# x4pred = data.frame(x = 1:nrow(datx))
# predx_bru1 = predict(bru_mod, x4pred, ~ exp(x + Intercept), n.samples=2000)
# daty = left_join(datx, predx_bru1)
# 
# ggplot() +
#   geom_point(data = daty, aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
#   geom_line(data = daty, aes(Year, median)) + 
#   geom_ribbon(data = daty, aes(Year, ymin=q0.025, ymax=q0.975), fill="skyblue4", alpha=0.25) +
#   theme_minimal() + 
#   #ggtitle(Hmisc::capitalize(daty$Host[1])) + 
#   #ylab("Number of viruses discovered") + 
#   ylab(expression(lambda)) +
#   xlab("Year") +
#   theme(plot.title=element_text(size=14, hjust=0.5),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11))
# 
# 
# # formula: linear effect of year + intercept
# datx = datx[ datx$Year >= 1970 ,]
# datx$x = 1:nrow(datx)
# form = Discovered ~ x + Intercept
# bru_mod = bru(form, datx, family = "poisson")
# summary(bru_mod)
# 
# # predict field for each year
# x4pred = data.frame(x = 1:nrow(datx))
# predx_bru1 = predict(bru_mod, x4pred, ~ exp(x + Intercept), n.samples=2000)
# daty = left_join(datx, predx_bru1)
# 
# ggplot() +
#   geom_point(data = daty, aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
#   geom_line(data = daty, aes(Year, median)) + 
#   geom_ribbon(data = daty, aes(Year, ymin=q0.025, ymax=q0.975), fill="skyblue4", alpha=0.25) +
#   theme_minimal() + 
#   ggtitle(Hmisc::capitalize(daty$Host[1])) + 
#   ylab("Number of viruses discovered") + 
#   xlab("Year") +
#   theme(plot.title=element_text(size=14, hjust=0.5),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11))
# 
# 
# 



# ======================= 2. Comparison of decay in total richness by 
