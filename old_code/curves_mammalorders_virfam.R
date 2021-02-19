


# ====================== Generate species viral discovery rate curves for mammals at the Order level =====================

# proof of concept/general trends
# fit curves across time epochs to examine changes in trends 
# how robust are these to adjusting for annual publication effort?

# root dir and dependencies
# dependencies and basedir
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "here")
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiscovery/code/pathogen_discovery/")

# domestic species to label
domestic = read.csv("./data/clover/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)

# associations data and no humans
clover = read.csv("./data/clover/Clover_v1.0_NBCIreconciled_20201218.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Host != "Homo sapiens") %>%
  dplyr::mutate(Domestic = ifelse(Host %in% domestic$Host, TRUE, FALSE)) %>%
  dplyr::filter(DetectionMethod != "Not specified")

# total viral family richness by order
tr = clover %>%
  group_by(HostOrder) %>%
  dplyr::summarise(VRichness = n_distinct(VirusFamily))

# publication effort by year at Order-level, with wild/domestic split
# this needs some thought and work if it's worth including
# pubs = read.csv("./output/host_effort/PubMed_Hosts_PubsByYear_Final.csv", stringsAsFactors = FALSE) %>%
#   dplyr::filter(Note == "") %>%
#   dplyr::select(1:3)
# pubs = pubs %>%
#   left_join(assoc[ , c("Host", "HostOrder", "Domestic") ]) %>%
#   dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0))
# pubs_order_all = pubsx %>%
#   dplyr::group_by(HostOrder, Year) %>%
#   dplyr::summarise(NumPubs = sum(NumPubs, na.rm=TRUE))
# pubs_order_dom = pubsx %>%
#   dplyr::group_by(HostOrder, Year, Domestic) %>%
#   dplyr::summarise(NumPubs = sum(NumPubs, na.rm=TRUE))


# 1. all species

# unique associations by order and year
dd = clover %>%
  dplyr::group_by(HostOrder, VirusFamily) %>%
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
  left_join(dd[ , c("HostOrder", "YearEarliest", "VirusFamily")]) %>%
  dplyr::group_by(HostOrder, YearEarliest) %>%
  dplyr::summarise(Discovered = sum(!is.na(VirusFamily)),
                   VirusFamily = paste(unique(VirusFamily), collapse=", ")) %>%
  left_join(dd[ !duplicated(dd$HostOrder), c("HostOrder", "VRichness") ]) %>%
  dplyr::arrange(desc(VRichness), HostOrder, YearEarliest) %>%
  dplyr::group_by(HostOrder) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered)) %>%
  dplyr::rename("Year" = YearEarliest)

ggplot(curves[ curves$Year <= 2010, ]) +
  geom_line(aes(Year, VirusCumulative)) +
  facet_wrap(~HostOrder, scales="free_y")

# combine with publication effort
# curves = left_join(curves, pubs_order_all) %>%
#   dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0))

# 2. split out by wild/domestic
  
# # unique associations by order, year
# ddw = clover %>%
#   dplyr::filter(!is.na(Year)) %>%
#   dplyr::group_by(HostOrder, Domestic, Virus) %>%
#   dplyr::summarise(Database = paste(unique(Database), collapse=", "),
#                    HostSpecies = paste(unique(Host), collapse=", "),
#                    NumRecords = length(Year),
#                    YearEarliest = min(Year, na.rm=TRUE),
#                    YearLatest = max(Year, na.rm=TRUE)) %>%
#   left_join(tr) %>%
#   dplyr::arrange(desc(VRichness), HostOrder, YearEarliest) 
# 
# # create cumulative curves of all years
# curvesw = expand.grid(unique(ddw$HostOrder), unique(ddw$Domestic), 1930:2016) %>%
#   dplyr::rename("HostOrder" = 1, "Domestic" = 2, "YearEarliest" = 3) %>%
#   left_join(ddw[ , c("HostOrder", "Domestic", "YearEarliest", "Virus")]) %>%
#   dplyr::group_by(HostOrder, Domestic, YearEarliest) %>%
#   dplyr::summarise(Discovered = sum(!is.na(Virus)),
#                    Virus = paste(unique(Virus), collapse=", ")) %>%
#   left_join(ddw[ !duplicated(ddw$HostOrder), c("HostOrder", "VRichness") ]) %>%
#   dplyr::arrange(desc(VRichness), HostOrder, Domestic, YearEarliest) %>%
#   dplyr::group_by(HostOrder, Domestic) %>%
#   dplyr::mutate(VirusCumulative = cumsum(Discovered)) %>%
#   dplyr::rename("Year" = YearEarliest)
# 
# # combine with publication effort
# # # need to debug this
# # curvesw = left_join(curves, pubs_order_dom) %>%
# #   dplyr::mutate(NumPubs = replace(NumPubs, is.na(NumPubs), 0))
# 
# ggplot(curvesw[ curvesw$Year <= 2010, ]) + 
#   geom_line(aes(Year, VirusCumulative, col=Domestic)) + 
#   facet_wrap(~HostOrder, scales="free_y")



# ================ for each Order, fit poisson model to discovery rates ===================

# fits Poisson model of discovery rates (novel virus counts per year) for a specified order with specified data
# equivalent to fitting inhomogenous 1D Poisson process but without smudging event times
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
  #summary(bru_mod)
  
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
result = lapply(unique(curves$HostOrder)[ 1:10 ], fitDiscoveryRateCurve, data=curves)

# extract estimates
fixed = do.call(rbind.data.frame, lapply(result, "[[", 1))
curve_preds = do.call(rbind.data.frame, lapply(result, "[[", 2))
write.csv(fixed, "./output/order_models/fixedeffects_allspecies_byorder_virusfamily_inlabrupois_20200106.csv", row.names=FALSE)
write.csv(curve_preds, "./output/order_models/curves_allspecies_byorder_virusfamily_inlabrupois_20200106.csv", row.names=FALSE)

# # 2. run models for wild species only
# result = lapply(unique(curvesw$HostOrder)[ 1:10 ], fitDiscoveryRateCurve, data=curvesw[ curvesw$Domestic == "Wild", ])
# 
# # extract estimates
# fixed = do.call(rbind.data.frame, lapply(result, "[[", 1))
# curve_preds = do.call(rbind.data.frame, lapply(result, "[[", 2))
# write.csv(fixed, "./output/order_models/fixedeffects_wild_byorder_inlabrupois.csv", row.names=FALSE)
# write.csv(curve_preds, "./output/order_models/curves_wild_byorder_inlabrupois.csv", row.names=FALSE)
# 






# ======================== visualise =====================

curve_preds$HostOrder = factor(curve_preds$HostOrder, levels=unique(curve_preds$HostOrder), ordered=TRUE)

# ggplot(curve_preds[ curve_preds$model == 1930, ]) +
#   geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
#   geom_line(aes(Year, median)) +
#   geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=HostOrder), alpha=0.25) +
#   theme_minimal() +
#   facet_wrap(~HostOrder) +
#   #ggtitle(Hmisc::capitalize(daty$Host[1])) +
#   ylab(expression(lambda)) +
#   xlab("Year") +
#   theme(plot.title=element_text(size=14, hjust=0.5),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11))
# 
# ggplot(curve_preds[ curve_preds$model == 1960, ]) +
#   geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
#   geom_line(aes(Year, median)) +
#   geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=HostOrder), alpha=0.25) +
#   theme_minimal() +
#   facet_wrap(~HostOrder) +
#   #ggtitle(Hmisc::capitalize(daty$Host[1])) +
#   ylab(expression(lambda)) +
#   xlab("Year") +
#   theme(plot.title=element_text(size=14, hjust=0.5),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11))
# 
# ggplot(curve_preds[ curve_preds$model == 2000, ]) +
#   geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
#   geom_line(aes(Year, median)) +
#   geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=HostOrder), alpha=0.25) +
#   theme_minimal() +
#   facet_wrap(~HostOrder) +
#   #ggtitle(Hmisc::capitalize(daty$Host[1])) +
#   ylab("Viral discovery rate (viruses/year)") +
#   xlab("Year") +
#   theme(plot.title=element_text(size=14, hjust=0.5),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=12),
#         axis.text = element_text(size=11))

curve_preds$model2 = paste(curve_preds$model, "-2010", sep="")
p1 = ggplot(curve_preds[ curve_preds$model %in% c(1930, 1960), ]) +
  geom_point(aes(Year, Discovered), size=1, col="grey70", alpha=0.8) +
  geom_line(aes(Year, median, group=factor(model2))) +
  geom_ribbon(aes(Year, ymin=q0.025, ymax=q0.975, fill=factor(model2)), alpha=0.3) +
  theme_minimal() +
  facet_wrap(~HostOrder, scales="free_y", nrow=2) +
  scale_fill_viridis_d( name="Time epoch", begin=0.2, end=0.7) +
  #ggtitle(Hmisc::capitalize(daty$Host[1])) +
  #ylab(expression(lambda)) +
  ylab("Virus family discovery rate (viral families/year)") +
  xlab("Year") +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=11),
        legend.title = element_text(size=12), 
        strip.text = element_text(size=13), 
        legend.position="bottom",
        legend.text = element_text(size=11), 
        axis.text = element_text(size=11))
ggsave(p1, file="./output/figures/Order_ViralDiscoveryRates_20201203.png", device="png", dpi=300, width=17, height=8, units="in", scale=0.95)

fixed$model2 = paste(fixed$model, "-2010", sep="")
#ggplot(fixed[ fixed$param %in% c("yearx") & fixed$model != "2000", ]) + 
p2 = ggplot(fixed[ fixed$param %in% c("yearx"), ]) + 
  geom_point(aes(HostOrder, median, col=model2, group=model2), position=position_dodge(width=0.5)) +
  geom_linerange(aes(HostOrder, ymin=q0.025, ymax=q0.975, col=model2, group=model2), position=position_dodge(width=0.5)) + 
  geom_hline(yintercept=0, lty=2) +
  theme_minimal() +
  ylab(expression(beta[year])) +
  scale_color_viridis_d( name="Time epoch", begin=0.2, end=0.7) +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=11), 
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11, angle=90))
ggsave(p2, file="./output/figures/Order_Betaestimates_20201203.png", device="png", dpi=300, width=8, height=5, units="in", scale=0.95)


# plot cumulative curves
curves2 = curves
curves2$HostOrder = factor(curves2$HostOrder, levels=unique(curves2$HostOrder), ordered=TRUE)
p3 =ggplot(curves2[ curves2$Year <= 2010, ]) + 
  geom_line(aes(Year, VirusCumulative, group=HostOrder), size=0.8, col="skyblue4") +
  facet_wrap(~HostOrder, scales="free_y") + 
  theme_minimal() + 
  ylab("Virus family richness") +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=11),
        legend.title = element_text(size=12), 
        strip.text = element_text(size=12.5), 
        legend.text = element_text(size=11), 
        axis.text = element_text(size=11))
ggsave(p3, file="./output/figures/Order_CumulativeCurves_20201203.png", device="png", dpi=300, width=12, height=8, units="in", scale=0.95)

  

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
