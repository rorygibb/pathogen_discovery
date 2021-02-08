

# ====================== Temporal changes in pathogen discovery over space and time =====================

# match host records to IUCN range maps in CLOVER

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiscovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "raster", "sf", "maptools")

# domestic species to label
domestic = read.csv("./code/pathogen_discovery/data/clover/domestic_status/HostLookup_Domestic.csv", stringsAsFactors = FALSE)

# associations data and no humans
clover = read.csv("./code/pathogen_discovery/data/clover/Clover_v1.0_NBCIreconciled_20201218.csv", stringsAsFactors = FALSE) %>%
  #dplyr::filter(Host != "Homo sapiens") %>%
  dplyr::mutate(Domestic = ifelse(Host %in% domestic$Host, TRUE, FALSE)) %>%
  dplyr::filter(DetectionMethod != "Not specified") %>%
  #dplyr::filter(Domestic == FALSE) %>%
  dplyr::filter(!is.na(Year)) %>%
  dplyr::filter(Year <= 2010)

# iucn ranges: 950 full matches
iucn = sf::st_read("./data/iucn_range/mammals_terrestrial/MAMMALS_TERRESTRIAL_ONLY.shp")
iucnf = sf::st_read("./data/iucn_range/mammals_freshwater/MAMMALS_FRESHWATER.shp")
iucn = rbind(iucn, iucnf)

# world polygons
data("wrld_simpl")



# ============= graph changes in discovery over time ================

dat = clover %>%
  group_by(DetectionMethod, Host, Virus) %>%
  dplyr::summarise(Year = min(Year), 
                   hv = paste(unique(Host), unique(Virus), collapse=", ")) %>%
  dplyr::arrange(Host, Virus) %>%
  dplyr::group_by(DetectionMethod, Year) %>%
  dplyr::summarise(Discovered = n_distinct(hv)) %>%
  dplyr::group_by(DetectionMethod) %>%
  dplyr::mutate(CumDiscovered = cumsum(Discovered))

# melt
dat = reshape2::melt(dat, id.vars=1:2)
dat$variable = as.vector(dat$variable)
dat$variable[ dat$variable == "Discovered" ] = "Novel host-virus associations"
dat$variable[ dat$variable == "CumDiscovered" ] = "Cumulative associations"
dat$variable = factor(dat$variable, levels = c("Novel host-virus associations", "Cumulative associations"), ordered=TRUE)
dat = rename(dat, "Method" = DetectionMethod)

p0 = ggplot(dat) + 
  geom_line(aes(Year, value, group=Method, col=Method), size=0.9) + 
  theme_classic() +
  #theme(legend.position=c(0.2, 0.9)) +
  ylab("Number of associations") +
  facet_wrap(~variable, nrow=2, scales="free_y") +
  #theme(legend.position="bottom") +
  scale_x_continuous(breaks=seq(1930, 2020, by=20), labels=seq(1930, 2020, by=20)) +
  scale_color_viridis_d(begin=0, end=0.85, option="viridis", direction=-1) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.background = element_blank(), 
        axis.title = element_text(size=13), 
        strip.text = element_text(size=13),
        #strip.text = element_blank(),
        axis.text = element_text(size=12))

ggsave(p0, file="./output/figures/VirusDiscovery_Trends.png", device="png", units="in", width=7, height=7, dpi=600, scale=0.9)


# full composite

p0 = ggplot(dat) + 
  geom_line(aes(Year, value, group=Method, col=Method), size=0.9) + 
  theme_classic() +
  #theme(legend.position=c(0.2, 0.9)) +
  ylab("Novel host-virus associations reported") +
  facet_wrap(~variable, nrow=2, scales="free_y") +
  theme(legend.position="none") +
  scale_x_continuous(breaks=seq(1930, 2020, by=20), labels=seq(1930, 2020, by=20)) +
  scale_color_viridis_d(begin=0, end=0.85, option="viridis", direction=-1) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.background = element_blank(), 
        axis.title = element_text(size=13), 
        strip.text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12))

# virus and host overall
dx = clover %>%
  group_by(Virus) %>%
  dplyr::summarise(Year = min(Year)) %>%
  dplyr::arrange(Virus) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(Discovered = n_distinct(Virus)) %>%
  dplyr::mutate(Cumulative = cumsum(Discovered),
                Type = "Viruses")
dy = clover %>%
  group_by(Host) %>%
  dplyr::summarise(Year = min(Year)) %>%
  dplyr::arrange(Host) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(Discovered = n_distinct(Host)) %>%
  dplyr::mutate(Cumulative = cumsum(Discovered),
                Type = "Hosts")
dat = rbind(dx, dy)

p1 = ggplot(dat) + 
  geom_line(aes(Year, Cumulative, group=Type), size=0.9, col="darkred") + 
  theme_classic() +
  #theme(legend.position=c(0.2, 0.9)) +
  ylab("Cumulative species") +
  facet_wrap(~Type, nrow=2, scales="free_y") +
  #theme(legend.position="bottom") +
  scale_x_continuous(breaks=seq(1930, 2020, by=20), labels=seq(1930, 2020, by=20)) +
  #scale_color_viridis_d(begin=0, end=0.85, option="viridis", direction=-1) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.background = element_blank(), 
        axis.title = element_text(size=13), 
        strip.text = element_text(size=13),
        #strip.text = element_blank(),
        axis.text = element_text(size=12))

p_composite = gridExtra::grid.arrange(p0, p1, nrow=2, heights=c(1, 0.5))
ggsave(p_composite, file="./output/figures/VirusDiscovery_Trends_composite.png", device="png", units="in", width=5, height=9, dpi=600, scale=0.8)


# ============= map unique host-virus combinations discovered over time =============

# template raster
tras = raster("./data/rasters/wc2.1_5m_elev/wc2.1_5m_elev.tif")
tras = aggregate(tras, fact=10)

# plot to generate maps of viral discoveries
virusDiscoveryMap = function(data, startyear, endyear, criteria="all"){
  
  if(criteria == "pcr_isolation"){
    uu = data[ data$DetectionMethod != "Antibodies", ]
  } else if(criteria == "isolation"){
    uu = data[ data$DetectionMethod == "Isolation/Observation", ]
  }  else{
    uu = data
  }
  
  # unique_combinations with first year
  uu = uu %>%
    group_by(Host, Virus) %>%
    dplyr::summarise(Year = min(Year, na.rm=TRUE)) %>%
    dplyr::filter(!is.na(Year))
  
  # all known pathogens by specified year and summarise by host (how many discovered)
  dx = uu[ uu$Year >= startyear & uu$Year <= endyear, ] %>%
    dplyr::group_by(Host) %>%
    dplyr::summarise(VirDisc = n_distinct(Virus))
  
  # combine with iucn
  iucnx = iucn %>%
    dplyr::filter(binomial %in% dx$Host) %>%
    #dplyr::summarise(geometry = sf::st_union(geometry)) %>%
    dplyr::left_join(dx, by=c("binomial" = "Host"))
  
  # rasterize
  rx = fasterize::fasterize(iucnx, tras, field="VirDisc", fun='sum')
  names(rx) = paste("UniqueHostVirusCombinations", startyear, endyear, sep="_")
  return(rx)
}


# ====================== partition style 1 ============================

# create maps for all criteria
d1 = virusDiscoveryMap(data=clover, startyear=1950, endyear=1970)
d2 = virusDiscoveryMap(data=clover, startyear=1971, endyear=1980)
d3 = virusDiscoveryMap(data=clover, startyear=1981, endyear=1990)
d4 = virusDiscoveryMap(data=clover, startyear=1991, endyear=2000)
d5 = virusDiscoveryMap(data=clover, startyear=2001, endyear=2005)
d6 = virusDiscoveryMap(data=clover, startyear=2006, endyear=2010)
discovery = stack(d1, d2, d3, d4, d5, d6)

# combine
rr = as.data.frame(discovery, xy=TRUE) %>%
  reshape2::melt(id.vars =1:2) %>%
  dplyr::mutate(startyear = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)),
                endyear = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 3)),
                facet = paste(startyear, endyear, sep="-"))
ww = st_as_sf(wrld_simpl)
rr$value[ rr$value == 0] = NA

plot_sero = ggplot() + 
  geom_raster(data = rr, aes(x, y, fill=value)) +
  geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25, size=0.25) +
  facet_wrap(~facet) + 
  scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
  theme_minimal() + xlab("") + ylab("") + 
  ggtitle("Novel host-virus associations reported (serology, PCR or isolation)") + 
  theme(plot.title=element_text(size=14, hjust=0.5), legend.title = element_blank())
ggsave(plot_sero, file="./output/figures/PathogenDiscovery_byyear_serology.png", device="png", units="in", width=15, height=8, dpi=600)


# create maps for isolation
d1 = virusDiscoveryMap(data=clover, startyear=1950, endyear=1970, criteria = "pcr_isolation")
d2 = virusDiscoveryMap(data=clover, startyear=1971, endyear=1980, criteria = "pcr_isolation")
d3 = virusDiscoveryMap(data=clover, startyear=1981, endyear=1990, criteria = "pcr_isolation")
d4 = virusDiscoveryMap(data=clover, startyear=1991, endyear=2000, criteria = "pcr_isolation")
d5 = virusDiscoveryMap(data=clover, startyear=2001, endyear=2005, criteria = "pcr_isolation")
d6 = virusDiscoveryMap(data=clover, startyear=2006, endyear=2010, criteria = "pcr_isolation")
discovery = stack(d1, d2, d3, d4, d5, d6)

# combine
rr = as.data.frame(discovery, xy=TRUE) %>%
  reshape2::melt(id.vars =1:2) %>%
  dplyr::mutate(startyear = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)),
                endyear = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 3)),
                facet = paste(startyear, endyear, sep="-"))
ww = st_as_sf(wrld_simpl)
rr$value[ rr$value == 0] = NA

plot_strict = ggplot() + 
  geom_raster(data = rr, aes(x, y, fill=value)) +
  geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25, size=0.25) +
  facet_wrap(~facet) + 
  scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
  theme_minimal() + xlab("") + ylab("") + 
  ggtitle("Novel host-virus associations reported (PCR or isolation)") + 
  theme(plot.title=element_text(size=16, hjust=0.5), legend.title = element_blank(),
        strip.text = element_text(size=13.5))
ggsave(plot_strict, file="./output/figures/PathogenDiscovery_byyear_strict.png", device="png", units="in", width=15, height=8, dpi=600)




# ====================== partition style 2: decades ============================

# create maps for all criteria
d1 = virusDiscoveryMap(data=clover, startyear=1930, endyear=1950)
d2 = virusDiscoveryMap(data=clover, startyear=1951, endyear=1970)
d3 = virusDiscoveryMap(data=clover, startyear=1971, endyear=1980)
d4 = virusDiscoveryMap(data=clover, startyear=1981, endyear=1990)
d5 = virusDiscoveryMap(data=clover, startyear=1991, endyear=2000)
d6 = virusDiscoveryMap(data=clover, startyear=2001, endyear=2010)
discovery = stack(d1, d2, d3, d4, d5, d6)

# combine
rr = as.data.frame(discovery, xy=TRUE) %>%
  reshape2::melt(id.vars =1:2) %>%
  dplyr::mutate(startyear = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)),
                endyear = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 3)),
                facet = paste(startyear, endyear, sep="-"))
ww = st_as_sf(wrld_simpl)
rr$value[ rr$value == 0] = NA

plot_sero = ggplot() + 
  geom_raster(data = rr, aes(x, y, fill=value)) +
  geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25, size=0.25) +
  facet_wrap(~facet) + 
  scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
  theme_minimal() + xlab("") + ylab("") + 
  ggtitle("Novel host-virus associations reported (serology, PCR or isolation)") + 
  theme(plot.title=element_text(size=16, hjust=0.5), legend.title = element_blank(),
        strip.text = element_text(size=15),
        legend.text = element_text(size=12))
ggsave(plot_sero, file="./output/figures/PathogenDiscovery_byyear_allcriteria_decades.png", device="png", units="in", width=14, height=7, dpi=600, scale=0.9)


# create maps for isolation
d1 = virusDiscoveryMap(data=clover, startyear=1930, endyear=1950, criteria = "pcr_isolation")
d2 = virusDiscoveryMap(data=clover, startyear=1951, endyear=1970, criteria = "pcr_isolation")
d3 = virusDiscoveryMap(data=clover, startyear=1971, endyear=1980, criteria = "pcr_isolation")
d4 = virusDiscoveryMap(data=clover, startyear=1981, endyear=1990, criteria = "pcr_isolation")
d5 = virusDiscoveryMap(data=clover, startyear=1991, endyear=2000, criteria = "pcr_isolation")
d6 = virusDiscoveryMap(data=clover, startyear=2001, endyear=2010, criteria = "pcr_isolation")
discovery = stack(d1, d2, d3, d4, d5, d6)

# combine
rr = as.data.frame(discovery, xy=TRUE) %>%
  reshape2::melt(id.vars =1:2) %>%
  dplyr::mutate(startyear = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)),
                endyear = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 3)),
                facet = paste(startyear, endyear, sep="-"))
ww = st_as_sf(wrld_simpl)
rr$value[ rr$value == 0] = NA

plot_strict = ggplot() + 
  geom_raster(data = rr, aes(x, y, fill=value)) +
  geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25, size=0.25) +
  facet_wrap(~facet) + 
  scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
  theme_minimal() + xlab("") + ylab("") + 
  ggtitle("Novel host-virus associations reported (PCR or isolation)") + 
  theme(plot.title=element_text(size=16, hjust=0.5), legend.title = element_blank(),
        strip.text = element_text(size=15),
        legend.text = element_text(size=12))
ggsave(plot_strict, file="./output/figures/PathogenDiscovery_byyear_strict_decades.png", device="png", units="in", width=14, height=7, dpi=600, scale=0.9)


plot_strict2 = ggplot() + 
  geom_raster(data = rr, aes(x, y, fill=value)) +
  geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25, size=0.25) +
  facet_wrap(~facet, ncol=2) + 
  scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
  theme_minimal() + xlab("") + ylab("") + 
  ggtitle("Novel host-virus associations reported (PCR or isolation)") + 
  theme(plot.title=element_text(size=16, hjust=0.5), legend.title = element_blank(),
        strip.text = element_text(size=15),
        legend.text = element_text(size=12))
ggsave(plot_strict2, file="./output/figures/PathogenDiscovery_byyear_strict_decades_vertical.png", device="png", units="in", width=10, height=9, dpi=600, scale=0.9)









### for filoviruses

dat = assoc[ assoc$VirusOrder == "Bunyavirales", ]

# create maps for isolation
d1 = virusDiscoveryMap(data = dat, startyear=1950, endyear=1980, criteria = "strict")
d2 = virusDiscoveryMap(data = dat, startyear=1980, endyear=1990, criteria = "strict")
d3 = virusDiscoveryMap(data = dat, startyear=1990, endyear=2000, criteria = "strict")
d4 = virusDiscoveryMap(data = dat, startyear=2000, endyear=2005, criteria = "strict")
d5 = virusDiscoveryMap(data = dat, startyear=2005, endyear=2010, criteria = "strict")
discovery = stack(d1, d2, d3, d4, d5)

# combine
rr = as.data.frame(discovery, xy=TRUE) %>%
  reshape2::melt(id.vars =1:2) %>%
  dplyr::mutate(startyear = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)),
                endyear = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 3)),
                facet = paste(startyear, endyear, sep="-"))
ww = st_as_sf(wrld_simpl)
rr$value[ rr$value == 0] = NA

pp = ggplot() + 
  geom_raster(data = rr, aes(x, y, fill=value)) +
  geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25, size=0.1) +
  facet_wrap(~facet) + 
  scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
  theme_minimal() + xlab("") + ylab("") + 
  ggtitle("Bunyavirus-host associations reported") + 
  theme(plot.title=element_text(size=14, hjust=0.5), legend.title = element_blank())








# for betacoronaviruses

hh = assoc[ assoc$VirusGenus == "Betacoronavirus", ] %>%
  dplyr::select(Host, Year) %>%
  distinct()
bc = read.csv("./code/pathogen_discovery/data/BetaCoV_hosts_Dec2020.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(1:2) %>%
  dplyr::rename("Host" = 1, "Year" = 2)
hh = rbind(hh, bc)

betacov_stack = stack()
for(yy in seq(1996, 2020, by=3)){
  
  hy = hh[ hh$Year <= yy, ]
  iucnx = iucn %>%
    dplyr::filter(binomial %in% hy$Host) %>%
    dplyr::mutate(dummy = 1)
  rx = fasterize::fasterize(iucnx, tras, field="dummy", fun='sum')
  names(rx) = yy
  betacov_stack = stack(betacov_stack, rx)
}

rr = as.data.frame(betacov_stack, xy=TRUE) %>%
  reshape2::melt(id.vars =1:2) %>%
  dplyr::mutate(year = substr(variable, 2, 30))

data("wrld_simpl")
ww = st_as_sf(wrld_simpl)

pb = ggplot() +
  geom_raster(data = rr, aes(x, y, fill=value)) +
  geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25) +
  facet_wrap(~year) +
  scale_fill_viridis_c(option="magma", direction=-1, na.value="white") +
  theme_minimal() + xlab("") + ylab("") +
  ggtitle("Betacoronavirus host richness") +
  theme(plot.title=element_text(size=14, hjust=0.5), legend.title = element_blank())
ggsave(pb, file="./output/figures/BetaCoV_hosts_byyear.png", device="png", units="in", width=15, height=12, dpi=600)



# # ============= for each time epoch starting in 1960 create rasters of geographical virus diversity =============
# 
# # template raster
# tras = raster("./data/rasters/wc2.1_5m_elev/wc2.1_5m_elev.tif")
# tras = aggregate(tras, fact=10)
# 
# # plot to generate maps of viral diversity
# # saves to output folder
# virusHostDiversityRasters = function(epoch, output_dir){
#   
#   # all known pathogens by specified year
#   dd = assoc[ assoc$Year <= epoch & !is.na(assoc$Year), ]
#   
#   # path_stack
#   path_stack = stack()
#   
#   # create raster of host richness per pathogen per epoch
#   for(path_name in unique(dd$Virus)){
#     
#     cat(paste(path_name, "...", sep=""))
#     dx = dd[ dd$Virus == path_name, ]
#     
#     iucnx = iucn %>%
#       dplyr::filter(binomial %in% dx$Host) %>%
#       #dplyr::summarise(geometry = sf::st_union(geometry)) %>%
#       dplyr::mutate(Virus = path_name,
#                     dummy = 1)
#     #return(iucnx)
#     
#     if(nrow(iucnx) == 0){ next }
#     
#     #plot(iucnx$geometry)
#     rx = fasterize::fasterize(iucnx, tras, field="dummy", fun='sum')
#     names(rx) = paste(path_name, "HostRichness", epoch, sep="_")
#     path_stack = stack(path_stack, rx)
#     
#   }
#   writeRaster(path_stack, filename=paste(output_dir, names(path_stack), ".tif", sep=""), bylayer=TRUE, format="GTiff", overwrite=TRUE)
#   return(path_stack)
# }
# 
# # create wildlife virus diversity maps across time epochs
# output_dir = "./output/virus_rasters/"
# vdiv_1960 = virusHostDiversityRasters(epoch = 1960, output_dir=output_dir)
# vdiv_1970 = virusHostDiversityRasters(epoch = 1970, output_dir=output_dir)
# vdiv_1980 = virusHostDiversityRasters(epoch = 1980, output_dir=output_dir)
# vdiv_1990 = virusHostDiversityRasters(epoch = 1990, output_dir=output_dir)
# vdiv_2000 = virusHostDiversityRasters(epoch = 2000, output_dir=output_dir)
# vdiv_2010 = virusHostDiversityRasters(epoch = 2010, output_dir=output_dir)
# 
# # combine into summary maps of total viral richness
# div_1960 = raster::calc(vdiv_1960, function(x) sum(x > 0 & !is.na(x))); names(div_1960) = "VRichness_1960"
# div_1970 = raster::calc(vdiv_1970, function(x) sum(x > 0 & !is.na(x))); names(div_1970) = "VRichness_1970"
# div_1980 = raster::calc(vdiv_1980, function(x) sum(x > 0 & !is.na(x))); names(div_1980) = "VRichness_1980"
# div_1990 = raster::calc(vdiv_1990, function(x) sum(x > 0 & !is.na(x))); names(div_1990) = "VRichness_1990"
# div_2000 = raster::calc(vdiv_2000, function(x) sum(x > 0 & !is.na(x))); names(div_2000) = "VRichness_2000"
# div_2010 = raster::calc(vdiv_2010, function(x) sum(x > 0 & !is.na(x))); names(div_2010) = "VRichness_2010"
# 
# ss = stack(div_1960, div_1970, div_1980, div_1990, div_2000, div_2010)
# rr = as.data.frame(ss, xy=TRUE) %>%
#   reshape2::melt(id.vars =1:2) %>%
#   dplyr::mutate(year = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)))
# data("wrld_simpl")
# ww = st_as_sf(wrld_simpl)
# rr$value[ rr$value == 0] = NA
# 
# pb = ggplot() + 
#   geom_raster(data = rr, aes(x, y, fill=value)) +
#   geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25) +
#   facet_wrap(~year) + 
#   scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
#   theme_minimal() + xlab("") + ylab("") + 
#   ggtitle("Total viral richness") + 
#   theme(plot.title=element_text(size=14, hjust=0.5), legend.title = element_blank())
# 
# # changes in decades
# d1 = div_1970 - div_1960
# d2 = div_1980 - div_1970
# d3 = div_1990 - div_1980
# d4 = div_2000 - div_1990
# d5 = div_2010 - div_2000
# changes = stack(d1, d2, d3, d4, d5)
# names(changes) = c("d1960-1970", "d1970-1980", "d1980-1990", "d1990-2000", "d2000-2010")
# 
# rr = as.data.frame(changes, xy=TRUE) %>%
#   reshape2::melt(id.vars =1:2) %>%
#   dplyr::mutate(variable = substr(variable, 2, 50))
# data("wrld_simpl")
# ww = st_as_sf(wrld_simpl)
# rr$value[ rr$value == 0] = NA
# 
# pb = ggplot() + 
#   geom_raster(data = rr, aes(x, y, fill=value)) +
#   geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25) +
#   facet_wrap(~variable) + 
#   scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
#   theme_minimal() + xlab("") + ylab("") + 
#   ggtitle("Viruses discovered") + 
#   theme(plot.title=element_text(size=14, hjust=0.5), legend.title = element_blank())



# ### BETACOVs
# 
# assocx = assoc[ assoc$VirusGenus == "Betacoronavirus", ]
# 
# virusHostDiversityRasters = function(assocx, years, output_dir){
#   
#   result = stack()
#   
#   for(i in years){
#     
#     print(i)
#     dx = assocx[ assocx$Year <= i & !is.na(assocx$Year), ]
#     iucnx = iucn %>%
#       dplyr::filter(binomial %in% dx$Host) %>%
#       #dplyr::summarise(geometry = sf::st_union(geometry)) %>%
#       dplyr::mutate(dummy = 1)
#     if(nrow(iucnx) == 0){ next }
#     
#     #plot(iucnx$geometry)
#     iucnx = as_Spatial(iucnx)
#     rx = raster::rasterize(iucnx, tras, field=iucnx$dummy, fun='sum')
#     names(rx) = paste("Betacoronavirus", "HostRichness", i, sep="_")
#     result = stack(result, rx)
#   }
#   
#   return(result)
# }
# 
# betacov_ras = virusHostDiversityRasters(assocx, years = seq(1990, 2015, by=3))
# 
# rr = as.data.frame(betacov_ras, xy=TRUE) %>%
#   reshape2::melt(id.vars =1:2) %>%
#   dplyr::mutate(year = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 3)))
# 
# library(maptools)
# data("wrld_simpl")
# ww = st_as_sf(wrld_simpl)
# ggplot() + 
#   #geom_raster(rr[ rr$year >= 1999, ], aes(x, y, fill=value)) + 
#   geom_sf(data=ww, aes(), col="grey50") +
#   facet_wrap(~year) + 
#   scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
#   coord_fixed() + 
#   theme_minimal()
# 
# pb = ggplot() + 
#   geom_raster(data = rr[ rr$year >= 1999, ], aes(x, y, fill=value)) +
#   geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25) +
#   facet_wrap(~year) + 
#   scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
#   theme_minimal() + xlab("") + ylab("") + 
#   ggtitle("Betacoronavirus host richness") + 
#   theme(plot.title=element_text(size=14, hjust=0.5), legend.title = element_blank())
# ggsave(pb, file="./output/figures/BetaCoV_hosts_byyear.png", device="png", units="in", width=15, height=8, dpi=600)
  















# create wildlife virus diversity maps across time epochs
output_dir = "./output/virus_rasters/"
vdiv_1960 = virusHostDiversityRasters(epoch = 1960, output_dir=output_dir)
vdiv_1970 = virusHostDiversityRasters(epoch = 1970, output_dir=output_dir)
vdiv_1980 = virusHostDiversityRasters(epoch = 1980, output_dir=output_dir)
vdiv_1990 = virusHostDiversityRasters(epoch = 1990, output_dir=output_dir)
vdiv_2000 = virusHostDiversityRasters(epoch = 2000, output_dir=output_dir)
vdiv_2010 = virusHostDiversityRasters(epoch = 2010, output_dir=output_dir)

# combine into summary maps of total viral richness
div_1960 = raster::calc(vdiv_1960, function(x) sum(x > 0 & !is.na(x))); names(div_1960) = "VRichness_1960"
div_1970 = raster::calc(vdiv_1970, function(x) sum(x > 0 & !is.na(x))); names(div_1970) = "VRichness_1970"
div_1980 = raster::calc(vdiv_1980, function(x) sum(x > 0 & !is.na(x))); names(div_1980) = "VRichness_1980"
div_1990 = raster::calc(vdiv_1990, function(x) sum(x > 0 & !is.na(x))); names(div_1990) = "VRichness_1990"
div_2000 = raster::calc(vdiv_2000, function(x) sum(x > 0 & !is.na(x))); names(div_2000) = "VRichness_2000"
div_2010 = raster::calc(vdiv_2010, function(x) sum(x > 0 & !is.na(x))); names(div_2010) = "VRichness_2010"

ss = stack(div_1960, div_1970, div_1980, div_1990, div_2000, div_2010)
rr = as.data.frame(ss, xy=TRUE) %>%
  reshape2::melt(id.vars =1:2) %>%
  dplyr::mutate(year = unlist(lapply(strsplit(as.vector(variable), "_"), "[", 2)))
data("wrld_simpl")
ww = st_as_sf(wrld_simpl)
rr$value[ rr$value == 0] = NA

pb = ggplot() + 
  geom_raster(data = rr, aes(x, y, fill=value)) +
  geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25) +
  facet_wrap(~year) + 
  scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
  theme_minimal() + xlab("") + ylab("") + 
  ggtitle("Total viral richness") + 
  theme(plot.title=element_text(size=14, hjust=0.5), legend.title = element_blank())

# changes in decades
d1 = div_1970 - div_1960
d2 = div_1980 - div_1970
d3 = div_1990 - div_1980
d4 = div_2000 - div_1990
d5 = div_2010 - div_2000
changes = stack(d1, d2, d3, d4, d5)
names(changes) = c("d1960-1970", "d1970-1980", "d1980-1990", "d1990-2000", "d2000-2010")

rr = as.data.frame(changes, xy=TRUE) %>%
  reshape2::melt(id.vars =1:2) %>%
  dplyr::mutate(variable = substr(variable, 2, 50))
data("wrld_simpl")
ww = st_as_sf(wrld_simpl)
rr$value[ rr$value == 0] = NA

pb = ggplot() + 
  geom_raster(data = rr, aes(x, y, fill=value)) +
  geom_sf(data=ww[ ww$NAME != "Antarctica", ], fill=NA, col="grey50", alpha=0.25) +
  facet_wrap(~variable) + 
  scale_fill_viridis_c(option="magma", direction=-1, na.value="white") + 
  theme_minimal() + xlab("") + ylab("") + 
  ggtitle("Viruses discovered") + 
  theme(plot.title=element_text(size=14, hjust=0.5), legend.title = element_blank())

