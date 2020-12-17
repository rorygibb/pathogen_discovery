

# ====================== Mapping temporal changes in pathogen discovery over space =====================

# match host records to IUCN range maps in CLOVER

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiscovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "raster", "sf", "maptools")

# associations data and no humans
assoc = read.csv("./data/clover_v1/Clover_v1.0_NBCIreconciled_20201211.csv", stringsAsFactors = FALSE) %>%
  #dplyr::filter(Database != "EID2") %>%
  dplyr::filter(YearType != "Nucleotide") %>%
  dplyr::filter(Host != "Homo sapiens") %>%
  dplyr::filter(DetectionMethod != "Not specified") 

# label domestics
load("./data/domestic/domestic_species.R")
domestic = Hmisc::capitalize(c(domestic, "canis familiaris", "bos frontalis", "bos grunniens", "bos grunniens mutus", "Bos taurus indicus", "Bos taurus primigenius"))
assoc$Domestic = ifelse(assoc$Host %in% domestic, "Domestic", "Wild")

# remove domestics
assoc = assoc[ assoc$Domestic == "Wild", ]
assoc = assoc[ assoc$Host_Original != "Canis lupus familiaris", ]

# iucn ranges: 950 full matches
iucn = sf::st_read("./data/iucn_range/mammals_terrestrial/MAMMALS_TERRESTRIAL_ONLY.shp")
iucnf = sf::st_read("./data/iucn_range/mammals_freshwater/MAMMALS_FRESHWATER.shp")
iucn = rbind(iucn, iucnf)

# world polygons
data("wrld_simpl")




# ============= for each time epoch starting in 1960 create rasters of geographical virus diversity =============

# template raster
tras = raster("./data/rasters/wc2.1_5m_elev/wc2.1_5m_elev.tif")
tras = aggregate(tras, fact=10)

# plot to generate maps of viral diversity
# saves to output folder
virusHostDiversityRasters = function(epoch, output_dir){
  
  # all known pathogens by specified year
  dd = assoc[ assoc$Year <= epoch & !is.na(assoc$Year), ]
  
  # path_stack
  path_stack = stack()
  
  # create raster of host richness per pathogen per epoch
  for(path_name in unique(dd$Virus)){
    
    cat(paste(path_name, "...", sep=""))
    dx = dd[ dd$Virus == path_name, ]
    
    iucnx = iucn %>%
      dplyr::filter(binomial %in% dx$Host) %>%
      #dplyr::summarise(geometry = sf::st_union(geometry)) %>%
      dplyr::mutate(Virus = path_name,
                    dummy = 1)
    #return(iucnx)
    
    if(nrow(iucnx) == 0){ next }
    
    #plot(iucnx$geometry)
    rx = fasterize::fasterize(iucnx, tras, field="dummy", fun='sum')
    names(rx) = paste(path_name, "HostRichness", epoch, sep="_")
    path_stack = stack(path_stack, rx)
    
  }
  writeRaster(path_stack, filename=paste(output_dir, names(path_stack), ".tif", sep=""), bylayer=TRUE, format="GTiff", overwrite=TRUE)
  return(path_stack)
}

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
  






# ============= map unique host-virus combinations discovered over time =============

# template raster
tras = raster("./data/rasters/wc2.1_5m_elev/wc2.1_5m_elev.tif")
tras = aggregate(tras, fact=10)

# plot to generate maps of viral discoveries
virusDiscoveryMap = function(startyear, endyear, criteria="all"){
  
  if(criteria == "strict"){
    uu = assoc[ assoc$DetectionMethod != "Antibodies", ]
  } else{
    uu = assoc
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

# create maps for serology
d1 = virusDiscoveryMap(startyear=1950, endyear=1970)
d2 = virusDiscoveryMap(startyear=1970, endyear=1980)
d3 = virusDiscoveryMap(startyear=1980, endyear=1990)
d4 = virusDiscoveryMap(startyear=1990, endyear=2000)
d5 = virusDiscoveryMap(startyear=2000, endyear=2005)
d6 = virusDiscoveryMap(startyear=2005, endyear=2010)
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
d1 = virusDiscoveryMap(startyear=1950, endyear=1970, criteria = "strict")
d2 = virusDiscoveryMap(startyear=1970, endyear=1980, criteria = "strict")
d3 = virusDiscoveryMap(startyear=1980, endyear=1990, criteria = "strict")
d4 = virusDiscoveryMap(startyear=1990, endyear=2000, criteria = "strict")
d5 = virusDiscoveryMap(startyear=2000, endyear=2005, criteria = "strict")
d6 = virusDiscoveryMap(startyear=2005, endyear=2010, criteria = "strict")
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
  theme(plot.title=element_text(size=14, hjust=0.5), legend.title = element_blank())
ggsave(plot_strict, file="./output/figures/PathogenDiscovery_byyear_strict.png", device="png", units="in", width=15, height=8, dpi=600)










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

