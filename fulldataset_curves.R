


# ====================== Generate species viral discovery curves for mammals using Shaw, GMPD2, HP3 ==================

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiscovery/code/pathogen_discovery/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2")

# associations data
assoc = read.csv("./data/clover_v1/Clover_reconciledassociations_v1_20201120.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Database != "EID2") %>%
  dplyr::filter(YearType != "Nucleotide")

# domestic
load("./data/domestic/domestic_species.R")
domestic = c(domestic, "canis familiaris", "bos frontalis")
assoc$Domestic = ifelse(assoc$Host_Harmonised %in% domestic, "Domestic", "Wild")

# total viral richness by species
tr = assoc %>%
  group_by(Host_Harmonised) %>%
  dplyr::summarise(VRichness = n_distinct(Pathogen_Harmonised))

# unique associations by year
dd = assoc %>%
  dplyr::group_by(Host_Harmonised, Pathogen_Harmonised) %>%
  dplyr::summarise(HostOrder = head(HostOrder, 1), 
                   HostFamily = head(HostFamily, 1), 
                   Database = paste(unique(Database), collapse=", "),
                   NumRecords = length(Year),
                   Year = min(Year, na.rm=TRUE),
                   Domestic = head(Domestic, 1)) %>%
  left_join(tr) %>%
  dplyr::arrange(desc(VRichness), Host_Harmonised, Year) %>%
  #dplyr::filter(Host_Harmonised != "homo sapiens") %>%
  dplyr::filter(Year != Inf) # temporary fix for pathogens with no year

# create cumulative curves of all years
curves = expand.grid(unique(dd$Host_Harmonised), 1930:2016) %>%
  dplyr::rename("Host_Harmonised" = 1, "Year" = 2) %>%
  left_join(dd[ , c("Host_Harmonised", "Year", "Pathogen_Harmonised")]) %>%
  dplyr::group_by(Host_Harmonised, Year) %>%
  dplyr::summarise(Discovered = sum(!is.na(Pathogen_Harmonised)),
                   Virus = paste(unique(Pathogen_Harmonised), collapse=", ")) %>%
  left_join(dd[ !duplicated(dd$Host_Harmonised), c("Host_Harmonised", "HostOrder", "HostFamily", "VRichness", "Domestic") ]) %>%
  dplyr::rename("Host" = Host_Harmonised) %>%
  dplyr::arrange(desc(VRichness), Host, Year) %>%
  dplyr::group_by(Host) %>%
  dplyr::mutate(VirusCumulative = cumsum(Discovered))
  

# some initial viz
ggplot(curves[ curves$Host %in% unique(curves$Host)[1:40], ]) + 
  geom_line(aes(Year, VirusCumulative, group=Host, col=Host), size=1)
ggplot(curves[ curves$Host %in% unique(curves$Host)[20:40], ]) + 
  geom_line(aes(Year, VirusCumulative, group=Host, col=Host), size=1)
ggplot(curves[ curves$Host %in% unique(curves$Host)[20:40], ]) + 
  geom_bar(aes(Year, Discovered), stat="identity") + 
  facet_wrap(~Host)

# mean viruses discovered per year per species over time
# interesting taper down from the mid-2000s: reflects shift towards GenBank and delay in publishing time?
ggplot(curves) + 
  geom_smooth(aes(Year, Discovered), method="gam", method.args=list(family="poisson"))
ggplot(curves) + 
  geom_smooth(aes(Year, Discovered, group=Domestic, col=Domestic), method="gam", method.args=list(family="poisson"))
ggplot(curves[ curves$Domestic == "Wild", ]) + 
  #geom_point(aes(Year, Discovered)) +
  geom_smooth(aes(Year, Discovered), method="gam", method.args=list(family="poisson"))

# mean cumulative viruses
ggplot(curves) + 
  geom_smooth(aes(Year, VirusCumulative), method="gam", method.args=list(family="poisson"))
ggplot(curves) + 
  geom_smooth(aes(Year, VirusCumulative, group=Domestic, col=Domestic), method="gam", method.args=list(family="poisson"))
ggplot(curves[ curves$Domestic == "Wild", ]) + 
  geom_smooth(aes(Year, VirusCumulative), method="gam", method.args=list(family="poisson"))





# =============== nonstationary PP ==============

library(NHPoisson)

# https://stats.stackexchange.com/questions/49012/how-to-estimate-poisson-process-using-r-or-how-to-use-nhpoisson-package


lambda=1/60 #1 event per minute
time.span=60*60*24 #24 hours, with time granularity one second

aux<-simNHP.fun(rep(lambda,time.span))

# fit PP
out<-fitPP.fun(posE=aux$posNH,n=time.span,start=list(b0=0)) # b0=0 is our guess at initial value for optimization, which is internally made with `nlminb` function


# nonstationary
time.span=60*60*24 #24 hours, with time granularity one second
all.seconds<-seq(1,time.span,length.out=time.span)
lambdas=0.05*exp(-0.0001*all.seconds) #we can't model a linear regression with NHPoisson. It must have the form with exp.
aux<-simNHP.fun(lambdas)

# fit model
out<-fitPP.fun(tind=TRUE,covariates=cbind(all.seconds),
               posE=aux$posNH,
               start=list(b0=0,b1=0),modSim=TRUE)

# fit for people
dat = curves[ curves$Host == "eidolon helvum", ]
ggplot(dat) + geom_line(aes(Year, VirusCumulative))

# create vector of all years rep'd by maximum number of possible events within a year
pois_dat = data.frame()
for(i in 1:nrow(dat)){
  dfx = data.frame(Year = rep(dat$Year[i], max(dat$Discovered)),
                   Disc = c(rep(1, dat$Discovered[i]), rep(0, max(dat$Discovered - dat$Discovered[i]))))
  pois_dat = rbind(pois_dat, dfx)  
}
pois_dat$ind = 1:nrow(pois_dat)

# vector of event locations
all.years = pois_dat$ind
all.disc = pois_dat$ind[ pois_dat$Disc == 1]

# fit model
out = fitPP.fun(tind=TRUE, covariates=cbind(all.years),
                posE=all.disc,
                start=list(b0=0,b1=0), modSim=TRUE)



# 
# 
# # ====================== Plot curves for species pathogen rarefaction ===========================
# 
# # root dir and dependencies
# # dependencies and basedir
# setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiversity_rarefaction/")
# pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2")
# 
# # full harmoniesd associations database (EID2, GMPD2, HP3)
# dd = read.csv("./output/data_processed/hostpathogen_harmonised/AllDatabases_Associations_Hosts_Harmonised_Oct2020.csv", stringsAsFactors = FALSE)
# dd = dd[ !is.na(dd$Database), ]
# dd = dd[ !is.na(dd$ParasiteType), ]
# dd = dd[ !is.na(dd$Host_Harmonised), ]
# 
# # publication effort by year
# pubs = read.csv("./output/data_processed/host_effort/PubMed_Hosts_PubsByYear_Final.csv", stringsAsFactors = FALSE) %>%
#   dplyr::mutate(Year = as.numeric(Year)) %>%
#   dplyr::filter(!is.na(Year)) %>%
#   dplyr::filter(!Note %in% c("No publications", "Lookup error"))
# pp = expand.grid(1950:2020, unique(pubs$Host)) %>%
#   rename("Year"=Var1, "Host_Harmonised"=Var2) %>%
#   left_join(pubs[ , c("Year", "Host", "NumPubs")], by=c("Year"="Year", "Host_Harmonised"="Host"))
# pp$NumPubs[ is.na(pp$NumPubs) ] = 0
# pp = pp %>%
#   dplyr::arrange(Host_Harmonised, Year) %>%
#   group_by(Host_Harmonised) %>%
#   dplyr::mutate(CumPubs = cumsum(NumPubs))
# 
# 
# # ggplot(pp[ pp$Host_Harmonised %in% sample(unique(pp$Host_Harmonised), 20, replace=F), ]) +
# #   geom_line(aes(Year, NumPubs)) +
# #   facet_wrap(~Host_Harmonised, scales="free_y") +
# #   theme_bw()
# # ggplot(pp[ pp$Host_Harmonised %in% sample(unique(pp$Host_Harmonised), 20, replace=F), ]) +
# #   geom_line(aes(Year, CumPubs)) +
# #   facet_wrap(~Host_Harmonised, scales="free_y") +
# #   theme_bw()
# 
# 
# # # create dataframe for everything
# # dd = expand.grid(1950:2020, unique(assoc$Host_Harmonised)) %>%
# #   rename("Year"=Var1, "Host_Harmonised"=Var2) %>%
# #   left_join(assoc[ !duplicated(assoc$Host_Harmonised), c("Host_Harmonised", "HostClass", "HostOrder", "HostFamily")], by=c("Host_Harmonised")) %>%
# #   left_join(assoc[ , c("Host_Harmonised", "Year", "PathogenName_Harmonised", "ParasiteType", "Database", "DetectionQuality", "HumanInfective_Any", "IsZoonotic") ], by=c("Host_Harmonised", "Year"))
# 
# ### mean curve
# 
# px = ggplot(vcurves[ vcurves$Year < 2017 & vcurves$HostClass == "Mammalia", ]) + 
#   geom_smooth(aes(Year, PathRich), method="gam", method.args=list(family="poisson")) +
#   #geom_point(aes(Year, PathRich, group=Host_Harmonised), alpha=0.2, size=0.1) +
#   theme_minimal() +
#   ylab("Mean pathogen richness (all species)")
# ggsave(px,file= "./output/figures/MeanDiscoveryCurve.png", device="png", units="in", width=6, height=5, scale=1)
# 
# vcurves$HostOrder[ vcurves$HostOrder == "Artiodactyla" ] = "Cetartiodactyla"
# vc2 = vcurves %>%
#   filter(Year < 2017,
#          HostClass == "Mammalia",
#          HostOrder %in% c("Rodentia", "Chiroptera", "Primates", "Cetartiodactyla", "Perissodactyla", "Carnivora", "Lagomorpha"))
# px_order = ggplot(vc2) + 
#   geom_smooth(aes(Year, PathRich, group=HostOrder, col=HostOrder), method="gam", method.args=list(family="poisson")) +
#   #geom_point(aes(Year, PathRich, group=Host_Harmonised), alpha=0.2, size=0.1) +
#   theme_minimal() +
#   ylab("Mean pathogen richness (all species)") + 
#   geom_hline(yintercept=0, lty=2)
# ggsave(px_order,file= "./output/figures/MeanDiscoveryCurve_byorder.png", device="png", units="in", width=6, height=5, scale=1)
# 
# px_order2 = ggplot(vc2[ vc2$TotalPathRich>9, ]) + 
#   geom_smooth(aes(Year, PathRich, group=HostOrder, col=HostOrder), method="gam", method.args=list(family="poisson")) +
#   #geom_point(aes(Year, PathRich, group=Host_Harmonised), alpha=0.2, size=0.1) +
#   theme_minimal() +
#   ylab("Mean pathogen richness (all species)") + 
#   geom_hline(yintercept=0, lty=2)

# ================ summarise by virus over time =====================



hostrange = dd[ dd$ParasiteType == "virus" & dd$HostClass =="Mammalia" & dd$Year > 1800 & dd$Year < 2020 & dd$Year != 1900, ] %>%
  group_by(Parasite) %>%
  dplyr::summarise(
    FirstYearReported = min(Year),
    HostRange = n_distinct(Host_Harmonised),
    OrderRange = n_distinct(HostOrder),
    FamilyRange = n_distinct(HostFamily)
    )

# virus first year of discovery
px = ggplot(hostrange) + 
  geom_point(aes(FirstYearReported, HostRange), size=2, alpha=0.5, col="skyblue4") + 
  geom_smooth(aes(FirstYearReported, HostRange), method="gam") +
  theme_classic()
ggsave(px,file= "./output/figures/HostRangeOverTime_Mammals.png", device="png", units="in", width=6, height=5, scale=1)

px = ggplot(hostrange) + 
  geom_point(aes(FirstYearReported, FamilyRange), size=2, alpha=0.5, col="skyblue4") + 
  geom_smooth(aes(FirstYearReported, FamilyRange), method="gam") +
  theme_classic() + ylab("Num Families infected")
ggsave(px,file= "./output/figures/HostRangeOverTime_Mammals_Family.png", device="png", units="in", width=6, height=5, scale=1)




# 
# dfx = curves[ curves$Host_Harmonised == "mastomys natalensis", ]
# a = ggplot(dfx) + geom_line(aes(Year, PathRich), col="coral2", size=0.5) + theme_bw()
# b = ggplot(dfx) + geom_line(aes(Year, NumPubs), col="blue", size=0.5) + theme_bw()
# 
# 
# 
# # viz
# ggplot(curves[ curves$TotalPathRich >= 40 & curves$Year <=2017 ,]) + 
#   geom_point(aes(Year, PathRich, size=NumRecords), alpha=0.3, col="skyblue4") +
#   geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
#   facet_wrap(~Host_Harmonised, scales="free_y") +
#   theme_minimal()
# 
# ggplot(curves[ curves$HostOrder == "Rodentia" & curves$TotalPathRich >= 5 & curves$Year <=2017 ,]) + 
#   geom_point(aes(Year, PathRich, size=NumRecords), alpha=0.3, col="skyblue4") +
#   geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
#   facet_wrap(~Host_Harmonised, scales="free_y") +
#   theme_minimal()
# 
# ggplot(vcurves[ vcurves$TotalPathRich >= 20 & vcurves$Year <=2017 ,]) + 
#   geom_point(aes(Year, PathRich, size=NumRecords), alpha=0.3, col="skyblue4") +
#   geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
#   facet_wrap(~Host_Harmonised, scales="free_y") +
#   theme_minimal()
# 
# ggplot(vcurves[ vcurves$TotalPathRich >= 20 & vcurves$Year <=2017 ,]) + 
#   geom_point(aes(Year, PathRich, size=NumRecords), alpha=0.3, col="skyblue4") +
#   geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
#   facet_wrap(~Host_Harmonised, scales="free_y") +
#   theme_minimal()
# 
# 
# ggplot(vcurves2[ vcurves2$TotalPathRich >= 20 & vcurves2$Year <=2017 ,]) + 
#   geom_point(aes(Year, PathRich, size=NumRecords), alpha=0.3, col="skyblue4") +
#   geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
#   facet_wrap(~Host_Harmonised, scales="free_y") +
#   theme_minimal()
# 
# ggplot(vcurves2[ vcurves2$TotalPathRich >= 20 & vcurves2$Year <=2017 ,]) + 
#   geom_point(aes(CumRecords, PathRich), alpha=0.3, col="skyblue4") +
#   #geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
#   facet_wrap(~Host_Harmonised, scales="free") +
#   theme_minimal()


ggplot(curves[ curves$HostFamily == "Muridae" & curves$Year <= 2017, ]) + 
  geom_line(aes(CumPubs, PathRich), col="blue", size=0.2) +
  geom_point(aes(CumPubs, PathRich), col="blue", size=0.4) +
  facet_wrap(~Host_Harmonised, scales="free") +
  theme_minimal()

ggplot(curves[ curves$TotalPathRich >= 30 & curves$Year <= 2017, ]) + 
  geom_point(aes(CumPubs, PathRich), col="blue", size=0.2) +
  facet_wrap(~Host_Harmonised, scales="free") +
  theme_minimal()
ggplot(vcurves[ vcurves$TotalPathRich >= 10 & vcurves$Year <= 2017, ]) + 
  geom_point(aes(CumPubs, PathRich), col="blue", size=0.2) +
  facet_wrap(~Host_Harmonised, scales="free") +
  theme_minimal()

plotOrder = function(order, curves_df, num_spp=25){
  
  cc = curves_df[ curves_df$HostOrder == order, ]
  host_order = cc[ !duplicated(cc$Host_Harmonised), ] %>%
    arrange(desc(TotalPathRich))
  host_order = host_order[ 1:num_spp, ]
  cc = cc[ cc$Host_Harmonised %in% host_order$Host_Harmonised, ]
  cc$Host_Harmonised = factor(cc$Host_Harmonised, levels=host_order$Host_Harmonised, ordered=TRUE)
  cc = cc[ cc$Year<=2010, ]
  
  ggplot(cc) + 
    #geom_line(aes(CumPubs, PathRich), col="blue", size=0.2) +
    geom_point(aes(CumPubs, PathRich), col="blue", size=1) +
    facet_wrap(~Host_Harmonised, scales="free") +
    theme_minimal()
}

plotOrder("Primates", curves, 36)
plotOrder("Rodentia", vcurves, 36)
plotOrder("Chiroptera", vcurves, 36)
plotOrder("Artiodactyla", vcurves, 36)



sx = Sys.time()
# search_term = paste("anas", "[TIAB]", "AND", "platyrhynchos", "[TIAB]", sep=" ")
# search_term = paste("mastomys", "[TIAB]", "AND", "natalensis", "[TIAB]", sep=" ")
search_term = paste("aepyceros", "[TIAB]", "AND", "melampus", "[TIAB]", sep=" ")

e = simpleError("test error")
search1 = tryCatch(EUtilsSummary(search_term, type="esearch", db="pubmed", datetype='pdat', mindate=1950, maxdate=2015), error=function(e) e)
search2 = tryCatch(EUtilsSummary(search_term, type="esearch", db="pubmed", datetype='pdat', mindate=1950, maxdate=2015, retmax=QueryCount(search1)), error=function(e) e)
yy = YearPubmed(EUtilsGet(search2))
ex = Sys.time()
ex-sx
hist(yy, 70)

search_term = paste("anas", "[TIAB]", "AND", "platyrhynchos", "[TIAB]", sep=" ")
e = simpleError("test error")
search = tryCatch(RISmed::EUtilsSummary(search_term, type="esearch", db="pubmed", datetype='pdat', mindate=1950, maxdate=2015), error=function(e) e)
QueryCount(search)


install.packages("easyPubMed")
library(easyPubMed)
searchx = "mastomys[TIAB] AND natalensis[TIAB] AND 1950[PDAT] : 2015[PDAT]"
ids = get_pubmed_ids(searchx)
batch_pubmed_download(pubmed_query_string = searchx, 
                      format = "xml", 
                      batch_size = 150,
                      dest_file_prefix = "easyPM_example")
ll = articles_to_list("easyPM_example01.txt")
article_to_df(ll[[1]])







# ===================== Summarise by species and year ======================

# can specify by pathogentype, database, detection quality
temporalRarefact = function(database="all", pathogen_type="all", detection_quality=0){
  
  # subset
  if(database=="all"){ 
    ddx = dd
  } else{
    ddx = dd[ dd$Database %in% database, ]
  }
  if(pathogen_type=="all"){
    ddx = ddx
  } else{
    ddx = ddx[ ddx$ParasiteType %in% pathogen_type, ]
  }
  ddx = ddx[ ddx$DetectionQuality >= detection_quality, ]
  
  # compile by species starting in 1950
  curve_calc = function(spp){
    cat(paste(spp, "...", sep=""))
    ddy = ddx[ ddx$Host_Harmonised == spp & !is.na(ddx$Host_Harmonised), ]
    res = data.frame()
    for(y in 1920:2019){
      resx = data.frame(Host_Harmonised = spp, HostClass = ddy$HostClass[1], HostOrder = ddy$HostOrder[1], HostFamily=ddy$HostFamily[1],
                        Year = y,
                        PathRich = n_distinct(ddy$PathogenName_Harmonised[ ddy$Year <= y ]), # num pathogens by this year
                        NumRecords = length(ddy$PathogenName_Harmonised[ ddy$Year == y])) # num records in this year
      res = rbind(res, resx)
    }
    res$CumRecords = cumsum(res$NumRecords)
    return(res)
  }
  spp_curves = do.call(rbind.data.frame, lapply(unique(ddx$Host_Harmonised), curve_calc))
  
  # add total pathogen richness
  total_pr = ddx %>%
    group_by(Host_Harmonised) %>%
    dplyr::summarise(TotalPathRich = n_distinct(PathogenName_Harmonised))
  result = left_join(spp_curves, total_pr, by="Host_Harmonised")
  return(result)
}

# run across all species
#curves = temporalRarefact(database="all", pathogen_type = "all", detection_quality=0)
#write.csv(curves, "./output/data_processed/curves/curves_allpathogens_alldb_detection0_Oct2020.csv", row.names=FALSE)
curves = read.csv("./output/data_processed/curves/curves_allpathogens_alldb_detection0_Oct2020.csv", stringsAsFactors = FALSE)

# viruses
# vcurves = temporalRarefact(database="all", pathogen_type = "virus", detection_quality=0)
# write.csv(vcurves, "./output/data_processed/curves/curves_viruses_alldb_detection0_Oct2020_2.csv", row.names=FALSE)
vcurves = read.csv("./output/data_processed/curves/curves_viruses_alldb_detection0_Oct2020_2.csv", stringsAsFactors=FALSE)

# viruses strict detection
#vcurves2 = temporalRarefact(database="all", pathogen_type = "virus", detection_quality=2)
#write.csv(vcurves2, "./output/data_processed/curves/curves_viruses_alldb_detection2_Oct2020.csv", row.names=FALSE)
vcurves2 = read.csv("./output/data_processed/curves/curves_viruses_alldb_detection2_Oct2020.csv", stringsAsFactors=FALSE)

# combine with pubmed publications per species
curves = left_join(curves[ curves$Host_Harmonised %in% pp$Host_Harmonised, ], pp, by=c("Year", "Host_Harmonised"))
vcurves = left_join(vcurves[ vcurves$Host_Harmonised %in% pp$Host_Harmonised, ], pp, by=c("Year", "Host_Harmonised"))
vcurves2 = left_join(vcurves2[ vcurves2$Host_Harmonised %in% pp$Host_Harmonised, ], pp, by=c("Year", "Host_Harmonised"))


# plot curves
plotCurveExamples = function(spp){
  dfx = curves[ curves$Host_Harmonised == spp & curves$Year <= 2015, ]
  a = ggplot(dfx) + geom_line(aes(Year, PathRich), col="coral2", size=0.5) + 
    ggtitle(spp) +
    theme_minimal() + 
    ylab("Path richness") +
    theme(axis.title.x=element_blank(), plot.title=element_text(hjust=0.5, size=16)) 
  b = ggplot(dfx) + geom_line(aes(Year, CumPubs), col="blue", size=0.5) + 
    theme_minimal() + 
    ylab("Cumulative pubs") +
    theme(axis.title.x=element_blank()) 
  c = ggplot(dfx) + geom_line(aes(Year, NumPubs), col="blue", size=0.5) + 
    ylab("Annual pubs") +
    theme_minimal() + 
    theme(axis.title.x=element_text(size=13)) 
  px = gridExtra::grid.arrange(a, b, c, nrow=3)
  return(px)
}

p1 = plotCurveExamples("felis silvestris")
p2 = plotCurveExamples("crocuta crocuta")
p3 = plotCurveExamples("mastomys natalensis")
p4 = plotCurveExamples("eidolon helvum")
p5 = plotCurveExamples("canis latrans")
p6 = plotCurveExamples("macaca mulatta")

ppp = gridExtra::grid.arrange(grobs=list(p1, p2, p3, p4, p5, p6), nrow=2, ncol=3, height=1.2, width=1)
ggsave(ppp, file="./publication_curve_examples.png", device="png", dpi=300, width=12, height=14, units="in", scale=0.8)


plotCurveExamples("myodes glareolus")


vcx = vcurves[ vcurves$HostClass == "Mammalia" & vcurves$Host_Harmonised != "homo sapiens", ] %>%
  group_by(Host_Harmonised) %>%
  arrange(desc(TotalPathRich))
vcx = vcx[ vcx$Host_Harmonised %in% unique(vcx$Host_Harmonised)[1:20], ]
vcx$Host_Harmonised = factor(vcx$Host_Harmonised, levels=unique(vcx$Host_Harmonised), ordered=TRUE)

p1 = ggplot(vcx[ vcx$Year <= 2016, ]) + 
  geom_point(aes(Year, PathRich), col="skyblue4", size=1.5) +
  theme_classic() +
  facet_wrap(~Host_Harmonised, scales="free_y") +
  xlab("Year") +
  ylab("Viral richness") + 
  ggtitle("Viral richness by year") +
  theme(plot.title = element_text(hjust=0.5, size=14))
p2 = ggplot(vcx[ vcx$Year <= 2016, ]) + 
  geom_point(aes(CumPubs, PathRich, col=Year), size=1.5) +
  theme_classic() +
  facet_wrap(~Host_Harmonised, scales="free") +
  xlab("Cumulative publications") +
  ylab("Viral richness") +
  ggtitle("Viral richness by cumulative publications") +
  theme(plot.title = element_text(hjust=0.5, size=14))
ggsave(p1, file="./output/figures/AllDatasets_ViralRichnessByYear.png", device="png", units="in", width=9, height=6, scale=1)
ggsave(p2, file="./output/figures/AllDatasets_ViralRichnessByCumPubs.png", device="png", units="in", width=9, height=6, scale=1)



