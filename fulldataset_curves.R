

# ====================== Plot curves for species pathogen rarefaction ===========================

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiversity_rarefaction/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2")

# full harmoniesd associations database (EID2, GMPD2, HP3)
dd = read.csv("./output/data_processed/hostpathogen_harmonised/AllDatabases_Associations_Hosts_Harmonised_Oct2020.csv", stringsAsFactors = FALSE)
dd = dd[ !is.na(dd$Database), ]
dd = dd[ !is.na(dd$ParasiteType), ]
dd = dd[ !is.na(dd$Host_Harmonised), ]

# publication effort by year
pubs = read.csv("./output/data_processed/host_effort/PubMed_Hosts_PubsByYear.csv", stringsAsFactors = FALSE)
pp = expand.grid(1950:2020, unique(dd$Host_Harmonised)) %>%
  rename("Year"=Var1, "Host_Harmonised"=Var2)


# # create dataframe for everything
# dd = expand.grid(1950:2020, unique(assoc$Host_Harmonised)) %>%
#   rename("Year"=Var1, "Host_Harmonised"=Var2) %>%
#   left_join(assoc[ !duplicated(assoc$Host_Harmonised), c("Host_Harmonised", "HostClass", "HostOrder", "HostFamily")], by=c("Host_Harmonised")) %>%
#   left_join(assoc[ , c("Host_Harmonised", "Year", "PathogenName_Harmonised", "ParasiteType", "Database", "DetectionQuality", "HumanInfective_Any", "IsZoonotic") ], by=c("Host_Harmonised", "Year"))
# 

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
    for(y in 1950:2019){
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
#vcurves = temporalRarefact(database="all", pathogen_type = "virus", detection_quality=0)
#write.csv(vcurves, "./output/data_processed/curves/curves_viruses_alldb_detection0_Oct2020.csv", row.names=FALSE)
vcurves = read.csv("./output/data_processed/curves/curves_viruses_alldb_detection0_Oct2020.csv", stringsAsFactors=FALSE)

# viruses strict detection
#vcurves2 = temporalRarefact(database="all", pathogen_type = "virus", detection_quality=2)
#write.csv(vcurves2, "./output/data_processed/curves/curves_viruses_alldb_detection2_Oct2020.csv", row.names=FALSE)
vcurves2 = read.csv("./output/data_processed/curves/curves_viruses_alldb_detection2_Oct2020.csv", stringsAsFactors=FALSE)

# combine with pubmed publications per species



# viz
ggplot(curves[ curves$TotalPathRich >= 40 & curves$Year <=2017 ,]) + 
  geom_point(aes(Year, PathRich, size=NumRecords), alpha=0.3, col="skyblue4") +
  geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
  facet_wrap(~Host_Harmonised, scales="free_y") +
  theme_minimal()

ggplot(curves[ curves$HostOrder == "Rodentia" & curves$TotalPathRich >= 5 & curves$Year <=2017 ,]) + 
  geom_point(aes(Year, PathRich, size=NumRecords), alpha=0.3, col="skyblue4") +
  geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
  facet_wrap(~Host_Harmonised, scales="free_y") +
  theme_minimal()

ggplot(vcurves[ vcurves$TotalPathRich >= 20 & vcurves$Year <=2017 ,]) + 
  geom_point(aes(Year, PathRich, size=NumRecords), alpha=0.3, col="skyblue4") +
  geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
  facet_wrap(~Host_Harmonised, scales="free_y") +
  theme_minimal()

ggplot(vcurves[ vcurves$TotalPathRich >= 20 & vcurves$Year <=2017 ,]) + 
  geom_point(aes(Year, PathRich, size=NumRecords), alpha=0.3, col="skyblue4") +
  geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
  facet_wrap(~Host_Harmonised, scales="free_y") +
  theme_minimal()


ggplot(vcurves2[ vcurves2$TotalPathRich >= 20 & vcurves2$Year <=2017 ,]) + 
  geom_point(aes(Year, PathRich, size=NumRecords), alpha=0.3, col="skyblue4") +
  geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
  facet_wrap(~Host_Harmonised, scales="free_y") +
  theme_minimal()

ggplot(vcurves2[ vcurves2$TotalPathRich >= 20 & vcurves2$Year <=2017 ,]) + 
  geom_point(aes(CumRecords, PathRich), alpha=0.3, col="skyblue4") +
  #geom_line(aes(Year, PathRich), col="coral2", size=0.5) +
  facet_wrap(~Host_Harmonised, scales="free") +
  theme_minimal()





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
