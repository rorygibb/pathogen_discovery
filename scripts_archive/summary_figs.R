

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
vir0 = vir

# subset to keep only years pre endyear
endyear = 2018
vir = vir[ vir$Year <= endyear, ]



# ================================

# summarise 

vir %>%
  dplyr::select(Host, HostOrder, Virus, Year, Domestic) %>%
  distinct() %>%
  group_by(Host, Virus) %>%
  summarize(Year = min(Year),
            Domestic = head(Domestic, 1),
            HostOrder = head(HostOrder, 1)) -> tx

hh = tx %>%
  dplyr::filter(Domestic == FALSE) %>%
  dplyr::group_by(HostOrder) %>%
  dplyr::summarise(HH = n_distinct(Host))

# num hosts, viruses, associations
n_distinct(vir$Host)
n_distinct(vir$Virus)
nrow(tx)

# virus status
unique(vir$Virus[ vir$PredictFlag == TRUE])
unique(vir$Virus[ vir$Database=="PREDICT"])

# ratified
nrow(vir %>% dplyr::filter(ICTVRatified==TRUE) %>% dplyr::select(Virus) %>% distinct())
nrow(vir %>% dplyr::filter(PredictFlag==TRUE) %>% dplyr::select(Virus) %>% distinct())


p1 = vir0 %>%
  dplyr::select(Host, Virus, Year) %>%
  distinct() %>%
  group_by(Host, Virus) %>%
  summarize(Year = min(Year)) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(Discovered = length(Year),
                   `Host species` = n_distinct(Host)) %>%
  ggplot() + 
  geom_line(aes(Year, Discovered), size=0.6, col="coral4") +
  geom_point(aes(Year, Discovered, size=`Host species`), alpha=0.5, col="darkred") +
  theme_classic() + 
  ylab("Novel host-virus associations") +
  xlab("Year") + 
  theme(axis.text = element_text(size=11),
        axis.title=element_text(size=12), 
        legend.position=c(0.09, 0.55),
        legend.title = element_text(size=10.5)) +
  scale_x_continuous(breaks=seq(1930, 2020, by=15), labels=seq(1930, 2020, by=15))
  
p0 = vir0 %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(NRecs = length(Year)) %>%
  ggplot() + 
  geom_line(aes(Year, NRecs), size=0.6, col="skyblue4") +
  geom_point(aes(Year, NRecs), size=1.5, alpha=0.5, col="grey20") +
  theme_classic() + 
  ylab("Total number of records") +
  xlab("Year") + 
  theme(axis.text = element_text(size=11),
        axis.title=element_text(size=12)) +
  scale_x_continuous(breaks=seq(1930, 2020, by=15), labels=seq(1930, 2020, by=15))

pc = gridExtra::grid.arrange(p0, p1)
ggsave(pc, file="./output/figures_2021/SI_OverallVIRIONTrends.jpeg", units="in", height=5.5, width=9.5, dpi=600)


vir0 %>%
  dplyr::select(Host, Virus, Database, Year) %>%
  distinct() %>%
  group_by(Host, Virus, Database) %>%
  summarize(Year = min(Year)) %>%
  dplyr::group_by(Database, Year) %>%
  dplyr::summarise(Discovered = length(Year)) %>%
  ggplot() + 
  geom_line(aes(Year, Discovered, col=Database), size=0.8) +
  theme_classic() + 
  ylab("Novel host-virus associations") +
  xlab("Year") + 
  theme(axis.text = element_text(size=11),
        axis.title=element_text(size=12)) +
  scale_x_continuous(breaks=seq(1930, 2020, by=15), labels=seq(1930, 2020, by=15))

px = vir0 %>%
  dplyr::select(Host, Virus, Database, Year) %>%
  distinct() %>%
  group_by(Host, Virus, Database) %>%
  summarize(Year = min(Year)) %>%
  dplyr::group_by(Database, Year) %>%
  dplyr::summarise(Discovered = length(Year)) %>%
  ggplot() + 
  geom_bar(aes(Year, Discovered, fill=Database), stat="identity") +
  theme_classic() + 
  ylab("Novel host-virus associations") +
  xlab("Year") + 
  theme(axis.text = element_text(size=11),
        axis.title=element_text(size=12)) +
  scale_x_continuous(breaks=seq(1930, 2020, by=15), labels=seq(1930, 2020, by=15))
ggsave(px, file="./output/figures_2021/SI_OverallVIRIONTrends_byDatabase.jpeg", units="in", height=2.8, width=8, dpi=600)
