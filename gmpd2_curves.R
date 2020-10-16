setwd("/Users/roryj/Documents/PhD/202008_PathogenDiversity/")
library(stringr)
library(dplyr)
library(ggplot2)

gg = read.csv("./data/GMPD2/GMPD2_datafiles/GMPD_main.csv", stringsAsFactors = FALSE)

# extract numerics (year)
years = regmatches(gg$Citation, gregexpr("[[:digit:]]+", gg$Citation))

# trim out second year marker
gg$Year = as.numeric(unlist(lapply(years, "[", 1)))
gg = gg[ !is.na(gg$Year), ]

# couple of refs w/o year
# unique(gg$Citation[ is.na(gg$Year) ])

# calculate species pathogen richness
prich = gg %>%
  group_by(HostCorrectedName) %>%
  dplyr::summarise(Order = head(HostOrder, 1), 
                   Family = head(HostFamily, 1),
                   ParasiteDiv = length(unique(ParasiteCorrectedName)),
                   VirDiv = length(unique(ParasiteCorrectedName[ ParType == "Virus" ])),
                   NonVirDiv = length(unique(ParasiteCorrectedName[ ParType != "Virus" ])))
prich = prich[ prich$HostCorrectedName != "no binomial name", ]

# temporal rarefaction function
tempRarefact = function(spp, partype){
  
  # calculate curve
  if(partype == "All"){  ggx = gg[ gg$HostCorrectedName == spp, ] }
  if(partype == "Virus"){  ggx = gg[ gg$HostCorrectedName == spp & gg$ParType == "Virus", ] }
  if(partype == "NonVirus"){  ggx = gg[ gg$HostCorrectedName == spp & gg$ParType != "Virus", ] }
  
  #resx = data.frame(Year = min(gg$Year, na.rm=T):max(gg$Year, na.rm=T), ParasiteDiv = NA)
  resx = data.frame(Year = min(ggx$Year, na.rm=T):max(ggx$Year, na.rm=T), ParasiteDiv = NA, NumPubs = NA)
  
  for(i in 1:nrow(resx)){
    ggxx = ggx[ ggx$Year <= resx$Year[i], ]
    resx$ParasiteDiv[i] = length(unique(ggxx$ParasiteCorrectedName))
    resx$NumPubs[i] = length(unique(ggx$Citation[ ggx$Year == resx$Year[i] ]))
  }
  
  # add metadata
  resx$Species = spp
  resx$Order = ggx$HostOrder[1]
  resx$Family = ggx$HostFamily[1]
  resx$TotalParasiteDiv = length(unique(ggx$ParasiteCorrectedName))
  resx
}

# run for all parasites
px = prich[ order(prich$ParasiteDiv, decreasing = TRUE), ]
result = do.call(rbind.data.frame, lapply(px$HostCorrectedName[ 1:36 ], tempRarefact, partype="All"))
result = do.call(rbind.data.frame, lapply(px$HostCorrectedName[ 37:(36*2) ], tempRarefact, partype = "All"))
result = do.call(rbind.data.frame, lapply(px$HostCorrectedName[ 1:100 ], tempRarefact, partype = "All"))

result$Species = factor(result$Species, levels=unique(result$Species), ordered=TRUE)
plot1 = ggplot(result, aes(Year, ParasiteDiv, fill=NumPubs, size=NumPubs)) + 
  geom_point(pch=21) +
  facet_wrap(~Species, scales="free") + 
  theme_classic() + 
  scale_fill_viridis_c(option="viridis", begin=0.1, end=0.95) +
  theme(strip.background = element_blank(),
        plot.title = element_text(size=15, hjust=0.5)) + 
  ggtitle("All parasite/pathogen richness (GMPD2)")

# viruses
px = prich[ order(prich$VirDiv, decreasing = TRUE), ]
result = do.call(rbind.data.frame, lapply(px$HostCorrectedName[ 1:36 ], tempRarefact, partype = "Virus"))
# result = do.call(rbind.data.frame, lapply(px$HostCorrectedName[ 37:(36*2) ], tempRarefact, partype = "Virus"))
# result = do.call(rbind.data.frame, lapply(px$HostCorrectedName[ 1:100 ], tempRarefact, partype = "Virus"))
result$Species = factor(result$Species, levels=unique(result$Species), ordered=TRUE)
plot2 = ggplot(result, aes(Year, ParasiteDiv, fill=NumPubs, size=NumPubs)) + 
  geom_point(pch=21) +
  facet_wrap(~Species, scales="free") + 
  theme_classic() + 
  scale_fill_viridis_c(option="viridis", begin=0.1, end=0.95) +
  theme(strip.background = element_blank(),
        plot.title = element_text(size=15, hjust=0.5)) + 
  ylab("VirDiv") +
  ggtitle("Virus richness (GMPD2)")

# viruses
px = prich[ order(prich$NonVirDiv, decreasing = TRUE), ]
result = do.call(rbind.data.frame, lapply(px$HostCorrectedName[ 1:36 ], tempRarefact, partype = "NonVirus"))
result$Species = factor(result$Species, levels=unique(result$Species), ordered=TRUE)
plot3 = ggplot(result, aes(Year, ParasiteDiv, fill=NumPubs, size=NumPubs)) + 
  geom_point(pch=21) +
  facet_wrap(~Species, scales="free") + 
  theme_classic() + 
  scale_fill_viridis_c(option="viridis", begin=0.1, end=0.95) +
  theme(strip.background = element_blank(),
        plot.title = element_text(size=15, hjust=0.5)) + 
  ylab("ParasiteDiv")  + 
  ggtitle("Non-Virus richness (GMPD2)")

# save
ggsave(plot1, file = "./output/figures/GMPD_ParasiteRich.png", device="png", width=12, height=10, units="in", dpi=600, scale=0.9)
ggsave(plot2, file = "./output/figures/GMPD_VirusRich.png", device="png", width=12, height=10, units="in", dpi=600, scale=0.9)
ggsave(plot3, file = "./output/figures/GMPD_NonVirusRich.png", device="png", width=12, height=10, units="in", dpi=600, scale=0.9)


# top 100
px = prich[ order(prich$ParasiteDiv, decreasing = TRUE), ]
result = do.call(rbind.data.frame, lapply(px$HostCorrectedName[ 1:100 ], tempRarefact, partype = "All"))
result$Species = factor(result$Species, levels=unique(result$Species), ordered=TRUE)
plot4 = ggplot(result, aes(Year, ParasiteDiv, fill=NumPubs, size=NumPubs)) + 
  geom_point(pch=21) +
  facet_wrap(~Species, scales="free") + 
  theme_classic() + 
  scale_fill_viridis_c(option="viridis", begin=0.1, end=0.95) +
  theme(strip.background = element_blank(),
        plot.title = element_text(size=15, hjust=0.5)) + 
  ggtitle("All parasite/pathogen richness (GMPD2)")
ggsave(plot4, file = "./output/figures/GMPD_ParasiteRich_100.png", device="png", width=22, height=18, units="in", dpi=600, scale=0.9)

px = prich[ order(prich$VirDiv, decreasing = TRUE), ]
result = do.call(rbind.data.frame, lapply(px$HostCorrectedName[ 1:100 ], tempRarefact, partype = "Virus"))
result$Species = factor(result$Species, levels=unique(result$Species), ordered=TRUE)
plot5 = ggplot(result, aes(Year, ParasiteDiv, fill=NumPubs, size=NumPubs)) + 
  geom_point(pch=21) +
  facet_wrap(~Species, scales="free") + 
  theme_classic() + 
  scale_fill_viridis_c(option="viridis", begin=0.1, end=0.95) +
  theme(strip.background = element_blank(),
        plot.title = element_text(size=15, hjust=0.5)) + 
  ylab("VirDiv") +
  ggtitle("Virus richness (GMPD2)")
ggsave(plot5, file = "./output/figures/GMPD_VirusRich_100.png", device="png", width=22, height=18, units="in", dpi=600, scale=0.9)

