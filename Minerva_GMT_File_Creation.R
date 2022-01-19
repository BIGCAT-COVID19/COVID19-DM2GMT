# Create COVID19 Disease Map GMT file
# Author: mkutmon
# Update: 2022-01-19
# Creates a GMT file with Entrez Gene identifiers from the COVID19 disease map (integrating maps from Minerva, WikiPathways and Reactome)

library(here)
library(rWikiPathways)
library(org.Hs.eg.db)

setwd(paste0(here(),"/DiseaseMapGMTFile"))

release.minerva <- "17.06.21"
version.wikipathways <- "19.01.22"
version.reactome <- "19.01.22"

# RETRIEVE PATHWAYS AND GENE IDENTIFIERS FROM COVID19 DISEASE MAP
source("https://git-r3lab.uni.lu/covid/models/-/raw/master/Integration/MINERVA_access/minerva_access_functions.R")

map <- "https://covid19map.elixir-luxembourg.org/minerva/api/"
map_components <- get_map_components(map)
mnv_refs <- get_components_annotations(map_components, "ENTREZ")


# FORMATTING OF DATA INTO GMT FORMAT
names(mnv_refs) <- paste0("MINERVA ", names(mnv_refs))

for (i in 1:length(mnv_refs)) {
  name <- names(mnv_refs[i])
  data <- data.frame(matrix(ncol = 1, nrow = 0))
  data[1,1] <- paste0(name)
  data[2,1] <- paste0("MINERVA_COVID19_",release.minerva)
  x <- 3
  for(j in 1:length(mnv_refs[[i]])) {
    val <- mnv_refs[[i]][[j]][1]
    if(!is.na(val)) {
      data[x,1] <- val
      x <- x+1
    }
  }
  if(dim(data)[1] > 2) {
    data_t <- t(data)
    write.table(data_t, paste0("output1/COVID19_DiseaseMap-20220119.gmt"),row.names = F, col.names = F, append = T, sep="\t", quote=F)
  }
}

# RETRIEVE PATHWAYS AND GENE IDENTIFIERS FROM WIKIPATHWAYS COVID PORTAL

wp.pathways <- rWikiPathways::getPathwayIdsByCurationTag("Curation:COVID19")
for(p in wp.pathways) {
  # exclude Reactome pathways (retrieved from Reactome directly)
  if(p != "WP5020" & p != "WP5021") {
    info <- rWikiPathways::getPathwayInfo(p)
    if(info$species == "Homo sapiens") {
      genes <- rWikiPathways::getXrefList(p,"L")
      data <- data.frame(matrix(ncol = 1, nrow = 0))
      data[1,1] <- paste0("WP_",info[3])
      data[2,1] <- paste0("WIKIPATHWAYS%",p,"%",version.wikipathways)
      x <- 3
      for(j in 1:length(genes)) {
        data[x,1] <- genes[j]
        x <- x+1
      }
      data_t <- t(data)
      write.table(data_t, paste0("output1/COVID19_DiseaseMap-20220119.gmt"),row.names = F, col.names = F, append = T, sep="\t", quote=F)
    }
  }
}

# RETRIEVE PATHWAYS AND GENE IDENTIFIERS FROM REACTOME

parent.pathways <- c("R-HSA-9694516","R-HSA-9678108")
#submaps <- c("R-HSA-9694614","R-HSA-9694676","R-HSA-9694682","R-HSA-9694635","R-HSA-9694322")

temp <- tempfile()
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",temp)
pathways <- read.table(unz(temp,"ReactomePathways.gmt"), sep="\t", header=FALSE,col.names = paste0("V",seq_len(4000)), quote="", fill=TRUE)
pathways <- pathways[pathways$V2 %in% parent.pathways,]
empty_columns <- colSums(is.na(pathways) | pathways == "") == nrow(pathways)
pathways <- pathways[, !empty_columns]

for (i in 1:nrow(pathways)) {
  name <- pathways[i,1]
  data <- data.frame(matrix(ncol = 1, nrow = 0))
  data[1,1] <- paste0("REACTOME_",name)
  data[2,1] <- paste0("REACTOME%",pathways[i,2],"%",version.reactome)
  
  genes <- t(pathways[i,3:ncol(pathways)])
  genes <- genes[genes[,1] != "",]
  ids <- select(org.Hs.eg.db, 
                 keys = as.vector(genes),
                 columns = c("ENTREZID", "SYMBOL"),
                 keytype = "SYMBOL")
  ids <- ids[!is.na(ids$ENTREZID),]
  x <- 3
  for(id in ids$ENTREZID) {
    data[x,1] <- id
    x <- x+1
  }
  if(dim(data)[1] > 2) {
    data_t <- t(data)
    write.table(data_t, paste0("output1/COVID19_DiseaseMap-20220119.gmt"),row.names = F, col.names = F, append = T, sep="\t", quote=F)
  }
}
