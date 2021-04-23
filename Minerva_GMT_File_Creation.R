# Create COVID19 Disease Map GMT file
# Author: mkutmon
# Update: 2021-04-16

library(here)
setwd(paste0(here(),"/DiseaseMapGMTFile"))


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
  data[1,1] <- paste0(name,"%MINERVA_April2021%",name,"%Homo sapiens")
  data[2,1] <- "MINERVA_COVID19"
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
    write.table(data_t, "COVID19_DiseaseMap_April2021.gmt",row.names = F, col.names = F, append = T, sep="\t", quote=F)
  }
}
