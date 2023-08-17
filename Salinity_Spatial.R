library(rgbif)
library(maps)
library(reshape2)
library(plyr)
library(data.table)
library(tidyverse)
library(fishtree)
library(raster)
library(ncdf4)
library(viridis)
library(foreach)
library(doParallel)
library(rgdal)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(lwgeom)

#--------- data ---------# 

#FishBase Data

fb.dat <- readRDS("data/fishbase_data.RDS")


#--------- data ---------# 
# Query GBIF
## Takes a long time, but reduced with dopar and foreach
## Change to T to run
if(F){
cl <- makePSOCKcluster(detectCores()-2) # leave two cores for background processes
# cl <- makePSOCKcluster(1) # if old PC use only 1 core
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers


system.time(
  gbif.dat <- foreach(i=rob.tree$tip.label,.verbose=T, .packages = c("data.table", "rgbif")) %dopar% {
    try(
      occ_search(scientificName = i, hasCoordinate=T, limit =1000,hasGeospatialIssue = F,basisOfRecord = "PRESERVED_SPECIMEN")$data,silent = T
    )
  }
)

try({gbif.i <- occ_search(scientificName = i, hasCoordinate=T, limit =1000);
message('\r', paste0(last,"/",max.sp, " species"), appendLF = FALSE)}
,silent = T)
if(!inherits(gbif.i,"try-error")){
  gbif.i <- NA
}
stopCluster(cl) # stop parallel backend


saveRDS(gbif.dat,"data/GBIF_data.RDS")
}


# Filter and save

gbif.dat <- readRDS("data/GBIF_data.RDS")


gbif.cols <-c("name","genus","species","acceptedScientificName","scientificName","key","decimalLatitude","decimalLongitude","speciesKey","depth","depthAccuracy","habitat","year","month","day")

g.l <- list() 
for(n in 1:length(gbif.dat)){
  d.n <- gbif.dat[[n]]
  
  if(length(d.n)>1){
    null.c <- sapply(gbif.cols,function(x) x %in% colnames(d.n))
    names.c <- names(null.c)[!null.c]
    d.nn <- data.table(d.n)
    d.nnn <- if(length(names.c)>0){d.nn[,c(names.c):=NA]}else(d.nn)
    g.l[[n]] <- d.nnn[,..gbif.cols]
  }else{NA
    
    nulls <- data.table()
    g.l[[n]] <- nulls[,c(gbif.cols):=NA]
  }
  
}


gbif <- do.call(rbind,g.l)

gbif.sum <- gbif[,.N,by=.(species)]


gbif.sum%>%
  ggplot(aes(x=N))+geom_histogram(binwidth = 50)

saveRDS(gbif,"data/GBIF_datatable.RDS")


gbif <- readRDS("data/GBIF_datatable.RDS")

#reduced for quick ploting
g.red <- gbif%>%
  group_by(scientificName)%>%
  dplyr::summarise(m.lat=mean(decimalLatitude),m.long=mean(decimalLongitude),n=length(decimalLongitude))


#--------- Maps ---------# 

# Establish maps for checking records.
#These were used extensively to check for issues in GBIF data and FishBase salinity categories

world <- map_data("world")

#using maps
world.map <- ggplot(world, aes(long, lat)) 

#using rnaturalearth
world2 <- rnaturalearth::ne_countries(scale = 'small', returnclass = c("sf"))

coast <- rnaturalearth::ne_coastline(scale = 'medium', returnclass = c("sf"))

world.map2 <- ggplot() +
  geom_sf(data = world2, size = .2, fill = "gray80", col = "gray90") +
  theme(panel.grid.major = element_line(color = gray(0.9), linetype = "dashed", size = 0.5))

coast.map <- ggplot() +
  geom_sf(data = coast, size = .3, col = "black") +
  theme(panel.grid.major = element_line(color = gray(0.9), linetype = "dashed", size = 0.5))

world.map+geom_polygon(aes(group = group), fill = "white", color = "gray40", size = .2) + geom_jitter(data = g.red,aes(m.long, m.lat), alpha=0.6, size = 1, color = "red") 



#--------- Spatial analysis ---------# 

# Load marine and estuarine shape files
# Using the function below, we can figure out if a record is from over the shelf break or in a pelagic realm. Also used to find TNC major habibat type with a different shape (`feow`)

#Spherical geometry (s2) switched off, otherwise polygon errors
sf_use_s2(FALSE)

#marine shapes
meow <- sf::read_sf("data/spatial/MEOW_PPOW/01_Data/WCMC-036-MEOW-PPOW-2007-2012.shp") 

#fw shapes
feow <- sf::read_sf("data/spatial/FEOW-TNC/FEOWv1_TNC.shp")

#simplify for size/time sake
meow2 <- meow %>% 
  sf::st_simplify(dTolerance = 01)

meow2 <- meow2 %>% 
  dplyr::group_by(PROVINC, REALM, TYPE) %>% 
  dplyr::summarise()

object.size(meow)
object.size(meow2)

#load estuaries shape (where from)
estuaries <- sf::read_sf("data/spatial/Estuaries2003/01_Data/14_001_UBC003_SAU_Estuaries2003_v2.shp")

#load countries data for inland/freshwater designation
countries <- ne_download(scale="large",type="countries",returnclass = "sf")



# Custom function to find over which polygon record resides

getRegionalInfo  <- function(lat1, long1, shape = meow) {
  #first, extract the co-ordinates (x,y - i.e., Longitude, Latitude)
  coords <- cbind(long1, lat1)
  
  df <-
    SpatialPointsDataFrame(coords, data.frame(value = paste(1:nrow(coords))))
  
  proj4string(df) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  
  df <- st_as_sf(df,
                 coords = c(x = "X", y = "Y"),
                 crs = "+proj=longlat +datum=WGS84 +no_defs")
  
  dsdat <- df %>% 
    st_join(shape)%>%
    mutate(dupl=duplicated(value))%>%
    filter(dupl==F)
  
  #due to simplification, records on shelf margin may be returned for both a MEOW and PEOW. Arbitrarily chose pelagic
  
  return(dsdat)
  
}




## Global marine species

#to make plotting faster
options(bitmapType = "cairo")

#map meows
meow.map <- world.map2+
  geom_sf(data = meow%>%filter(TYPE=="MEOW"),fill="gray10",size = .2, alpha=1,inherit.aes = F)


#data table to speed up calculations

fb.gbif <- fb.dat%>%
  full_join(tibble(copy(gbif))%>%dplyr::select(species,decimalLatitude,decimalLongitude)%>%dplyr::rename(Species=species))%>%
  mutate(rec=1:nrow(.))%>%
  as.data.table()


#ID records as neritic with meow ecoregion or with biome from PPOW to find pelagics. Also pass through and inland based on shape data polygons
mht <-  fb.gbif[!is.na(decimalLatitude),{getRegionalInfo(lat1=decimalLatitude,long1 = decimalLongitude,shape=meow%>%dplyr::select(TYPE))%>%mutate(Species=Species,mar.fresh2=mar.fresh2,rec=rec)}]%>%  
  st_join(countries%>%dplyr::select(CONTINENT))%>% #is record inland
  st_join(estuaries%>%dplyr::select(LABEL))%>% #is record estuarine
  st_join(feow%>%dplyr::select(MHT_TXT)) #is record fw

#save data
saveRDS(list(mht=mht,map=meow.map,countries=countries),"data/mht_data.RDS")


