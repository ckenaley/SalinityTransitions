library(tidyverse)
library(randomForest)
library(data.table)
library(sf)
library(units)
data(levitus, package="ocedata")
attach(levitus)
library(ggOceanMaps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(lwgeom)
library(CoordinateCleaner)
library(parallel)
library(doParallel)
library(worms)
library(h2o)
library(fishtree)
library(xtable)

#--------- data and functions ---------# 

### function for mapping


world2 <- rnaturalearth::ne_countries(scale = 'small', returnclass = c("sf"))
world.map2 <- ggplot() +
  geom_sf(data = world2, size = .2, fill = "gray80", col = "gray90") +
  theme(panel.grid.major = element_line(color = gray(0.9), linetype = "dashed", size = 0.5))



map_rec <- function(r = NULL) {
  if (!is.numeric(r))
    r <- mht2[Species == r]$rec
  
  c <- mht2[rec %in% r, geometry] %>% st_coordinates()
  
  ln <- range(c[, 1])
  lt <- range(c[, 2])
  box <- 10
  world.map2 +
    stat_sf_coordinates(
      data = mht2 %>% dplyr::filter(rec %in% r),
      aes(geometry = geometry, col = sal),
      col = "red"
    ) + xlim(range(ln - box, ln + box)) + ylim(range(lt - box, lt + box))
  
}

phy <- fishtree_phylogeny()

phy$tip.label <- gsub("_"," ",phy$tip.label)
phy <- keep.tip(phy,d$Species)

#cities data from
#https://geo.nyu.edu/catalog/stanford-yk247bg4748


#fw data from http://data.freshwaterbiodiversity.eu/
#paper https://www.nature.com/articles/sdata2017141

fw <- read_delim("data/fw.csv",delim=";")

colnames(fw) <- colnames(fw) %>% gsub("\\.","_",.) %>% gsub("\\d_","",.)

fw <- fw %>% 
  mutate_at(vars(Fishbase_Valid_Species_Name),function(x) gsub("\\."," ",x)) %>% 
  dplyr::rename(Species=Fishbase_Valid_Species_Name,fw_status=Occurrence_Status) %>% 
  dplyr::select(Species,fw_status) %>% 
  group_by(Species) %>% 
  filter(fw_status=="valid"&!duplicated(Species)) %>% 
  arrange(Species)

#data from GBIF analysis in script "XXXX"

mht.dat <- readRDS("data/mht_data.RDS") 
mht <- mht.dat$mht
meow.map <- mht.dat$map
countries <- mht.dat$countries

#salinity
sss <- expand_grid(lat=levitus$latitude,lon=levitus$longitude) 
sss$sal <- levitus$SSS %>% matrix(ncol=1) %>% unlist()


#Spherical geometry (s2) switched off, otherwise polygon errors
sf_use_s2(FALSE)

mht2 <- copy(mht)

#remove duplicated position records
mht2[,dupl:=duplicated(geometry),by=.(Species)]
mht2 <- copy(mht2[dupl==F,])

mht2[,lon:=st_coordinates(geometry)[,1] ]
mht2[,lat:=st_coordinates(geometry)[,2] ]


#flagging coords according to tests (outliers, at institutions, etc.)
mht2$flagged <- clean_coordinates(x = mht2,
                           lon = "lon",
                           lat = "lat",
                           species = "Species",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                     "zeros","outliers"),
                          value="flagged") # most test are on by default


mht2 <- copy(mht2[flagged==T,])

#for salinity data
mht2[,lon:=st_coordinates(geometry)[,1] %>% plyr::round_any(0.5)]
mht2[,lat:=st_coordinates(geometry)[,2] %>% plyr::round_any(0.5)]


#cities
cities<- sf::read_sf("data/spatial/stanford_cities/ne_10m_urban_areas_landscan.shp") 

cities2 <- cities  %>% 
  dplyr::group_by(name_conve) %>% 
  dplyr::summarise() %>% 
  sf::st_simplify(dTolerance = 0.1)

cities_poly <- sf:::as_Spatial(cities)

city_cent <- cities %>% group_by(name_conve) %>% 
  st_centroid


dist_city <- function(x,y){
  d <- st_distance(x %>% st_centroid(),y)
  d2 <- dist2land2(x,bin = F,shp=cities_poly[which.min(d[1,]),])
  return(d2)
}

countries2 <- countries %>% 
  dplyr::group_by(CONTINENT) %>% 
  dplyr::summarise() %>% 
  sf::st_simplify(dTolerance = 1)

countries2 <- st_cast(countries, "MULTILINESTRING")

nec <- ne_countries("large")

dist2land2 <- function(x,bin=FALSE,shp=nec){
  xy <- x %>% unlist
  xy <- data.frame(lon=xy[1],lat=xy[2])
  colnames(xy) <- c("lon","lat")
  res <- ggOceanMaps::dist2land(xy,shapefile = shp,verbose=F,binary=bin)$ldist 
  return(res)
}


#--------- spatial analysis  ---------# 

#use shapes to estimate mar, fresh, or estuary
#the really consuming parts, took ~15h on 8 cores

cl <- makePSOCKcluster(detectCores()-2) # leave two cores for background processes
# cl <- makePSOCKcluster(1) # if old PC use only 1 core
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers

mht_par <- function(x=mht2,ii=i){
  x2 <- copy(x[Species==ii,])
  x2[,dtco:=min(st_distance(geometry,countries2)),by=.(rec)]
  
  #distance to city, uses centroid of cities to find closet city shape, dist computed as distance from urban area
  x2[,dtc:=dist_city(x=geometry,y=city_cent),by=.(rec)]
  
  #distance to coastline
  x2[,dts:=dist2land2(geometry),by=.(rec)]
  return(x2)
}

# change to T to run
if(F){
    mht.l <- foreach(i=mht2$Species %>% unique,.verbose=T, .packages = c("data.table", "sf","ggOceanMaps"),.errorhandling = "pass") %dopar% {
  #distance to coast
  mht_par(ii=i)
}

stopCluster(cl) # stop parallel backend

        
        
saveRDS(mht.l,"data/mht_par_results.RDS")
}

mht.l <- readRDS("data/mht_par_results.RDS") 

#remove null results
probs <- (1:length(mht.l))[sapply(mht.l,function(x) !"data.table" %in% class(x))]#just last NA for species
mht2 <- do.call(rbind,mht.l[-probs])

#add salinity
mht2 <- mht2 %>% 
  left_join(sss)

#is vector on land
mht2[,onland:=ifelse(dts==0,T,F),by=.(rec)]


#change to T to run
if(F){
#worms_l <- list()
sp <- mht2$Species %>%  unique
sp_n <- seq(1,mht2$Species %>% unique %>% length,50)
l=2
 for(i in l:(-length(sp_n)-1)){
  s <- sp_n[i-1]
  e <- sp_n[i]-1
  print(sp[s:e])
  sp.i <- sp[s:e][!sp[s:e]%in%bad]
  #sp.i <- sp[10951:11001][!sp[s:e]%in%bad]
  worms_l[[i]] <- wormsbymatchnames(taxon_names = sp.i ,verbose =F, ids = FALSE,marine_only = "false",sleep_btw_chunks_in_sec = 0.01)
  l=i
 }

worms_rec <- rbindlist(worms_l) %>% 
  dplyr::rename(Species=scientificname) %>% 
  select(Species,isMarine:match_type)

worms_rec <- worms_rec %>% 
  mutate(worms_marfresh=ifelse(isMarine==1&isFreshwater==0&isBrackish==0,"Marine",NA)) %>% 
  mutate(worms_marfresh=ifelse(isMarine==1&isFreshwater==1|isBrackish==1,"Euryhaline",worms_marfresh)) %>% 
  mutate(worms_marfresh=ifelse(isMarine==0&isFreshwater==1&isBrackish==0,"Fresh",worms_marfresh))
                        
saveRDS(worms_rec,"data/WORMS_records.RDS")
}

worms_rec <- readRDS("data/WORMS_records.RDS")


#is record marine
mht2[,marine:=FALSE]
mht2[,marine:=ifelse(!is.na(TYPE)|!is.na(sal)&sal>30,TRUE,marine),by=.(rec)]
mht2[,marine:=ifelse(onland==T,F,marine),by=.(rec)]

#is record in an estuary
mht2[,estuary:=FALSE]
mht2[,estuary:=ifelse(!is.na(LABEL),TRUE,estuary)] #in estuary polygon


#is record likely from a market
mht2[,market:=ifelse(dtc==0,TRUE,FALSE),by=.(rec)]

#is record freshwater
mht2[,freshwater:=FALSE]
mht2[,freshwater:=ifelse(!is.na(sal),FALSE,freshwater)] #if salinity is not na, then false
mht2[,freshwater:=ifelse(!is.na(TYPE),FALSE,freshwater)] #if in a meow, then false
mht2[,freshwater:=ifelse(onland==T&dtco>=set_units(5e3,m),TRUE,freshwater)] #if beyond 5km of coast, then true
mht2[,freshwater:=ifelse(!grepl("coastal|No Data|oceanic islands",MHT_TXT)&!is.na(MHT_TXT),TRUE,freshwater)] #if in feow and not a coastal or oceanic record
mht2[,freshwater:=ifelse(onland==T&is.na(LABEL),TRUE,freshwater)] #if on land and not in and estuary polygon then true

mht3<- mht2 %>% 
  filter(dupl==F) %>% 
  filter(market==F) %>% 
  dplyr::select(Species,mar.fresh2,onland,estuary,marine,freshwater,rec,dtco) %>%
  group_by(Species) %>% 
  summarize_at(vars(estuary:freshwater),function(x) sum(as.numeric(x),na.rm=T)/n())
 
 
dtco <- NULL
dtco<- mht2 %>% 
  mutate(sss=ifelse(estuary==TRUE,15,sal)) %>% 
  mutate(sss=ifelse(freshwater==TRUE,0,sss)) %>% 
  mutate(sss=ifelse(onland==F&dtco>set_units(3e3,m),35,sss)) %>%
  group_by(Species) %>% 
  mutate(dtco=ifelse(onland==T,dtco*set_units(-1,m),dtco)) %>% 
  dplyr::summarize(m_dtco=(median(dtco)/set_units(1000,m) %>% as.numeric()),
            sal=mean(sss,na.rm = T) %>% log,
            mar.fresh2=first(mar.fresh2)) %>% 
  mutate(sal=ifelse(is.infinite(sal),0,sal)) %>% 
  dplyr::filter(!is.na(sal)) 


mht4 <- mht3 %>%       
  mutate_at(vars(estuary:freshwater),function(x) asin(sqrt(x))) %>% 
  left_join(dtco) %>%
 left_join(fw) %>% 
  mutate(fw_status=ifelse(!is.na(fw_status),1,0)) %>% 
  left_join(worms_rec %>% dplyr::select(Species,isMarine,isFreshwater)) %>%
  na.omit



#--------- ML analysis  ---------# 

#first designation of major halohabitat

rf_dat <-  mht4 %>% dplyr::select(-Species)
rf_dat$mar.fresh2 <- as.factor(rf_dat$mar.fresh2)
rf <-randomForest(mar.fresh2~.,data= rf_dat, ntree=10000,mtry=2) 

#add predicitons
rf_dat$rf_pred<-mht4$rf_pred<- rf$predicted


#GBM analysis in h2O

#initialize
h2o.init(nthreads = 2) 

h2o_data <- mht4%>% dplyr::select(estuary:isFreshwater,-rf_pred)

# convert split train and test data frames to h2o objects
mht.h2o <- h2o.splitFrame(as.h2o(h2o_data),.8)

test <- mht.h2o[[1]]
train <- mht.h2o[[2]]

y <- "mar.fresh2"

# Identify the predictor columns
x <- setdiff(names(train), y)

# Convert response to factor
train$mar.fresh2 <- as.factor(train$mar.fresh2)
test$mar.fresh2 <- as.factor(test$mar.fresh2)

# Run model
model <- h2o.gbm(
  x = x,
  y = y,
  training_frame = train,
  ntrees = 10000,
  learn_rate = 0.1,
  sample_rate = 1.0,
  max_depth = 5,
  col_sample_rate_per_tree = 0.8,
  seed = 1,
  keep_cross_validation_predictions = T
  
)


perf <- h2o.performance(model, train)
print(perf)

# Make predictions and add to data

pred <- h2o.predict(model,newdat=as.h2o(mht4))
mht4$h2o_pred <- pred %>% as.tibble() %>% pull(predict)

# covert to matrices
res <- mht4 %>% 
  dplyr::rename(FB_salinity=mar.fresh2) %>% 
  dplyr::select(Species,FB_salinity,rf_pred,h2o_pred) %>% mutate_at(vars(FB_salinity:h2o_pred),function(x) as.numeric(as.factor(x))) %>% 
  pivot_longer(FB_salinity:h2o_pred)

m.l <- list()
for(i in res$Species){ 
ii <- res %>% dplyr::filter(Species==i) %>% dplyr::select(-Species)

m <- matrix(ncol=length(ii$name),nrow=length(ii$name))
rownames(m) <- ii$name
colnames(m) <- ii$name
for(x in 1:length(ii$name)){
  for(y in 1:length(ii$name)){
    s <- ii[x,2]==ii[y,2] 
    m[y,x] <- s%>% as.numeric()
  }
}
m.l[[paste0(i)]] <- m
}


#save
#saveRDS(mht4,"data/ML_data.RDS")

res_matrix <- Reduce(`+`, m.l)/length(mht4$Species %>% unique())

saveRDS(list(rf=rf,h2o=model,matrices=m.l,res_matrix=res_matrix),"data/ML_results.RDS")


#now designate diadromy

#load data and predictions saved for table in supplementary material Table S1

HT_dat <-readRDS("data/Final_Halotolerance.RDS")%>% select(Species,Final_HT)

# Corush (2019) data
diad<- read_csv("data/corush_states.csv") %>% 
  pivot_longer(Fresh:Diadromous) %>% 
  dplyr::rename(state=name,Species=value) %>% 
  mutate(Species=gsub("_"," ",Species),
  Final_HT=ifelse(is.na(Final_HT),FB_HT,Final_HT)) %>% 
  na.omit

# span of distance from coast
span_dtco <- mht2 %>% 
  group_by(Species) %>% 
  mutate(dtco=ifelse(onland==T,dtco*set_units(-1,m),dtco)) %>% 
  dplyr::summarize(span_dtco=(diff(range(dtco))/set_units(1000,m) %>% as.numeric()))


mht5 <- mht4 %>% 
  left_join(span_dtco) %>%
  left_join(HT_dat) %>% 
  full_join(diad) %>% 
  mutate(Final_HT=ifelse(!is.na(state)&state=="Diadromous"&Final_HT=="Euryhaline",paste0(Final_HT,"-",state),Final_HT)) %>% 
 mutate(dup=duplicated(Species)) %>% 
  filter(dup==F) %>% 
  select(-dup) %>% 
  filter(Species %in% phy$tip.label)

mht5 %>% view


#random forest 

rf_dat2 <-  mht5 %>% dplyr::select(-Species,-state,-mar.fresh2) 
rf_dat2$Final_HT <- as.factor(rf_dat2$Final_HT)
rf2 <-randomForest(Final_HT~.,data= rf_dat2, ntree=10000,mtry=2) 


rf_dat2$rf_pred2<-mht5$rf_pred2<- rf2$predicted

mht5 %>% filter(rf_pred2!=Final_HT) %>% view()


#GBM h2o analysis

#initialize
h2o.init(nthreads = 2) 

h2o_data2 <- mht5%>% dplyr::select(estuary:span_dtco,Final_HT,-rf_pred2,-state,-mar.fresh2)
# convert, split, train and test data frames to h2o objects

#split data
mht.h2o2 <- h2o.splitFrame(as.h2o(h2o_data2),.8)

test2 <- mht.h2o2[[1]]
train2 <- mht.h2o2[[2]]
print(dim(train))

print(dim(test))
y <- "Final_HT"

# Identify the predictor columns
x <- setdiff(names(train2), y)

# Convert response to factor
train2$Final_HT <- as.factor(train2$Final_HT)
test2$Final_HT <- as.factor(test2$Final_HT)

model2 <- h2o.gbm(
  x = x,
  y = y,
  training_frame = train2,
  ntrees = 10000,
  learn_rate = 0.1,
  sample_rate = 1.0,
  max_depth = 5,
  col_sample_rate_per_tree = 0.8,
  seed = 1,
  keep_cross_validation_predictions = T
  
)


perf2 <- h2o.performance(model2, train2)
print(perf2)


pred2 <- h2o.predict(model2,newdat=as.h2o(mht5))
mht5$h2o_pred <- pred2 %>% as.tibble() %>% pull(predict)

mht5 %>% filter(Final_HT!=h2o_pred) %>% view

res2 <- mht5 %>% 
  dplyr::select(Species,Final_HT,rf_pred2,h2o_pred) %>% 
  mutate_at(vars(Final_HT:h2o_pred),function(x) as.numeric(as.factor(x))) %>% 
  pivot_longer(Final_HT:h2o_pred)

m.l2 <- list()
for(i in res2$Species){ 
  ii <- res2 %>% dplyr::filter(Species==i) %>% dplyr::select(-Species)
  
  m <- matrix(ncol=length(ii$name),nrow=length(ii$name))
  rownames(m) <- ii$name
  colnames(m) <- ii$name
  for(x in 1:length(ii$name)){
    for(y in 1:length(ii$name)){
      s <- ii[x,2]==ii[y,2] 
      m[y,x] <- s%>% as.numeric()
    }
  }
  m.l2[[paste0(i)]] <- m
}


#save
res_matrix2 <- Reduce(`+`, m.l2)/length(mht5$Species %>% unique())

saveRDS(list(rf=rf2,h2o=model2,matrices=m.l,res_matrix=res_matrix2),"data/Diad_ML_results.RDS")


#---------ML results reported in paper ---------# 

ML_HT_res <- readRDS("data/ML_results.RDS")$res_matrix
Diad_HT_res <- readRDS("data/Diad_ML_results.RDS")$res_matrix

Final_class <- readRDS("data/Diad_Final_Halotolerance.RDS")

#number of species in each state

Final_class %>% 
  group_by(Final_HTb) %>% 
  count %>% 
  ungroup %>% 
  mutate(per=n/sum(n))

#number of times we had to depend on literature.
Final_class %>% 
  group_by(FB_HT!=rf_pred|FB_HT!=h2o_pred,FB_HT) %>% 
  count %>% 
  arrange(FB_HT)

Final_class %>% 
  mutate(HT_agree=Final_HT==FB_HT) %>% 
  group_by(HT_agree,FB_HT,Final_HT) %>%
  count %>% 
  group_by(FB_HT) %>% 
  mutate(per=n/sum(n)) %>% 
  arrange(FB_HT) %>% 
  filter(HT_agree==F) %>% 
  mutate(FB_wrong=sum(per))



#which Orders are the culprits in FB
ords <- fishtree_taxonomy() %>% filter(rank=="order") %>% 
 filter(!grepl("Incertae sedis",name))
  
ords_sp <- lapply(ords$name, function(x) tibble(species=fishtree_taxonomy(rank=x)[[1]]$species,order=x)) %>% do.call(rbind,.)



Final_class %>% 
  mutate(HT_agree=Final_HT==FB_HT) %>% 
  left_join(ords_sp %>% rename(Species=species)) %>% 
  group_by(order,FB_HT,Final_HT,HT_agree) %>% 
  count %>% 
  ungroup() %>% 
  left_join(
    Final_class %>% 
      left_join(ords_sp %>% rename(Species=species)) %>% 
      group_by(order,FB_HT) %>% 
      summarize(n_species=n())
  ) %>% 
  filter(HT_agree==F) %>% 
  mutate(FB_wrong=n/n_species) %>% 
  arrange(-n) 

Final_class %>% 
  mutate(HT_agree=Final_HT==FB_HT) %>% 
  left_join(ords_sp %>% rename(Species=species)) %>% 
  group_by(order,FB_HT,Final_HT,HT_agree,verification) %>% 
  count %>% 
  filter(order=="Cypriniformes")
  

Final_class %>% 
  mutate(HT_agree=Final_HT==FB_HT) %>% 
  left_join(ords_sp %>% rename(Species=species)) %>% 
  group_by(order,FB_HT,Final_HT,HT_agree,verification) %>% 
  filter(order=="Cypriniformes") %>% 
  view


#1.6% Originally marine, but not

#Diad

Final_class%>% 
  filter(Final_HTb!=h2o_pred_diad|Final_HTb!=rf_pred_diad) %>% filter(rf_pred_diad=="Euryhaline-Diadromous"|Final_HTb=="Euryhaline-Diadromous"|h2o_pred_diad=="Euryhaline-Diadromous") 

#number of times we had to depend on literature.
Final_class %>% 
  group_by(Final_HTb==rf_pred_diad&Final_HTb==h2o_pred_diad,Orig_Diad) %>% 
  count %>% 
  filter(Orig_Diad=="Diadromous") %>%
  ungroup() %>% 
  mutate(per=n/sum(n))



Final_class %>% 
  mutate(Diad_agree=Final_Diad==Orig_Diad) %>% 
  group_by(Diad_agree) %>%
  filter(Orig_Diad=="Diadromous") %>% 
  count %>% 
  ungroup() %>% 
  mutate(per=n/sum(n))



#12.8% Originally diad, but not

Final_class %>% 
  mutate(Diad_agree=Final_Diad==Orig_Diad) %>% 
  group_by(Diad_agree) %>%
  filter(Orig_Diad=="Nondiadromous") %>% 
  count %>% 
  ungroup() %>% 
  mutate(per=n/sum(n))

#0.1% Originally nondiad, but diad

Final_class %>% 
filter(Final_HTb=="Euryhaline-Diadromous"|Orig_Diad=="Diadromous",Final_Diad!=Orig_Diad) %>% 
  view



  
