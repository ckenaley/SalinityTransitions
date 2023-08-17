library(ape)
library(tidyverse)
library(fishtree)
library(phytools)
library(parallel)
library(foreach)
library(doParallel)
library(corHMM)
library(ggbeeswarm)
library(ggtree)



#----- Data -----#
dat <- readRDS("data/Diad_Final_Halotolerance.RDS") %>% select(Species,Final_HTb)

HT <- dat$Final_HTb
names(HT) <- dat$Species

phy <- fishtree_phylogeny()

phy$tip.label <- gsub("_"," ",phy$tip.label)
phy <- keep.tip(phy,dat$Species)


states <- sort(unique(HT))

#----- Fit Models with freshwater at root -----#

fitD_ARD <- fitDiscrete(phy,HT,"ARD",ncores=6,pi=c(0,0,1,0))
fitD_SYM <- fitDiscrete(phy,HT,"SYM",ncores=6,pi=c(0,0,1,0))
fitD_ER <- fitDiscrete(phy,HT,"ER",ncores=6,pi=c(0,0,1,0))

lapply(list(fitD_ARD,fitD_SYM,fitD_ER),AIC)

#very similar results with simple ape::ace (no limitation in that prior wasn't/can't be set)

fit_ace <- ace(x = HT,phy = phy,type="discrete",model="ARD")
Q_ace <- matrix(fit_ace$rates[fit_ace$index.matrix],ncol=4,nrow=4,dimnames = list(sort(unique(HT)),sort(unique(HT))))


## ordered single rate, no M <-> F
model<-matrix(c(0,1,1,1,
                1,0,1,0,
                1,1,0,0,
                1,1,0,0
                ),4,4,byrow=TRUE,
              dimnames=list(states,states))
fitD_SYM_ord<-fitDiscrete(phy,HT,model=model,ncores=6,pi=c(0,0,1,0))

## ordered all rates diff, no M <-> F, no ED->M

model<-matrix(c(0,1,2,3,
                4,0,5,0,
                6,7,0,0,
                8,9,0,0),4,4,byrow=TRUE,
              dimnames=list(states,states))
fitD_ARDor<-fitDiscrete(phy,HT,model=model,ncores=6,pi=c(0,0,1,0))


lapply(list(fitD_ARD,fitD_SYM,fitD_ER,fitD_ARDor), AIC)

Q_ARD <- geiger::as.Qmatrix.gfit(fitD_ARD)
Q_SYM <- geiger::as.Qmatrix.gfit(fitD_SYM)
Q_ER <- geiger::as.Qmatrix.gfit(fitD_ER)

#saveRDS(list(fitD_ARD,fitD_SYM,fitD_ER,fitD_ARDor,fitOrdered1),"data/MK_fits.RDS")

fits <- readRDS("data/MK_fits.RDS")
#following https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12423

#set root prior to fresh
pi<-setNames(c(0,0,1,0),sort(unique(HT)))



#----- Run 1000 simmaps with ARD model in parallel-----#

start_time <- Sys.time()

cl <- parallel::makeCluster(spec=detectCores()-1)
doParallel::registerDoParallel(cl)
simmap_l <- foreach(i=1:100,.packages = c("phytools")) %dopar% {
  phytools::make.simmap(phy,HT,nsim=10,model = "ARD", pi=pi,Q=Q_ARD)
}

parallel::stopCluster(cl)

end_time <- Sys.time()

simmap_trees <- simmap_l%>% do.call(c,.)
res <- list(times=c(start_time,end_time),res=simmap_trees,dat=dat,phy=phy)

dat_name <- paste0("data/simmaps_ARD_",Sys.Date(),".RDS")

saveRDS(res,dat_name)


#SYM, JIC one would like it.
start_time <- Sys.time()
#https://privefl.github.io/blog/a-guide-to-parallelism-in-r/

cl <- parallel::makeCluster(spec=detectCores()-1)
doParallel::registerDoParallel(cl)
simmap_l <- foreach(i=1:100,.packages = c("phytools")) %dopar% {
  phytools::make.simmap(phy,HT,nsim=10,model = "SYM", pi=pi,Q=Q_SYM)
}

parallel::stopCluster(cl)

end_time <- Sys.time()

simmap_trees <- simmap_l%>% do.call(c,.)
res <- list(times=c(start_time,end_time),res=simmap_trees,dat=dat$d,phy=dat$phy)

dat_name <- paste0("data/simmaps_SYM_",Sys.Date(),".RDS")

saveRDS(res,dat_name)



