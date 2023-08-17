#----- Load packages -----#
library(ape)  
library(phytools)  
library(diversitree)
library(scales)  
library(HDInterval)  
library(coda)  
library(tidyverse)

dat <- readRDS("data/Diad_Final_Halotolerance.RDS") %>% select(Species,Final_HTb)


HT_sse <- dat$Final_HTb
HT_sse <- recode(HT_sse,Euryhaline=1,`Euryhaline-Diadromous` =2,Fresh=3, Marine=4)
names(HT_sse) <- dat$Species
phy <- fishtree_phylogeny()

phy$tip.label <- gsub("_"," ",phy$tip.label)
phy <- keep.tip(phy,dat$Species)

phy_ul <- force.ultrametric(phy)

states <- sort(unique(dat$Final_HTb))


samp.f <- c(length(dat$Final_HTb[dat$Final_HTb=="Euryhaline"])/4000,
length(dat$Final_HTb[dat$Final_HTb=="Euryhaline-Diadromous"])/223,
length(dat$Final_HTb[dat$Final_HTb=="Marine"])/16764,
length(dat$Final_HTb[dat$Final_HTb=="Fresh"])/15170)


HT.musse=make.musse(phy_ul,HT_sse,length(unique(HT_sse)), sampling.f=samp.f)

HT.start=starting.point.musse(phy_ul,length(unique(HT_sse)))
  
#----- Fitting the model -----#
HT.fituncons=find.mle(HT.musse,HT.start,control=list(maxit=100000))
  
#constrain so ED->M=0, F<->M=0 , ext and speciation vary with each state
HT.cons1 <- constrain(HT.musse, q24 ~ 0, 
                             q34 ~ 0, q43 ~ 0)
HT.cons1.fit <- find.mle(HT.cons1,  HT.start[argnames(HT.cons1)],control=list(maxit=100000))
  
  #constrain so that M <->F must, but go through only Eury.
  
HT.cons2 <- constrain(HT.musse, q23 ~ 0,q24 ~ 0, 
                        q32 ~ 0, q34 ~ 0,
                        q42~0,q43~0)
HT.cons2.fit <- find.mle(HT.cons2,  HT.start[argnames(HT.cons2)],control=list(maxit=100000))
  
HT.cons2 <- constrain(HT.musse, q23 ~ 0,q24 ~ 0, 
                        q32 ~ 0, q34 ~ 0,
                        q42~0,q43~0)
  
  
saveRDS(list(HT.fituncons,HT.cons1.fit,HT.cons2.fit),"data/Salinity_MuSSE_models.RDS")

lapply(list(HT.fituncons,HT.cons1.fit,HT.cons2.fit),AIC)
  #unconstrained model fits best
  
  #----- MCMC -----#
  
prior <- make.prior.exponential(1/2)
prelim <- diversitree::mcmc(HT.musse, coef(HT.fituncons, full=TRUE), nsteps=200, prior=prior, w=1, print.every=10)
  
w <- diff(sapply(prelim[2:(ncol(prelim)-1)], quantile, c(0.05, 0.95)))
samples <- diversitree::mcmc(HT.musse, coef(HT.fituncons, full=TRUE), nsteps=2000, prior = prior, w=w, print.every=50)
  
  
  
  #----- Load sampels from HPC run -----#
  
  samples <- readRDS("data/Salinity_MuSSE_res.RDS")$samples

  samples <- samples[201:2000,]
  
  sam <- na.omit(samples[,c(2:9,22)])
  
  samp_final <- sam
  


#----- Extracting rates -----#

#Diversification rates
dr_Eury <- samp_final$lambda1-samp_final$mu1
dr_EuryDiad <- samp_final$lambda2-samp_final$mu2
dr_Fresh <- samp_final$lambda3-samp_final$mu3
dr_Marine <- samp_final$lambda4-samp_final$mu4

dr_all <- cbind(dr_Eury,dr_EuryDiad,dr_Fresh,dr_Marine)

#Speciation rates
sr_Eury <- samp_final$lambda1
sr_EuryDiad <- samp_final$lambda2
sr_Fresh <- samp_final$lambda3
sr_Marine <- samp_final$lambda4


sr_all <- cbind(sr_Eury,sr_EuryDiad,sr_Fresh,sr_Marine)

#Extinction rates
xr_Eury <- samp_final$mu1
xr_EuryDiad <- samp_final$mu2
xr_Fresh <- samp_final$mu3
xr_Marine <- samp_final$mu4

xr_all <- cbind(xr_Eury,xr_EuryDiad,xr_Fresh,xr_Marine)


#----- Plotting results -----#
mode <- function(s) {
  d <- density(s)
  d$x[which.max(d$y)]
}

xr_all_tab <- xr_all %>% data.frame %>% mutate(param="Ext. rate")
sr_all_tab <- sr_all %>% data.frame %>% mutate(param="Spec. rate")
dr_all_tab <- dr_all %>% data.frame %>% mutate(param="Net. div. rate")
colnames(dr_all_tab) <- colnames(sr_all_tab) <- colnames(xr_all_tab) <- c(sort(unique(dat$Final_HTb)),"param")


dr_all_tab <-  rbind(xr_all_tab,sr_all_tab,dr_all_tab)%>% pivot_longer(Euryhaline:Marine) 


ul <-dr_all_tab %>% 
  group_by(name,param) %>% 
  summarise(lower=hdi(value)["lower"],
            upper=hdi(value)["upper"])

hdi_tab <- dr_all_tab %>% 
  group_by(name,param) %>% 
  summarize(dt_x=density(value)$x,
            dt_y=density(value)$y)
 

hdi_tab2<- hdi_tab %>% 
  left_join(ul) %>% 
  summarize(x=c(min(lower), dt_x[dt_x>=min(lower) & dt_x<=max(upper)],max(upper)),
            y=c(0, dt_y[dt_x>=min(lower) & dt_x<=max(upper)],0))

hdi_tab3<- dr_all_tab %>%
  group_by(name,param) %>% 
  left_join(ul) %>% 
  summarize(min=min(lower),max=max(upper),mode=mode(value)) %>% 
  mutate(y=seq(-5,-20,-5)[order(mode)])

library(viridis)
library(cowplot)
cols <- setNames(pal_aaas("default")(4),states)




p_xr <- hdi_tab2 %>% filter(param=="Ext. rate")%>% ggplot(aes(x,y,fill=name))+geom_polygon(alpha=0.8)+geom_line(data=hdi_tab  %>% filter(param=="Ext. rate"),aes(dt_x,dt_y),linewidth=0.2)+geom_segment(data=hdi_tab3%>% filter(param=="Ext. rate"),aes(x=min,xend=max,y=y,yend=y,col=name),linewidth=1)+theme_classic(10)+geom_point(data=hdi_tab3%>% filter(param=="Ext. rate"),aes(x=mode,y=y,col=name),size=1.25)+theme(axis.line.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),legend.position = "none")+ylab("")+xlab("Extinction Rate")+scale_x_continuous(breaks=seq(0,.2,.075))+xlim(c(-0.01,0.25))+scale_fill_manual(values=cols)+scale_colour_manual(values=cols)+labs(colour = "State",fill="State") 

#print(p_xr)


p_sr <- hdi_tab2 %>% filter(param=="Spec. rate")%>% ggplot(aes(x,y,fill=name))+geom_polygon(alpha=0.8)+geom_line(data=hdi_tab  %>% filter(param=="Spec. rate"),aes(dt_x,dt_y),linewidth=0.2)+geom_segment(data=hdi_tab3%>% filter(param=="Spec. rate"),aes(x=min,xend=max,y=y,yend=y,col=name),linewidth=1)+theme_classic(10)+geom_point(data=hdi_tab3%>% filter(param=="Spec. rate"),aes(x=mode,y=y,col=name),size=1.25)+theme(axis.line.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),legend.position = "none")+ylab("")+xlab("Speciation Rate")+scale_x_continuous(breaks=seq(0,.2,.075))+xlim(c(0.05,0.41))+scale_fill_manual(values=cols)+scale_colour_manual(values=cols)+labs(colour = "State",fill="State") 

#print(p_sr)

p_div <- hdi_tab2 %>% filter(param=="Net. div. rate") %>% ggplot(aes(x,y,fill=name))+geom_polygon(alpha=0.8)+geom_line(data=hdi_tab%>% filter(param=="Net. div. rate"),aes(dt_x,dt_y),linewidth=0.2)+geom_segment(data=hdi_tab3%>% filter(param=="Net. div. rate"),aes(x=min,xend=max,y=y,yend=y,col=name),linewidth=1)+theme_classic(10)+geom_point(data=hdi_tab3%>% filter(param=="Net. div. rate"),aes(x=mode,y=y,col=name),size=1.25)+theme(axis.line.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),legend.position = c(0.6,0.6),legend.text = element_text(size=10),legend.key.size = unit(0.5, 'cm'))+ylab("")+xlab("Net Diversification Rate")+scale_x_continuous(breaks=seq(0,.2,.075))+xlim(c(0,0.35))+scale_fill_manual(values=cols)+scale_colour_manual(values=cols)+labs(colour = "State",fill="State") 
#print(p_div)


lower_p <- plot_grid(p_sr,p_xr,nrow=1,labels=c("B","C"))

pdf("figures/Div_rates.pdf",h=7,w=5.5)
p_all <- plot_grid(p_div,lower_p,nrow=2,rel_heights = c(1.3,1),labels=c("A",NA))
print(p_all)
dev.off()




