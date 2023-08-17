library(ggbeeswarm)
library(ggtree)
library(phytext2)
library(ape)
library(ggsci)
library(phytools) 
library(circlize) 
library(scales) 
library(diversitree)  
library(fishtree)
source("helper_functions.R")

#----- Data -----#

dat <- readRDS("data/Diad_Final_Halotolerance.RDS") %>% select(Species,Final_HTb)

phy <- fishtree_phylogeny()
phy$tip.label <- gsub("_"," ",phy$tip.label)
phy <- keep.tip(phy,dat$Species)

#load simmaps from best-fit ARD model
res_ARD <- readRDS("data/simmaps_ARD_2023-06-09.RDS")


#get times for each transitions using custom function
trans_times <- get_trans_times(res_ARD$res)



#summary of changes
class(res_ARD$res) <- c("multiSimmap","multiPhylo")

#trans matrix
trans_ARD <- lapply(res_ARD$res,function(x) describe.simmap(x)$Tr)

#node states, values with highest PP
trans_nodes <- summary(res_ARD$res)

node_states <- trans_nodes$ace %>% data.frame %>% 
  mutate(node=rownames(trans_nodes$ace) %>% as.numeric) %>% 
  na.omit() %>% 
  pivot_longer(Euryhaline:Marine) %>% 
  group_by(node) %>% 
  summarize(state=name[which.max(value)])


#----- make gg tree and clade labels -----#
phy_gg <- ggtree(phy,ladderize = F,col="#7a7777",alpha=0.2)
phy_gg_fan <- ggtree(phy,ladderize = F,layout = "fan",open.angle = 180,col="#7a7777",alpha=0.2)


near_names <- c("Holostei",
                "Osteoglossomorpha",
                "Elopocephalai",#=Elopomorpha
                "Neoteleostei",#=Euteleosteo+Lepidogalaxias
                "Otocephala",
                "Otomorpha",
                "Ostariophysi",
                "Otophysi",
                "Euteleosteomorpha",
                "Acanthomorphata",
                "Percomorphaceae")

ranks <- fishtree_taxonomy() %>% filter(name%in%near_names) %>% select(rank) %>% unique() %>% mutate(order=1:n())

#nodes and names for labels (requires some filtering from fishtree data)
tx <- fishtree_taxonomy() %>% filter(name%in%near_names) %>% 
  filter(!grepl("Incert",name)) %>% 
  pull(name) %>% 
  lapply(., function(x) tibble(species=fishtree_taxonomy(rank=x)[[1]]$species,name=x,rank=fishtree_taxonomy() %>% filter(name==x) %>% pull(rank))) %>% 
  do.call(rbind,.) %>% 
  filter(species%in%res_ARD$res[[1]]$tip.label) %>% 
  group_by(name,rank) %>% 
  filter(!name%in%c("Amiiformes","Lepidogalaxii","Lepidogalaxiiformes","Pholidichthyiformes","Stylephoriformes")) %>% 
  summarize(node=getMRCA(res_ARD$res[[1]],species),n=mean(grep(species,res_ARD$res[[1]]$tip.label)))

clade_aes <- tibble(rank=unique(tx$rank)) %>% 
  left_join(ranks)



tx <- tx %>% left_join(clade_aes) %>% arrange(order) 

#make blocks of clades
tx$block <-c(1,2,3,5,4,5,4,5,5)


tx <- tx %>% 
  group_by(block) %>% 
  mutate(col=rep(c("gray55","gray70","gray85","gray100"),20)[1:n()]) %>% 
  group_by(block) %>% 
  arrange(order,.by_group = T) %>% 
  mutate(order=1:n(),
         size=seq(40,10,length.out=n()),
         offset=seq(120000,5000,length.out=n()))

clade_labs <- lapply(1:nrow(tx),function(x) geom_cladelabel(node=tx[x,]$node,label=tx[x,]$name,angle=0,color=c(tx[x,]$col,"black"),extend = F,barsize =tx[x,]$size,offset = tx[x,]$offset,align = T,fontsize=6))



td <- dat %>% 
  rename(state=Final_HTb,label=Species)


#fix "."
node_states <- node_states%>% mutate(state=gsub("\\.","-",state))



anc_map <- phy2 +geom_point(aes(col=state.x),alpha=0.5,size=1.5)+geom_tippoint(aes(col=state.y),alpha=0.5,size=0.5)+scale_color_aaas(na.translate = F)+theme(legend.position=c(0.55,0.6),plot.margin=grid::unit(c(-5,0,-30,0), "cm"),legend.text = element_text(size=30),title = element_blank())+guides(colour = guide_legend(override.aes = list(alpha = 1,size = 10)))


anc_map+clade_labs

#save a rough draft of Figure for editing
ggsave("figures/anc_map_ARD.pdf", width = 50, height = 30, units = "cm", limitsize = FALSE)



#----- Figures 2-4, S1-3 -----#

#examples
exs <- c(17713,18989,16159,12960,11594,9759)
exs <- setNames(exs,c("Catfishes","Clupeiformes","Stomiatii","cottoids","Carangiomorpha","Ovalentaria"))



# Figures 4, S1 Ovalentaria+(carngimorpharia+anabantaria)

#look a tree to see where important clades are
oval_phy <- ggtree(clade_simms[[1]],ladderize = F,col="#7a7777",alpha=0.2)


#marine clades to collpases
coll_n <- c(4204,4163,3684,2821)
coll_n <- setNames(coll_n,c("Carangimorpha","Mugiliformes+Ambassidae","Blenniiformes+Pomacentridae+Others","Belonoidei"))

#fw clades to highlight
fw_oval <- c(4603,3053,2900,2374,2913)
fw_oval <- setNames(fw_oval,c("Anabantaria","Cichliformes","Andrianichthyidae","Cyprinodontiformes","Atheriniformes"))

oval_simms <- lapply(res_ARD$res,function(x) extract.clade.simmap(x,exs[6]))

class(oval_simms) <- c("multiSimmap","multiPhylo")

oval_trans_times <- get_trans_times(oval_simms)

oval_trans <- lapply(oval_simms,get_trans_loc)

oval_trans <- lapply(1:length(oval_trans),function(x) oval_trans[[x]] %>% mutate(simm=x)) %>% do.call(rbind,.)


#fresh chord diagrams

oval_trans_1 <- simmap_trans(lapply(oval_simms,function(x) extract.clade.simmap(x,fw_oval[1])))

oval_trans_2 <- simmap_trans(lapply(oval_simms,function(x) extract.clade.simmap(x,fw_oval[2])))

oval_trans_3 <- simmap_trans(lapply(oval_simms,function(x) extract.clade.simmap(x,fw_oval[3])))

oval_trans_4 <- simmap_trans(lapply(oval_simms,function(x) extract.clade.simmap(x,fw_oval[4])))

oval_trans_5 <- simmap_trans(lapply(oval_simms,function(x) extract.clade.simmap(x,fw_oval[5])))

pdf("figures/oval_mar_chord.pdf")
par(mfrow = c(4,1),cex = 0.75)
quick_chord(oval_mar_trans_1$trans,10)
title(names(coll_n)[1])
quick_chord(oval_mar_trans_2$trans,2)
title(names(coll_n)[2])
quick_chord(oval_mar_trans_3$trans,5)
title(names(coll_n)[3])
quick_chord(oval_mar_trans_4$trans,5)
title(names(coll_n)[4])
dev.off()

#marine chord diagrams

oval_mar_trans_1 <- simmap_trans(lapply(oval_simms,function(x) extract.clade.simmap(x,coll_n[1])))

oval_mar_trans_2 <- simmap_trans(lapply(oval_simms,function(x) extract.clade.simmap(x,coll_n[2])))

oval_mar_trans_3 <- simmap_trans(lapply(oval_simms,function(x) extract.clade.simmap(x,coll_n[3])))

oval_mar_trans_4 <- simmap_trans(lapply(oval_simms,function(x) extract.clade.simmap(x,coll_n[4])))


pdf("figures/oval_mar_chord.pdf")
par(mfrow = c(5,1),cex = 0.75)
quick_chord(oval_trans_1$trans,0.5)
title(names(fw_oval)[1])
quick_chord(oval_trans_2$trans,2)
title(names(fw_oval)[2])
quick_chord(oval_trans_3$trans,1)
title(names(fw_oval)[3])
quick_chord(oval_trans_4$trans,5)
title(names(fw_oval)[4])
quick_chord(oval_trans_5$trans,5)
title(names(fw_oval)[5])
dev.off()

#fresh highlight
oval_clade_fresh <- gg_simmap_clade(simms=res_ARD$res,clade=exs[6],coll=coll_n,highlight = fw_oval,trans_hl = "->Fresh")

#marine higlight
oval_clade_marine <- gg_simmap_clade(simms=res_ARD$res,clade=exs[6],coll=fw_oval,highlight = coll_n,trans_hl = "->Marine")


oval_clade_marine$plot<- oval_clade_marine$plot %<+% td 


oval_clade_fresh$plot$data$state.y <- gsub("\\.","-",oval_clade_fresh$plot$data$state.y )
oval_clade_fresh$plot$data$state.x <- gsub("\\.","-",oval_clade_fresh$plot$data$state.x )

other_nodes <-lapply(c(fw_oval,coll_n),function(x) nodepath(oval_simms[[1]],2368,x)) %>% unlist() %>% unique()
other_nodes2 <- lapply(c(fw_oval,coll_n),function(x) getDescendants(oval_simms[[1]],x)[1:2]) %>% unlist() %>% unique()
other_nodes3 <- lapply(c(other_nodes2),function(x) getDescendants(oval_simms[[1]],x)[1:2]) %>% unlist() %>% unique()
other_nodes4 <- lapply(c(other_nodes3),function(x) getDescendants(oval_simms[[1]],x)[1:2]) %>% unlist() %>% unique()

#function to move up tree and get node states and times
up_tree <- function(node,tree,h=0.3){
  d <- getDescendants(tree,node)
  nh <- lapply(d,function(x) nodeheight(tree,x)) %>% unlist
  nh2 <- (nh-min(nh))/diff(range(nh))
  d2 <- d[which(nh2<=h)]
  return(d2)
}

mar_nodes <-up_tree(2368,oval_simms[[1]],0.4)


#node data
nd_mar <- oval_clade_marine$plot$data %>% filter(node%in% c(mar_nodes)) 
nd_fw <- oval_clade_fresh$plot$data %>% filter(node%in% c(other_nodes,other_nodes2,other_nodes3,other_nodes4,fw_oval,coll_n)) 

nd_eury <- oval_clade_fresh$plot$data %>% filter(clade=="Cyprinodontiformes" & state.x=="Euryhaline") 



oval_clade_marine$plot+scale_color_aaas(na.translate = F)+theme(legend.position = c(0.1,0.3),legend.text = element_text(size=6),legend.title = element_text(size=8))+geom_tippoint(aes(col=state.y),size=1,shape=15)+scale_shape(na.translate = F)+scale_size(range=c(0.5,3),guide = F)+geom_point(data=nd_mar,aes(x=x,y=y,col=state.x),size=2,shape=15,inherit.aes = F)+scale_x_continuous(breaks = t_breaks,label = t_labs)+theme_tree2()+geom_point(data=nd_mar,aes(x=x,y=y,col=state.x),size=2,shape=15,inherit.aes = F)

ggsave("figures/oval_marine_clade.pdf")

oval_clade_fresh$plot<- oval_clade_fresh$plot %<+% td 


t_breaks <- max(oval_clade_fresh$plot$data$x,na.rm = T)-seq(0,round_any(max(oval_clade_fresh$plot$data$x,na.rm = T),20),20)
t_labs <- seq(0,round_any(max(oval_clade_fresh$plot$data$x,na.rm = T),20),20)

oval_clade_fresh$plot+scale_color_aaas(na.translate = F)+theme_tree(plot.margin=margin(15, 15, 15, 15))+theme(legend.position = c(0.1,0.7),legend.text = element_text(size=6),legend.title = element_text(size=8))+geom_tippoint(aes(col=state.y),size=1,shape=15)+scale_shape(na.translate = F)+scale_size(range=c(0.5,3),guide = F)+geom_point(data=nd_fw,aes(x=x,y=y,col=state.x),size=2,shape=15,inherit.aes = F)+geom_point(data=nd_eury,aes(x=x,y=y,col=state.x),size=2,shape=15,inherit.aes = F)+scale_x_continuous(breaks = t_breaks,label = t_labs)+theme_tree2()

ggsave("figure/oval_fw_clade.pdf")




#violin plots

bee_plots <- lapply(1:4,function(x) bee_simms(oval_simms,fw_oval[x],n=2))

bee_plots[[5]] <-  bee_simms(oval_simms,fw_oval[5],n=6)
bee_data <- lapply(1:5,function(x) bee_plots[[x]]["trans"][[1]]) %>% do.call(rbind,.)


bee_data$clade <- factor(bee_data$clade, levels = unique(bee_data$clade)[c(1,2,5,3,4)])
bee_data$trans <- gsub("(.*)EuryhalineDiadromous(.*)","\\1Eury.-Diad\\2",bee_data$trans)
bee_data$trans <- gsub("(.*)Euryhaline(.*)","\\1Eury.\\2",bee_data$trans)
bee_data %>%
  ggplot( aes(x = trans, y =t))  +
  geom_quasirandom(alpha=0.2,dodge.width=1,size=0.5)+theme_classic(20)+xlab("Transition")+ylab("MYA")+theme(axis.text.x = element_text(angle=75,vjust = 1,hjust = 1,size=15))+ylim(c(0,max(clade_trans_tt$t)))+facet_wrap(.~clade,ncol=1)

ggsave("figures/bee_plots.pdf",height=11,width=4,limitsize = F)



#mar violin plots

mar_bee_plots <- lapply(1:4,function(x) bee_simms(oval_simms,coll_n[x],n=5))
mar_bee_data <- lapply(1:4,function(x) mar_bee_plots[[x]]["trans"][[1]]) %>% do.call(rbind,.)
mar_bee_data$clade <- as.factor(mar_bee_data$clade)
mar_bee_data$clade <- factor(mar_bee_data$clade, levels = levels(mar_bee_data$clade)[c(3,4,2,1)])
mar_bee_data$trans <- gsub("(.*)EuryhalineDiadromous(.*)","\\1Eury.-Diad\\2",mar_bee_data$trans)
mar_bee_data$trans <- gsub("(.*)Euryhaline(.*)","\\1Eury.\\2",mar_bee_data$trans)
mar_bee_data %>%
  ggplot( aes(x = trans, y =t))  +
  geom_quasirandom(alpha=0.2,dodge.width=0.5,size=0.5)+theme_classic(20)+xlab("Transition")+ylab("MYA")+theme(axis.text.x = element_text(angle=75,vjust = 1,hjust = 1,size=15))+ylim(c(0,max(clade_trans_tt$t)))+facet_wrap(.~clade,ncol=1)

ggsave("figures/oval_mar_bee_plots.pdf",height=11,width=4,limitsize = F)



#Figure S2 Catfish 

#estract clade
cat_phy <- extract.clade.simmap(res_ARD$res[[1]],exs[1])

#extract simmaps
cat_simms <- lapply(res_ARD$res,function(x) extract.clade.simmap(x,exs[1]))
class(cat_simms) <- c("multiSimmap","multiPhylo")


cat_phy_gg <- ggtree(cat_phy,ladderize = F,col="#7a7777",alpha=0.2)

cat_phy_gg <- cat_phy_gg%<+% nd %<+% td


#collapse fresh clade
coll_n <- c(219)
coll_n <- setNames(coll_n,c('Ictaluroidea'))

#ID marine clade
cat_mar <- c(130)
cat_mar <- setNames(cat_mar,c("Arioidea"))


cat_trans_1 <- simmap_trans(lapply(cat_simms,function(x) extract.clade.simmap(x,130)))

cat_trans_1$mean_trans

pdf("figures/cat_chord.pdf")
par(cex = 2.5)
quick_chord(cat_trans_1$trans,5)
title("Arioidea")
dev.off()

#custom functions to plot simmmaps
cat_clade <- gg_simmap_clade(simms=res_ARD$res,clade=exs[1],coll=coll_n,highlight = cat_mar,trans_hl = "->Marine")

  
cat_clade$plot$data$state.y <- gsub("\\.","-",cat_clade$plot$data$state.y )
cat_clade$plot$data$state.x <- gsub("\\.","-",cat_clade$plot$data$state.x )

#add some some node values
cat_nodes <-lapply(c(cat_mar,coll_n),function(x) nodepath(cat_phy,129,x)) %>% unlist() %>% unique()
cat_nodes2 <- lapply(c(cat_mar),function(x) getDescendants(cat_phy,x)) %>% unlist() %>% unique()

nd_mar <- cat_clade$plot$data %>% filter(node%in% c(cat_nodes,cat_nodes2,coll_n,cat_mar)) 
nd_fw <- cat_clade$plot$data %>% filter(node%in%  c(cat_nodes,coll_n,cat_mar))... 

cat_clade$plot$data$state.y <- factor(cat_clade$plot$data$state.y,levels=c(cat_clade$plot$data$state.y %>% unique,"Eury.-Diad."))


#for plotting time scale
t_breaks <- max(cat_clade$plot$data$x,na.rm = T)-seq(0,round_any(max(cat_clade$plot$data$x,na.rm = T),10),10)
t_labs <- seq(0,round_any(max(cat_clade$plot$data$x,na.rm = T),10),10)

cat_clade$plot+scale_color_aaas(na.translate = F)+theme_tree(plot.margin=margin(15, 15, 15, 15))+theme(legend.position = c(0.05,0.1),legend.text = element_text(size=8),legend.title = element_text(size=10))+geom_tippoint(aes(col=state.y),size=1,shape=15)+scale_shape(na.translate = F)+scale_size(range=c(0.5,3),guide = F)+geom_point(data=nd_mar,aes(x=x,y=y,col=state.x),size=2,shape=15,inherit.aes = F)+geom_point(data=nd_fw,aes(x=x,y=y,col=state.x),size=2,shape=15,inherit.aes = F)+scale_x_continuous(breaks = t_breaks,label = t_labs)+theme_tree2()

ggsave("figures/cat_clade.pdf")

#cat fish violin plots
cat_bee_plots <- lapply(1,function(x) bee_simms(cat_simms,cat_mar[x],n=5,inc=c(inc[3:5],"Marine->Euryhaline")))
cat_bee_data <- lapply(1,function(x) cat_bee_plots [[x]]["trans"][[1]]) %>% do.call(rbind,.)


cat_bee_data$trans <- gsub("(.*)EuryhalineDiadromous(.*)","\\1Eury.-Diad\\2",cat_bee_data$trans)
cat_bee_data$trans <- gsub("(.*)Euryhaline(.*)","\\1Eury.\\2",cat_bee_data$trans)
cat_bee_data %>%
  ggplot( aes(x = trans, y =t))  +
  geom_quasirandom(alpha=0.2,dodge.width=1,size=0.5)+theme_classic(20)+xlab("Transition")+ylab("MYA")+theme(axis.text.x = element_text(angle=75,vjust = 1,hjust = 1,size=15))+ylim(c(0,40))+facet_wrap(.~clade,ncol=1)

ggsave("figures/cat_bee_plots.pdf",height=5,width=3,limitsize = F)



#Figures S3 clupeiformes


clupe_simms <- lapply(res_ARD$res,function(x) extract.clade.simmap(x,exs[2]))

class(clupe_simms) <- c("multiSimmap","multiPhylo")

clup_phy <- extract.clade.simmap(res_ARD$res[[1]],exs[2])
clup_phy_gg <- ggtree(clup_phy,ladderize = F,col="#7a7777",alpha=0.2)


clup_phy_gg <- clup_phy_gg%<+% nd %<+% td

#mar clades to collapse
coll_n <- c(161,281,218,241,202)

coll_n <- setNames(coll_n,c('Marine Dorosomatidae',"Engraulidae 1","Chirocentridae+ Dussumieriidae + Pristigasteridae","Engraulidae 2","Alosidae"))

#fw clades to highlight
clup_fw <- c(268,196)
clup_fw <- setNames(clup_fw,c("FW Engraulini","FW Dorosomatidae"))



clupe_trans_1 <- simmap_trans(lapply(clupe_simms,function(x) extract.clade.simmap(x,268)))

clupe_trans_2 <- simmap_trans(lapply(clupe_simms,function(x) extract.clade.simmap(x,196)))





pdf("figures/clupe_chord.pdf")
par(mfrow = c(2,1),oma=c(0, 5, 5, 5),cex = 1.25)
quick_chord(clupe_trans_1$trans,1)
title("FW Engraulini")
quick_chord(clupe_trans_2$trans,0.25)
title("FW Dorosomatidae")
dev.off()


clup_clade <- gg_simmap_clade(simms=res_ARD$res,clade=exs[2],coll=coll_n,highlight = clup_fw,trans_hl = "->Fresh")

clup_clade$plot <- clup_clade$plot %<+% td

clup_clade$plot$data$state.y <- gsub("\\.","-",clup_clade$plot$data$state.y )
clup_clade$plot$data$state.x <- gsub("\\.","-",clup_clade$plot$data$state.x )

clupe_nodes <-lapply(c(clup_fw,coll_n,297),function(x) nodepath(clup_phy,151,x)) %>% unlist() %>% unique()

clupe_nodes2 <- lapply(c(clup_fw,coll_n),function(x) getDescendants(clup_phy,x)[1:2]) %>% unlist() %>% unique()

clupe_nodes3 <- lapply(c(other_nodes2),function(x) getDescendants(clup_phy,x)[1:2]) %>% unlist() %>% unique()

clupe_nodes4 <- getDescendants(clup_phy,268)

other_nodes5 <- lapply(c(),function(x) getDescendants(oval_simms[[1]],x)) %>% unlist() %>% unique()

#node data 
nd_mar <- clup_clade$plot$data %>% filter(node%in% c(clupe_nodes,clupe_nodes2,coll_n,fw_oval)) 
nd_fw <- clup_clade$plot$data %>% filter(node%in% c(clupe_nodes,clupe_nodes2,clupe_nodes3,clupe_nodes4,clup_fw,coll_n)) 


t_breaks <- max(clup_clade$plot$data$x,na.rm = T)-seq(0,round_any(max(clup_clade$plot$data$x,na.rm = T),40),40)
t_labs <- seq(0,round_any(max(clup_clade$plot$data$x,na.rm = T),40),40)


clup_clade$plot+scale_color_aaas(na.translate = F)+theme_tree(plot.margin=margin(15, 15, 15, 15))+theme(legend.position = c(0.1,0.3),legend.text = element_text(size=8),legend.title = element_text(size=10))+geom_tippoint(aes(col=state.y),size=1,shape=15)+scale_shape(na.translate = F)+scale_size(range=c(0.5,3),guide = F)+geom_point(data=nd_mar,aes(x=x,y=y,col=state.x),size=2,shape=15,inherit.aes = F)+geom_point(data=nd_fw,aes(x=x,y=y,col=state.x),size=2,shape=15,inherit.aes = F)+scale_x_continuous(breaks = t_breaks,label = t_labs)+theme_tree2()


ggsave("figures/clupe_fw.pdf")


#violin plots

clup_bee_plots <- lapply(1:2,function(x) bee_simms(clupe_simms,clup_fw[x],n=2))
clup_bee_data <- lapply(1:2,function(x) clup_bee_plots [[x]]["trans"][[1]]) %>% do.call(rbind,.)
clup_bee_data$clade <- factor(clup_bee_data$clade, levels = unique(clup_bee_data$clade)[c(1,2)])
clup_bee_data$trans <- gsub("(.*)EuryhalineDiadromous(.*)","\\1Eury.-Diad\\2",clup_bee_data$trans)
clup_bee_data$trans <- gsub("(.*)Euryhaline(.*)","\\1Eury.\\2",clup_bee_data$trans)

clup_bee_data %>%
  ggplot( aes(x = trans, y =t))  +
  geom_quasirandom(alpha=0.2,dodge.width=1,size=0.5)+theme_classic(20)+xlab("Transition")+ylab("MYA")+theme(axis.text.x = element_text(angle=75,vjust = 1,hjust = 1,size=15))+ylim(c(0,40))+facet_wrap(.~clade,ncol=1)

ggsave("figures/clup_bee_plots .pdf",height=8,width=3,limitsize = F)


#----- Chord diagrams (Figure 2)-----#
# summarize transitions on multiSimmap object
simmap_trans <- function (trees) 
{
  if(!"multiSimmap" %in% class(trees)) class(trees) <- "multiSimmap"
  
  trans <- lapply(trees,function(x) describe.simmap(x)$Tr)
  
  mean_trans <- Reduce('+', trans)/length(trans)
  
  states <- sort(unique(getStates(trees[[1]])))
  return(list(trans=trans,mean_trans=mean_trans,states=states))
}


trans_ARD <- simmap_trans(res_ARD$res)

states <- trans_ARD$states

cols <- setNames(pal_aaas("default")(4),states)
group = structure(gsub("\\d", "", states), names = states)


pdf("chord_ARD2.pdf")
if(T){
  par(mfrow = c(2,2),mar=c(0.5, 0.5, 0.5, 0.5),oma=c(0, 10, 5, 10),cex = 0.75)
  
  #Eury
  arr.col = data.frame(rep(states[1],3), states[-1],
                       c("black", "black", "black"))
  
  
  chordDiagram(trans_ARD$mean_trans, group=group,row.col = c(alpha(cols[1], 0.9),alpha(cols[2], 0.1),alpha(cols[3], 0.1),alpha(cols[4], 0.1)), annotationTrackHeight =  c(0.25, .1),
               column.col = cols, grid.col = cols, scale = F,
               directional = 1, direction.type = "arrows",link.arr.col = arr.col,link.arr.length = 0.1,link.arr.type = "big.arrow",h.ratio = 0.2,big.gap = 3,annotationTrack = c("name", "grid"))
  title("Euryhaline")
  
  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.axis(h = "top", labels.cex = 0.7, major.at=seq(0,sum(trans_ARD$mean_trans),100),
                sector.index = sector.name, track.index = 2,minor.ticks=1)
  }, bg.border = NA)
  
  
  #Eury-Diad
  arr.col = data.frame(rep(states[2],2), states[-c(2,4)],
                       rep(cols[2],2))
  
  chordDiagram(trans_ARD$mean_trans, group=group,row.col = c(alpha(cols[1], 0.1),alpha(cols[2], 0.9),alpha(cols[3], 0.1),alpha(cols[4], 0.1)), annotationTrackHeight =  c(0.25, .1), 
               column.col = cols, grid.col = cols, scale = F,
               directional = 1, direction.type = "arrows",link.arr.col = arr.col,link.arr.length = 0.1,link.arr.type = "big.arrow",h.ratio = 0.2,big.gap = 3,annotationTrack = c("name", "grid"))
  
  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.axis(h = "top", labels.cex = 0.7, major.at=seq(0,sum(trans_ARD$mean_trans),100),
                sector.index = sector.name, track.index = 2,minor.ticks=1)
  }, bg.border = NA)
  
  title("Euryhaline-Diadromous")
  
  #Fresh
  arr.col = data.frame(rep(states[3],2), states[-c(3,4)],
                       rep(cols[3],2))
  
  chordDiagram(trans_ARD$mean_trans, group=group,row.col = c(alpha(cols[1], 0.1),alpha(cols[2], 0.1),alpha(cols[3], 0.9),alpha(cols[4], 0.1)), annotationTrackHeight =  c(0.25, .1), 
               column.col = cols, grid.col = cols, scale = F,
               directional = 1, direction.type = "arrows",link.arr.col = arr.col,link.arr.length = 0.1,link.arr.type = "big.arrow",h.ratio = 0.2,big.gap = 3,annotationTrack = c("name", "grid"))
  
  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.axis(h = "top", labels.cex = 0.7, major.at=seq(0,sum(trans_ARD$mean_trans),100),
                sector.index = sector.name, track.index = 2,minor.ticks=1)
  }, bg.border = NA)
  
  title("Fresh")
  
  #Marine
  arr.col = data.frame(rep(states[4],2), states[-c(3,4)],"black")
  #rep(cols[4],2))
  
  
  chordDiagram(trans_ARD$mean_trans, group=group,row.col = c(alpha(cols[1], 0.1),alpha(cols[2], 0.1),alpha(cols[3], 0.1),alpha(cols[4], 0.9)), annotationTrackHeight = c(0.25, .1), 
               column.col = cols, grid.col = cols, scale = F,
               directional = 1, direction.type = "arrows",link.arr.col = arr.col,link.arr.length = 0.1,link.arr.type = "big.arrow",  h.ratio = 0.2,big.gap = 3,annotationTrack = c("name", "grid"))
  
  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.axis(h = "top", labels.cex = 0.7, major.at=seq(0,sum(trans_ARD$mean_trans),100),
                sector.index = sector.name, track.index = 2,minor.ticks=1)
  }, bg.border = NA)
  
  
  title("Marine")
}
dev.off()