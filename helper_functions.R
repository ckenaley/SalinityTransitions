
### retrieve and summarize transitions on multiSimmap object
simmap_trans <- function (simm) 
{
  if(!"multiSimmap" %in% class(simm)) class(simm) <- "multiSimmap"
 
  trans <- lapply(simm,function(x) describe.simmap(x)$Tr)
  
  trans <- lapply(trans,add_miss)
  mean_trans <- Reduce('+', trans)/length(trans)
  
  states <- sort(unique(getStates(simm[[1]])))
  return(list(trans=trans,mean_trans=mean_trans,states=states))
}


#get transition times for multipSimmap object
get_trans_times <- function(simm){
  plot(simm[[1]],plot=F)
trans_time<- lapply(1:length(simm),function(x) {n <- x; 
maxH <- nodeHeights(simm[[n]]) %>% max;
markChanges(simm[[n]],plot=FALSE)  %>% data.frame %>% mutate(tree=n,maxH=maxH) %>% rownames_to_column(var="trans")}) %>% do.call(rbind,.) %>% mutate(trans=gsub("(\\w+)\\.\\.(\\w+)\\.*\\d*","\\1->\\2",trans)) %>% 
  mutate(trans=gsub("\\d+|\\.","",trans))
dev.off()
return(trans_time)
}

#get transition locations on tree

get_trans_loc <- function (tree,plot=F,obj_xx=obj_xx) 
{
  if(plot) plot(tree,plot=F)
  states <- sort(unique(getStates(tree)))
  nc <- sapply(tree$maps, length) - 1
  ii <- which(nc > 0)
  nc <- nc[ii]
  
  trans_l <- list()
  xx <- yy <- vector()
  for (i in 1:length(ii)) {
    for (j in 1:nc[i]) {
      ss <- names(tree$maps[[ii[i]]])[j + 1] #state
      mm <- tree$edge[ii[i], 1] #mother
      dd <- tree$edge[ii[i], 2] #daughter
      x <- rep(obj$xx[mm] + cumsum(tree$maps[[ii[i]]])[j], 
               2)
      trans <-  paste(names(tree$maps[[ii[i]]])[j:(j + 1)], collapse = "->")
      
      trans_l[[paste0(i,"-",j)]] <- tibble(trans=trans,x2=x[1],parent=mm,node=dd)
    }
  }
  trans <- do.call(rbind,trans_l)
}

#get earliest transition of some time after node, takes df of cols trans, node, x2 (node depth), see get_trans_loc()

get_early_trans <- function(phy,trans,to,n){
  trans %>% 
    filter(node %in%getDescendants(phy,n)) %>% 
    filter(grepl(to,trans)) %>% 
    arrange(x2)
}


#get number of transitions

get_n <- function(x){
  nm <- rownames(x)
  x2 <- x %>% data.frame()
  colnames(x2) <- rownames(x) 
  from <- lapply(1:length(nm), function(x) rep(nm[x],length(nm))) %>% unlist
  
  x3 <- x2 %>% pivot_longer(cols = 1:length(nm)) %>% mutate(from=from) %>% 
    group_by_all() %>% 
    mutate(from_to=paste(from,name,sep="->")) %>% 
    filter(name!=from) %>% ungroup %>% select(from_to,value)
  
  return(x3)
}

#get transiton from multisimmap object
get_trans <- function(x){
  capture.output(plot(x,plot=F))
  c <- invisible(markChanges(x,plot=F) %>% data.frame)
  c <- c %>% mutate(trans=rownames(c)) %>% mutate(trans=str_replace(trans,"\\.\\.","->"))%>% mutate(trans=str_replace(trans,"\\.\\d+",""))
  rownames(c) <- NULL
  return(c)
}

#get states at nodes for multisimmap, return as data.frame.
get_states <- function(x){
  s <- getStates(x)
  s <-data.frame(node=names(s) %>% as.numeric(),state=s)
  rownames(s) <- NULL
  return(s)
}


#plot simmap transitions of interest over clades, highlighting some, collapsing others
gg_simmap_clade <- function(simms,clade,coll,highlight,trans_hl="->Fresh"){
  #extract the clade
  class(simms) <- c("multiSimmap","multiPhylo")
  clade_simms <- lapply(simms,function(x) extract.clade.simmap(x,clade))
  
  class(clade_simms) <- c("multiSimmap","multiPhylo")
  
  
  clade_trans_times <- get_trans_times(clade_simms)
  
  plot(clade_simms[[1]],plot=F)
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  obj_xx <- obj$xx
  
  
  clade_trans <- lapply(clade_simms,get_trans_loc) %>% do.call(rbind,.)
  
  clade_trans_sum <- clade_trans %>% 
    group_by(trans,node,parent) %>% 
    summarize(x2=median(x2)) 
  
  
  coll_n <- coll
  
  trans_hl2 <- clade_trans$trans[grep(trans_hl,clade_trans$trans)] %>% unique
  
  if(!is.null(highlight)){
  hl_trans <- lapply(highlight,function(x) get_early_trans(phy=clade_simms[[1]],trans=clade_trans_sum,to=trans_hl,n=x))
  
 
  hl_paths <- lapply(1:length(hl_trans), function(x) 
    hl_trans[[x]] %>% 
      group_by(node,parent,trans) %>% 
      summarize(path=nodepath(clade_simms[[1]],highlight[x],node)) %>% mutate(clade=names(highlight)[x])
  ) %>% do.call(rbind,.)
  
  clade_trans2 <- clade_trans%>% 
    filter(node%in% (hl_paths %>% pull(path)),trans%in% trans_hl2) %>% 
    left_join(hl_paths %>% group_by(parent,node) %>% summarize(clade=clade[1])) 
  
  clade_trans_sum2 <- clade_trans2 %>% 
    group_by(trans,node,parent,clade) %>% 
    summarise(n_trans=n(),med_x=median(x2))
  
  }else{
    clade_trans2 <- clade_trans %>% filter(trans%in%trans_hl2)
    clade_trans_sum2 <-   clade_trans2 %>%group_by(trans,node,parent) %>% 
      summarise(n_trans=n(),med_x=median(x2))
  }
  clade_ace <- summary(clade_simms)$ace
  
  clade_nodes <- clade_ace%>% data.frame %>% 
    mutate(node=rownames(clade_ace) %>% as.numeric) %>% 
    na.omit() %>% 
    pivot_longer(Euryhaline:Marine) %>% 
    group_by(node) %>% 
    summarize(state=name[which.max(value)])

 
  clade_phy <- ggtree(clade_simms[[1]],ladderize = F,col="#7a7777",alpha=0.2)
  clade_phy <-  clade_phy %<+%  clade_trans2%<+%  clade_trans_sum2 %<+% clade_nodes
  
  #estract the subclades
  
  
  if(!is.null(coll_n)){
    p_col <- scaleClade(clade_phy,coll_n[1], .4,vertical_only = T)  
    
    for(i in coll_n[-1]){
      
      p_col <- scaleClade(p_col, i, .4,vertical_only = T)   
      
    }
    
    for(i in coll_n){
      
      p_col <- collapse(p_col,i, 'mixed', fill="gray90")    
      
    }
  
    clad_l <- list()
    for(i in 1:length(coll_n)){
      clad_l[[i]] <-  geom_cladelabel(label = names(coll_n[i]),node=coll_n[i],align=T)
    }
  }else{
    clad_l <- NULL
  p_col <- clade_phy}
  
  if(!is.null(highlight)){
    clad_l2 <- list()
    for(i in 1:length(highlight)){
      clad_l2[[i]] <-  geom_cladelabel(label = names(highlight[i]),node=highlight[i],align=T)
    }
    
    
  }else{
    clad_l2 <- NULL
    
  }
    #p_col$data$trans <- gsub("Euryhaline-Diadromous","EuryD", p_col$data$trans)
    #p_col$data$trans <- gsub("Euryhaline","Eury", p_col$data$trans)
    p_col_trans <- p_col+clad_l+clad_l2

    
  return(list(plot=p_col_trans,trans_data=clade_trans2,trans_sum=clade_trans_sum2))
  
  
}


add_miss <- function(x,ord=c("Euryhaline","Euryhaline-Diadromous","Fresh","Marine")){
  not_in <- setdiff(ord,rownames(x))
  not_n <- length(not_in)
  
  new_mat <- matrix(0,ncol=4,nrow = 4)
  new_mat[1:nrow(x),1:nrow(x)] <- x
  
  rownames(new_mat) <- colnames(new_mat) <- c(rownames(x),not_in)
  new_mat <- new_mat[ord,ord]
  return(new_mat)
}


bee_simms <- function(simms,clade=fw_oval[2],n=3,inc=trans_inc){
 
  clade_simms <- lapply(simms,function(x) extract.clade.simmap(x,clade))
  class(clade_simms) <- c("multiSimmap","multiPhylo")
  
  
  clade_trans_times <- get_trans_times(clade_simms)
  
  clade_trans_tt <- clade_trans_times %>% 
    group_by(tree) %>% 
    mutate(t=maxH-x) %>% 
    dplyr::arrange(t*-1,.by_group=TRUE) %>%
    group_by(trans,tree) %>% 
    dplyr::mutate(trans_n=1:n()) %>% 
    group_by(trans,tree) %>% 
    dplyr::mutate(trans_sum=cumsum(trans_n) ,trans_prop=trans_sum/max(trans_sum)) 
  
  
  clade_mean <- clade_trans_tt %>%
    group_by(tree,trans) %>% 
    dplyr::mutate(n=1:n()) %>%
    dplyr::arrange(tree,trans) %>% 
    dplyr::mutate(sum_n=sum(n)) %>% 
    group_by(trans) %>%
    dplyr::summarize(m_n=mean(sum_n),m_t=mean(t)) 
  
  clade_d <- clade_trans_tt %>%
    dplyr::mutate(keep=1,clade=names(clade)) %>% 
    dplyr::left_join(clade_mean) %>% 
    filter(m_n>n)
  
  dum <- tibble(trans=inc) %>% group_by(trans) %>% dplyr::summarize(t=runif(min=0,max=max(clade_trans_tt$t),1000))
  
  
  
  p <-clade_d %>% 
    ggplot( aes(x = trans, y =t))  +
    geom_quasirandom(alpha=0.2,dodge.width=1)+theme_classic(15)+xlab("Transition")+ylab("MYA")+theme(axis.text.x = element_text(angle=70,vjust = 1,hjust = 1))+ylim(c(0,max(clade_trans_tt$t)))+geom_quasirandom(data=dum,dodge.width=1,alpha=0)
  
  
  print(p)
  return(list(p=p,trans=clade_d))
}

#quickly draw a chor diagram

quick_chord <- function(trans,tick=1){

  mean_trans <-   mean_trans <- Reduce('+', trans)/length(trans)
  states <- rownames(mean_trans)
  cols <- setNames(pal_aaas("default")(4),states)
  
  arr.col <- data.frame(rep(states[1],3), states[-1],
                        c("black", "black", "black"))
  
  
  chordDiagram(mean_trans,row.col = c(alpha(cols[1], 0.6),alpha(cols[2], 0.6),alpha(cols[3], 0.6),alpha(cols[4], 0.6)), annotationTrackHeight = c(0.15, .1), column.col = cols, grid.col = cols, scale = F, directional = 1, direction.type = "arrows",link.arr.col = arr.col,link.arr.length = 0.1,link.arr.type = "big.arrow",  h.ratio = 0.1,big.gap = 3,annotationTrack = c("name", "grid"))
  
  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.axis(h = "top", labels.cex = 0.5, major.at=seq(0,sum(mean_trans),tick),
                sector.index = sector.name, track.index = 2,minor.ticks=tick/2)
  }, bg.border = NA)
  
  
}
