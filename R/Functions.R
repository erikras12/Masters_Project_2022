deroman <- function(x){ x %>% str_remove(.,"Gy_chr") %>% 
    ifelse(. %in% c("Un","M"),., as.roman(.) %>% 
             as.integer() %>% as.character() %>% str_pad(.,2,"left",0))
}

get_pos <- function(CHROM,genome){
  tibble(CHROM = CHROM, 
         POS = sample(x=1:genome$length[genome$chrom == CHROM],size=1))
}

makeManhattanPlots_mult_trt <- function(DMSfile, annotFile, GYgynogff, mycols=c("grey50","grey50","darkred","darkred"), 
         mytitle = "Manhattan plot of DMS"){
  #GA_genome.fa.sizes.txt is a file with chromosome sizes and names
  genome <- GYgynogff %>%
    mutate(chrom_nr=chrom %>% deroman(),
           chrom_order=factor(chrom_nr) %>% as.numeric()) %>% 
    arrange(chrom_order) %>%
    mutate(gstart=lag(length,default=0) %>% cumsum(),
           gend=gstart+length, 
           type=LETTERS[2-(chrom_order%%2)],
           gmid=(gstart+gend)/2)
  
  #genome without M re-type:
  genome2=genome[genome$chrom_nr!="M",] %>%
    mutate(type=rep(c("A","B"),length(length)/2))
  
  region=as.factor(ifelse(annotFile$prom==1,"promoter",
                          ifelse(annotFile$exon==1,"exon",
                                 ifelse(annotFile$intron==1, "intron","intergenic"))))
  
  mydata = tibble(chrom=DMSfile$chr,
                  pos=DMSfile$pos,
                  meth.diff=DMSfile$avg.meth,
                  region=region)
  
  # table(DMSfile$chr)## check that chrXIX and chrUN are well removed!!
  
  # join DMS and genomic position
  data = left_join(mydata, genome2) %>% 
    mutate(gpos = pos + gstart,abs(meth.diff)<15,"not significant","significant")
  # all signif
  
  #plot only significant DMS:
  ggplot()+
    geom_rect(data=genome2,aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf,fill=type), alpha=.2)+
    geom_point(data=data[abs(data$meth.diff)>0,],
               aes(x=gpos,y=meth.diff,col=region,shape=region),fill="white", size = 2)+
    scale_color_manual(values = mycols)+
    scale_shape_manual(values=c(21,21,21,21))+
    scale_fill_manual(values=c(A=rgb(.9,.9,.9),B=NA),guide="none")+
    scale_x_continuous(breaks=genome2$gmid,labels=genome2$chrom %>% str_remove(.,"Gy_chr"),
                       position = "top",expand = c(0,0))+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.line=element_blank(),
          axis.title = element_blank(),
          strip.placement = "outside")+
    ggtitle(mytitle)
}

makeManhattanPlots <- function(DMSfile, annotFile, GYgynogff, mycols=c("grey50","grey50","darkred","darkred"), 
                               mytitle = "Manhattan plot of DMS"){
  #GA_genome.fa.sizes.txt is a file with chromosome sizes and names
  genome <- GYgynogff %>%
    mutate(chrom_nr=chrom %>% deroman(),
           chrom_order=factor(chrom_nr) %>% as.numeric()) %>% 
    arrange(chrom_order) %>%
    mutate(gstart=lag(length,default=0) %>% cumsum(),
           gend=gstart+length, 
           type=LETTERS[2-(chrom_order%%2)],
           gmid=(gstart+gend)/2)
  
  #genome without M re-type:
  genome2=genome[genome$chrom_nr!="M",] %>%
    mutate(type=rep(c("A","B"),length(length)/2))
  
  region=as.factor(ifelse(annotFile$prom==1,"promoter",
                          ifelse(annotFile$exon==1,"exon",
                                 ifelse(annotFile$intron==1, "intron","intergenic"))))
  
  mydata = tibble(chrom=DMSfile$chr,
                  pos=DMSfile$start,
                  meth.diff=DMSfile$meth.diff,
                  qval=DMSfile$qvalue,
                  region=region)
  
  # table(DMSfile$chr)## check that chrXIX and chrUN are well removed!!
  
  # join DMS and genomic position
  data = left_join(mydata, genome2) %>% 
    mutate(gpos=pos+gstart,significance= ifelse(abs(qval>0.0125) | abs(meth.diff)<15,"not significant","significant"))
  
  table(data$significance) # all signif
  
  #plot only significant DMS:
  ggplot()+
    geom_rect(data=genome2,aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf,fill=type), alpha=.2)+
    geom_point(data=data[abs(data$meth.diff)>15 & data$significance=="significant",],
               aes(x=gpos,y=meth.diff,col=region,shape=region),fill="white", size = 2)+
    scale_color_manual(values = mycols)+
    scale_shape_manual(values=c(21,21,21,21))+
    scale_fill_manual(values=c(A=rgb(.9,.9,.9),B=NA),guide="none")+
    scale_x_continuous(breaks=genome2$gmid,labels=genome2$chrom %>% str_remove(.,"Gy_chr"),
                       position = "top",expand = c(0,0))+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.line=element_blank(),
          axis.title = element_blank(),
          strip.placement = "outside")+
    ggtitle(mytitle)
}
