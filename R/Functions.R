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

subset_unite<-function(myDiff_obj,unite_obj){
  pos_of_interest_trt <- paste(myDiff_obj[["chr"]],myDiff_obj[["start"]])
  
  #identify no. pf CpGs
  unite_obj$start

  length(myDiff_obj$chr)
  #table(unite_obj$chr == myDiff_obj$chr & 
  #       unite_obj$start == myDiff_obj$start)
  
  #see if there are any duplicates
  table(duplicated(paste(unite_obj$chr,unite_obj$start)))
  
  #check wether number of True's == no of CpGs identified above
  table(paste(unite_obj$chr, unite_obj$start) %in% paste(myDiff_obj$chr, myDiff_obj$start))
  
  #add columns to unite object and my diff object (by combining chrom and pos)
  unite_obj$pos <- paste(unite_obj$chr, unite_obj$start)
  myDiff_obj$pos <- paste(myDiff_obj$chr, myDiff_obj$start)
  #make object containing CpGs which have matching pos in both unite and myDiff 
  rows_to_keep<-which(unite_obj$pos %in% myDiff_obj$pos)
  head(rows_to_keep)
  
  #subset the unite
  methyl_unite_rows_selected_trt<-methylKit::select(x = unite_obj,i = rows_to_keep)
}

percMeth_for_trt<-function(methobj,trt){
  #create object which will subset by treatment
  obj2<-methylKit::reorganize(methobj, treatment = myMetaData$trtG1G2_NUM[myMetaData$trtG1G2_NUM %in% trt],
                              sample.ids = myMetaData$SampleID[myMetaData$trtG1G2_NUM %in% trt])
  
  #get perc meth from previous object
  perc_obj2<-percMethylation(obj2)
  #calcuate average of each CpG for each sample with the specific treatment
  row_mean_trt_1 <-rowMeans(perc_obj2,na.rm = T)
}

percMeth_for_trt_sex<-function(methobj,trt){
  obj2<-methylKit::reorganize(methobj, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% trt],
                              sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% trt])
  
  perc_obj2<-percMethylation(obj2)
  row_mean_trt_1 <-rowMeans(perc_obj2,na.rm = T)
  
}

