suppressMessages(
  {
    library(motifscanR)
    library(Biostrings)
    library(GenomicRanges)
    library(extrafont)
    library(ggrepel)
    library(ggplot2)    
    }
  )

Differential_motif_enrichment_analysis<-function(input,genome,outdir,
                                                 adjusted_p_value_cutoff=0.1,
                                                 Mval_cutoff=0,top_n=20,
                                                 cutoff.matrix.loc='./',cutoff.matrix.name=NULL){
  analysis_result<-read.table(input,sep='\t',header=T)
  if(sum(c(analysis_result$padj<adjusted_p_value_cutoff)&c(analysis_result$Mval<Mval_cutoff))==0|
     sum(c(analysis_result$padj<adjusted_p_value_cutoff)&c(analysis_result$Mval>Mval_cutoff))==0){
    png(gettextf('%s/error_information.png',outdir),family='Arial')
    plot(1,1,col='white')
    text(1,1,'Adjusted p-value cutoff and M value are too stringent')
    dev.off()
    print('Adjusted p-value cutoff and M value are too stringent')
    return(0)
  }

  temp1<-analysis_result[c(analysis_result$padj<adjusted_p_value_cutoff)&c(analysis_result$Mval< -Mval_cutoff),]
  #filtering peaks with size less than 100bp
  temp1<-temp1[(temp1$end-temp1$start)>=100,]
  print(gettextf('#upregulated peaks in condition1: %s',dim(temp1)[[1]]))
  if(length(temp1)==0){
    png(gettextf('%s/error_information.png',outdir),family='Arial')
    plot(1,1,col='white')
    text(1,1,'No peaks left after filtered by peak size in condition1')
    dev.off()
    print('No peaks left after filtered by peak size in condition1')
    return(0)    
  }
  condition1_specific_peaks<-GRanges(seqnames=c(as.character(temp1$chrom)),
                                     ranges=IRanges(start=c(temp1$start),
                                                    end=c(temp1$end)))

  PFMatrixList_list<-c(getJasparMotifs(species="Homo sapiens",collection="CORE"),
                      #getJasparMotifs(species="Homo sapiens",collection="CNE"),
                      getJasparMotifs(species="Homo sapiens",collection="UNVALIDATED"),
                      getJasparMotifs(species="Mus musculus",collection="CORE"),
                      getJasparMotifs(species="Mus musculus",collection="UNVALIDATED"),
                      getJasparMotifs(species="Mus musculus",collection="PBM_HOMEO"))

  total_tf_names <- unlist(lapply(PFMatrixList_list, function(pwm) pwm@name))
  tf_names_collect <- lapply(unique(total_tf_names), function(x) names(total_tf_names)[which(total_tf_names %in% x)])
  unique_names <- unlist(lapply(tf_names_collect, function(x) x[1]))
  PFMatrixList <- PFMatrixList_list[unique_names]

  motif_pos1<-motifScan(PFMatrixList,
                        condition1_specific_peaks,
                        genome=genome,
                        out='matches',
                        cutoff.matrix.loc=cutoff.matrix.loc,
                        cutoff.matrix.name=cutoff.matrix.name)

  message("condition1 motifScan finished...\n")

  temp2<-analysis_result[c(analysis_result$padj<adjusted_p_value_cutoff)&c(analysis_result$Mval>Mval_cutoff),]
  temp2<-temp2[(temp2$end-temp2$start)>=100,]
  print(gettextf('#upregulated peaks in condition2: %s',dim(temp2)[[1]]))
  if(length(temp2)==0){
    png(gettextf('%s/error_information.png',outdir),family='Arial')
    plot(1,1,col='white')
    text(1,1,'No peaks left after filtered by peak size in condition2')
    dev.off()
    print('No peaks left after filtered by peak size in condition2')
    return(0)    
  }  
  condition2_specific_peaks<-GRanges(seqnames=c(as.character(temp2$chrom)),
                                     ranges=IRanges(start=c(temp2$start),
                                                    end=c(temp2$end)))

  motif_pos2<-motifScan(PFMatrixList,
                        condition2_specific_peaks,
                        genome=genome,
                        out='matches',
                        cutoff.matrix.loc=cutoff.matrix.loc,
                        cutoff.matrix.name=cutoff.matrix.name)

  message("condition2 motifScan finished...\n")

  test_results <- motifEnrichment(motif_pos2,motif_pos1)
  test_results<-na.omit(test_results)

  test_results[test_results$Fold_change==Inf,'Fold_change']<-max(test_results[test_results$Fold_change!=Inf,]$Fold_change)
  test_results$Annotation<-'3-Other'
  rows<-order(test_results$Corrected_P_value)[c(1:top_n)]
  red<-rownames(test_results)[c(test_results[['Fold_change']]>1)&
                              c(test_results[['Corrected_P_value']]<0.01)]
  blue<-rownames(test_results)[c(test_results[['Fold_change']]<1)&
                              c(test_results[['Corrected_P_value']]<0.01)]      
  condition1<-strsplit(colnames(analysis_result)[4],'.mean')[[1]]
  condition2<-strsplit(colnames(analysis_result)[5],'.mean')[[1]]
  test_results[red,'Annotation']<-gettextf('1-Enriched in %s (n=%d)',condition2,length(red))
  test_results[blue,'Annotation']<-gettextf('2-Enriched in %s (n=%d)',condition1,length(blue))

  color_list<-c()
  if(length(red)>0){
    color_list<-c(color_list,'red')
  }

  if(length(blue)>0){
    color_list<-c(color_list,'blue')
  }

  color_list<-c(color_list,'grey')

  g.plot<-ggplot(test_results, aes(x = log2(Fold_change), y = -log10(Corrected_P_value))) +
  geom_point(aes(color = Annotation)) +
  scale_color_manual(values = color_list ) +
  theme_bw(base_size = 40) + theme(legend.position = "right")+
  geom_text_repel(
      data = subset(test_results, rownames(test_results) %in% rownames(test_results[rows,])),
      aes(label = Motif),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      max.overlaps=100
    ) +
  xlab(expression(paste(Log[2],'(fold change)')))+
  ylab(expression(paste(Log[10],'(corrected p-value)')))

  ggsave(g.plot,file=gettextf('%s/Motif_enrichment_plot.pdf',outdir),width=18,height=12)
  ggsave(g.plot,file=gettextf('%s/Motif_enrichment_plot.png',outdir),width=18,height=12)

  colnames(test_results) <- c('Motif',
                              paste0('Num_', condition1), paste0('Num_', condition2),
                              'Fold_change', 'Enriched_P_value',
                              'Depleted_P_value', 'Corrected_P_value')

  write.table(test_results,gettextf('%s/Differential_motif_enrichment_analysis_%s_vs_%s.txt',
                                    outdir,condition1,condition2),
              sep='\t',quote=F,row.names=F)

}

argv<-commandArgs(TRUE)

help_doc="
Usage: Rscript Differential_motif_enrichment_analysis.R [options]
Options:
    --input=CHARACTER                    Input file name, output from Differential_analysis.R.
    --species=CHARACTER                  Genome version of input data, such as hg19, hg38, mm9, mm10.
    --cutoff_matrix_loc=CHARACTER        Path to motif cutoff matrix.
    --outdir=CHARACTER                   Output directory.
    --adjusted_p_value_cutoff=FLOAT      Adjusted p value used for identifing significant HVRs [Default: 0.1].
    --log2_fold_change=INTEGER           Log 2 fold change used for identifing significant differential peaks [Default: 0].
    --top_n=INTEGER                      Top ranked enriched motif showed in the plot [Default: 20].
    -h,--help                            Show this help message and exit.
"

valid_keys<-c('--input','--species','--cutoff_matrix_loc','--outdir','--adjusted_p_value_cutoff', '--top_n',
              '--log2_fold_change')

if(length(argv)==0){
  stop("At least one argument must be supplied",call.=TRUE)
}else if(length(argv)==1){
  if(argv[1]=='--help'|argv[1]=='-h'){
    cat(help_doc)
  }else{
    stop("Invalid argument",call.=TRUE)
  }
}else{
  argv.list<-list()
  argv.list[['--adjusted_p_value_cutoff']]<-0.1
  argv.list[['--top_n']]<-20
  argv.list[['--log2_fold_change']]<-0
  for(i in argv){
    temp<-strsplit(i,'=')[[1]]
    if(length(temp)==2){
      key<-temp[1]
      value<-temp[2]
      if(key=='--top_n'){
          argv.list[[key]]<-as.integer(value)
      }else if(key=='--adjusted_p_value_cutoff'){
          argv.list[[key]]<-as.numeric(value)
      }else if(key=='--log2_fold_change'){
          argv.list[[key]]<-as.numeric(value)
      }else{
          argv.list[[key]]<-value
      }
    }
  }
  if(sum(names(argv.list) %in% valid_keys)==length(argv.list)){
    print(argv.list)
    input<-argv.list[['--input']]
    species<-argv.list[['--species']]
    cutoff.matrix.loc<-argv.list[['--cutoff_matrix_loc']]
    outdir<-argv.list[['--outdir']]
    adjusted_p_value_cutoff<-argv.list[['--adjusted_p_value_cutoff']]
    log2_fold_change<-argv.list[['--log2_fold_change']]
    top_n<-argv.list[['--top_n']]
    if(species=='hg19'){
      suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
      genome<-BSgenome.Hsapiens.UCSC.hg19
      cutoff.matrix.name<-'BSgenome_hg19_humanANDmouse_motifset'
    }else if(species=='hg38'){
      suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
      genome<-BSgenome.Hsapiens.UCSC.hg38
      cutoff.matrix.name<-'BSgenome_hg38_humanANDmouse_motifset'
    }else if(species=='mm9'){
      suppressMessages(library(BSgenome.Mmusculus.UCSC.mm9))
      genome<-BSgenome.Mmusculus.UCSC.mm9
      cutoff.matrix.name<-'BSgenome_mm9_humanANDmouse_motifset'
    }else if(species=='mm10'){
      suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
      genome<-BSgenome.Mmusculus.UCSC.mm10
      cutoff.matrix.name<-'BSgenome_mm10_humanANDmouse_motifset'
    }else{
      stop("Invalid species",call.=TRUE)
    }
    Differential_motif_enrichment_analysis(input,genome,outdir,
                                           adjusted_p_value_cutoff=adjusted_p_value_cutoff,
                                           Mval_cutoff=log2_fold_change,top_n=top_n,
                                           cutoff.matrix.loc=cutoff.matrix.loc,cutoff.matrix.name=cutoff.matrix.name)
  }else{
    stop("Invalid argument",call.=TRUE)
  }
}
