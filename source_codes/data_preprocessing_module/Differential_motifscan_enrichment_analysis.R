suppressMessages(
  {
    library(motifscanR)
    library(Biostrings)
    library(GenomicRanges)
    library(extrafont)
    }
  )

Differential_motif_enrichment<-function(input_dir,genome,outdir,
                                        adjusted_p_value_cutoff=0.1,
                                        Mval_cutoff=0,top_n=20){
  for(input in Sys.glob(gettextf('%s/Differential_analysis*.txt',input_dir))){
    temp<-unlist(strsplit(input,"/"))
    temp1<-unlist(strsplit(temp[length(temp)],"\\."))
    temp2<-unlist(strsplit(temp1[1],"Differential_analysis_"))

    Differential_motif_enrichment_analysis(input,genome,outdir,
                                           adjusted_p_value_cutoff=adjusted_p_value_cutoff,
                                           Mval_cutoff=0,top_n=top_n,name=temp2[2])    
  }
}

Differential_motif_enrichment_analysis<-function(input,genome,outdir,
                                                 adjusted_p_value_cutoff=0.1,
                                                 Mval_cutoff=0,top_n=20,name=''){
  analysis_result<-read.table(input,sep='\t',header=T)
  if(sum(c(analysis_result$padj<adjusted_p_value_cutoff)&c(analysis_result$Mval<Mval_cutoff))==0|
     sum(c(analysis_result$padj<adjusted_p_value_cutoff)&c(analysis_result$Mval>Mval_cutoff))==0){
    print(input)
    png(gettextf('%s/%s_motif_enrichment_plot.png',outdir,name),family='Arial')
    plot(1,1)
    text(1,1,'Adjusted p-value cutoff and M value are too stringent')
    dev.off()
    print('Adjusted p-value cutoff and M value are too stringent')
    return(0)
  }

  temp1<-analysis_result[c(analysis_result$padj<adjusted_p_value_cutoff)&c(analysis_result$Mval<Mval_cutoff),]
  temp1<-temp1[(temp1$end-temp1$start)>=100,]
  condition1_specific_peaks<-GRanges(seqnames=c(as.character(temp1$chrom)),
                                     ranges=IRanges(start=c(temp1$start),
                                                    end=c(temp1$end)))
  PFMatrixList <- getJasparMotifs(species = genome@metadata$organism)

  motif_pos1<-motifScan(PFMatrixList,
                        condition1_specific_peaks,
                        genome=genome,
                        out='matches',
                        thread=6)

  message("condition1 motifScan finished...\n")

  temp2<-analysis_result[c(analysis_result$padj<adjusted_p_value_cutoff)&c(analysis_result$Mval>Mval_cutoff),]
  temp2<-temp2[(temp2$end-temp2$start)>=100,]  
  condition2_specific_peaks<-GRanges(seqnames=c(as.character(temp2$chrom)),
                                     ranges=IRanges(start=c(temp2$start),
                                                    end=c(temp2$end)))

  motif_pos2<-motifScan(PFMatrixList,
                        condition2_specific_peaks,
                        genome=genome,
                        out='matches',
                        thread=6)

  message("condition2 motifScan finished...\n")

  test_results <- motifEnrichment(motif_pos1, motif_pos2)
  test_results<-na.omit(test_results)

  test_results[test_results$Fold_change==Inf,'Fold_change']<-max(test_results[test_results$Fold_change!=Inf,]$Fold_change)

  png(gettextf('%s/%s_motif_enrichment_plot.png',outdir,name),family='Arial')

  par(mar=c(5,6,4,2),cex.main=1,font.main=1)
  plot(log2(test_results$Fold_change),-log10(test_results$Corrected_P_value),
       xlab=expression(paste(Log[2],'(fold change)')),
       ylab=expression(paste(Log[10],'(p-value)')),
       cex=1.8,cex.lab=1.8,cex.axis=1.8)
  for(i in order(test_results$Corrected_P_value)[c(1:top_n)]){
    if(log2(test_results$Fold_change[i])>0){
      points(log2(test_results$Fold_change[i]),
             -log10(test_results$Corrected_P_value[i]),col='red',cex=1.8)
      text(log2(test_results$Fold_change[i]),
           -log10(test_results$Corrected_P_value[i]),
           test_results$motif[i],col='red')
    }else{
      points(log2(test_results$Fold_change[i]),
             -log10(test_results$Corrected_P_value[i]),col='blue',cex=1.8)
      text(log2(test_results$Fold_change[i]),
           -log10(test_results$Corrected_P_value[i]),
           test_results$motif[i],col='blue')
    }
  }
  condition1<-strsplit(colnames(analysis_result)[4],'.mean')[[1]]
  condition2<-strsplit(colnames(analysis_result)[5],'.mean')[[1]]
  legend("top",c(gettextf('Enriched in %s',condition1),
                 gettextf('Enriched in %s',condition2)),
         inset=0.0025,pch=c(21),cex=1.5,col=c('red','blue'),
         xpd=T)

  dev.off()
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
    --outdir=CHARACTER                   Output directory
    --adjusted_p_value_cutoff=FLOAT      Adjusted p value used for identifing significant HVRs [Default: 0.1]
    --top_n=INTEGER                      Top ranked enriched motif showed in the plot [Default: 20].
    -h,--help                            Show this help message and exit
"

valid_keys<-c('--input','--species','--outdir','--adjusted_p_value_cutoff', '--top_n')

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
  for(i in argv){
    temp<-strsplit(i,'=')[[1]]
    if(length(temp)==2){
      key<-temp[1]
      value<-temp[2]
      if(key=='--top_n' | key=='--adjusted_p_value_cutoff'){
        argv.list[[key]]<-as.numeric(value)
      }else{
        argv.list[[key]]<-value
      }
    }
  }
  if(sum(names(argv.list) %in% valid_keys)==length(argv.list)){
    print(argv.list)
    input_dir<-argv.list[['--input']]
    species<-argv.list[['--species']]
    outdir<-argv.list[['--outdir']]
    adjusted_p_value_cutoff<-argv.list[['--adjusted_p_value_cutoff']]
    top_n<-argv.list[['--top_n']]
    if(species=='hg19'){
      suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
      genome<-BSgenome.Hsapiens.UCSC.hg19
    }else if(species=='hg38'){
      suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
      genome<-BSgenome.Hsapiens.UCSC.hg38
    }else if(species=='mm9'){
      suppressMessages(library(BSgenome.Mmusculus.UCSC.mm9))
      genome<-BSgenome.Mmusculus.UCSC.mm9
    }else if(species=='mm10'){
      suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
      genome<-BSgenome.Mmusculus.UCSC.mm10
    }else{
      stop("Invalid species",call.=TRUE)
    }
    Differential_motif_enrichment(input_dir,genome,outdir,
                                  adjusted_p_value_cutoff=adjusted_p_value_cutoff,
                                  Mval_cutoff=0,top_n=top_n)
  }else{
    stop("Invalid argument",call.=TRUE)
  }
}
