suppressMessages(
  {
    library(motifscanR)
    library(Biostrings)
    library(GenomicRanges)
    library(TFBSTools)
    library(extrafont)
    library(pcaMethods)
    library(RColorBrewer)
    library(ggrepel)
    library(ggplot2)
    }
  )

color_map<-function(xs,n_bins=10,cmap='RdYlBu'){
    vmin<-min(xs)
    vmax<-max(xs)
    values<-seq(vmin,vmax,(vmax-vmin)/n_bins)
    palette<-rev(brewer.pal(n_bins,cmap))
    color_map_df<-list()
    color_map_df[['color']]<-palette
    color_map_df[['value1']]<-values[1:n_bins]
    color_map_df[['value2']]<-values[2:(n_bins+1)]
    color_map_df<-as.data.frame(color_map_df)
    return(as.character(lapply(xs,function(x){if(x>=vmax){return(palette[length(palette)])}else if(x<=vmin){
        return(palette[1])}else{return(as.character(color_map_df$color[c(color_map_df$value1<=x)&c(color_map_df$value2>x)]))}})))
}

TF_activity_score_scatter_plot<-function(pcs,xs,motif,outdir){
    colors<-color_map(xs)
    png(gettextf('%s/%s_TF_activity_scatter_plot.png',outdir,motif),family='Arial')
    par(xpd=T,mar=par()$mar+c(0,0.5,0,7),cex.main=1,font.main=1)
    plot(scores(pcs)[,'PC1'],scores(pcs)[,'PC2'],
         bg=colors,
         col='black',pch=21,
         xlab=gettextf('PC 1 (%.1f%%)',R2cum(pcs)['PC1']*100),
         ylab=gettextf('PC 2 (%.1f%%)',(R2cum(pcs)['PC2']-R2cum(pcs)['PC1'])*100),
         main=motif,
         cex.lab=1.5,cex=2,cex.axis=1.5,cex.main=1.8)

    legend('topright',rev(c('Low','','Medium','','High')),
            inset=c(-0.42,0),
            pch=c(15),
            cex=1.5,
            col=rev(color_map(c(0,0.25,0.5,0.75,1))),ncol=1,
            xpd=T)
    dev.off()
}

Differential_TF_activity_analysis<-function(input,genome,
                                            interested_motifs,metadata,interested_variable,
                                            outdir,adjusted_p_value_cutoff=0.1,Mval_cutoff=0,top_n=20){
    input_files<-lapply(input,function(x){read.table(x,sep='\t',header=T)})
    flags<-unlist(lapply(input_files,function(x){sum(x[,ncol(x)]<adjusted_p_value_cutoff)!=0}))
    if(all(flags==F)){
        stop('Adjusted p-value cutoff is too small')
    }
    num_of_samples<-0
    for(i in colnames(input_files[[1]])){
        if(grepl('.read_cnt',i)){
            num_of_samples<-num_of_samples+1
        }
    }

    temp<-do.call(rbind,lapply(input_files[flags],function(x){
        x[x[,ncol(x)]<adjusted_p_value_cutoff,
          colnames(x)[1:(4+num_of_samples-1)]]
    }))

    print(dim(temp))




    zscores_matrix<-scale(t(temp[,c(4:(4+num_of_samples-1))]))
    print(dim(zscores_matrix))
    pcs<-pca(zscores_matrix,nPcs=2)
    print(pcs)
    return(1)



}

argv<-commandArgs(TRUE)

help_doc="
Usage: Rscript Differential_TF_activity_analysis.R [options]
Options:
    --input=CHARACTER                    Input file names, output from Hypervariable_analysis.R.
                                         Files separated by comma(e.g. ../Proximal_hypervariable_analysis.txt,../Distal_hypervariable_analysis.txt).
    --species=CHARACTER                  Genome version of input data, such as hg19, hg38, mm9, mm10.
    --metadata=CHARACTER                 Metadata for each sample (i.e. Grouping variables).
    --interested_variable=CHARACTER      One of the columns of Metadata.
    --interested_motifs=CHARACTERs       Show activity scores of these TF motifs in the plot, e.g. TP63.
    --outdir=CHARACTER                   Output directory.
    --adjusted_p_value_cutoff=DOUBLE     Adjusted p value used for identifing significant HVRs [Default: 0.1].
    --top_n=INTEGER                      Top ranked enriched motif showed in the plot. [Default: 20].
    -h,--help                            Show this help message and exit.
"

valid_keys<-c('--input','--species','--metadata','--interested_variable',
              '--interested_motifs','--outdir','--adjusted_p_value_cutoff','--top_n')

if(length(argv)==0){
    stop("At least one argument must be supplied")
}else if(length(argv)==1){
    if(argv[1]=='--help'|argv[1]=='-h'){
        cat(help_doc)
    }else{
        stop("Invalid argument")
    }
}else{
    argv.list<-list()
    argv.list[['--top_n']]<-20
    argv.list[['--adjusted_p_value_cutoff']]<-0.1
    for(i in argv){
        temp<-strsplit(i,'=')[[1]]
        if(length(temp)==2){
            key<-temp[1]
            value<-temp[2]
            if(key=='--adjusted_p_value_cutoff'){
                argv.list[[key]]<-as.numeric(value)
            }else if(key=='--interested_motifs'){
                argv.list[[key]]<-as.character(strsplit(value,',')[[1]])
            }else{
              argv.list[[key]]<-value
            }
        }
    }
    if(sum(names(argv.list) %in% valid_keys)==length(argv.list)){
        print(argv.list)
        input<-strsplit(argv.list[['--input']],',')[[1]]
        if(length(input)<1){
            stop('Require at least one input file')
        }
        species<-argv.list[['--species']]
        metadata<-argv.list[['--metadata']]
        interested_variable<-as.character(argv.list[['--interested_variable']])
        interested_motifs<-argv.list[['--interested_motifs']]
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
            stop("Invalid species")
        }
        Differential_TF_activity_analysis(input,
                                          genome,
                                          interested_motifs,
                                          metadata,
                                          interested_variable,
                                          outdir,
                                          adjusted_p_value_cutoff=adjusted_p_value_cutoff,top_n=top_n)
    }else{
        stop("Invalid argument")
    }
}
