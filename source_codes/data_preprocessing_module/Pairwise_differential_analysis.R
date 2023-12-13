suppressMessages(
  {
    library(MAnorm2)
    library(extrafont)
    library(scales)
    }
)

Pairwise_differential_analysis<-function(input,metadata,interest_variable,outdir,
                                          filtered_chromosomes=c('chrX','chrY','chrM'),fdr_cutoff=0.1){
    meta_data<-read.table(metadata,sep=',',header=T)
    classes<-as.character(unique(meta_data[[interest_variable]]))
    for(i in c(1:(length(classes)-1))){
        for(j in c((i+1):length(classes))){
            condition1<-classes[i]
            condition2<-classes[j]
            print(paste0('Differential analysis: ',condition1,'_vs_',condition2))
            Differential_analysis(input,
                                  metadata,
                                  interest_variable,
                                  condition1,condition2,
                                  outdir,
                                  filtered_chromosomes=filtered_chromosomes,
                                  fdr_cutoff=adjusted_p_value_cutoff)
        }
    }
}

Differential_analysis<-function(input,metadata,interest_variable,
                                condition1,condition2,outdir,
                                filtered_chromosomes=c('chrX','chrY','chrM'),fdr_cutoff=0.1){
    raw_reads_count<-read.table(input,header=T,sep='\t')
    raw_reads_count<-raw_reads_count[!c(raw_reads_count$chrom %in% filtered_chromosomes),]
    num_of_samples<-0
    for(i in colnames(raw_reads_count)){
        if(grepl('.read_cnt',i)){
            num_of_samples<-num_of_samples+1
        }
    }

    reads_count<-colnames(raw_reads_count)[c(4:(4+num_of_samples-1))]
    occupancy<-colnames(raw_reads_count)[c((4+num_of_samples):(4+num_of_samples+num_of_samples-1))]

    meta_data<-read.table(metadata,sep=',',header=T)
    bioConds_list<-list()
    temp_reads_count<-reads_count[meta_data[[interest_variable]]==condition1]
    temp_occupancy<-occupancy[meta_data[[interest_variable]]==condition1]
    norm<-normalize(raw_reads_count,
                    temp_reads_count,
                    temp_occupancy)
    bioConds_list[[condition1]]<-bioCond(norm[temp_reads_count],
                                         norm[temp_occupancy],name=condition1)

    temp_reads_count<-reads_count[meta_data[[interest_variable]]==condition2]
    temp_occupancy<-occupancy[meta_data[[interest_variable]]==condition2]
    norm<-normalize(norm,
                    temp_reads_count,
                    temp_occupancy)
    bioConds_list[[condition2]]<-bioCond(norm[temp_reads_count],
                                         norm[temp_occupancy],name=condition2)

    bioConds_list<-normBioCond(bioConds_list)

    bioConds_list<-fitMeanVarCurve(bioConds_list,method='local',occupy.only=TRUE)
    png(gettextf('%s/%s_vs_%s_MVC_plot.png',outdir,condition1,condition2),family='Arial')

    par(mar=c(5,6,4,2),cex.main=1,font.main=1)
    means<-c(bioConds_list[[condition1]]$sample.mean,bioConds_list[[condition2]]$sample.mean)
    log10vars<-log10(c(bioConds_list[[condition1]]$sample.var/bioConds_list[[condition1]]$fit.info$ratio.var,
                       bioConds_list[[condition2]]$sample.var/bioConds_list[[condition2]]$fit.info$ratio.var))
    smoothScatter(means,log10vars,
                  xlab=expression(paste('Observed mean ',log[2],'(read count)')),
                  ylab=expression(paste(Log[10],'(observed variance)')),
                  main=gettextf('d0=%.2f',bioConds_list[[condition1]]$fit.info$df.prior),cex.main=2,cex.lab=2,cex.axis=1.8)
    min_mean<-min(means)
    max_mean<-max(means)
    xs<-seq(min_mean,to=max_mean,length.out=1000)
    ys<-log10(bioConds_list[[condition1]]$fit.info$predict(xs))
    lines(xs,ys,col='red',lwd=3)

    dev.off()

    results<-diffTest(bioConds_list[[condition1]],bioConds_list[[condition2]])

    png(gettextf('%s/%s_vs_%s_MA_plot.png',outdir,condition1,condition2),family='Arial')
    par(mar=c(5,6,4,2),cex.main=1,font.main=1)
    A_vals<-rowMeans(results[c(1,2)])
    M_vals<-results[['Mval']]
    smoothScatter(A_vals,M_vals,
                  xlab='A value',ylab='M value',
                  cex.main=2,cex.lab=2,cex.axis=1.8)

    flag1<-c(results$padj<fdr_cutoff)&c(results$Mval<0)
    points(A_vals[flag1],M_vals[flag1],col='blue')

    flag2<-c(results$padj<fdr_cutoff)&c(results$Mval>0)
    points(A_vals[flag2],M_vals[flag2],col='red')

    legend("topright",c(gettextf('%s up (n=%s)',condition2,sum(flag2)),
                        gettextf('%s up (n=%s)',condition1,sum(flag1))),
           inset=0.0025,pch=c(21),cex=1.2,col=c('red','blue'),
           xpd=T)
    dev.off()

    write.table(cbind(raw_reads_count[c('chrom','start','end')],results),
                gettextf('%s/Differential_analysis_%s_vs_%s.txt',
                         outdir,condition1,condition2),
                sep='\t',quote=F,row.names=F)

}

argv<-commandArgs(TRUE)

help_doc="
Usage: Rscript Differential_analysis.R [options]
Options:
    --input=CHARACTER                    Input file names, output from profile_bins.
    --metadata=CHARACTER                 Metadata for each sample.
    --interested_variable=CHARACTER      One of the columns of Metadata.
    --outdir=CHARACTER                   Output directory.
    --filtered_chromosomes=CHARACTERs    Remove these chromosomes before analysis [Default: chrX,chrY,chrM].
    --adjusted_p_value_cutoff=DOUBLE     Adjusted p value used for identifing significant differential peaks [Default: 0.1].
    -h,--help                            Show this help message and exit.
"

valid_keys<-c('--input','--metadata','--interested_variable','--outdir',
              '--filtered_chromosomes','--adjusted_p_value_cutoff')

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
    argv.list[['--filtered_chromosomes']]<-c('chrX','chrY','chrM')
    argv.list[['--adjusted_p_value_cutoff']]<-0.1
    for(i in argv){
        temp<-strsplit(i,'=')[[1]]
        if(length(temp)==2){
            key<-temp[1]
            value<-temp[2]
            if(key=='--filtered_chromosomes'){
                argv.list[[key]]<-as.character(strsplit(value,',')[[1]])
            }else if(key=='--adjusted_p_value_cutoff'){
                argv.list[[key]]<-as.numeric(value)
            }else{
                argv.list[[key]]<-value
            }
        }
    }
    if(sum(names(argv.list) %in% valid_keys)==length(argv.list)){
        print(argv.list)
        input<-argv.list[['--input']]
        metadata<-argv.list[['--metadata']]
        interest_variable<-argv.list[['--interested_variable']]
        outdir<-argv.list[['--outdir']]
        filtered_chromosomes<-argv.list[['--filtered_chromosomes']]
        adjusted_p_value_cutoff<-argv.list[['--adjusted_p_value_cutoff']]

        Pairwise_differential_analysis(input,metadata,interest_variable,outdir,
                                       filtered_chromosomes=filtered_chromosomes,
                                       fdr_cutoff=adjusted_p_value_cutoff)

    }else{
        stop("Invalid argument")
    }
}
