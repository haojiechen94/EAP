suppressMessages(
  {
    library(MAnorm2)
    library(extrafont)
    library(pcaMethods)
    library(scales)
    library(Rtsne)
    }
)

#calculate single end p-value
single.end.pvalue<-function(biocond){
    df<-attr(biocond,"df")
    p.value<-pf(biocond$fold.change,df[1],df[2],lower.tail=F)
    return(p.value)
}

#plot MVC and highlight HVRs
plot_MVC_and_HVRs<-function(MVC_fitting,outdir,type,flag){
    png(gettextf('%s/%s_MVC_and_HVRs_plot.png',outdir,type),family='Arial')
    par(mar=c(5,6,4,2),cex.main=1,font.main=1)
    smoothScatter(MVC_fitting[[1]]$sample.mean,log10(MVC_fitting[[1]]$sample.var),
                  xlab=expression(paste('Observed mean ',log[2],'(read count)')),
                  ylab=expression(paste(Log[10],'(observed variance)')),
                  main=gettextf('%s regions',type),cex.main=2,cex.lab=2,cex.axis=1.8)
    xs<-seq(min(MVC_fitting[[1]]$sample.mean),max(MVC_fitting[[1]]$sample.mean),0.1)
    lines(xs,log10(MVC_fitting[[1]]$fit.info$predict(xs)),lwd=3,col='red')
    points(MVC_fitting[[1]]$sample.mean[flag],log10(MVC_fitting[[1]]$sample.var)[flag],col='red')
    legend("bottomleft",gettextf('Significant HVRs (n=%s)',sum(flag)),
           inset=0.0025,
           pch=c(21),
           cex=1.5,
           col='red',
           xpd=T)
    dev.off()

    pdf(gettextf('%s/%s_MVC_and_HVRs_plot.pdf',outdir,type))
    par(mar=c(5,6,4,2),cex.main=1,font.main=1)
    smoothScatter(MVC_fitting[[1]]$sample.mean,log10(MVC_fitting[[1]]$sample.var),
                  xlab=expression(paste('Observed mean ',log[2],'(read count)')),
                  ylab=expression(paste(Log[10],'(observed variance)')),
                  main=gettextf('%s regions',type),cex.main=2,cex.lab=2,cex.axis=1.8)
    xs<-seq(min(MVC_fitting[[1]]$sample.mean),max(MVC_fitting[[1]]$sample.mean),0.1)
    lines(xs,log10(MVC_fitting[[1]]$fit.info$predict(xs)),lwd=3,col='red')
    points(MVC_fitting[[1]]$sample.mean[flag],log10(MVC_fitting[[1]]$sample.var)[flag],col='red')
    legend("bottomleft",gettextf('Significant HVRs (n=%s)',sum(flag)),
           inset=0.0025,
           pch=c(21),
           cex=1.5,
           col='red',
           xpd=T)
    dev.off()
}

#PCA analysis
PCA_analysis<-function(zscore_matrix,color_by,outdir){
    color_by<-as.character(color_by)
    unique_labels<-unique(color_by)
    color_list<-hue_pal()(length(unique_labels))
    names<-c()
    values<-c()
    for(i in c(1:length(unique_labels))){
        names<-c(names,unique_labels[i])
        values<-c(values,color_list[[i]])
    }
    color_map<-setNames(values,names)
    colors<-c()
    for(i in color_by){
        colors<-c(colors,color_map[[i]])
    }

    pcs<-pca(zscore_matrix,nPcs=2,method="nipals",cv="q2")

    png(gettextf('%s/PCA_scatter_plot.png',outdir),family='Arial')
    par(xpd=T,mar=par()$mar+c(0,0.5,0,7),cex.main=1,font.main=1)
    plot(round(scores(pcs)[,'PC1'],3),
         round(scores(pcs)[,'PC2'],3),
         bg=colors,
         col='black',pch=21,
         xlab=gettextf('PC 1 (%.1f%%)',pcs@R2[1]*100),
         ylab=gettextf('PC 2 (%.1f%%)',pcs@R2[2]*100),
         cex.lab=1.5,cex=2,cex.axis=1.5,cex.main=1.8)
    legend('topright',names,
           inset=c(-0.42,0),
           pch=c(19),
           cex=1.3,
           col=values,ncol=1,
           xpd=T)
    dev.off()

    pdf(gettextf('%s/PCA_scatter_plot.pdf',outdir))
    par(xpd=T,mar=par()$mar+c(0,0.5,0,7),cex.main=1,font.main=1)
    plot(round(scores(pcs)[,'PC1'],3),
         round(scores(pcs)[,'PC2'],3),
         bg=colors,
         col='black',pch=21,
         xlab=gettextf('PC 1 (%.1f%%)',pcs@R2[1]*100),
         ylab=gettextf('PC 2 (%.1f%%)',pcs@R2[2]*100),
         cex.lab=1.5,cex=2,cex.axis=1.5,cex.main=1.8)
    legend('topright',names,
           inset=c(-0.42,0),
           pch=c(19),
           cex=1.3,
           col=values,ncol=1,
           xpd=T)
    dev.off()

}

#TSNE visualization
TSNE_visualization<-function(zscore_matrix,top_number_of_PCs,perplexity,color_by,outdir){
    
    color_by<-as.character(color_by)
    unique_labels<-unique(color_by)
    color_list<-hue_pal()(length(unique_labels))
    names<-c()
    values<-c()
    for(i in c(1:length(unique_labels))){
        names<-c(names,unique_labels[i])
        values<-c(values,color_list[[i]])
    }
    color_map<-setNames(values,names)
    colors<-c()
    for(i in color_by){
        colors<-c(colors,color_map[[i]])
    }
    
    set.seed(1234)
    tsne_out<-Rtsne(zscore_matrix,dims=2,perplexity=perplexity,initial_dims=top_number_of_PCs)

    png(gettextf('%s/TSNE_plot.png',outdir),family='Arial')
    par(xpd=T,mar=par()$mar+c(0,0,0,7),cex.main=1,font.main=1)
    plot(tsne_out$Y[,1],tsne_out$Y[,2],
        col='black',pch=21,bg=colors,
        cex=1.5,cex.main=1.5,
        xlab='t-SNE dimension 1',ylab='t-SNE dimension 2',cex.lab=1.5,cex.axis=1.5)

    legend('topright',names,
            inset=c(-0.42,0),
            pch=c(19),
            cex=1.3,
            col=values,ncol=1,
            xpd=T)
    dev.off()

    pdf(gettextf('%s/TSNE_plot.pdf',outdir))
    par(xpd=T,mar=par()$mar+c(0,0,0,7),cex.main=1,font.main=1)
    plot(tsne_out$Y[,1],tsne_out$Y[,2],
        col='black',pch=21,bg=colors,
        cex=1.5,cex.main=1.5,
        xlab='t-SNE dimension 1',ylab='t-SNE dimension 2',cex.lab=1.5,cex.axis=1.5)

    legend('topright',names,
            inset=c(-0.42,0),
            pch=c(19),
            cex=1.3,
            col=values,ncol=1,
            xpd=T)
    dev.off()    

}

#HVRs analysis
Hypervariable_analysis<-function(input_proximal,input_distal,metadata,
                                 categorical_variable,outdir,
                                 top_number_of_PCs=2,
                                 perplexity=0,
                                 filtered_chromosomes=c('chrX','chrY','chrM'),fdr_cutoff=0.1){
    method<-'loc'
    proixmal_raw_reads_count<-read.table(input_proximal,header=T,sep='\t')
    proixmal_raw_reads_count<-proixmal_raw_reads_count[!c(proixmal_raw_reads_count$chrom %in% filtered_chromosomes),]
    distal_raw_reads_count<-read.table(input_distal,header=T,sep='\t')
    distal_raw_reads_count<-distal_raw_reads_count[!c(distal_raw_reads_count$chrom %in% filtered_chromosomes),]
    num_of_samples<-0
    for(i in colnames(proixmal_raw_reads_count)){
        if(grepl('.read_cnt',i)){
            num_of_samples<-num_of_samples+1
        }
    }
    if(perplexity==0){
         perplexity<-floor((num_of_samples-1)/3)
    }

    reads_count<-colnames(proixmal_raw_reads_count)[c(4:(4+num_of_samples-1))]
    occupancy<-colnames(proixmal_raw_reads_count)[c((4+num_of_samples):(4+num_of_samples+num_of_samples-1))]

    proixmal_normalized_data<-MAnorm2::normalize(proixmal_raw_reads_count,
                                                 reads_count,occupancy,offset=0.5,baseline='pseudo-reference')
    proixmal_normalized_data<-proixmal_normalized_data[rowMeans(proixmal_normalized_data[reads_count])>0.5,]
    proixmal_biocond<-bioCond(proixmal_normalized_data[reads_count],
                              proixmal_normalized_data[occupancy],name='Proximal')

    distal_normalized_data<-MAnorm2::normalize(distal_raw_reads_count,
                                               reads_count,occupancy,offset=0.5,baseline='pseudo-reference')
    distal_normalized_data<-distal_normalized_data[rowMeans(distal_normalized_data[reads_count])>0.5,]
    distal_biocond<-bioCond(distal_normalized_data[reads_count],
                            distal_normalized_data[occupancy],name='Distal')

    proximal_conds_list<-list(proixmal_biocond)
    proximal_MVC_fitting<-fitMeanVarCurve(proximal_conds_list,method=method,occupy.only=F,args.lp=list(nn=1.0))
    proixmal_biocond<-estParamHyperChIP(proximal_MVC_fitting[[1]])
    proixmal_varTest<-varTestBioCond(proixmal_biocond)
    proximal_p_values<-single.end.pvalue(proixmal_varTest)
    p_value<-proximal_p_values
    proximal_fdrs<-p.adjust(proximal_p_values,method='fdr')
    fdr<-proximal_fdrs
    proximal_result<-cbind(proixmal_normalized_data[c('chrom','start','end')],
                           proixmal_normalized_data[reads_count],
                           proixmal_normalized_data[occupancy],
                           proixmal_varTest[,c('observed.mean','observed.var','prior.var','fold.change')],
                           p_value,fdr)
    flag<-proximal_result$fdr<fdr_cutoff
    plot_MVC_and_HVRs(proximal_MVC_fitting,outdir,'Proximal',flag)

    distal_conds_list<-list(distal_biocond)
    distal_MVC_fitting<-fitMeanVarCurve(distal_conds_list,method=method,occupy.only=F,args.lp=list(nn=1.0))
    distal_biocond<-estParamHyperChIP(distal_MVC_fitting[[1]])
    distal_varTest<-varTestBioCond(distal_biocond)
    distal_p_values<-single.end.pvalue(distal_varTest)
    p_value<-distal_p_values
    distal_fdrs<-p.adjust(distal_p_values,method='fdr')
    fdr<-distal_fdrs
    distal_result<-cbind(distal_normalized_data[c('chrom','start','end')],
                         distal_normalized_data[reads_count],
                         distal_normalized_data[occupancy],
                         distal_varTest[,c('observed.mean','observed.var','prior.var','fold.change')],
                         p_value,fdr)
    flag<-distal_result$fdr<fdr_cutoff
    plot_MVC_and_HVRs(distal_MVC_fitting,outdir,'Distal',flag)

    proximal_zscore_matrix<-scale(t(proximal_result[proximal_result$fdr<fdr_cutoff,reads_count]))
    distal_zscore_matrix<-scale(t(distal_result[distal_result$fdr<fdr_cutoff,reads_count]))
    print('#proximal HVRs')
    print(sum(proximal_result$fdr<fdr_cutoff))
    print('#distal HVRs')
    print(sum(distal_result$fdr<fdr_cutoff))
    zscore_matrix<-cbind(proximal_zscore_matrix,distal_zscore_matrix)
    zscore_matrix<-round(zscore_matrix,5)

    meta_data<-read.table(metadata,sep=',',header=T)
    color_by<-meta_data[[categorical_variable]]
    PCA_analysis(zscore_matrix,color_by,outdir)
    TSNE_visualization(zscore_matrix,top_number_of_PCs,perplexity,color_by,outdir)
    write.table(proximal_result,gettextf('%s/Proximal_hypervariable_analysis.txt',outdir),sep='\t',quote=F,row.names=F,qmethod="double")
    write.table(distal_result,gettextf('%s/Distal_hypervariable_analysis.txt',outdir),sep='\t',quote=F,row.names=F,qmethod="double")
    save(proximal_result,file=gettextf('%s/Proximal_hypervariable_analysis.RData',outdir))
    save(distal_result,file=gettextf('%s/Distal_hypervariable_analysis.RData',outdir))

}

argv<-commandArgs(TRUE)

help_doc="
Usage: Rscript Hypervariable_analysis.R [options]
Options:
    --input=CHARACTER                    Input file names, output from profile_bins and separate_proximal_and_distal_peak_regions.py.
                                         Files separated by comma(i.e. ',', e.g. ../proximal_peaks.txt,../distal_peaks.txt).
    --metadata=CHARACTER                 Metadata for each sample.
    --categorical_variable=CHARACTER     One of columns in Metadata, assign colors to samples based on this categorical variable.
    --top_number_of_PCs=INTEGER          Top number of principal components used for TSNE dimesion reduction [Default: 2].
    --perplexity=INTEGER                 This parameter control how many nearest neighbours are taken into account when constructing the embedding in the low-dimensional space [Default: floor((num_of_samples-1)/3)].
    --outdir=CHARACTER                   Output directory.
    --filtered_chromosomes=CHARACTERs    Remove these chromosomes before analysis [Default: chrX,chrY,chrM].
    --adjusted_p_value_cutoff=DOUBLE     Adjusted p value used for identifing significant HVRs [Default: 0.1].
    -h,--help                            Show this help message and exit.
"

valid_keys<-c('--input','--metadata','--outdir','--categorical_variable',
              '--filtered_chromosomes','--adjusted_p_value_cutoff','--top_number_of_PCs','--perplexity')

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
    argv.list[['--top_number_of_PCs']]<-2
    argv.list[['--perplexity']]<-0
    for(i in argv){
        temp<-strsplit(i,'=')[[1]]
        if(length(temp)==2){
            key<-temp[1]
            value<-temp[2]
            if(key=='--filtered_chromosomes'){
                argv.list[[key]]<-as.character(strsplit(value,',')[[1]])
            }else if(key=='--adjusted_p_value_cutoff'){
                argv.list[[key]]<-as.numeric(value)
            }else if(key=='--top_number_of_PCs'){
              argv.list[[key]]<-as.integer(value)
            }else if(key=='--perplexity'){
              argv.list[[key]]<-as.integer(value)
            }else{
                argv.list[[key]]<-value
            }
        }
    }
    if(sum(names(argv.list) %in% valid_keys)==length(argv.list)){
        print(argv.list)
        temp<-strsplit(argv.list[['--input']],',')
        input_proximal<-temp[[1]][1]
        input_distal<-temp[[1]][2]
        outdir<-argv.list[['--outdir']]
        categorical_variable<-argv.list[['--categorical_variable']]
        metadata<-argv.list[['--metadata']]
        filtered_chromosomes<-argv.list[['--filtered_chromosomes']]
        adjusted_p_value_cutoff<-argv.list[['--adjusted_p_value_cutoff']]
        top_number_of_PCs<-argv.list[['--top_number_of_PCs']]
        perplexity<-argv.list[['--perplexity']]
        Hypervariable_analysis(input_proximal,
                               input_distal,
                               metadata,
                               categorical_variable,
                               outdir,
                               filtered_chromosomes=filtered_chromosomes,
                               fdr_cutoff=adjusted_p_value_cutoff,
                               top_number_of_PCs=top_number_of_PCs,
                               perplexity=perplexity)
    }else{
        stop("Invalid argument")
    }
}
