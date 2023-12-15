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
    library(Rtsne)
    library(scales)
    library(ComplexHeatmap)    
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


PCA_visualization<-function(pcs,color_by,outdir){
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

    png(gettextf('%s/PCA_scatter_plot_1.png',outdir),family='Arial')
    par(xpd=T,mar=par()$mar+c(0,0.5,0,7),cex.main=1,font.main=1)
    plot(scores(pcs)[,'PC1'],scores(pcs)[,'PC2'],
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

    pdf(gettextf('%s/PCA_scatter_plot_1.pdf',outdir))
    par(xpd=T,mar=par()$mar+c(0,0.5,0,7),cex.main=1,font.main=1)
    plot(scores(pcs)[,'PC1'],scores(pcs)[,'PC2'],
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
    write.table(scores(pcs)[,c('PC1','PC2')],gettextf('%s/PCA_plot.txt',outdir),sep='\t',quote=FALSE)

}

TF_activity_score_scatter_plot<-function(pcs,xs,motif,outdir){
    colors<-color_map(xs)
    png(gettextf('%s/%s_TF_activity_scatter_plot.png',outdir,motif),family='Arial')
    par(xpd=T,mar=par()$mar+c(0,0.5,0,7),cex.main=1,font.main=1)
    plot(scores(pcs)[,'PC1'],scores(pcs)[,'PC2'],
         bg=colors,
         col='black',pch=21,
         xlab=gettextf('PC 1 (%.1f%%)',pcs@R2[1]*100),
         ylab=gettextf('PC 2 (%.1f%%)',pcs@R2[2]*100),
         main=motif,
         cex.lab=1.5,cex=2,cex.axis=1.5,cex.main=1.8)

    legend('topright',rev(c('Low','','Medium','','High')),
            inset=c(-0.42,0),
            pch=c(15),
            cex=1.5,
            col=rev(color_map(c(0,0.25,0.5,0.75,1))),ncol=1,
            xpd=T)
    dev.off()

    pdf(gettextf('%s/%s_TF_activity_scatter_plot.pdf',outdir,motif))
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


#TSNE visualization
TF_activity_score_TSNE_plot<-function(tsne_out,xs,motif,outdir,perplexity){  
    colors<-color_map(xs)

    png(gettextf('%s/%s_TF_activity_score_TSNE_plot.png',outdir,motif),family='Arial')
    par(xpd=T,mar=par()$mar+c(0,0,0,7),cex.main=1,font.main=1)
    plot(tsne_out$Y[,1],tsne_out$Y[,2],
         bg=colors,
         col='black',pch=21,
         xlab='t-SNE dimension 1',ylab='t-SNE dimension 2',
         main=motif,
         cex.lab=1.5,cex=2,cex.axis=1.5,cex.main=1.8)

    legend('topright',rev(c('Low','','Medium','','High')),
            inset=c(-0.42,0),
            pch=c(15),
            cex=1.5,
            col=rev(color_map(c(0,0.25,0.5,0.75,1))),ncol=1,
            xpd=T)
    dev.off()

    pdf(gettextf('%s/%s_TF_activity_score_TSNE_plot.pdf',outdir,motif))
    par(xpd=T,mar=par()$mar+c(0,0,0,7),cex.main=1,font.main=1)
    plot(tsne_out$Y[,1],tsne_out$Y[,2],
         bg=colors,
         col='black',pch=21,
         xlab='t-SNE dimension 1',ylab='t-SNE dimension 2',
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


#TSNE visualization
TSNE_visualization<-function(original_zscores_matrix,top_number_of_PCs,perplexity,color_by,outdir){
    
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
    tsne_out<-Rtsne(original_zscores_matrix,dims=2,perplexity=perplexity,initial_dims=top_number_of_PCs)

    png(gettextf('%s/TSNE_plot_1.png',outdir),family='Arial')
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

    pdf(gettextf('%s/TSNE_plot_1.pdf',outdir))
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
    write.table(tsne_out$Y,gettextf('%s/TSNE_plot1.txt',outdir),sep='\t',quote=FALSE)
    return(tsne_out)
}

#heatmap
Top_ranked_TF_activity_heatmap<-function(results_df,temp_df,metadata,top_n,interested_variable,num_of_samples,outdir){
    columns<-c()
    selected_TF_motifs<-c()
    labels<-c()
    for(class in colnames(results_df)[c(1:(length(colnames(results_df))-1))]){
        columns<-c(columns,colnames(temp_df)[c(1:(length(colnames(temp_df))-1))][metadata[[interested_variable]]==class])
        labels<-c(labels,metadata[[interested_variable]][metadata[[interested_variable]]==class])
        for(motif in results_df[rev(order(results_df[[class]]))[c(1:top_n)],][['motif']]){
            if(! motif %in% selected_TF_motifs){
                selected_TF_motifs<-c(selected_TF_motifs,motif)
            }
        }
    }

    TF_scores<-list()
    rownames<-c()
    for(motif in selected_TF_motifs){
        TF_scores[[motif]]<-colMeans(temp_df[temp_df$motif == motif, columns])
        rownames<-c(rownames,motif)
    }
    TF_scores_df<-t(as.data.frame(TF_scores))

    color_by<-metadata[[interested_variable]]
    color_by<-as.character(color_by)
    unique_labels<-unique(color_by)
    color_list<-hue_pal()(length(unique_labels))
    names<-c()
    values<-c()
    for(i in c(1:length(unique_labels))){
        names<-c(names,unique_labels[i])
        values<-c(values,color_list[[i]])
    }
    color_map1<-setNames(values,names)

    heatmap_annotation<-HeatmapAnnotation(Annotation=labels,col=list(Annotation=color_map1))
    png(gettextf('%s/TF_activity_score_heatmap.png',outdir),
                 width= 30 * num_of_samples,height= 20 * top_n * (length(colnames(results_df))-1),family='Arial')
    draw(Heatmap(
        TF_scores_df, 
        cluster_rows=F,
        cluster_columns=F,
        top_annotation=heatmap_annotation,
        row_labels=rownames,
        heatmap_legend_param=list(title ='TF activity score')))
    dev.off()

    pdf(gettextf('%s/TF_activity_score_heatmap.pdf',outdir),
                 width= floor(0.5 * num_of_samples),height= floor(0.25 * top_n * (length(colnames(results_df))-1)))
    draw(Heatmap(
        TF_scores_df, 
        cluster_rows=F,
        cluster_columns=F,
        top_annotation=heatmap_annotation,
        row_labels=rownames,
        heatmap_legend_param=list(title ='TF activity score')))
    dev.off()    


}

Differential_TF_activity_analysis<-function(input,genome,
                                            interested_motifs,metadata,interested_variable,outdir,
                                            top_number_of_PCs=2,perplexity=0,
                                            adjusted_p_value_cutoff=0.1,Mval_cutoff=0,top_n=20,
                                            cutoff.matrix.loc='./',cutoff.matrix.name=NULL){
    load(input[1])
    load(input[2])
    print('#proximal HVRs')
    print(sum(proximal_result$fdr<adjusted_p_value_cutoff))
    print('#distal HVRs')
    print(sum(distal_result$fdr<adjusted_p_value_cutoff))

    num_of_samples<-0
    for(i in colnames(proximal_result)){
        if(grepl('.read_cnt',i)){
            num_of_samples<-num_of_samples+1
        }
    }
    if(perplexity==0){
         perplexity<-floor((num_of_samples-1)/3)
    }        

    proximal_zscore_matrix<-scale(t(proximal_result[proximal_result$fdr<adjusted_p_value_cutoff,
                                                    c(4:(4+num_of_samples-1))]))
    distal_zscore_matrix<-scale(t(distal_result[distal_result$fdr<adjusted_p_value_cutoff,
                                                c(4:(4+num_of_samples-1))]))
    zscore_matrix<-cbind(proximal_zscore_matrix,distal_zscore_matrix)

    if(dim(zscore_matrix)[[2]]==0){
        png(gettextf('%s/error_information.png',outdir),family='Arial')
        plot(1,1,col='white')
        text(1,1,'Adjusted p-value cutoff are too stringent')
        dev.off()
        print('Adjusted p-value cutoff are too stringent')
        return(0)
    }
    original_zscores_matrix<-round(zscore_matrix,5)


    temp<-rbind(proximal_result[proximal_result$fdr<adjusted_p_value_cutoff,],
                                distal_result[distal_result$fdr<adjusted_p_value_cutoff,])
    #filtering peaks with size less than 100bp
    temp<-temp[(temp$end-temp$start)>=100,]
    if(length(temp)==0){
        png(gettextf('%s/error_information.png',outdir),family='Arial')
        plot(1,1,col='white')
        text(1,1,'No peaks left after filtered by peak size')
        dev.off()
        print('No peaks left after filtered by peak size')
        return(0)    
    }

    hyper_variable_peaks<-GRanges(seqnames=c(as.character(temp$chrom)),
                                  ranges=IRanges(start=c(temp$start),
                                                 end=c(temp$end)))

    
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

    motif_pos<-motifScan(PFMatrixList,
                         hyper_variable_peaks,
                         genome=genome,
                         out='matches',
                         cutoff.matrix.loc=cutoff.matrix.loc,
                         cutoff.matrix.name=cutoff.matrix.name)
    message("motifScan finished...")
    match_pos<-as.matrix(motif_pos)
    colnames(match_pos)<-name(PFMatrixList)
    match_pos<-match_pos[,colSums(match_pos)>0]

    zscores_matrix<-scale(t(temp[,c(4:(4+num_of_samples-1))]))
    zscores_matrix<-zscores_matrix[,colSums(is.na(zscores_matrix))==0]

    means<-rowMeans(zscores_matrix)
    temp1<-list()
    for(sample in rownames(zscores_matrix)){
        temp2<-c()
        for(motif in colnames(match_pos)){
            activity<-mean(zscores_matrix[sample,match_pos[,motif]])-means[[sample]]
            temp2<-c(temp2,activity)
        }
        temp1[[sample]]<-temp2
    }
    temp_df<-as.data.frame(temp1)
    meta_data<-read.table(metadata,sep=',',header=T,stringsAsFactors=F)


    results<-list()
    for(class1 in unique(meta_data[[interested_variable]])){
        t.values<-c()
        interested_samples<-meta_data[[interested_variable]]==class1
        other_samples<-meta_data[[interested_variable]]!=class1
        for(i in c(1:dim(temp_df)[1])){
            a<-t.test(temp_df[i,interested_samples],temp_df[i,other_samples],alternative='greater')
            t.values<-c(t.values,a$statistic[['t']])
        }
        results[[class1]]<-t.values
    }
    results_df<-as.data.frame(results)

    for(class in colnames(results_df)){
        data<-list()
        data[['t.statistic']]<-results_df[[class]][rev(order(results_df[[class]]))]
        data[['Rank']]<-c(1:length(results_df[[class]]))
        data<-as.data.frame(data)

        ggplot(data[c(1:top_n),],aes(Rank,t.statistic,label=colnames(match_pos)[rev(order(results_df[[class]]))][c(1:top_n)]))+geom_text_repel(size = 5,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"),max.overlaps=100) +geom_point(color='red')+theme_classic(base_size=16)+ggtitle(gettextf("%s:",class))
        ggsave(gettextf('%s/%s_TF_activity_rank_plot.png',outdir,class))
        ggsave(gettextf('%s/%s_TF_activity_rank_plot.pdf',outdir,class))
    }

    results_df[['motif']]<-colnames(match_pos)
    write.table(results_df,gettextf('%s/Differential_TF_activity_analysis.txt',outdir),sep='\t',quote=F,row.names=F)

    pcs<-pca(original_zscores_matrix,nPcs=top_number_of_PCs,method="nipals",cv="q2")
    color_by<-meta_data[[interested_variable]]
    PCA_visualization(pcs,color_by,outdir)
    tsne_out<-TSNE_visualization(original_zscores_matrix,top_number_of_PCs,perplexity,color_by,outdir)   

    for(motif in interested_motifs){
        if(motif %in% colnames(match_pos)){
            xs<-as.vector(colMeans(temp_df[colnames(match_pos)==motif,]))
            TF_activity_score_scatter_plot(pcs,xs,motif,outdir)
            TF_activity_score_TSNE_plot(tsne_out,xs,motif,outdir,perplexity)
        }else{
            print(gettextf('Motif: %s not available.',motif))
            png(gettextf('%s/%s_error_information.png',outdir,motif))
            plot(1,1,col='white')
            text(1,1,gettextf('Motif: %s not available.',motif))
            dev.off()
        }
    }

    temp_df[['motif']]<-colnames(match_pos)
    write.table(temp_df,gettextf('%s/TF_activity_scores.txt',outdir),sep='\t',quote=F,row.names=F)

    Top_ranked_TF_activity_heatmap(results_df,temp_df,meta_data,top_n,interested_variable,num_of_samples,outdir)

}

argv<-commandArgs(TRUE)

help_doc="
Usage: Rscript Differential_TF_activity_analysis.R [options]
Options:
    --input=CHARACTER                    Input file names, output from Hypervariable_analysis.R.
                                         Files separated by comma(e.g. ../Proximal_hypervariable_analysis.RData,../Distal_hypervariable_analysis.RData).
    --species=CHARACTER                  Genome version of input data, such as hg19, hg38, mm9, mm10.
    --metadata=CHARACTER                 Metadata for each sample (i.e. Grouping variables).
    --interested_variable=CHARACTER      One of the columns of Metadata.
    --interested_motifs=CHARACTERs       Show activity scores of these TF motifs in the plot, e.g. TP63.
    --top_number_of_PCs=INTEGER          Top number of principal components used for TSNE dimesion reduction [Default: 2].
    --perplexity=INTEGER                 This parameter control how many nearest neighbours are taken into account when constructing the embedding in the low-dimensional space [Default: floor((num_of_samples-1)/3)].    
    --cutoff_matrix_loc=CHARACTER        Path to motif cutoff matrix.    
    --outdir=CHARACTER                   Output directory.
    --adjusted_p_value_cutoff=DOUBLE     Adjusted p value used for identifing significant HVRs [Default: 0.1].
    --top_n=INTEGER                      Top ranked enriched motif showed in the plot. [Default: 20].
    -h,--help                            Show this help message and exit.
"

valid_keys<-c('--input','--species','--metadata','--interested_variable',
              '--interested_motifs','--outdir','--adjusted_p_value_cutoff','--top_n',
              '--cutoff_matrix_loc','--top_number_of_PCs','--perplexity')

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
    argv.list[['--top_number_of_PCs']]<-2
    argv.list[['--perplexity']]<-0
    for(i in argv){
        temp<-strsplit(i,'=')[[1]]
        if(length(temp)==2){
            key<-temp[1]
            value<-temp[2]
            if(key=='--adjusted_p_value_cutoff'){
                argv.list[[key]]<-as.numeric(value)
            }else if(key=='--interested_motifs'){
                argv.list[[key]]<-as.character(strsplit(value,',')[[1]])
            }else if(key=='--top_number_of_PCs'){
              argv.list[[key]]<-as.integer(value)
            }else if(key=='--perplexity'){
              argv.list[[key]]<-as.integer(value)
            }else if(key=='--top_n'){
              argv.list[[key]]<-as.integer(value)
            }else{
                argv.list[[key]]<-value
            }
        }
    }
    if(sum(names(argv.list) %in% valid_keys)==length(argv.list)){
        print(argv.list)
        temp<-strsplit(argv.list[['--input']],',')
        input<-temp[[1]]

        species<-argv.list[['--species']]
        metadata<-argv.list[['--metadata']]
        interested_variable<-as.character(argv.list[['--interested_variable']])
        interested_motifs<-argv.list[['--interested_motifs']]
        cutoff.matrix.loc<-argv.list[['--cutoff_matrix_loc']]
        outdir<-argv.list[['--outdir']]
        adjusted_p_value_cutoff<-argv.list[['--adjusted_p_value_cutoff']]
        top_n<-argv.list[['--top_n']]
        top_number_of_PCs<-argv.list[['--top_number_of_PCs']]
        perplexity<-argv.list[['--perplexity']]        
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
            stop("Invalid species")
        }
        Differential_TF_activity_analysis(input,
                                          genome,
                                          interested_motifs,
                                          metadata,
                                          interested_variable,
                                          outdir,
                                          adjusted_p_value_cutoff=adjusted_p_value_cutoff,top_n=top_n,
                                          cutoff.matrix.loc=cutoff.matrix.loc,
                                          cutoff.matrix.name=cutoff.matrix.name,
                                          top_number_of_PCs=top_number_of_PCs,
                                          perplexity=perplexity)
    }else{
        stop("Invalid argument")
    }
}
