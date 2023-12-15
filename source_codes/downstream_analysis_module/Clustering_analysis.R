suppressMessages(
  {
    library(extrafont)
    library(scales)
    library(dendextend)
    library(pcaMethods)
    library(ComplexHeatmap)
    }
  )


calculate_cutoff<-function(R2,outdir){
    PCs<-unlist(lapply(c(length(R2):1),function(x){gettextf('PC%s',x)}))
    a_list<-as.vector(R2)
    a_list<-sort(a_list)
    slope<-(max(a_list)-min(a_list))/length(a_list)
    rank<-floor(optimize(number_of_points_below_line,lower=1,upper=length(a_list),
                         my_list=a_list,slope=slope)$minimum)
    cutoff<-a_list[rank]
    png(gettextf('%s/Explained_variance_ratio_in_each_PC.png',outdir),family='Arial')
    par(mar=c(5,6,4,2),cex.main=1,font.main=1)
    plot(1:length(a_list),a_list,type="l",lwd=3,
         ylab='Explained variance ratio',xlab='Principal components',
         cex.lab=1.5,xaxt='n',cex.axis=1.2,
         main=gettextf('%s: %.4f',PCs[rank],a_list[rank]))
    intercept<-cutoff-(slope*rank)
    abline(v=rank,h=cutoff,lty=2,col=8)
    points(rank,cutoff,pch=16,cex=0.9,col=2)
    abline(coef=c(intercept,slope),col=2)
    axis(side=1,
         at=c(1:length(a_list)),
         labels=PCs,cex.axis=1.2)  
    dev.off()

    pdf(gettextf('%s/Explained_variance_ratio_in_each_PC.pdf',outdir))
    par(mar=c(5,6,4,2),cex.main=1,font.main=1)
    plot(1:length(a_list),a_list,type="l",lwd=3,
         ylab='Explained variance ratio',xlab='Principal components',
         cex.lab=1.5,xaxt='n',cex.axis=1.2,
         main=gettextf('%s: %.4f',PCs[rank],a_list[rank]))
    intercept<-cutoff-(slope*rank)
    abline(v=rank,h=cutoff,lty=2,col=8)
    points(rank,cutoff,pch=16,cex=0.9,col=2)
    abline(coef=c(intercept,slope),col=2)
    axis(side=1,
         at=c(1:length(a_list)),
         labels=PCs,cex.axis=1.2)  
    dev.off()
    print('the best number of top ranked PCs:')
    print(c(length(R2):1)[rank])
    return(c(length(R2):1)[rank])

}


number_of_points_below_line<-function(my_list,slope,x){
    value<-my_list[x]
    intercept<-value-(slope*x)
    ranks<- 1:length(my_list)
    return(sum(my_list<=(ranks*slope+intercept)))
}


Clustering_analysis<-function(input,metadata,
                              number_of_cluster,interested_variable,outdir,
                              name='Test',adjusted_p_value_cutoff=0.1,top_number_of_PCs=1){

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

    zscores_matrix<-round(zscore_matrix,5)

    meta_data<-read.table(metadata,sep=',',header=T)
    color_by<-meta_data[[interested_variable]]
    color_by<-as.character(color_by)
    unique_labels<-unique(color_by)
    color_list<-hue_pal()(length(unique_labels))
    names<-c()
    values<-c()
    for(i in c(1:length(unique_labels))){
        names<-c(names,unique_labels[i])
        values<-c(values,color_list[[i]])
    }
    color_map_annotation<-setNames(values,names)
    
    color_list<-hue_pal()(number_of_cluster)
    pcs<-pca(zscores_matrix,nPcs=ceiling(dim(zscores_matrix)[1]*3/4),
        method="nipals",cv="q2"
        )
    best_number_of_top_ranked_PCs<-calculate_cutoff(pcs@R2,outdir)
    if(top_number_of_PCs==0){
        top_number_of_PCs<-best_number_of_top_ranked_PCs
    }
    print(top_number_of_PCs)
    clusters<-hclust(dist(round(scores(pcs)[,c(1:top_number_of_PCs)],3)),method="ward.D")
    #save(zscores_matrix,file=gettextf('%s/zscores_matrix.RData',outdir))
    clusters<-color_branches(clusters,k=number_of_cluster,col=color_list)

    meta_data[['cluster']]<-as.character(lapply(cutree(clusters,k=number_of_cluster),function(x){return(paste0('C',x))}))
    write.table(meta_data,gettextf('%s/%s_clustering_result.txt',outdir,name),sep=',',quote=F,row.names=F)

    #color_by<-meta_data[['cluster']]
    #PCA_analysis(zscores_matrix,color_by,outdir)

    color_by<-meta_data[['cluster']]
    color_by<-as.character(color_by)
    unique_labels<-unique(color_by)
    color_list<-hue_pal()(length(unique_labels))
    names<-c()
    values<-c()
    for(i in c(1:length(unique_labels))){
        names<-c(names,unique_labels[i])
        values<-c(values,color_list[[i]])
    }
    color_map_cluster<-setNames(values,names)

    heatmap_annotation<-HeatmapAnnotation(Cluster=meta_data[['cluster']],
                                          Annotation=meta_data[[interested_variable]],
                                          col=list(Cluster=color_map_cluster,
                                                   Annotation=color_map_annotation))

    png(gettextf('%s/%s_clustering_plot.png',outdir,name),family='Arial')
    draw(Heatmap(t(zscores_matrix), 
                 cluster_rows=T,
                 cluster_columns=clusters,
                 heatmap_legend_param=list(title='value'),
                 top_annotation=heatmap_annotation,
                 show_row_names=F))
    dev.off()

    pdf(gettextf('%s/%s_clustering_plot.pdf',outdir,name))
    draw(Heatmap(t(zscores_matrix), 
                 cluster_rows=T,
                 cluster_columns=clusters,
                 heatmap_legend_param=list(title='value'),
                 top_annotation=heatmap_annotation,
                 show_row_names=F))
    dev.off()


}

argv<-commandArgs(TRUE)

help_doc="
Usage: Rscript Clustering_analysis.R [options]
Options:
    --input=CHARACTER                    Input file names, output from Hypervariable_analysis.R.
                                         Files separated by comma(e.g. ../Proximal_hypervariable_analysis.txt,../Distal_hypervariable_analysis.txt)
    --metadata=CHARACTER                 Metadata for each sample.
    --top_number_of_PCs=INTEGER          Top number of principal components used for samples clustering [Default: 0, automaticly choose the best number of top ranked PCs].    
    --number_of_cluster=INTEGER          Classified these samples into k cluster [Default: 2].
    --interested_variable=CHARACTER      One of the columns of Metadata.
    --outdir=CHARACTER                   Output directory.
    --name=CHARACTER                     Name of output file (e.g. {output_file_name}_clustering_reuslt.txt, {output_file_name}_clustering_plot.png) [Default: Test].
    --adjusted_p_value_cutoff=DOUBLE     Adjusted p value used for identifing significant HVRs [Default: 0.1].
    -h,--help                            Show this help message and exit.
"

valid_keys<-c('--input','--metadata','--number_of_cluster',
              '--outdir','--name','--adjusted_p_value_cutoff',
              '--interested_variable','--top_number_of_PCs')

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
    argv.list[['--adjusted_p_value_cutoff']]<-0.1
    argv.list[['--number_of_cluster']]<-2
    argv.list[['--top_number_of_PCs']]<-0
    argv.list[['--name']]<-'Test'
    for(i in argv){
        temp<-strsplit(i,'=')[[1]]
        if(length(temp)==2){
            key<-temp[1]
            value<-temp[2]
            if(key=='--adjusted_p_value_cutoff' | key=='--number_of_cluster' | key=='--top_number_of_PCs'){
                argv.list[[key]]<-as.numeric(value)
            }else{
                argv.list[[key]]<-value
            }
        }
    }
    if(sum(names(argv.list) %in% valid_keys)==length(argv.list)){
        print(argv.list)
        temp<-strsplit(argv.list[['--input']],',')
        input<-temp[[1]]
        metadata<-argv.list[['--metadata']]
        top_number_of_PCs<-argv.list[['--top_number_of_PCs']]
        number_of_cluster<-argv.list[['--number_of_cluster']]
        outdir<-argv.list[['--outdir']]
        name<-argv.list[['--name']]
        adjusted_p_value_cutoff<-argv.list[['--adjusted_p_value_cutoff']]
        interested_variable<-argv.list[['--interested_variable']]
        Clustering_analysis(input,
                            metadata,
                            number_of_cluster,
                            interested_variable,
                            outdir,
                            name=name,
                            adjusted_p_value_cutoff=adjusted_p_value_cutoff,
                            top_number_of_PCs=top_number_of_PCs)
    }else{
        stop("Invalid argument")
    }
}
