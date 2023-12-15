suppressMessages(
  {
    library(extrafont)
    library(scales)
    }
  )


get_signature_genes_list<-function(signature_genes){
    signature_genes_df<-read.table(signature_genes,sep='\t',stringsAsFactors=F)
    signature_genes_list<-as.character(signature_genes_df[,c(3:ncol(signature_genes_df))])
    signature_genes_set_name<-signature_genes_df[,1]
    return(list(signature_genes_set_name,signature_genes_list))
}

get_genes_to_peaks_links<-function(peaks_to_genes_links){
    proximal_peaks_to_genes_links_df<-read.table(peaks_to_genes_links,sep='\t',stringsAsFactors=F,header=T)
    overlapping_genes<-c()
    chroms<-c()
    starts<-c()
    ends<-c()
    for(i in c(1:dim(proximal_peaks_to_genes_links_df)[1])){
        genes<-proximal_peaks_to_genes_links_df[i,'overlapping_genes']
        if(genes!=''){
            chrom<-proximal_peaks_to_genes_links_df[i,'chrom']
            start<-proximal_peaks_to_genes_links_df[i,'start']
            end<-proximal_peaks_to_genes_links_df[i,'end']
            for(gene in unlist(strsplit(genes,','))){
                overlapping_genes<-c(overlapping_genes,gene)
                chroms<-c(chroms,chrom)
                starts<-c(starts,start)
                ends<-c(ends,end)
            }
        }
    }
    genes_to_bins_links_df<-data.frame(overlapping_gene=overlapping_genes,chrom=chroms,start=starts,end=ends,
                                       stringsAsFactors=F)
    return(genes_to_bins_links_df)
}

link_genes_to_peaks<-function(signature_genes_list,genes_to_bins_links_df,genomic_positions){
    rows<-c()
    for(gene in signature_genes_list){
        for(i in rownames(genes_to_bins_links_df[genes_to_bins_links_df['overlapping_gene']==gene,])){
            chrom<-genes_to_bins_links_df[i,'chrom']
            start<-genes_to_bins_links_df[i,'start']
            end<-genes_to_bins_links_df[i,'end']
            temp<-genomic_positions[c(genomic_positions['chrom']==chrom) & 
                                    c(genomic_positions['start']==start) & 
                                    c(genomic_positions['end']==end),]
            if(dim(temp)[1]!=0){
                rows<-c(rows,rownames(temp))
            }
        }    
    }
    return(rows)
}

Signature_genes_score_estimation<-function(input,peaks_to_genes_links,signature_genes,metadata,interested_variable,outdir){
    load(input)

    num_of_samples<-0
    for(i in colnames(proximal_result)){
        if(grepl('.read_cnt',i)){
            num_of_samples<-num_of_samples+1
        }
    }   

    genomic_positions<-proximal_result[,c(1:3)]
    proximal_zscore_matrix<-scale(t(proximal_result[,c(4:(4+num_of_samples-1))]))

    result<-get_signature_genes_list(signature_genes)
    signature_genes_set_name<-result[[1]]
    signature_genes_list<-result[[2]]

    genes_to_bins_links_df<-get_genes_to_peaks_links(peaks_to_genes_links)

    rows<-link_genes_to_peaks(signature_genes_list,genes_to_bins_links_df,genomic_positions)

    signature_scores<-rowMeans(proximal_zscore_matrix[,rows])-rowMeans(proximal_zscore_matrix)
    meta_data<-read.table(metadata,sep=',',header=T,stringsAsFactors=F)

    signature_scores_df<-cbind(signature_scores,meta_data[interested_variable])
    colnames(signature_scores_df)<-c('signature_scores','interested_variable')

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

    png(gettextf('%s/%s.png',outdir,signature_genes_set_name),family='Arial')
    par(mar=c(5,6,4,2),cex.main=1,font.main=1)
    boxplot(signature_scores~interested_variable,signature_scores_df,
            xlab=interested_variable,ylab='Signature score',
            main=signature_genes_set_name,cex.main=1.5,
            cex.lab=2,cex.axis=1.8,col=color_map_annotation)

    abline(h=0,lty=2,col='red')
    dev.off()

    pdf(gettextf('%s/%s.pdf',outdir,signature_genes_set_name))
    par(mar=c(5,6,4,2),cex.main=1,font.main=1)
    boxplot(signature_scores~interested_variable,signature_scores_df,
            xlab=interested_variable,ylab='Signature score',
            main=signature_genes_set_name,cex.main=1.5,
            cex.lab=2,cex.axis=1.8,col=color_map_annotation)

    abline(h=0,lty=2,col='red')
    dev.off()


}


argv<-commandArgs(TRUE)

help_doc="
Usage: Rscript Signature_genes_score_estimation.R [options]
Options:
    --input=CHARACTER                    Input file name, output from Hypervariable_analysis.R.(e.g. ../Proximal_hypervariable_analysis.txt)
    --peaks_to_genes_links=CHARACTER     Peaks to genes annotation, output from link_peaks_to_genes.py in Data preprocessing module.(e.g. ../step6_reads_counting/proximal_peaks_to_genes_links)
    --signature_genes=CHARACTER          Signature genes list in GMT format.
    --metadata=CHARACTER                 Metadata for each sample.
    --interested_variable=CHARACTER      One of the columns of Metadata.
    --outdir=CHARACTER                   Output directory.
    -h,--help                            Show this help message and exit.
"

valid_keys<-c('--input','--peaks_to_genes_links','--signature_genes','--metadata',
              '--interested_variable','--outdir')

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
    for(i in argv){
        temp<-strsplit(i,'=')[[1]]
        if(length(temp)==2){
            key<-temp[1]
            value<-temp[2]
            argv.list[[key]]<-value
        }
    }
    if(sum(names(argv.list) %in% valid_keys)==length(argv.list)){
        print(argv.list)
        input<-argv.list[['--input']]
        peaks_to_genes_links<-argv.list[['--peaks_to_genes_links']]
        signature_genes<-argv.list[['--signature_genes']]
        metadata<-argv.list[['--metadata']]        
        interested_variable<-argv.list[['--interested_variable']]
        outdir<-argv.list[['--outdir']]

        Signature_genes_score_estimation(input,peaks_to_genes_links,signature_genes,metadata,interested_variable,outdir)


    }else{
        stop("Invalid argument")
    }
}