# peaks_annotation.py
# 2022-03-29
# Haojie Chen

"""
peaks_annotation.py --input=path --ref=reference_genome_file 
                        [--name=<str>] [--distance=2000] 
                        [--outdir=output_directory]
                        [--pathname=<str>]

--input=<str>          Peaks file in bed format.

--ref=<str>            Path to genome annotation file, for example hg19_refGene.gtf. 
                       File must be in GTF format. The priority in genomic annotation:
                       promoter > exon > intron > intergenic.

[--name=<str>]         Name of output file(s).
                       Default: Test

[--distance=<int>]     Distance used to define promoter regions.
                       Default: 2000

[--outdir=<str>]       Output directory for the processed result.
                       Default: current directory

[--pathname=<str>]     Find all pathname matching this pattern. Using relative pathname (e.g. /usr/src/*.bed) 
                       to find any matching files in a specified directory.

--help/-h              print this page.                       

"""

from sys import argv, stderr, stdin, stdout
from getopt import getopt
import copy
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import glob

def get_genomic_features(path,distance,protein_coding_genes=True):
    genomic_features={'promoter':{},'exon':{},'gene_body':{},'5UTR':{},'3UTR':{}}
    with open(path) as infile:
        for line in infile:
            temp=line.strip().split('\t')
            chrom,start,end,strand,Type=[temp[0],int(temp[3]),int(temp[4]),temp[6],temp[2]]
            dic={}
            for i in temp[8].strip(';').split(';'):
                k,v=i.strip().split(' ')
                v=v.replace('"','')
                dic[k]=v
            if protein_coding_genes:
                if 'NR' in dic['transcript_id']:
                    continue
            if Type=='transcript':
                if chrom in genomic_features['promoter']:
                    genomic_features['gene_body'][chrom].append([start,end,dic['gene_name']])
                    if strand=='+':
                        genomic_features['promoter'][chrom].append([start-distance,start+distance,dic['gene_name']])
                    else:
                        genomic_features['promoter'][chrom].append([end-distance,end+distance,dic['gene_name']])
                else:
                    genomic_features['gene_body'][chrom]=[[start,end,dic['gene_name']]]
                    if strand=='+':
                        genomic_features['promoter'][chrom]=[[start-distance,start+distance,dic['gene_name']]]
                    else:
                        genomic_features['promoter'][chrom]=[[end-distance,end+distance,dic['gene_name']]]
            elif Type=='exon':
                if chrom in genomic_features['exon']:
                    genomic_features['exon'][chrom].append([start,end,dic['gene_name']+':'+dic['exon_number']])
                else:
                    genomic_features['exon'][chrom]=[[start,end,dic['gene_name']+':'+dic['exon_number']]]
            elif Type=='5UTR':
                if chrom in genomic_features['5UTR']:
                    genomic_features['5UTR'][chrom].append([start,end,dic['gene_name']+':'+dic['exon_number']])
                else:
                    genomic_features['5UTR'][chrom]=[[start,end,dic['gene_name']+':'+dic['exon_number']]]
            elif Type=='3UTR':
                if chrom in genomic_features['3UTR']:
                    genomic_features['3UTR'][chrom].append([start,end,dic['gene_name']+':'+dic['exon_number']])
                else:
                    genomic_features['3UTR'][chrom]=[[start,end,dic['gene_name']+':'+dic['exon_number']]]                          

    for Type in genomic_features:
        for chrom in genomic_features[Type]:
            genomic_features[Type][chrom]=sorted(genomic_features[Type][chrom],key=lambda x:x[0])

    return genomic_features


def get_peaks(path):
    peaks_dic={}
    with open(path) as infile:
        for line in infile:
            temp=line.strip().split('\t')
            chrom,start,end,peak_id=[temp[0],int(temp[1]),int(temp[2]),temp[3]]
            if chrom in peaks_dic:
                peaks_dic[chrom].append([start,end])
            else:
                peaks_dic[chrom]=[[start,end]]
    for chrom in peaks_dic:
        peaks_dic[chrom]=sorted(peaks_dic[chrom],key=lambda x:x[0])
    return peaks_dic


def overlap(s1,e1,s2,e2):
    if min(e1,e2)-max(s1,s2)>=0:
        return 1
    else:
        return -1


def annotate_it(peaks_dic,genomic_features_dic):
    result=[]
    other_peaks={}
    for chrom in peaks_dic:
        index=0
        for peak_start,peak_end in peaks_dic[chrom]:
            temp=[]
            temp_index=[]
            if chrom in genomic_features_dic:
                while index<len(genomic_features_dic[chrom]):
                    start,end,name=genomic_features_dic[chrom][index]
                    if overlap(peak_start,peak_end,start,end)==1:
                        temp.append(':'.join([name,chrom,str(start),str(end)]))
                        index+=1
                        temp_index.append(index)
                    else:
                        if peak_start>end:
                            index+=1
                        elif peak_end<start:
                            break
            if len(temp)>0:
                index=temp_index[0]
                result.append([chrom,peak_start,peak_end,list(set(temp))])
            else:
                if chrom in other_peaks:
                    other_peaks[chrom].append([peak_start,peak_end])
                else:
                    other_peaks[chrom]=[[peak_start,peak_end]]

    return result,other_peaks


def annotate_peaks(peaks_dic,genomic_features_dic):
    promoter_result,other_peaks=annotate_it(peaks_dic,genomic_features_dic['promoter'])
    exon_result,other_peaks=annotate_it(other_peaks,genomic_features_dic['exon'])
    intron_result,other_peaks=annotate_it(other_peaks,genomic_features_dic['gene_body'])
    intergenic_result=[]
    for chrom in other_peaks:
        for peak_start,peak_end in other_peaks[chrom]:
            intergenic_result.append([chrom,peak_start,peak_end,[]])
    annotate_peaks_result={
        'promoter':promoter_result,
        'exon':exon_result,
        'intron':intron_result,
        'intergenic':intergenic_result,
    }
    return annotate_peaks_result


def write_out(result,path):
    with open(path,'w') as outfile:
        for i in result:
            outfile.write('\t'.join([i[0],str(i[1]),str(i[2]),','.join(i[3])])+'\n')


def write_out_genomic_features(genomic_features_dic,out_dir):
    for Type in genomic_features_dic:
        with open('%s/%s.bed'%(out_dir,Type),'w') as outfile:
            for chrom in genomic_features_dic[Type]:
                for start,end,name in genomic_features_dic[Type][chrom]:
                    outfile.write('\t'.join([chrom,str(start),str(end),name])+'\n')


def write_out_distributions(distributions,labels,out_dir,name):
    with open('%s/%s_distributions.txt'%(out_dir,name),'w') as outfile:
        outfile.write('\t'.join(['Sample_ID']+['promoter','exon','intron','intergenic'])+'\n')
        for i,j in zip(labels,distributions):
            outfile.write('\t'.join([i]+[str(k) for k in j])+'\n')


def plot_it(distributions,out_dir,name,labels=['Test']):
    plt.switch_backend('agg')
    plt.figure(figsize=(6,len(distributions)))
    for i,j in enumerate(distributions):
        plt.barh([i]*len(j),
                 [k*1.0/sum(j) for k in j],
                 left=[sum(j[:k])*1.0/sum(j) for k in range(len(j))],
                 color=['salmon','yellowgreen','skyblue','orange'],edgecolor='black')
    plt.xlabel('Fraction',size=20)
    plt.xticks([0.0,0.2,0.4,0.6,0.8,1.0],size=20)
    plt.yticks([i for i in range(len(labels))],labels,size=20)
    patches=[plt.scatter([],[],marker='s',color='salmon',edgecolor='black',label='promoter',linewidths=1),
             plt.scatter([],[],marker='s',color='yellowgreen',edgecolor='black',label='exon',linewidths=1),
             plt.scatter([],[],marker='s',color='skyblue',edgecolor='black',label='intron',linewidths=1),
             plt.scatter([],[],marker='s',color='orange',edgecolor='black',label='intergenic',linewidths=1)]
    plt.legend(bbox_to_anchor=(1,1),handles=patches,fontsize=20,markerscale=3)
    plt.savefig('%s/%s_distribution_in_genome.pdf'%(out_dir,name),bbox_inches='tight')
    plt.savefig('%s/%s_distribution_in_genome.png'%(out_dir,name),bbox_inches='tight')
    plt.close('all')


def main():
    path=''
    out_dir=False
    ref_path=''
    distance=2000
    name='Test'
    pathname=False

    try:
        opts,args=getopt(argv[1:],'h',['input=','name=','outdir=','ref=','distance=','pathname=','help'])
        for i,j in opts:   
            if i=="-h" or i=="--help":
                stdout.write(__doc__)
                exit(0)
            elif i=='--input':
                path=j
            elif i=='--outdir':
                out_dir=j
            elif i=='--name':
                name=j
            elif i=='--ref':
                ref_path=j
            elif i=='--pathname':
                pathname=j
            elif i=='--distance':
                try:
                    distance=int(j)
                    assert distance>0
                except (ValueError, AssertionError):
                    raise Exception("Invalid argument for %s: %s"%(i,j))
            else:
                raise Exception("Internal errors occur when parsing command line arguments.")
    except Exception as e:
        stderr.write("%s\n" % e)
        stderr.write("Type 'python bed_coverage_heatmap.py --help' for more information.\n")
        exit(1)

    if not out_dir:
        out_dir=os.getcwd()

    start_alaysis=time.time()

    if pathname:
        print(pathname,out_dir,ref_path,distance)
        genomic_features_dic=get_genomic_features(ref_path,distance)
        labels=[]
        distributions=[]
        for path in glob.glob(pathname):
            peaks_dic=get_peaks(path)
            ID=path.split('/')[-1].split('_peaks.bed')[0]
            labels.append(ID)
            result=annotate_peaks(peaks_dic,genomic_features_dic)
            distribution=[]
            print(ID)
            for Type in ['promoter','exon','intron','intergenic']:
                print(Type,': ',len(result[Type]))
                distribution.append(len(result[Type]))
                write_out(result[Type],'%s/%s_%s_annotations.bed'%(out_dir,ID,Type)) 
            distributions.append(distribution)
        plot_it(distributions,out_dir,name,labels=labels)
        write_out_distributions(distributions,labels,out_dir,name)
        alaysis_finished=time.time()
        print('Time: ',alaysis_finished-start_alaysis)
    else:
        print(path,out_dir,ref_path,distance)
        peaks_dic=get_peaks(path)
        genomic_features_dic=get_genomic_features(ref_path,distance)
        result=annotate_peaks(peaks_dic,genomic_features_dic)
        distribution=[]
        for Type in ['promoter','exon','intron','intergenic']:
            print(Type,': ',len(result[Type]))
            distribution.append(len(result[Type]))
            write_out(result[Type],'%s/%s_%s_annotations.bed'%(out_dir,name,Type))
        distributions=[distribution]
        labels=[name]
        write_out_distributions(distributions,labels,out_dir,name)
        plot_it(distributions,out_dir,name)
        alaysis_finished=time.time()
        print('Time: ',alaysis_finished-start_alaysis)


if __name__ == '__main__':
    main()





