# functional_annotation_of_differential_peaks.py
# 2022-04-20
# Haojie Chen

"""
functional_annotation_of_differential_peaks.py --input_dir=directory --annotation_file=path_to_file
                                               --gene_sets=path_to_functional_annotation_library
                                               [--adjusted_p_value_cutoff=0.1] [--p_value_cutoff=0.05] 
                                               [--outdir=current_directory]


--input_dir=<str>                   The direcoty of output files from differential analysis. 
                                    (e.g. /usr/src/step7_differential_analysis/)

--annotation_file=<str>             Path to genome annotation file, for example hg19_refGene.gtf. 
                                    File must be in GTF format.

--gene_sets=<str>                   Gene set library for enrichment analyis. Each gene set is described 
                                    by a name, a description, and the genes in the gene set. 
                                    (i.e. File must be in GMT format.)

[--adjusted_p_value_cutoff=0.1]     Adjusted p value used for identifing significant differential peaks.
                                    Default: 0.1

[--p_value_cutoff=0.05]             P value used for identifying significantly over-represented gene sets.
                                    Default: 0.05

[--Type=overlapping_genes]          Using overlapping_genes or proximal_genes or nearest_genes as input.
                                    Default: overlapping_genes                                    

[--outdir=<str>]                    Output directory for the processed result.
                                    Default: current directory

--help/-h                           print this page.
"""

import glob
import pandas as pd
from sys import argv, stderr, stdin, stdout
from getopt import getopt
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import statsmodels.stats.multitest
import scipy.stats

def get_genomic_features(path,promoter_distance,protein_coding_genes=True):
    genomic_features={'TSS':{},'gene':{},'name':{}}
    all_genes_list=set()
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
                all_genes_list.add(dic['gene_name'])
                if chrom in genomic_features['gene']:
                    if strand=='+':
                        genomic_features['gene'][chrom].append([start-promoter_distance,end,dic['gene_name']])
                        genomic_features['TSS'][chrom].append([start,dic['gene_name']])
                    else:
                        genomic_features['gene'][chrom].append([start,end+promoter_distance,dic['gene_name']])
                        genomic_features['TSS'][chrom].append([end,dic['gene_name']])
                else:
                    if strand=='+':
                        genomic_features['gene'][chrom]=[[start-promoter_distance,end,dic['gene_name']]]
                        genomic_features['TSS'][chrom]=[[start,dic['gene_name']]]
                    else:
                        genomic_features['gene'][chrom]=[[start,end+promoter_distance,dic['gene_name']]]
                        genomic_features['TSS'][chrom]=[[end,dic['gene_name']]]

    for Type in genomic_features:
        for chrom in genomic_features[Type]:
            genomic_features[Type][chrom]=sorted(genomic_features[Type][chrom],key=lambda x:x[0])
    
    genomic_features['name']={}
    for chrom in genomic_features['TSS']:
        genomic_features['name'][chrom]=[i[1] for i in genomic_features['TSS'][chrom]]
        genomic_features['TSS'][chrom]=[i[0] for i in genomic_features['TSS'][chrom]]

    return genomic_features,all_genes_list

def get_peaks(path,adjusted_p_value_cutoff=0.1):
    differential_analysis_result=pd.read_csv(path,sep='\t')
    temp=path.split('/')[-1].split('.')[0].split('_')
    condition1=temp[2]
    condition2=temp[4]
    peaks_dic={condition1:{},condition2:{}}
    for i in differential_analysis_result[differential_analysis_result['padj']<adjusted_p_value_cutoff].index:
        chrom,start,end,Mval=differential_analysis_result.loc[i,['chrom','start','end','Mval']]
        if Mval<0:         
            if chrom in peaks_dic[condition1]:
                peaks_dic[condition1][chrom].append([start,end])
            else:
                peaks_dic[condition1][chrom]=[[start,end]]
        else:
            if chrom in peaks_dic[condition2]:
                peaks_dic[condition2][chrom].append([start,end])
            else:
                peaks_dic[condition2][chrom]=[[start,end]]            

    for condition in peaks_dic:
        for chrom in peaks_dic[condition]:
            peaks_dic[condition][chrom]=sorted(peaks_dic[condition][chrom],key=lambda x:x[0])
    return peaks_dic

def overlap(s1,e1,s2,e2):
    if min(e1,e2)-max(s1,s2)>=0:
        return 1
    else:
        return -1

def get_overlapping_genes(peaks_dic,genomic_features_dic):
    result=[]
    for chrom in peaks_dic:
        index=0
        for peak_start,peak_end in peaks_dic[chrom]:
            temp=[]
            temp_index=[]
            if chrom in genomic_features_dic:
                while index<len(genomic_features_dic[chrom]):
                    start,end,name=genomic_features_dic[chrom][index]
                    if overlap(peak_start,peak_end,start,end)==1:
                        temp.append(name)
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

    return result

def bisect(sorted_arr,target,right=True):
    start=0
    end=len(sorted_arr)
    if right:
        while end>start:
            mid=(start+end)//2
            if target<sorted_arr[mid]:
                end=mid
            else:
                start=mid+1
    else:
        while end>start:
            mid=(start+end)//2
            if target>sorted_arr[mid]:
                start=mid+1
            else:
                end=mid
    return start

def get_proximal_genes_and_nearest_genes(peaks_dic,genomic_features_dic,proximal_distance):
    proximal_result=[]
    nearest_result=[]
    for chrom in peaks_dic:
        for peak_start,peak_end in peaks_dic[chrom]:
            temp1=[]
            temp2=[]
            if chrom in genomic_features_dic['TSS']:
                left=bisect(genomic_features_dic['TSS'][chrom],peak_start-proximal_distance,right=False)
                right=bisect(genomic_features_dic['TSS'][chrom],peak_end+proximal_distance,right=True)
                if right-left>0:
                    temp1=genomic_features_dic['name'][chrom][left:right]
                    mid=(peak_start+peak_end)*1.0/2
                    distances=[abs(mid-genomic_features_dic['TSS'][chrom][i]) for i in range(left,right)]
                    min_dis=min(distances)
                    index=distances.index(min_dis)
                    temp2=[temp1[index]]
            proximal_result.append([chrom,peak_start,peak_end,list(set(temp1))])
            nearest_result.append([chrom,peak_start,peak_end,temp2])

    return proximal_result,nearest_result

def link_peaks_to_genes(input_dir,ref_path,promoter_distance,proximal_distance):
    genomic_features_dic,all_genes_list=get_genomic_features(ref_path,promoter_distance)
    result_dic={}
    for path in glob.glob(input_dir+'/Differential_analysis*.txt'):
        peaks_dic=get_peaks(path)
        analysis_ID=path.split('/')[-1].split('.')[0].split('Differential_analysis_')[1]
        result_dic[analysis_ID]={}
        for condition in peaks_dic:
            overlapping_genes=get_overlapping_genes(peaks_dic[condition],genomic_features_dic['gene'])
            proximal_genes,nearest_genes=get_proximal_genes_and_nearest_genes(peaks_dic[condition],genomic_features_dic,proximal_distance)
            result=[]
            for i,j,k in zip(overlapping_genes,proximal_genes,nearest_genes):
                result.append(i[:-1]+[','.join(i[-1]),','.join(j[-1]),','.join(k[-1])])
            temp_df=pd.DataFrame(result,columns=['chrom','start','end','overlapping_genes','proximal_genes','nearest_genes'])
            result_dic[analysis_ID][condition]=temp_df
    return result_dic,all_genes_list

def get_gene_sets(path):
    gene_sets_dic={}
    with open(path) as infile:
        for line in infile:
            temp=line.strip().split('\t')
            gene_sets_dic[temp[0]]=set(temp[2:])
    return gene_sets_dic

def KEGG_enrichment(gene_list,gene_sets_dic,all_genes_list,p_value_cutoff=0.05):
    result=[]
    for pathway in gene_sets_dic:
        a=gene_list&gene_sets_dic[pathway]
        b=all_genes_list&gene_sets_dic[pathway]
        c=gene_list-gene_sets_dic[pathway]
        d=all_genes_list-gene_sets_dic[pathway]
        odd_ratio,p_value=scipy.stats.fisher_exact([[len(a),len(b)],
                                                    [len(c),len(d)]],
                                                   alternative='greater')
        result.append([pathway,a,odd_ratio,p_value])

    result1=[]
    for i,j in zip(result,statsmodels.stats.multitest.multipletests([i[-1] for i in result],method='fdr_bh')[1]):
        result1.append(i+[j])

    df=pd.DataFrame(result1,columns=['pathway','overlapping_genes','odd_ratio','p_value','p_adjust'])
    return df[df['p_value']<p_value_cutoff].sort_values(by='p_value')

def bar_plot(df1,df2,condition1,condition2,out_dir,name,top_number=10):
    plt.switch_backend('agg')
    plt.rc('xtick',labelsize=20)
    plt.figure(figsize=(6,14))
    gs=mpl.gridspec.GridSpec(12,6)
    ax=plt.subplot(gs[0:5,:])    
    if len(df1)<=top_number:
        plt.barh([i for i in range(len(df1))],
                 [-np.log10(i) for i in df1['p_value']],
                 color='salmon',edgecolor='black')
        plt.yticks([i for i in range(len(df1))],df1['pathway'],size=20)
        plt.xlabel('-Log10(p-value)',size=20)
        plt.ylim(10,-1)
        plt.title(condition1,size=20)
    else:
        plt.barh([i for i in range(top_number)],
                 [-np.log10(i) for i in df1.loc[df1.index[:top_number],'p_value']],
                 color='salmon',edgecolor='black')
        plt.yticks([i for i in range(top_number)],df1.loc[df1.index[:top_number],'pathway'],size=20)
        plt.xlabel('-Log10(p-value)',size=20)
        plt.ylim(top_number,-1)
        plt.title(condition1,size=20)

    ax=plt.subplot(gs[6:11,:])  
    if len(df2)<=top_number:
        plt.barh([i for i in range(len(df2))],
                 [-np.log10(i) for i in df2['p_value']],
                 color='yellowgreen',edgecolor='black')
        plt.yticks([i for i in range(len(df2))],df2['pathway'],size=20)
        plt.xlabel('-Log10(p-value)',size=20)
        plt.ylim(10,-1)
        plt.title(condition2,size=20)
    else:
        plt.barh([i for i in range(top_number)],
                 [-np.log10(i) for i in df2.loc[df2.index[:top_number],'p_value']],
                 color='yellowgreen',edgecolor='black')
        plt.yticks([i for i in range(top_number)],df2.loc[df2.index[:top_number],'pathway'],size=20)
        plt.xlabel('-Log10(p-value)',size=20)
        plt.ylim(top_number,-1)
        plt.title(condition2,size=20)

    plt.savefig('%s/%s_pathway_enrichment.pdf'%(out_dir,name),bbox_inches='tight')
    plt.savefig('%s/%s_pathway_enrichment.png'%(out_dir,name),bbox_inches='tight')
    plt.close('all')

def main():
    Type='overlapping_genes'
    p_value_cutoff=0.05
    adjusted_p_value_cutoff=0.1
    promoter_distance=2000
    proximal_distance=50000

    try:
        opts,args=getopt(argv[1:],'h',['input_dir=','annotation_file=','gene_sets=','adjusted_p_value_cutoff=','p_value_cutoff=','Type=','outdir=','help'])
        for i,j in opts:   
            if i=="-h" or i=="--help":
                stdout.write(__doc__)
                exit(0)
            elif i=='--input_dir':
                input_dir=j
            elif i=='--annotation_file':
                annotation_file=j
            elif i=='--gene_sets':
                gene_sets=j
            elif i=='--Type':
                Type=j
            elif i=='--outdir':
                out_dir=j
            elif i=='--p_value_cutoff':
                try:
                    p_value_cutoff=float(j)
                    assert p_value_cutoff>=0 and p_value_cutoff<=1
                except (ValueError, AssertionError):
                    raise Exception("Invalid argument for %s: %s"%(i,j))
            elif i=='--adjusted_p_value_cutoff':
                try:
                    adjusted_p_value_cutoff=float(j)
                    assert adjusted_p_value_cutoff>=0 and adjusted_p_value_cutoff<=1
                except (ValueError, AssertionError):
                    raise Exception("Invalid argument for %s: %s"%(i,j))                    
            else:
                raise Exception("Internal errors occur when parsing command line arguments.")
    except Exception as e:
        stderr.write("%s\n" % e)
        stderr.write("Type 'python gene_set_enrichment_analysis.py --help' for more information.\n")
        exit(1)

    if not out_dir:
        out_dir=os.getcwd()
    start_alaysis=time.time()    

    result_dic,all_genes_list=link_peaks_to_genes(input_dir,annotation_file,promoter_distance,proximal_distance)
    gene_sets_dic=get_gene_sets(gene_sets)
    for analysis_ID in result_dic:
        enrichment_result_list=[]
        conditions=[]
        for condition in result_dic[analysis_ID]:
            conditions.append(condition)
            gene_list=set()
            result_dic[analysis_ID][condition].to_csv('%s/%s_%s_specific_peaks_associated_genes.txt'%(out_dir,analysis_ID,condition),sep='\t',index=False)
            for i in result_dic[analysis_ID][condition][Type]:
                if type(i)==str:
                    for j in i.split(','):
                        gene_list.add(j)
            enrichment_result_df=KEGG_enrichment(gene_list,gene_sets_dic,all_genes_list,p_value_cutoff=p_value_cutoff)
            enrichment_result_df.to_csv('%s/%s_functional_enrichment_based_on_%s_specific_peaks_associated_genes.txt'%(out_dir,analysis_ID,condition),sep='\t',index=False)
            enrichment_result_list.append(enrichment_result_df)
        bar_plot(enrichment_result_list[0],enrichment_result_list[1],conditions[0],conditions[1],out_dir,analysis_ID,top_number=10)

if __name__ == '__main__':
    main()
