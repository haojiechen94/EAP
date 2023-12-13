# gene_set_enrichment_analysis.py
# 2022-04-11
# Haojie Chen

"""
gene_set_enrichment_analysis.py --input=path --database=reference_genome_file 
                                [--p_value_cutoff=0.05] [--pathname=<str>]                        

--input=<str>                       The output file from link_peaks_to_genes.py. (e.g. Test_peaks_to_genes_links.txt)

--database=<str>                    Gene set library for enrichment analyis. Each gene set is described by a name, 
                                    a description, and the genes in the gene set. (i.e. File must be in GMT format.)

--ref=<str>                         Path to genome annotation file, for example hg19_refGene.gtf. 
                                    File must be in GTF format.

[--name=<str>]                      Name of output file(s).
                                    Default: Test

[--Type=overlapping_genes]          Using overlapping_genes or proximal_genes or nearest_genes as input.
                                    Default: overlapping_genes

[--p_value_cutoff=0.05]             P value used for identifying significantly over-represented gene sets.
                                    Default: 0.05

[--outdir=<str>]                    Output directory for the processed result.
                                    Default: current directory

[--pathname=<str>]                  Find all pathname matching this pattern. Using relative pathname 
                                    (e.g. /usr/src/*peaks_to_genes_links.txt) to find any matching files 
                                    in a specified directory.

--help/-h                           print this page.

Description: This tool takes peaks associated genes as input (TXT format) and performs Gene Ontology 
             (Biological Processes) and KEGG Pathway Enrichment Analysis of the peaks associated genes.

"""

from sys import argv, stderr, stdin, stdout
from getopt import getopt
import copy
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.stats.multitest
import time
import glob
import pandas as pd
import scipy.stats


def get_genomic_features(path,protein_coding_genes=True):
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
    return all_genes_list

def get_gene_list(path,Type='overlapping_genes'):
    df=pd.read_csv(path,sep='\t')
    gene_list=set()
    for i in df.index:
        if type(df.loc[i,Type])==str:
            for j in df.loc[i,Type].split(','):
                gene_list.add(j)
    return gene_list

def get_gene_sets(path):
    gene_sets_dic={}
    with open(path) as infile:
        for line in infile:
            temp=line.strip().split('\t')
            gene_sets_dic[temp[0]]=set(temp[2:])
    return gene_sets_dic


#+--------------------------+-----------------------------+----------------------------+
#|                          | Genes in interest gene list | Genes not in the gene list |
#+--------------------------+-----------------------------+----------------------------+
#| Genes in the pathway     |            a                |            b               |
#+--------------------------+-----------------------------+----------------------------+
#| Genes not in the pathway |            c                |            d               |
#+--------------------------+-----------------------------+----------------------------+

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

def bar_plot(df,out_dir,name,top_number=10):
    if len(df)<=top_number:
        plt.switch_backend('agg')
        plt.rc('xtick',labelsize=20)
        plt.figure(figsize=(6,len(df)))
        plt.barh([i for i in range(len(df))],
                 [-np.log10(i) for i in df['p_value']],
                 color='salmon',edgecolor='black')
        plt.yticks([i for i in range(len(df))],df['pathway'],size=20)
        plt.xlabel('-Log10(p-value)',size=20)
        plt.ylim(len(df),-1)
    else:
        plt.switch_backend('agg')
        plt.rc('xtick',labelsize=20)
        plt.figure(figsize=(6,top_number))
        plt.barh([i for i in range(top_number)],
                 [-np.log10(i) for i in df.loc[df.index[:top_number],'p_value']],
                 color='salmon',edgecolor='black')
        plt.yticks([i for i in range(len(df))],df.loc[df.index[:top_number],'pathway'],size=20)
        plt.xlabel('-Log10(p-value)',size=20)
        plt.ylim(top_number,-1)        
    plt.savefig('%s/%s_pathway_enrichment.pdf'%(out_dir,name),bbox_inches='tight')
    plt.savefig('%s/%s_pathway_enrichment.png'%(out_dir,name),bbox_inches='tight')
    plt.close('all')

def bubble_heatmap(dic,out_dir,name,p_value_cutoff=0.05,top_number=10):
    ID_list=list(dic.keys())
    selected_pathways=[]
    for ID in ID_list:
        df=dic[ID]
        if sum(df['p_value']<p_value_cutoff)<=top_number:
            for i in df[df['p_value']<p_value_cutoff]['pathway']:
                if i not in selected_pathways:
                    selected_pathways.append(i)
        else:
            for i in df.loc[df.index[:top_number],'pathway']:
                if i not in selected_pathways:
                    selected_pathways.append(i)
    print(selected_pathways)
    data=[]
    for pathway in selected_pathways:
        temp=[]
        for ID in ID_list:
            df=dic[ID]
            p_value=df.loc[df['pathway']==pathway,'p_value'].values[0]
            temp.append(-np.log10(p_value))
        data.append(temp)

    plt.switch_backend('agg')
    plt.rc('xtick',labelsize=20)
    plt.rc('ytick',labelsize=20)
    plt.figure(figsize=(len(ID_list),len(selected_pathways)))
    for y,i in enumerate(data):
        for x,j in enumerate(i):
            plt.scatter(x,y,s=j*50,c='red',edgecolor='black')
    plt.xticks([i for i in range(len(ID_list))],ID_list,rotation=90)
    plt.yticks([i for i in range(len(selected_pathways))],selected_pathways)
    patches=[plt.scatter([],[],marker='o',s=-np.log10(0.01)*50,color='white',edgecolor='black',label='0.01'),
             plt.scatter([],[],marker='o',s=-np.log10(0.001)*50,color='white',edgecolor='black',label='0.001'),
             plt.scatter([],[],marker='o',s=-np.log10(0.0001)*50,color='white',edgecolor='black',label='0.0001'),
             plt.scatter([],[],marker='o',s=-np.log10(0.00001)*50,color='white',edgecolor='black',label='0.00001')]
    plt.legend(bbox_to_anchor=(1,1),handles=patches,fontsize=20,markerscale=1,title='P-value',title_fontsize=20)    
    plt.savefig('%s/%s_pathway_enrichment.pdf'%(out_dir,name),bbox_inches='tight')
    plt.savefig('%s/%s_pathway_enrichment.png'%(out_dir,name),bbox_inches='tight')
    plt.close('all')    

def main():
    path=''
    out_dir=False
    ref_path=''
    database=''
    name='Test'
    Type='overlapping_genes'
    pathname=False
    p_value_cutoff=0.05

    try:
        opts,args=getopt(argv[1:],'h',['input=','database=','ref=','name=','Type=','p_value_cutoff=',
                                       'outdir=','pathname=','help'])
        for i,j in opts:   
            if i=="-h" or i=="--help":
                stdout.write(__doc__)
                exit(0)
            elif i=='--input':
                path=j
            elif i=='--database':
                database=j
            elif i=='--ref':
                ref_path=j
            elif i=='--name':
                name=j
            elif i=='--Type':
                Type=j
            elif i=='--outdir':
                outdir=j
            elif i=='--pathname':
                pathname=j
            elif i=='--p_value_cutoff':
                try:
                    p_value_cutoff=float(j)
                    assert p_value_cutoff>=0 and p_value_cutoff<=1
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

    if pathname:
        print(pathname,database,ref_path,Type,out_dir)        
        all_genes_list=get_genomic_features(ref_path,protein_coding_genes=True)
        gene_sets_dic=get_gene_sets(database)
        dic={}
        for path in glob.glob(pathname):
            ID=path.split('/')[-1].split('_')[0]
            gene_list=get_gene_list(path)
            result=KEGG_enrichment(gene_list,gene_sets_dic,all_genes_list,
                                   p_value_cutoff=1)
            result[result['p_value']<p_value_cutoff].to_csv('%s/%s_pathway_enrichment_analysis.txt'%(out_dir,ID),
                                                            sep='\t',index=False)
            dic[ID]=result

        bubble_heatmap(dic,out_dir,name,p_value_cutoff=0.05,top_number=10)            
        alaysis_finished=time.time()
        print('Time: ',alaysis_finished-start_alaysis)                     
    else:
        print(path,database,ref_path,Type,out_dir)
        all_genes_list=get_genomic_features(ref_path,protein_coding_genes=True)
        gene_list=get_gene_list(path)
        gene_sets_dic=get_gene_sets(database)
        result=KEGG_enrichment(gene_list,gene_sets_dic,all_genes_list,
                               p_value_cutoff=p_value_cutoff)
        bar_plot(result,out_dir,name)
        result.to_csv('%s/%s_pathway_enrichment_analysis.txt'%(out_dir,name),sep='\t',index=False)
        alaysis_finished=time.time()
        print('Time: ',alaysis_finished-start_alaysis)

if __name__ == '__main__':
    main()



