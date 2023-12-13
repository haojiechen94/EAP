# link_peaks_to_genes.py
# 2022-03-30
# Haojie Chen

"""
link_peaks_to_genes.py --input=path --ref=reference_genome_file 
                        [--name=<str>] [--promoter_distance=2000]
                        [--proximal_distance=2000] 
                        [--outdir=output_directory]
                        [--pathname=<str>]

--input=<str>                   Peaks file in bed format.

--ref=<str>                     Path to genome annotation file, for example hg19_refGene.gtf. 
                                File must be in GTF format.

[--name=<str>]                  Name of output file(s).
                                Default: Test

[--promoter_distance=<int>]     Distance used to define promoter regions.
                                Default: 2000

[--proximal_distance=<int>]     Distance used to identify proximal genes.
                                Default: 50000

[--outdir=<str>]                Output directory for the processed result.
                                Default: current directory

[--pathname=<str>]              Find all pathname matching this pattern. Using relative pathname (e.g. /usr/src/*.bed) 
                                to find any matching files in a specified directory.

--help/-h                       print this page.

Description: This tool takes the coordinates of peaks as input (BED format) and identifies genes assocaited with
             these peaks, such as proximal genes (genes with TSS located within 50kb from the boundary of the 
             peak), overlapping genes (gene's gene body and promoter overlaped with the peak region) and nearest 
             gene (the closest gene among proximal genes).
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

def get_genomic_features(path,promoter_distance,protein_coding_genes=True):
    genomic_features={'TSS':{},'gene':{},'name':{}}
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

    return genomic_features

def get_peaks(path):
    peaks_dic={}
    with open(path) as infile:
        for line in infile:
            temp=line.strip().split('\t')
            chrom,start,end=[temp[0],int(temp[1]),int(temp[2])]
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
                        temp_index.append(index)
                        index+=1                        
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

def link_peaks_to_genes(path,ref_path,promoter_distance,proximal_distance):
    peaks_dic=get_peaks(path)
    genomic_features_dic=get_genomic_features(ref_path,promoter_distance)
    overlapping_genes=get_overlapping_genes(peaks_dic,genomic_features_dic['gene'])
    proximal_genes,nearest_genes=get_proximal_genes_and_nearest_genes(peaks_dic,genomic_features_dic,proximal_distance)
    result=[]
    for i,j,k in zip(overlapping_genes,proximal_genes,nearest_genes):
        result.append(i+[j[-1],k[-1]])
    return result

def link_peaks_to_genes2(pathname,ref_path,promoter_distance,proximal_distance,out_dir):
    genomic_features_dic=get_genomic_features(ref_path,promoter_distance)
    for path in glob.glob(pathname):
        peaks_dic=get_peaks(path)
        ID=path.split('/')[-1].split('_')[0]
        overlapping_genes=get_overlapping_genes(peaks_dic,genomic_features_dic['gene'])
        proximal_genes,nearest_genes=get_proximal_genes_and_nearest_genes(peaks_dic,genomic_features_dic,proximal_distance)
        result=[]
        for i,j,k in zip(overlapping_genes,proximal_genes,nearest_genes):
            result.append(i+[j[-1],k[-1]])
        write_out(result,out_dir,ID)
    return result

def write_out(result,out_dir,name):
    with open('%s/%s_peaks_to_genes_links.txt'%(out_dir,name),'w') as outfile:
        outfile.write('\t'.join(['chrom','start','end','overlapping_genes','proximal_genes','nearest_genes'])+'\n')
        for i in result:
            outfile.write('\t'.join([i[0],str(i[1]),str(i[2]),','.join(i[3]),','.join(i[4]),','.join(i[5])])+'\n')

def main():
    path=''
    out_dir=False
    ref_path=''
    promoter_distance=2000
    proximal_distance=50000
    name='Test'
    pathname=False

    try:
        opts,args=getopt(argv[1:],'h',['input=','name=','outdir=','ref=','promoter_distance=','proximal_distance=','pathname=','help'])
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
            elif i=='--promoter_distance':
                try:
                    promoter_distance=int(j)
                    assert promoter_distance>0
                except (ValueError, AssertionError):
                    raise Exception("Invalid argument for %s: %s"%(i,j))
            elif i=='--proximal_distance':
                try:
                    proximal_distance=int(j)
                    assert proximal_distance>0
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
        link_peaks_to_genes2(pathname,ref_path,promoter_distance,proximal_distance,out_dir)
        alaysis_finished=time.time()
        print('Time: ',alaysis_finished-start_alaysis)
    else:
        result=link_peaks_to_genes(path,ref_path,promoter_distance,proximal_distance)
        write_out(result,out_dir,name)
        alaysis_finished=time.time()
        print('Time: ',alaysis_finished-start_alaysis)

if __name__ == '__main__':
    main()



