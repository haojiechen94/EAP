# bed_coverage_heatmap.py
# 2022-03-15
# Haojie Chen

"""
bed_coverage_heatmap.py --input=path --ref=reference_genome_file 
                        [--shift_size=100] [--distance=2000] 
                        [--outdir=output_directory]
                        [--matrix=<boolean>]

--input=<str>          Input mapped reads file in bed format.

--ref=<str>            Path to genome annotation file, for example hg19_refGene.gtf. 
                       File must be in GTF format.

[--name=<str>]         Name of output file(s).
                       Default: Test

[--shift_size=<int>]   Distance used to shift the reads toward binding site.
                       Default: 100


[--distance=<int>]     Distance used to extend from TSS and TES.
                       Default: 2000

[--outdir=<str>]       Output directory for the processed result.
                       Default: current directory

[--matrix]             Output the reads densitiy matrix.
                       Default: OFF

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


def get_genomic_regions(path,distance):
    genomic_regions_dic={}
    with open(path) as infile:
        for line in infile:
            temp=line.strip().split('\t')
            if temp[2]=='transcript':
                chrom,strand,TSS,TES=[temp[0],temp[6],temp[3],temp[4]]
                if chrom not in genomic_regions_dic:
                    genomic_regions_dic[chrom]=[]
                if strand=='+':
                    genomic_regions_dic[chrom].append([int(TSS),int(TES)])
                else:
                    genomic_regions_dic[chrom].append([int(TES),int(TSS)])
    return genomic_regions_dic

def get_mapped_reads(path,shift_size):
    reads_dic={}
    total_number_of_reads=0
    with open(path) as infile:
        for line in infile:
            total_number_of_reads+=1
            temp=line.strip().split('\t')
            chrom,start,end,name,score,strand=temp[:6]
            start=int(start)
            end=int(end)
            if chrom in reads_dic:
                if strand=='+':
                    reads_dic[chrom].append(start+shift_size)
                else:
                    reads_dic[chrom].append(end-shift_size)
            else:
                if strand=='+':
                    reads_dic[chrom]=[start+shift_size]
                else:
                    reads_dic[chrom]=[end-shift_size]
    for chrom in reads_dic:
        reads_dic[chrom]=sorted(reads_dic[chrom])
    return reads_dic,total_number_of_reads

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

def couting_reads(path,ref_path,distance,shift_size):
    print('Reading genes annotation gtf file...')
    genomic_regions_dic=get_genomic_regions(ref_path,distance)
    print('Reading mapped reads bed file...')
    reads_dic,total_number_of_reads=get_mapped_reads(path,shift_size)
    reads_densities=[]
    for chrom in genomic_regions_dic:
        for start,end in genomic_regions_dic[chrom]:
            reads_densities_gene=[]
            if end>start:
                temp=list(np.linspace(start-distance,start,51))[:-1]+list(np.linspace(start,end,101))+list(np.linspace(end,end+distance,51))[1:]
                for i,j in zip(temp[:-1],temp[1:]):
                    num=bisect(reads_dic[chrom],j,right=True)-bisect(reads_dic[chrom],i,right=False) if chrom in reads_dic else 0
                    length=j-i
                    RPKM=num*1.0/length/total_number_of_reads*(10**9)
                    reads_densities_gene.append(RPKM)
                if sum(reads_densities_gene)>0:
                    reads_densities.append(reads_densities_gene)
            else:
                temp=list(np.linspace(start+distance,start,51))[:-1]+list(np.linspace(start,end,101))+list(np.linspace(end,end-distance,51))[1:]
                for i,j in zip(temp[:-1],temp[1:]):
                    num=bisect(reads_dic[chrom],i,right=True)-bisect(reads_dic[chrom],j,right=False) if chrom in reads_dic else 0
                    length=i-j
                    RPKM=num*1.0/length/total_number_of_reads*(10**9)
                    reads_densities_gene.append(RPKM)
                if sum(reads_densities_gene)>0:
                    reads_densities.append(smoothing(reads_densities_gene)) 
    reads_densities=sorted(reads_densities,key=lambda x:np.mean(x),reverse=True)

    return reads_densities

def smoothing(a_list,step=5):
    if len(a_list)<5:
        return [np.mean(a_list)]*len(a_list)
    else:
        smoothing_list=[]
        b_list=a_list+[a_list[-1]]*step
        i=0
        while i<len(a_list):
            smoothing_list.append(np.mean(b_list[i:i+step]))
            i+=1
        return smoothing_list

def plot_heatmap(matrix,out_dir,distance,name):
    plt.switch_backend('agg')
    plt.figure(figsize=(6,15))
    gs=mpl.gridspec.GridSpec(15,10,wspace=1,hspace=0.1)
    ax=plt.subplot(gs[0:4,0:9])
    matrix=np.array(matrix)
    ys=np.mean(matrix,axis=0)
    plt.plot([i for i in range(len(ys))],ys)
    plt.ylabel('RPKM')
    plt.xlim(0,len(ys)-1)
    plt.xticks([],[])
    ys=ys[ys!=0]
    TSS_enrichment_score=max(ys)/min(ys)    
    plt.title('%s (TSS enrichment score:%.3f)'%(name,TSS_enrichment_score))


    ax=plt.subplot(gs[4:,0:9])
    vmin=np.quantile(matrix,1.0/10)
    vmax=np.quantile(matrix,9.0/10)
    print('Min: ',vmin,'Max: ',vmax)
    plt.imshow(matrix,cmap='Reds',aspect='auto',vmin=vmin,vmax=vmax)
    plt.xticks([0,50,150,200],['-%.1fkb'%(distance*1.0/1000),'TSS','TES','+%.1fkb'%(distance*1.0/1000)])
    plt.yticks([],[])


    ax=plt.subplot(gs[4:,9])
    norm=mpl.colors.Normalize(vmin=vmin,
                              vmax=vmax)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=plt.get_cmap('Reds'),norm=norm)
    cb1.set_label('RPKM',size=10)
    cb1.ax.tick_params(labelsize=8)

    plt.savefig('%s/%s_heatmap.pdf'%(out_dir,name),bbox_inches='tight')
    plt.savefig('%s/%s_heatmap.png'%(out_dir,name),bbox_inches='tight')
    plt.close('all')

def write_out(matrix,out_dir,name):
    with open('%s/%s_heatmap.txt'%(out_dir,name),'w') as outfile:
        for i in matrix:
            outfile.write('\t'.join([str(j) for j in i])+'\n')

def main():
    path=''
    out_dir=False
    ref_path=''
    distance=2000
    shift_size=100
    out=False
    name='Test'

    try:
        opts,args=getopt(argv[1:],'h',['input=','name=','outdir=','ref=','distance=','shift_size=','matrix','help'])
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
            elif i=='--matrix':
                out=True               
            elif i=='--shift_size':                
                try:
                    shift_size=int(j)
                    assert shift_size>0
                except (ValueError, AssertionError):
                    raise Exception("Invalid argument for %s: %s"%(i,j))                 
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
    print(path,out_dir,ref_path,distance,shift_size)
    matrix=couting_reads(path,ref_path,distance,shift_size)
    plot_heatmap(matrix,out_dir,distance,name)
    alaysis_finished=time.time()
    print('Time: ',alaysis_finished-start_alaysis)
    if out:
        write_out(matrix,out_dir,name)



if __name__ == '__main__':
    main()

