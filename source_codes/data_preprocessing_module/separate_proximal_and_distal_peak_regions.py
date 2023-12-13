# separate_proximal_and_distal_regions.py
# 2021-08-20
# Haojie Chen

"""
separate_proximal_and_distal_regions.py --input=path [--outdir=output_directory]
                                        --ref=reference_genome_file [--distance=2000]
                                        [--bed]

--input=<str>          Input file from the output of profile_bins.
[--outdir=<str>]       Output directory for the processed result.
                       Default: current directory
--ref=<str>            Path to genome annotation file, for example hg19_refGene.gtf. 
                       File must be in GTF format.
[--distance=<int>]     Distance used to defined proximal regions. Transcription start site 
                       (TSS) +/- distance defined as proximal regions. Distance must greater
                       than 0.
                       Default: 2000
[--bed]                Writing out the proximal and distal regions in BED format.
                       Default: OFF

--help/-h              print this page.                       

"""
from sys import argv, stderr, stdin, stdout
from getopt import getopt
import copy
import pandas as pd
import os


def overlap(read1,read2):
    if min(read1[1],read2[1])-max(read1[0],read2[0])>=0:
        return 1
    else:
        return -1

#merge regions
def merge(arr):
    sort_arr=copy.deepcopy(sorted(arr))
    temp=[]
    for i in sort_arr:
        if temp and temp[-1][1]>=i[0]-1:
            temp[-1][1]=max(temp[-1][1],i[1])
        else:
            temp.append(i)
    return temp

#create promoter regions
def get_proximal_regions(path,distance):
    dic={}
    with open(path) as infile:
        for line in infile:
            temp=line.strip().split('\t')
            if temp[2]=='transcript':
                chrom,strand,TSS,TES=[temp[0],temp[6],temp[3],temp[4]]
                if chrom not in dic:
                    dic[chrom]=[]
                if strand=='+':
                    dic[chrom].append([int(TSS)-distance,int(TSS)+distance])
                else:
                    dic[chrom].append([int(TES)-distance,int(TES)+distance])
    dic1={chrom:merge(dic[chrom]) for chrom in dic}
    return dic1

#separate proximal and distal_bins
def separate_proximal_and_distal_bins(path,out_dir,ref_path,distance=2000,peak_bed=False):  
    raw_reads_count_df=pd.read_csv(path,sep='\t')
    promoter_dic=get_proximal_regions(ref_path,distance)
    chrom_list=[chrom for chrom in list(set(raw_reads_count_df['chrom'])) if chrom in promoter_dic]
    tumor_bins={}
    for chrom in chrom_list:
        tumor_bins[chrom]={}
        tumor_bins[chrom]['bins']=[]
        tumor_bins[chrom]['id']=[]

    for i in raw_reads_count_df.index:
        chrom=raw_reads_count_df.loc[i,'chrom']
        start=raw_reads_count_df.loc[i,'start']
        end=raw_reads_count_df.loc[i,'end']
        if chrom in chrom_list:
            tumor_bins[chrom]['bins'].append([start,end])
            tumor_bins[chrom]['id'].append(i)

    distal_bins={}
    proximal_bins={}
    for chrom in chrom_list:
        distal_bins[chrom]={}
        distal_bins[chrom]['bins']=[]
        distal_bins[chrom]['id']=[]
        proximal_bins[chrom]={}
        proximal_bins[chrom]['bins']=[]
        proximal_bins[chrom]['id']=[]

    for chrom in tumor_bins:
        temp_bins=tumor_bins[chrom]['bins']
        IDs=tumor_bins[chrom]['id']
        promoters=promoter_dic[chrom]
        index_of_bins=0
        index_of_promoter=0
        if len(temp_bins)==0 or len(promoters)==0:
            distal_bins[chrom]['bins']=temp_bins
            distal_bins[chrom]['id']=IDs
            continue
        current_bins=temp_bins[index_of_bins]
        current_ID=IDs[index_of_bins]
        current_promoter=promoters[index_of_promoter]
        while True:
            if overlap(current_bins,current_promoter)==1:
                proximal_bins[chrom]['bins'].append(current_bins)
                proximal_bins[chrom]['id'].append(current_ID)
                if current_bins[1]<=current_promoter[1]:
                    index_of_bins+=1
                    if index_of_bins>len(temp_bins)-1:
                        break
                    elif index_of_promoter>len(promoters)-1:
                        index_of_bins-=1
                        break
                    current_bins=temp_bins[index_of_bins]
                    current_ID=IDs[index_of_bins]
                else:
                    index_of_bins+=1
                    index_of_promoter+=1
                    if index_of_bins>len(temp_bins)-1:
                        break
                    elif index_of_promoter>len(promoters)-1:
                        index_of_bins-=1
                        break
                    current_bins=temp_bins[index_of_bins]
                    current_ID=IDs[index_of_bins]
                    current_promoter=promoters[index_of_promoter]
            else:
                if current_bins[1]<=current_promoter[1]:
                    distal_bins[chrom]['bins'].append(current_bins)
                    distal_bins[chrom]['id'].append(current_ID)
                    index_of_bins+=1
                    if index_of_bins>len(temp_bins)-1:
                        break
                    elif index_of_promoter>len(promoters)-1:
                        index_of_bins-=1
                        break
                    current_bins=temp_bins[index_of_bins]
                    current_ID=IDs[index_of_bins]
                else:
                    index_of_promoter+=1
                    if index_of_bins>len(temp_bins)-1:
                        break
                    elif index_of_promoter>len(promoters)-1:
                        index_of_bins-=1
                        break
                    current_promoter=promoters[index_of_promoter]
            if index_of_bins>len(temp_bins)-1:
                break
            elif index_of_promoter>len(promoters)-1:
                index_of_bins-=1
                break
        if index_of_bins<len(temp_bins)-1:
            for i in range(index_of_bins+1,len(temp_bins)):
                distal_bins[chrom]['bins'].append(temp_bins[i])
                distal_bins[chrom]['id'].append(IDs[i])

    stdout.write("#tumor bin: %d\n#distal bin: %d \n#proximal bin: %d\n"%(sum([len(tumor_bins[i]['id']) for i in tumor_bins]),
                                                                          sum([len(distal_bins[i]['id']) for i in distal_bins]), 
                                                                          sum([len(proximal_bins[i]['id']) for i in proximal_bins])))

    distal_peaks_ids=[]
    for i in distal_bins:
        distal_peaks_ids.extend(distal_bins[i]['id'])
    proximal_peaks_ids=[]
    for i in proximal_bins:
        proximal_peaks_ids.extend(proximal_bins[i]['id'])
        
    raw_reads_count_df.loc[proximal_peaks_ids,:].to_csv('%s/proximal_peak_regions_%dbp.txt'%(out_dir,distance),sep='\t',index=False)    
    raw_reads_count_df.loc[distal_peaks_ids,:].to_csv('%s/distal_peak_regions_%dbp.txt'%(out_dir,distance),sep='\t',index=False)
    if peak_bed:
        raw_reads_count_df.loc[proximal_peaks_ids,['chrom','start','end']].to_csv('%s/proximal_peak_regions_%dbp.bed'%(out_dir,distance),sep='\t',index=False,header=False)
        raw_reads_count_df.loc[distal_peaks_ids,['chrom','start','end']].to_csv('%s/distal_peak_regions_%dbp.bed'%(out_dir,distance),sep='\t',index=False,header=False)
        raw_reads_count_df.loc[:,['chrom','start','end']].to_csv('%s/All_peak_regions.bed'%(out_dir),sep='\t',index=False,header=False)            


def main():
    path=''
    out_dir=False
    ref_path=''
    distance=2000
    peak_bed=False
    try:
        opts,args=getopt(argv[1:],'h',['input=','outdir=','ref=','distance=','bed','help'])
        if len(args)!=0:
            if len(args)>1:
                raise Exception("Only one file required: %r" % args)
            path=args[0]
        for i,j in opts:   
            if i=="-h" or i=="--help":
                stdout.write(__doc__)
                exit(0)
            elif i=='--input':
                path=j
            elif i=='--outdir':
                out_dir=j
            elif i=='--ref':
                ref_path=j
            elif i=='--bed':
                peak_bed=True
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
        stderr.write("Type 'python separate_proximal_and_distal_regions.py --help' for more information.\n")
        exit(1)


    if not out_dir:
        out_dir=os.getcwd()

    separate_proximal_and_distal_bins(path,out_dir,ref_path,distance=distance,peak_bed=peak_bed)



if __name__ == '__main__':
    main()
