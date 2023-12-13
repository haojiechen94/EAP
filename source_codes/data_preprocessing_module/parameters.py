# parameters.py
# 2022-04-18
# Haojie Chen

"""
parameters.py --peaks=pathname --summits=pathname --reads==pathname --black_list=path --metadata=path 
              [--typical_bin_size=2000] [--outdir=current_directory]  [--keep_peaks=All]

--peaks=<str>                     Peaks in BED format identified by macs. Find all pathname matching this pattern. 
                                  Using relative pathname (e.g. /usr/src/*_peaks.bed) to find any matching files in a 
                                  specified directory.

--summits=<str>                   Summits of peaks in BED format identified by macs. Find all pathname matching 
                                  this pattern. Using relative pathname (e.g. /usr/src/*_summits.bed) to find 
                                  any matching files in a specified directory.

--reads=<str>                     Mapped reads in BED format. Find all pathname matching this pattern. 
                                  Using relative pathname (e.g. /usr/src/*_drm.bed) to find any matching files in a 
                                  specified directory.

--black_list=<str>                Any reference bin overlapping some region in the list is removed from the output.

--metadata=<str>                  Metadata for the input samples, including sample names, input files and 
                                  clinical information.

[--typical_bin_size=<int>]        After merging peak regions from all the samples, the resulting regions with 
                                  a size smaller than or matching <int> bps are directly taken as bins. Each 
                                  of the others is divided into non-overlapping bins of <int> bps (except
                                  the bins at the edge of merged peaks).
                                  Default: 2000

[--outdir=<str>]                  Output directory for the parameters file.
                                  Default: current directory

[--keep_peaks=<int>]              If set, for each peak file, only <int> peaks will be used. In this case, peaks 
                                  are sorted by the 5th column of each peak file, and the <int> ones with the 
                                  greatest scores are used.
                                  Default: Keep all peaks

[--sequencing_type=<str>]         ChIP or ATAC
                                  Default: ATAC

--help/-h                         print this page.                       

"""

from sys import argv, stderr, stdin, stdout
from getopt import getopt
import os
import glob
import pandas as pd



def get_peaks_file(pathname):
    peaks_dic={}
    for path in glob.glob(pathname):
        ID=path.split('/')[-1].split('_peaks.bed')[0]
        peaks_dic[ID]=path
    return peaks_dic

def get_summits_file(pathname):
    summits_dic={}
    for path in glob.glob(pathname):
        ID=path.split('/')[-1].split('_summits.bed')[0]
        summits_dic[ID]=path
    return summits_dic

def get_reads_file(pathname,sequencing_type):
    reads_dic={}
    if sequencing_type=='ATAC':
        for path in glob.glob(pathname):
            ID=path.split('/')[-1].split('_ss1_drm.bed')[0]
            reads_dic[ID]=path
    elif sequencing_type=='ChIP':
         for path in glob.glob(pathname):
            ID=path.split('/')[-1].split('_treatment_ss100_drm.bed')[0]
            reads_dic[ID]=path
    else:
        print(sequencing_type,'Unknown sequencing type')
        exit(1)

    return reads_dic

def create_parameters_file(metadata,peaks,summits,reads,black_list,typical_bin_size,out_dir,sequencing_type,keep_peaks_number=False):
    shiftsize=1 if sequencing_type=='ATAC' else 100
    peaks_dic=get_peaks_file(peaks)
    summits_dic=get_summits_file(summits)
    reads_dic=get_reads_file(reads,sequencing_type)
    metadata_df=pd.read_csv(metadata,sep=',')
    parameters='labs='+','.join(metadata_df['sample_name'].tolist())+'\n'
    parameters=parameters+'reads='+','.join([reads_dic[i] for i in metadata_df['sample_name'].tolist()])+'\n'
    parameters=parameters+'peaks='+','.join([peaks_dic[i] for i in metadata_df['sample_name'].tolist()])+'\n'
    parameters=parameters+'summits='+','.join([summits_dic[i] for i in metadata_df['sample_name'].tolist()])+'\n'
    parameters=parameters+'typical-bin-size=%s\n'%(typical_bin_size)
    parameters=parameters+'keep-dup=all\n'
    parameters=parameters+'shiftsize=%s\n'%(shiftsize)
    parameters=parameters+'filter=%s\n'%(black_list)
    if keep_peaks_number:
        parameters=parameters+'--keep-peaks=%d\n'%(keep_peaks_number)

    with open(out_dir+'/parameters.txt','w') as outfile:
        outfile.write(parameters)

def main():

    out_dir=False
    keep_peaks=False
    peaks=''
    summits=''
    reads=''
    black_list=''
    metadata=''
    genome=''
    typical_bin_size=2000
    sequencing_type='ATAC'

    try:
        opts,args=getopt(argv[1:],'h',['peaks=','summits=','reads=','black_list=','metadata=','typical_bin_size=','outdir=','keep_peaks=','sequencing_type=','help'])
        for i,j in opts:   
            if i=="-h" or i=="--help":
                stdout.write(__doc__)
                exit(0)
            elif i=='--peaks':
                peaks=j
            elif i=='--summits':
                summits=j
            elif i=='--reads':
                reads=j
            elif i=='--black_list':
                black_list=j
            elif i=='--metadata':
                metadata=j
            elif i=='--outdir':
                out_dir=j
            elif i=='--sequencing_type':
                sequencing_type=j                
            elif i=='--keep_peaks':
                try:
                    keep_peaks=int(j)
                    assert keep_peaks>0
                except (ValueError, AssertionError):
                    raise Exception("Invalid argument for %s: %s"%(i,j))
            elif i=='--typical_bin_size':
                try:
                    typical_bin_size=int(j)
                    assert typical_bin_size>0
                except (ValueError, AssertionError):
                    raise Exception("Invalid argument for %s: %s"%(i,j))                      
            else:
                raise Exception("Internal errors occur when parsing command line arguments.")
    except Exception as e:
        stderr.write("%s\n" % e)
        stderr.write("Type 'python parameters.py --help' for more information.\n")
        exit(1)

    if not out_dir:
        out_dir=os.getcwd()

    create_parameters_file(metadata,peaks,summits,reads,black_list,typical_bin_size,out_dir,sequencing_type,keep_peaks_number=False)

if __name__ == '__main__':
    main()






