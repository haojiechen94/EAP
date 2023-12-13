# generate_stats.py
# 2022-04-27
# Haojie Chen

"""
create_report.py --input_dir=folder --name=<str> --sequencing_type=ATAC/ChIPPE/ChIPSE [--outdir=output_directory]

--input_dir=<str>              Folder of the output.
--name=<str>                   Name for the analysis.
--sequencing_type=<str>        ATAC or ChIPPE or ChIPSE
[--outdir=<str>]               Output directory for the processed result.
                               Default: current directory

--help/-h                      print this page.   

Required packages: bs4, pylatex, cairosvg, fitz, using pip to install these packages.

[--cover_image=<str>]          Path to the cover image (*.png) 
[--logo_image=<str>]           Path to the logo image (*.png)
[--workflow_chart=<str>]       Path to the workflow chart (*.png)

"""


import zipfile
import glob
import os
from sys import argv, stderr, stdin, stdout
from getopt import getopt
from pylatex import Document,PageStyle,Head,Foot,MiniPage,StandAloneGraphic,MultiColumn,Tabu,LongTabu,LargeText,MediumText,LineBreak,NewPage, Tabularx,TextColor,simple_page_number,Section,Subsection,SmallText
from pylatex.utils import bold,NoEscape
from datetime import date
import re
from bs4 import BeautifulSoup
import cairosvg
import pandas as pd

def get_QC_png(sequencing_type,input_dir,out_dir):
    png_dic={}
    out_data=[]
    if sequencing_type=='ATAC':
        for path in glob.glob(input_dir+'/step1_fastqc_and_trim_galore/fastqc/*/*/*.zip'):
            temp=path.split('/')
            sample_name=temp[-3]
            R1_or_R2=temp[-2]
            if sample_name not in png_dic:
                png_dic[sample_name]={}
                png_dic[sample_name]['before']={}
                png_dic[sample_name]['after']={}
            zip_file=zipfile.ZipFile(path)
            zip_list=zip_file.namelist()
            for f in zip_list:
                if 'per_base_quality.png' in f:
                    zip_file.extract(f,input_dir+'/step1_fastqc_and_trim_galore/fastqc/%s/%s/'%(sample_name,R1_or_R2))
                    png_dic[sample_name]['before'][R1_or_R2]=input_dir+'/step1_fastqc_and_trim_galore/fastqc/%s/%s/%s'%(sample_name,R1_or_R2,f)
            zip_file.close()
        for path in glob.glob(input_dir+'/step1_fastqc_and_trim_galore/trim_galore/*/*.zip'):
            temp=path.split('/')
            sample_name=temp[-2]
            R1_or_R2='R1' if temp[-1].split('_')[-2]=='1' else 'R2'
            zip_file=zipfile.ZipFile(path)
            zip_list=zip_file.namelist()
            for f in zip_list:
                if 'per_base_quality.png' in f:
                    zip_file.extract(f,input_dir+'/step1_fastqc_and_trim_galore/trim_galore/%s/'%(sample_name))
                    png_dic[sample_name]['after'][R1_or_R2]=input_dir+'/step1_fastqc_and_trim_galore/trim_galore/%s/%s'%(sample_name,f)
            zip_file.close()
        for sample_name in png_dic:
            for before_or_after in png_dic[sample_name]:
                for R1_or_R2 in png_dic[sample_name][before_or_after]:
                    out_data.append([sample_name,before_or_after,R1_or_R2,png_dic[sample_name][before_or_after][R1_or_R2]])
    elif sequencing_type=='ChIPPE':
        for path in glob.glob(input_dir+'/step1_fastqc_and_trim_galore/fastqc/*/*/*/*.zip'):
            temp=path.split('/')
            sample_name=temp[-4]
            treatment_or_control=temp[-3]
            R1_or_R2=temp[-2]
            if sample_name not in png_dic:
                png_dic[sample_name]={}
                png_dic[sample_name]['before']={}
                png_dic[sample_name]['before']['treatment']={}
                png_dic[sample_name]['before']['control']={}
                png_dic[sample_name]['after']={}
                png_dic[sample_name]['after']['treatment']={}
                png_dic[sample_name]['after']['control']={}
            zip_file=zipfile.ZipFile(path)
            zip_list=zip_file.namelist()
            for f in zip_list:
                if 'per_base_quality.png' in f:
                    zip_file.extract(f,input_dir+'/step1_fastqc_and_trim_galore/fastqc/%s/%s/%s/'%(sample_name,treatment_or_control,R1_or_R2))
                    png_dic[sample_name]['before'][treatment_or_control][R1_or_R2]=input_dir+'/step1_fastqc_and_trim_galore/fastqc/%s/%s/%s/%s'%(sample_name,treatment_or_control,R1_or_R2,f)
            zip_file.close()
        for path in glob.glob(input_dir+'/step1_fastqc_and_trim_galore/trim_galore/*/*/*.zip'):
            temp=path.split('/')
            sample_name=temp[-3]
            treatment_or_control=temp[-2]
            R1_or_R2='R1' if temp[-1].split('_')[-2]=='1' else 'R2'            
            zip_file=zipfile.ZipFile(path)
            zip_list=zip_file.namelist()
            for f in zip_list:
                if 'per_base_quality.png' in f:
                    zip_file.extract(f,input_dir+'/step1_fastqc_and_trim_galore/trim_galore/%s/%s/'%(sample_name,treatment_or_control))                 
                    png_dic[sample_name]['after'][treatment_or_control][R1_or_R2]=input_dir+'/step1_fastqc_and_trim_galore/trim_galore/%s/%s/%s'%(sample_name,treatment_or_control,f)
            zip_file.close()
        for sample_name in png_dic:
            for before_or_after in png_dic[sample_name]:
                for treatment_or_control in png_dic[sample_name][before_or_after]:
                    for R1_or_R2 in png_dic[sample_name][before_or_after][treatment_or_control]:
                        out_data.append([sample_name,before_or_after,treatment_or_control,R1_or_R2,png_dic[sample_name][before_or_after][treatment_or_control][R1_or_R2]])            
    elif sequencing_type=='ChIPSE':
        for path in glob.glob(input_dir+'/step1_fastqc_and_trim_galore/fastqc/*/*/*.zip'):
            temp=path.split('/')
            sample_name=temp[-3]
            treatment_or_control=temp[-2]
            if sample_name not in png_dic:
                png_dic[sample_name]={}
                png_dic[sample_name]['before']={}
                png_dic[sample_name]['after']={}
            zip_file=zipfile.ZipFile(path)
            zip_list=zip_file.namelist()
            for f in zip_list:
                if 'per_base_quality.png' in f:
                    zip_file.extract(f,input_dir+'/step1_fastqc_and_trim_galore/fastqc/%s/%s/'%(sample_name,treatment_or_control))
                    png_dic[sample_name]['before'][treatment_or_control]=input_dir+'/step1_fastqc_and_trim_galore/fastqc/%s/%s/%s'%(sample_name,treatment_or_control,f)
            zip_file.close()
        for path in glob.glob(input_dir+'/step1_fastqc_and_trim_galore/trim_galore/*/*/*.zip'):
            temp=path.split('/')
            sample_name=temp[-3]
            treatment_or_control=temp[-2]
            zip_file=zipfile.ZipFile(path)
            zip_list=zip_file.namelist()
            for f in zip_list:
                if 'per_base_quality.png' in f:
                    zip_file.extract(f,input_dir+'/step1_fastqc_and_trim_galore/trim_galore/%s/%s/'%(sample_name,treatment_or_control))
                    png_dic[sample_name]['after'][treatment_or_control]=input_dir+'/step1_fastqc_and_trim_galore/trim_galore/%s/%s/%s'%(sample_name,treatment_or_control,f)
            zip_file.close()
        for sample_name in png_dic:
            for before_or_after in png_dic[sample_name]:
                for treatment_or_control in png_dic[sample_name][before_or_after]:
                    out_data.append([sample_name,before_or_after,treatment_or_control,png_dic[sample_name][before_or_after][treatment_or_control]])   
    with open(out_dir+'/Quality_control_for_bases.txt','w') as outfile:
        for i in out_data:
            outfile.write('\t'.join(i)+'\n')
    return png_dic

def get_mapping_statistics(sequencing_type,input_dir,out_dir):
    mapping_statistics_dic={}
    if sequencing_type=='ATAC':
        for path in glob.glob(input_dir+'/step2_mapping/*/*.mapstats'):
            temp=path.split('/')
            sample_name=temp[-2]
            mapping_statistics_dic[sample_name]={}
            with open(path) as infile:
                for line in infile:
                    if '# reads processed:' in line:
                        mapping_statistics_dic[sample_name]['Total_number_of_reads']=int(re.search(r'\d+',line.strip()).group())
                    if '# reads with at least one reported alignment:' in line:
                        mapping_statistics_dic[sample_name]['Number_of_uiquely_mapped_reads']=int(re.search(r'\d+',line.strip()).group())    
        for path in glob.glob(input_dir+'/step2_mapping/*/*.dupstats'):
            temp=path.split('/')
            sample_name=temp[-2]
            with open(path) as infile:
                for line in infile:
                    if 'Unique read pairs:' in line:
                        mapping_statistics_dic[sample_name]['Number_of_deduplicated_reads']=int(re.search(r'\d+',line.strip()).group())
    elif sequencing_type=='ChIPPE' or sequencing_type=='ChIPSE':
        for path in glob.glob(input_dir+'/step2_mapping/*/*/*.mapstats'):
            temp=path.split('/')
            sample_name=temp[-3]
            treatment_or_control=temp[-2]
            if sample_name not in mapping_statistics_dic:
                mapping_statistics_dic[sample_name]={}
                mapping_statistics_dic[sample_name]['treatment']={}
                mapping_statistics_dic[sample_name]['control']={}
            with open(path) as infile:
                for line in infile:
                    if '# reads processed:' in line:
                        mapping_statistics_dic[sample_name][treatment_or_control]['Total_number_of_reads']=int(re.search(r'\d+',line.strip()).group())
                    if '# reads with at least one reported alignment:' in line:
                        mapping_statistics_dic[sample_name][treatment_or_control]['Number_of_uiquely_mapped_reads']=int(re.search(r'\d+',line.strip()).group())
        for path in glob.glob(input_dir+'/step2_mapping/*/*/*.dupstats'):              
            temp=path.split('/')
            sample_name=temp[-3]
            treatment_or_control=temp[-2]
            with open(path) as infile:
                for line in infile:
                    if 'Unique read' in line:
                        mapping_statistics_dic[sample_name][treatment_or_control]['Number_of_deduplicated_reads']=int(re.search(r'\d+',line.strip()).group())
    print(mapping_statistics_dic)
    with open(out_dir+'/Quality_control_for_mapping1.txt','w') as outfile:
        outfile.write('\t'.join(['Sample ID','Total number of reads (pairs)','Number of uiquely mapped reads (pairs)','Number of deduplicated reads (pairs)'])+'\n')
        for sample_name in mapping_statistics_dic:
            if sequencing_type=='ATAC':
                temp=[sample_name,mapping_statistics_dic[sample_name]['Total_number_of_reads'],mapping_statistics_dic[sample_name]['Number_of_uiquely_mapped_reads'],mapping_statistics_dic[sample_name]['Number_of_deduplicated_reads']]
                outfile.write('\t'.join([str(i) for i in temp])+'\n')    
            if sequencing_type=='ChIPPE' or sequencing_type=='ChIPSE':
                temp=[sample_name+':treatment',mapping_statistics_dic[sample_name]['treatment']['Total_number_of_reads'],mapping_statistics_dic[sample_name]['treatment']['Number_of_uiquely_mapped_reads'],mapping_statistics_dic[sample_name]['treatment']['Number_of_deduplicated_reads']]
                outfile.write('\t'.join([str(i) for i in temp])+'\n')
                temp=[sample_name+':control',mapping_statistics_dic[sample_name]['control']['Total_number_of_reads'],mapping_statistics_dic[sample_name]['control']['Number_of_uiquely_mapped_reads'],mapping_statistics_dic[sample_name]['control']['Number_of_deduplicated_reads']]
                outfile.write('\t'.join([str(i) for i in temp])+'\n')                            
    return mapping_statistics_dic

def get_peaks_calling_statistics(input_dir,out_dir):
    peaks_calling_statistics_dic={}
    for path in glob.glob(input_dir+'/step3_peaks_calling/*/*_peaks.bed'):
        sample_name=path.split('/')[-2]
        peaks_calling_statistics_dic[sample_name]={}
        number=0
        with open(path) as infile:
            for line in infile:
                number+=1
        peaks_calling_statistics_dic[sample_name]['Number_of_peaks']=number
    with open(input_dir+'/step6_reads_counting/NA_profile_bins_log.txt') as infile:
        for line in infile:
            if 'labs=' in line:
                labs=line.strip().split('=')[-1].split(',')
    with open(input_dir+'/step6_reads_counting/NA_profile_bins_log.txt') as infile:
        peaks=[]
        this=False
        for line in infile:
            if this:
                temp=line.strip().split('\t')
                peaks.append([int(temp[1]),int(temp[2])])
            if 'peak_file\tpeaks_after_filtering\tpeaks_after_merging' in line:
                this=True
            if len(peaks)==len(labs):
                break
        for i,j in zip(labs,peaks):
            peaks_calling_statistics_dic[i]['Number_of_peaks_after_filtering']=j[0]
            peaks_calling_statistics_dic[i]['Number_of_peaks_after_merging']=j[1]
    with open(input_dir+'/step6_reads_counting/NA_profile_bins_log.txt') as infile: 
        bins=[]       
        this=False
        for line in infile:
            if this:
                temp=line.strip().split('\t')
                bins.append([int(temp[1]),int(temp[2]),float(temp[3][:-1])])
            if "read_file\treads_after_filtering\treads_within_bins\twithin_ratio" in line:
                this=True
            if len(bins)==len(labs):
                break
        for i,j in zip(labs,bins):
            peaks_calling_statistics_dic[i]['Number_of_reads_after_filtering']=j[0]
            peaks_calling_statistics_dic[i]['Number_of_reads_within_peaks']=j[1]
            peaks_calling_statistics_dic[i]['Ratio_of_reads_within_peaks']=j[2] 
    with open(out_dir+'/Quality_control_for_peaks_calling1.txt','w') as outfile:
        outfile.write('\t'.join(['Sample ID','Number of peaks','Number of peaks after filtering','Number of peaks after merging','Number_of reads after filtering','Number of reads_within_peaks','Ratio of reads within peaks'])+'\n')
        for sample_name in peaks_calling_statistics_dic:
            temp=[sample_name,peaks_calling_statistics_dic[sample_name]['Number_of_peaks'],peaks_calling_statistics_dic[sample_name]['Number_of_peaks_after_filtering'],peaks_calling_statistics_dic[sample_name]['Number_of_peaks_after_merging'],peaks_calling_statistics_dic[sample_name]['Number_of_reads_after_filtering'],peaks_calling_statistics_dic[sample_name]['Number_of_reads_within_peaks'],peaks_calling_statistics_dic[sample_name]['Ratio_of_reads_within_peaks']]
            outfile.write('\t'.join([str(i) for i in temp])+'\n')

    return peaks_calling_statistics_dic  

def get_known_result(path):
    prefix='/'.join(path.split('/')[:-1]+['knownResults/'])

    html=open(path,'r',encoding='utf-8')
    bsObj=BeautifulSoup(html,features='lxml')
    temp=bsObj.find('table').findAll('tr')[0].findAll('td')    

    columns=[temp[0].text.strip().replace('\n',''),temp[2].text.strip().replace('\n',''),temp[3].text.strip().replace('\n',''),
             temp[4].text.strip().replace('\n',''),temp[5].text.strip().replace('\n',''),temp[6].text.strip().replace('\n',''),
             temp[7].text.strip().replace('\n',''),temp[8].text.strip().replace('\n',''),temp[9].text.strip().replace('\n',''),
             temp[11].text.strip().replace('\n','')]

    data=[]
    for i in bsObj.find('table').findAll('tr')[1:]:
        temp=i.findAll('td')
        data.append([
            temp[0].text.strip().replace('\n',''),
            temp[2].text.strip().replace('\n',''),
            temp[3].text.strip().replace('\n',''),
            temp[4].text.strip().replace('\n',''),
            temp[5].text.strip().replace('\n',''),
            temp[6].text.strip().replace('\n',''),
            temp[7].text.strip().replace('\n',''),
            temp[8].text.strip().replace('\n',''),
            temp[9].text.strip().replace('\n',''),
            prefix+temp[11].find('a')['href'].split('/')[-1],
        ])
    html.close()
    return pd.DataFrame(data,columns=columns)

def get_homer_result(path):
    prefix='/'.join(path.split('/')[:-1]+['homerResults/'])
    html=open(path,'r',encoding='utf-8')
    bsObj=BeautifulSoup(html,features='lxml')
    temp=bsObj.find('table').findAll('tr')[0].findAll('td')
    columns=[temp[0].text.strip().replace('\n',''),temp[2].text.strip().replace('\n',''),temp[3].text.strip().replace('\n',''),
             temp[4].text.strip().replace('\n',''),temp[5].text.strip().replace('\n',''),temp[6].text.strip().replace('\n',''),
             temp[7].text.strip().replace('\n',''),temp[8].text.strip().replace('\n','')]

    data=[]
    for i in bsObj.find('table').findAll('tr')[1:]:
        temp=i.findAll('td')
        data.append([
            temp[0].text.strip().replace('\n',''),
            temp[2].text.strip().replace('\n',''),
            temp[3].text.strip().replace('\n',''),
            temp[4].text.strip().replace('\n',''),
            temp[5].text.strip().replace('\n',''),
            temp[6].text.strip().replace('\n',''),
            temp[7].text.split('More Information')[0],
            prefix+temp[8].find('a')['href'].split('/')[-1].split('.')[0]+'.logo.svg'
        ])
    html.close()
    return pd.DataFrame(data,columns=columns)

def create_report(name,sequencing_type,input_dir,out_dir,logo_image_path='/bdp-picb/nfspv/dev-rsgeno--rsgeno/references/logo_and_cover_image/logo.png',cover_image_path='/bdp-picb/nfspv/dev-rsgeno--rsgeno/references/logo_and_cover_image/cover.png',workflow_chart_path='/picb/rsgeno/chenhj/reports/out_dir/reports/v3/ChIP/workflow_chart.png'):
    geometry_options={
        "head":"40pt",
        "margin":"0.5in",
        "bottom":"0.6in",
        "includeheadfoot":True
    }
    doc=Document(geometry_options=geometry_options)

    #Generating first page style
    first_page=PageStyle("firstpage")

    #Header logo image
    with first_page.create(Head("L")) as header_left:
        with header_left.create(MiniPage(width=NoEscape(r"0.49\textwidth"),pos='c')) as logo_wrapper:
            logo_file=logo_image_path
            logo_wrapper.append(StandAloneGraphic(image_options="width=60px",filename=logo_file))

    #Add document title
    with first_page.create(Head("R")) as right_header:
        with right_header.create(MiniPage(width=NoEscape(r"0.49\textwidth"),pos='c',align='r')) as title_wrapper:
            title_wrapper.append(LargeText(bold("Epigenome Analysis Platform")))
            title_wrapper.append(LineBreak())
            title_wrapper.append(LineBreak())
            title_wrapper.append(MediumText(bold("Report")))


    # Add footer
    with first_page.create(Foot("C")) as footer:
        with footer.create(Tabularx(
                "X X X X",
                width_argument=NoEscape(r"\textwidth"))) as footer_table:
            footer_table.add_hline(color="black")
            document_details = MiniPage(width=NoEscape(r"0.25\textwidth"),
                                        pos='t', align='r')
            document_details.append(simple_page_number())
            footer_table.add_row(['', '', '', document_details])


    #Add logo image
    doc.preamble.append(first_page)
    doc.change_document_style("firstpage")
    doc.add_color(name="lightgray",model="gray",description="0.80")

    # Add cover image
    with doc.create(LongTabu("X[c] X[3c] X[c]",row_height=1.5)) as data_table:
        cheque_file=cover_image_path
        cheque=StandAloneGraphic(cheque_file, image_options="width=250px")
        [data_table.add_empty_row() for i in range(8)]
        data_table.add_row(['',cheque,''])

    #Add name for analysis
    with doc.create(LongTabu("X[1c] X[4c] X[1c]",row_height=2)) as data_table:
        [data_table.add_empty_row() for i in range(2)]
        data_table.add_row(['','Analysis for %s'%(name),''],mapper=[bold,LargeText])  
        [data_table.add_empty_row() for i in range(2)]
        data_table.add_row(['','Date: %s'%(date.today().strftime('%d-%B-%Y')),''],mapper=[LargeText])
    doc.append(NewPage())

    doc.append(NoEscape(r'\tableofcontents'))
    doc.append(NewPage())

    #Introduction
    with doc.create(Section('Introduction')):
        doc.append('In this analyis, we used FastQC for reads bases quality control and Trim-galore for cutting adapters and low quality bases in reads. Resulting ChIP/ATAC-seq reads were aligned to the user specified reference genome by using Bowtie. Then PCR duplicates were removing based on the genomic position. The remaining reads were used for peaks calling by using MACS. Finally, the count table was generated by using MAnorm2-utils. Downstream analyses were performed based on this count table (e.g. differential analysis, differential TF motif enrichment analysis and differential functional annotation). Users can further perform interactive analysis and customize their analysis by selecting an appropriate paramter (e.g. adjusted p-value cutoff, number of clusters) based on their context and requirements (Tools implemented in EAP).')
        with doc.create(LongTabu("X[5c]",row_height=2)) as data_table:
            cheque_file=workflow_chart_path
            cheque=StandAloneGraphic(cheque_file, image_options="width=400px") 
            data_table.add_row([cheque])     
            data_table.add_row([bold('Workflow chart')])

    doc.append(NewPage())

    #Quality control for bases
    with doc.create(Section('Quality control for bases')):
        doc.append('Here we showed the bases quality before and after trimming. Each plot shows the Q-scores in different positions on sequencing reads. Q-score represents base call quality, defined by the following equation: Q-score=-10log10(e) where e is the estimated probability of the base call being wrong. Good quality (Illumina) is generally Q-score>28. Concerning quality (Illumina) is Q-score<20.')
        with doc.create(Subsection('The figures below show the quality control of bases before and after trimming.')):            
            png_dic=get_QC_png(sequencing_type,input_dir,out_dir)
            with doc.create(LongTabu("X[2l] X[2l]",row_height=1.5)) as data_table:
                [data_table.add_empty_row() for i in range(2)]
                if sequencing_type=='ATAC':
                    for sample_name in png_dic:
                        data_table.add_row([sample_name,''])
                        data_table.add_row(['R1','R2'])
                        data_table.add_row(['Before trimming',''])
                        data_table.add_row([StandAloneGraphic(png_dic[sample_name]['before']['R1'],image_options="width=180px"),StandAloneGraphic(png_dic[sample_name]['before']['R2'],image_options="width=180px")])
                        data_table.add_row(['After trimming',''])
                        data_table.add_row([StandAloneGraphic(png_dic[sample_name]['after']['R1'],image_options="width=180px"),StandAloneGraphic(png_dic[sample_name]['after']['R2'],image_options="width=180px")])
                elif sequencing_type=='ChIPPE':
                    for sample_name in png_dic:
                        data_table.add_row([sample_name,''])
                        data_table.add_row(['Treatment',''])
                        data_table.add_row(['R1','R2'])
                        data_table.add_row(['Before trimming',''])
                        data_table.add_row([StandAloneGraphic(png_dic[sample_name]['before']['treatment']['R1'],image_options="width=180px"),StandAloneGraphic(png_dic[sample_name]['before']['treatment']['R2'],image_options="width=180px")])
                        data_table.add_row(['After trimming',''])
                        data_table.add_row([StandAloneGraphic(png_dic[sample_name]['after']['treatment']['R1'],image_options="width=180px"),StandAloneGraphic(png_dic[sample_name]['after']['treatment']['R2'],image_options="width=180px")])
                        data_table.add_row(['Control',''])
                        data_table.add_row(['R1','R2'])
                        data_table.add_row(['Before trimming',''])
                        data_table.add_row([StandAloneGraphic(png_dic[sample_name]['before']['control']['R1'],image_options="width=180px"),StandAloneGraphic(png_dic[sample_name]['before']['control']['R2'],image_options="width=180px")])
                        data_table.add_row(['After trimming',''])
                        data_table.add_row([StandAloneGraphic(png_dic[sample_name]['after']['control']['R1'],image_options="width=180px"),StandAloneGraphic(png_dic[sample_name]['after']['control']['R2'],image_options="width=180px")])                        
                elif sequencing_type=='ChIPSE':
                    for sample_name in png_dic:
                        print(sample_name)
                        data_table.add_row([sample_name,''])
                        data_table.add_row(['Treatment','Control'])
                        data_table.add_row(['Before trimming',''])
                        data_table.add_row([StandAloneGraphic(png_dic[sample_name]['before']['treatment'],image_options="width=180px"),StandAloneGraphic(png_dic[sample_name]['before']['control'],image_options="width=180px")])
                        data_table.add_row(['After trimming',''])
                        data_table.add_row([StandAloneGraphic(png_dic[sample_name]['after']['treatment'],image_options="width=180px"),StandAloneGraphic(png_dic[sample_name]['after']['control'],image_options="width=180px")])                    

    doc.append(NewPage())
    #Quality control for reads mapping
    with doc.create(Section('Quality control for reads mapping')):
        doc.append('Here we showed the mapping statistics in this analysis. Subsection 3.1 showed total number of reads (pairs), number of uniquely mapped reads (pairs) (i.e. a read that maps to a signle position in reference genome) and number of deduplicated reads (pairs) (only keep at most one read at each genomic position). Subsection 3.2 showed the reads distribution across gene body (from TSS-2kb to TES+2kb). This coverage is calculated as the number of reads in each bin (normalized by Reads Per Kilobase per Million mapped reads, RPKM), where bins are a list of consecutive windows for counting reads.')
        with doc.create(Subsection('The table below shows the number of reads keeped after reads mapping and removing PCR duplicates.')):

            with doc.create(LongTabu("X[2.5l] X[2l] X[2l] X[2l]",
                                     row_height=1.5)) as data_table:
                data_table.add_hline()
                data_table.add_row(["Sample_ID",
                                    "Total number of reads (pairs)",
                                    "Number of uniquely mapped reads (pairs)",
                                    "Number of deduplicated reads (pairs)"],
                                   mapper=bold,color="lightgray")
                data_table.add_hline()

                mapping_statistics_dic=get_mapping_statistics(sequencing_type,input_dir,out_dir)

                index=0
                for sample_name in mapping_statistics_dic:
                    if index%2==0:
                        color='white'
                    else:
                        color='lightgray'
                    if sequencing_type=='ATAC':
                        data_table.add_row([sample_name,
                                            '%d'%(mapping_statistics_dic[sample_name]['Total_number_of_reads']),
                                            '%d (%.2f%%)'%(mapping_statistics_dic[sample_name]['Number_of_uiquely_mapped_reads'],
                                                               mapping_statistics_dic[sample_name]['Number_of_uiquely_mapped_reads']/mapping_statistics_dic[sample_name]['Total_number_of_reads']*100),
                                            '%d (%.2f%%)'%(mapping_statistics_dic[sample_name]['Number_of_deduplicated_reads'],
                                                               mapping_statistics_dic[sample_name]['Number_of_deduplicated_reads']/mapping_statistics_dic[sample_name]['Total_number_of_reads']*100)],color=color)
                    elif sequencing_type=='ChIPPE' or sequencing_type=='ChIPSE':
                        data_table.add_row([sample_name+':treatment',
                                            '%d'%(mapping_statistics_dic[sample_name]['treatment']['Total_number_of_reads']),
                                            '%d (%.2f%%)'%(mapping_statistics_dic[sample_name]['treatment']['Number_of_uiquely_mapped_reads'],
                                                               mapping_statistics_dic[sample_name]['treatment']['Number_of_uiquely_mapped_reads']/mapping_statistics_dic[sample_name]['treatment']['Total_number_of_reads']*100),
                                            '%d (%.2f%%)'%(mapping_statistics_dic[sample_name]['treatment']['Number_of_deduplicated_reads'],
                                                               mapping_statistics_dic[sample_name]['treatment']['Number_of_deduplicated_reads']/mapping_statistics_dic[sample_name]['treatment']['Total_number_of_reads']*100)],color=color)
                        data_table.add_row([sample_name+':control',
                                            '%d'%(mapping_statistics_dic[sample_name]['control']['Total_number_of_reads']),
                                            '%d (%.2f%%)'%(mapping_statistics_dic[sample_name]['control']['Number_of_uiquely_mapped_reads'],
                                                               mapping_statistics_dic[sample_name]['control']['Number_of_uiquely_mapped_reads']/mapping_statistics_dic[sample_name]['control']['Total_number_of_reads']*100),
                                            '%d (%.2f%%)'%(mapping_statistics_dic[sample_name]['control']['Number_of_deduplicated_reads'],
                                                               mapping_statistics_dic[sample_name]['control']['Number_of_deduplicated_reads']/mapping_statistics_dic[sample_name]['control']['Total_number_of_reads']*100)],color=color)                        
                    index+=1

                data_table.add_hline()
        doc.append(NewPage())
        with doc.create(Subsection('The figures below show the reads distribution across gene body.')):
            if sequencing_type=='ATAC':
                paths=glob.glob(input_dir+'/step2_mapping/*/*_heatmap.png')
            elif sequencing_type=='ChIPPE' or sequencing_type=='ChIPSE':
                paths=glob.glob(input_dir+'/step2_mapping/*/*/*_heatmap.png')
            with open(out_dir+'/Quality_control_for_mapping2.txt','w') as outfile:
                for path in paths:
                    outfile.write(path+'\n')

            with doc.create(LongTabu("X[2l] X[2l] X[2l] X[2l]",row_height=1.5)) as data_table:
                [data_table.add_empty_row() for i in range(2)]
                index=0
                step=4
                while index+step<=len(paths):
                    data_table.add_row([StandAloneGraphic(path,image_options="width=100px") for path in paths[index:index+step]])
                    index+=step
                if index<len(paths):
                    data_table.add_row([StandAloneGraphic(path,image_options="width=100px") for path in paths[index:len(paths)]]+['']*(index+step-len(paths)))

    doc.append(NewPage())
    #Quality control for peak calling
    with doc.create(Section('Quality control for peaks calling')):
        doc.append('In this section, we showed the statistics of peaks calling analysis. The table showed number of peaks (i.e. number of peaks identified by MACS), Number of peaks after filtering (i.e. number of peaks after removing those in blacklisted regions), Number of peaks after merging (i.e. number of peaks after merging peaks from each sample), Number of reads after filtering (i.e. number of reads after removing those overlapped with blacklisted regions), Number of reads within peaks (i.e. number of reads that fall into a peak) and Reads within peaks ratio (defined as the fraction of reads that fall into a peak and is often used as a measure of ChIP/ATAC-seq quality.).')
        with doc.create(Subsection('The table below shows the number of peaks and reads within peaks ratio in each sample.')):
            peaks_calling_statistics_dic=get_peaks_calling_statistics(input_dir,out_dir)
            with doc.create(LongTabu("X[3l] X[2l] X[2l] X[2l] X[2l] X[2l] X[2l]",
                                     row_height=1.5)) as data_table:
                data_table.add_hline()
                data_table.add_row(["Sample_ID",
                                    "Number of peaks",
                                    "Number of peaks after filtering",
                                    "Number of peaks after merging",
                                    "Number of reads after filtering",
                                    "Number of reads within peaks",
                                    "Reads within peaks ratio"],
                                   mapper=bold,color="lightgray"
                                   )
                data_table.add_hline()
                index=0
                for sample_name in peaks_calling_statistics_dic:
                    if index%2==0:
                        color='white'
                    else:
                        color='lightgray'
                    data_table.add_row([sample_name,'%d'%peaks_calling_statistics_dic[sample_name]["Number_of_peaks"],'%d'%peaks_calling_statistics_dic[sample_name]["Number_of_peaks_after_filtering"],'%d'%peaks_calling_statistics_dic[sample_name]["Number_of_peaks_after_merging"],'%d'%peaks_calling_statistics_dic[sample_name]["Number_of_reads_after_filtering"],'%d'%peaks_calling_statistics_dic[sample_name]["Number_of_reads_within_peaks"],'%.2f%%'%peaks_calling_statistics_dic[sample_name]["Ratio_of_reads_within_peaks"]],color=color)
                data_table.add_hline()

    doc.append(NewPage())

    #Peaks annotations               
    with doc.create(Section('Genomic annotations of peaks')):
        doc.append('In this section, peaks were annotated by relative location relationship with genomic features (Priority: promoter > exon > intron > intergenic).')
        with doc.create(Subsection('The figure below show the proportion of peaks assigned to different genomic features in each sample.')):
             with doc.create(LongTabu("X[2l]",row_height=1.5)) as data_table:
                [data_table.add_empty_row() for i in range(2)]
                data_table.add_row([StandAloneGraphic(input_dir+'/step5_peaks_annotations/peaks_annotation/Test_distribution_in_genome.png',image_options="width=280px")])
                with open(out_dir+'/Quality_control_for_peaks_calling2.txt','w') as outfile:
                    outfile.write(input_dir+'/step5_peaks_annotations/peaks_annotation/Test_distribution_in_genome.png')
    doc.append(NewPage())

    #Motif enrichment
    with doc.create(Section('Motif enrichment')):
        doc.append('In this section, we showed the motif enrichment results in each sample, including Known motif enrichment and de novo motif enrichment. Two motif enrichment analysis were both used to detecting motifs that are more enriched in ChIP/ATAC-seq peaks than random genomic regions.')
        with open(out_dir+'/Motif_enrichment_known_motif.txt','w') as outfile:
            for path in glob.glob(input_dir+'/step4_motif_enrichment/*/knownResults.html'):
                outfile.write(path+'\n')
        with open(out_dir+'/Motif_enrichment_de_novo_motif.txt','w') as outfile:
            for path in glob.glob(input_dir+'/step4_motif_enrichment/*/homerResults.html'):
                outfile.write(path+'\n')

        with doc.create(Subsection('Known motif enrichment')):
            doc.append('The table below shows the known motif enrichment reuslt using Homer.\n')
            for path in glob.glob(input_dir+'/step4_motif_enrichment/*/knownResults.html'):
                sample_name=path.split('/')[-2]
                known_result_df=get_known_result(path)     
                doc.append(LargeText(sample_name))
                with doc.create(LongTabu("X[3l] X[0.6l] X[3l] X[1l] X[1.5l] X[0.8l] X[0.8l] X[0.8l] X[0.8l] X[0.8l]",row_height=1.5)) as data_table:
                    [data_table.add_empty_row() for i in range(2)]
                    data_table.add_hline()
                    data_table.add_row(['Motif']+known_result_df.columns.tolist()[:-1],mapper=bold,color="lightgray")
                    data_table.add_hline() 
                    for i in known_result_df.index[:10]:
                        cairosvg.svg2png(url=known_result_df.loc[i,'SVG'],write_to=known_result_df.loc[i,'SVG'].split('.')[0]+'.png')
                        if i%2==0:                    
                            data_table.add_row([StandAloneGraphic(known_result_df.loc[i,'SVG'].split('.')[0]+'.png',image_options="width=80px")]+['%s'%j for j in known_result_df.iloc[i,:-1]],color="white")
                        else:
                            data_table.add_row([StandAloneGraphic(known_result_df.loc[i,'SVG'].split('.')[0]+'.png',image_options="width=80px")]+['%s'%j for j in known_result_df.iloc[i,:-1]],color="lightgray")
                    data_table.add_hline()
                    doc.append(NewPage())

        with doc.create(Subsection('de novo motif enrichment')):
            doc.append('The table below shows the de novo motif enrichment reuslt using homer.(* denotes possible false positive)\n')
            for path in glob.glob(input_dir+'/step4_motif_enrichment/*/homerResults.html'):
                sample_name=path.split('/')[-2]
                homer_result_df=get_homer_result(path)     
                doc.append(LargeText(sample_name))          
                with doc.create(LongTabu("X[3l] X[1l] X[1l] X[1l] X[1l] X[1l] X[1l] X[3l]",
                                         row_height=1.5)) as data_table:
                    [data_table.add_empty_row() for i in range(2)]
                    data_table.add_hline()
                    data_table.add_row(['Motif']+homer_result_df.columns.tolist()[:-1],
                                       mapper=[bold,SmallText,SmallText,SmallText,SmallText],color="lightgray")
                    data_table.add_hline() 
                    for i in homer_result_df.index[:10]:
                        cairosvg.svg2png(url=homer_result_df.loc[i,'Motif File'],write_to=homer_result_df.loc[i,'Motif File'].split('.')[0]+'.png')
                        if i%2==0:                    
                            data_table.add_row([StandAloneGraphic(homer_result_df.loc[i,'Motif File'].split('.')[0]+'.png',image_options="width=80px")]+[SmallText(SmallText('%s'%j)) for j in homer_result_df.iloc[i,:-1]],
                                                   color="white")
                        else:
                            data_table.add_row([StandAloneGraphic(homer_result_df.loc[i,'Motif File'].split('.')[0]+'.png',image_options="width=80px")]+['%s'%j for j in homer_result_df.iloc[i,:-1]],
                                                   color="lightgray") 
                    data_table.add_hline() 
                    doc.append(NewPage())               

    #Motif enrichment
    with doc.create(Section('Differential analysis')): 
        doc.append('Here we first performed differential analysis based on the user specified variable of interest, then we performed differential TF motif enrichment based on those differential enriched/accassible peaks (DEPs/DAPs), identifying TF motifs that wre enriched in one set of peaks regions relative to another set. Finally, we performed functional annotation to the nearest genes of DEPs/DAPs to identify enriched biological themes in different conditions.')
        doc.append('The figures below show the differential analysis and functional annotations for different comparisons.')
        png_dic={}
        for path in glob.glob(input_dir+'/step7_differential_analysis/*MA_plot.png'):
            analysis_ID=path.split('/')[-1].split('_MA_plot')[0]
            png_dic[analysis_ID]={}
            png_dic[analysis_ID]['MA_plot']=path
        for path in glob.glob(input_dir+'/step7_differential_analysis/*MVC_plot.png'):
            analysis_ID=path.split('/')[-1].split('_MVC_plot')[0]
            png_dic[analysis_ID]['MVC_plot']=path

        for path in glob.glob(input_dir+'/step8_functional_enrichment/*motif_enrichment_plot.png'):
            analysis_ID=path.split('/')[-1].split('_motif_enrichment_plot')[0]
            png_dic[analysis_ID]['motif_enrichment_plot']=path

        for path in glob.glob(input_dir+'/step8_functional_enrichment/*pathway_enrichment.png'):
            analysis_ID=path.split('/')[-1].split('_pathway_enrichment')[0]
            png_dic[analysis_ID]['pathway_enrichment']=path            

        with doc.create(LongTabu("X[2l] X[2l] X[2l]",row_height=1.5)) as data_table:
            [data_table.add_empty_row() for i in range(2)]
            for analysis_ID in png_dic:
                data_table.add_row([analysis_ID,'',''],mapper=[LargeText])
                data_table.add_row(['MVC plot','MA plot','Differential motif enrichment'])
                data_table.add_row([StandAloneGraphic(png_dic[analysis_ID]['MVC_plot'],image_options="width=100px"),StandAloneGraphic(png_dic[analysis_ID]['MA_plot'],image_options="width=100px"),StandAloneGraphic(png_dic[analysis_ID]['motif_enrichment_plot'],image_options="width=100px")])
        with open(out_dir+'/Differential_analysis.txt','w') as outfile:
            outfile.write('\t'.join(['Analysis ID','MVC plot','MA plot','Differential motif enrichment','Differential functional enrichment'])+'\n')
            for analysis_ID in png_dic:
                outfile.write('\t'.join([analysis_ID,png_dic[analysis_ID]['MVC_plot'],png_dic[analysis_ID]['MA_plot'],png_dic[analysis_ID]['motif_enrichment_plot'],png_dic[analysis_ID]['pathway_enrichment']])+'\n')

        with doc.create(LongTabu("X[2l]",row_height=1.5)) as data_table:
            [data_table.add_empty_row() for i in range(2)]
            for analysis_ID in png_dic:
                data_table.add_row([analysis_ID],mapper=[LargeText])
                data_table.add_row(['Functional enrichment'])
                data_table.add_row([StandAloneGraphic(png_dic[analysis_ID]['pathway_enrichment'],image_options="width=300px")])
    doc.append(NewPage())
    with doc.create(Section('References')): 
        doc.append('[1] Andrews, S. (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data.\n')
        doc.append('[2] Felix Krueger, Frankie James, Phil Ewels, Ebrahim Afyounian, & Benjamin Schuster-Boeckler. (2021). FelixKrueger/TrimGalore: v0.6.7 - DOI via Zenodo (0.6.7). Zenodo. https://doi.org/10.5281/zenodo.5127899.\n')
        doc.append('[3] Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol 10:R25.\n')
        doc.append('[4] Zhang, Y., Liu, T., Meyer, C.A. et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol 9, R137 (2008). https://doi.org/10.1186/gb-2008-9-9-r137.\n')
        doc.append('[5] Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589.\n')
        doc.append('[6] Tu, S., Li, M., Chen, H., Tan, F., Xu, J., Waxman, D. J., Zhang, Y., & Shao, Z. (2021). MAnorm2 for quantitatively comparing groups of ChIP-seq samples. Genome research, 31(1), 131–145. https://doi.org/10.1101/gr.262675.120.\n')            
    doc.generate_pdf(os.path.join(out_dir,'report'),clean_tex=False)


def main():
    out_dir=False

    try:
        opts,args=getopt(argv[1:],'h',['input_dir=','name=','sequencing_type=','outdir=','cover_image=','logo_image=','workflow_chart=','help'])
        for i,j in opts:   
            if i=="-h" or i=="--help":
                stdout.write(__doc__)
                exit(0)
            elif i=='--input_dir':
                input_dir=j
            elif i=='--name':
                name=j
            elif i=='--sequencing_type':
                if j in ['ATAC','ChIPPE','ChIPSE']:
                    sequencing_type=j
                else:
                    raise Exception("--sequencing_type parameter only takes ATAC or ChIPPE or ChIPSE.")
            elif i=='--outdir':
                out_dir=j
            elif i=='--cover_image':
                cover_image=j
            elif i=='--logo_image':
                logo_image=j
            elif i=='--workflow_chart':
                workflow_chart=j
            else:
                raise Exception("Internal errors occur when parsing command line arguments.")
    except Exception as e:
        stderr.write("%s\n" % e)
        stderr.write("Type 'python create_report.py --help' for more information.\n")
        exit(1)

    if not out_dir:
        out_dir=os.getcwd()
      
    print(input_dir,name,sequencing_type,out_dir)
    create_report(name,sequencing_type,input_dir,out_dir,logo_image_path=logo_image,cover_image_path=cover_image,workflow_chart_path=workflow_chart)


if __name__ == '__main__':
    main()
