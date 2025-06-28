# generate_stats.py
# 2021-10-12
# Haojie Chen

"""
generate_stats.py --input=folder --name=<str> --tech=ChIP/ATAC [--output_file_name=<str>] [--outdir=output_directory]

--input=<str>                  Folder of the output.
--name=<str>                   Name for the analysis.
--tech=ChIP/ATAC               ChIP/ATAC.
[--output_file_name=<str>]     Output file name.
                               Default: temp
[--outdir=<str>]               Output directory for the processed result.
                               Default: current directory

--help/-h                      print this page.   

Required packages: bs4, pylatex, cairosvg, fitz, using pip to install these packages.

--cover_image=<str>            Path to the cover image (*.png) 
--logo_image=<str>             Path to the logo image (*.png)

"""



import zipfile
import os
import glob
import re
from bs4 import BeautifulSoup
import pandas as pd
from sys import argv, stderr, stdin, stdout
from getopt import getopt
from pylatex import Document,PageStyle,Head,Foot,MiniPage,StandAloneGraphic,MultiColumn,Tabu,LongTabu,LargeText,MediumText,LineBreak,NewPage, Tabularx,TextColor,simple_page_number,Section,Subsection,SmallText
from pylatex.utils import bold,NoEscape
from datetime import date
import cairosvg
import fitz

def get_known_result(path,Motif_7_folder_abs):
    prefix=glob.glob(Motif_7_folder_abs+'/*/knownResults/')[0]
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


def get_homer_result(path,Motif_7_folder_abs):
    prefix=glob.glob(Motif_7_folder_abs+'/*/homerResults/')[0]
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

def pdf2png(Annotation_6_folder_abs):
    png_paths_dic={'Peak_distribution_over_chromosomes':[],
                   'Peak_distribution_over_genome_location':[],
                   'Peak_function_annotation':[]}
    zoom_x,zoom_y,rotation_angle=[5,5,0]
    for pdf_path in glob.glob(os.path.join(Annotation_6_folder_abs,'*/*.pdf')):
        name=pdf_path.split('/')[-1].split('.')[0]
        pdf=fitz.open(pdf_path)
        for pg in range(0,pdf.pageCount):
            page=pdf[pg]
            trans=fitz.Matrix(zoom_x,zoom_y).preRotate(rotation_angle)
            pm=page.getPixmap(matrix=trans,alpha=False)
            png_path=os.path.join('/'.join(pdf_path.split('/')[:-1]),'_'.join([name.replace(' ','_'),str(pg),".png"]))
            pm.writePNG(png_path)
            if 'Peak_distribution_over_chromosomes' in png_path:
                png_paths_dic['Peak_distribution_over_chromosomes'].append(png_path)
            if 'Peak_distribution_over_genome_location' in png_path:
                png_paths_dic['Peak_distribution_over_genome_location'].append(png_path)
            if 'Peak_function_annotation' in png_path:
                png_paths_dic['Peak_function_annotation'].append(png_path)
    return png_paths_dic

def generate(name='ATAC-seq',
             folder='/media/chenhaojie/Data/epigenome_platform/chip/',
             cover_image_path='/home/chenhaojie/Pictures/cover.png',
             logo_image_path='/home/chenhaojie/Pictures/1.png',
             out_dir='/media/chenhaojie/Data/epigenome_platform',
             output_file_name='report',
             tech='ChIP'):
    #1-QC
    QC_1_folder_abs=folder+'/1-QC/'
    QC_1_image_paths=[]
    for path in glob.glob(os.path.join(QC_1_folder_abs,'*_fastqc.zip')):
        zip_file=zipfile.ZipFile(path)
        zip_list=zip_file.namelist()
        for f in zip_list:
            if 'per_base_quality.png' in f:
                zip_file.extract(f,QC_1_folder_abs)
                QC_1_image_paths.append(os.path.join(QC_1_folder_abs,f))
        zip_file.close()

    #2-Trim
    Trim_2_folder_abs=folder+'/2-Trim/'
    Trim_2_image_paths=[]
    for path in glob.glob(os.path.join(Trim_2_folder_abs,'*_fastqc.zip')):
        zip_file=zipfile.ZipFile(path)
        zip_list=zip_file.namelist()
        for f in zip_list:
            if 'per_base_quality.png' in f:
                zip_file.extract(f,Trim_2_folder_abs)
                Trim_2_image_paths.append(os.path.join(Trim_2_folder_abs,f))
        zip_file.close()

    #3-Alignment
    dic={}
    Alignment_3_folder_abs=folder+'/3-Alignment/'
    if tech=='ChIP':
        for path in glob.glob(os.path.join(Alignment_3_folder_abs,'*.stat')):
            key=''
            if path.split('/')[-1].split('_')[-1]=='sort.stat':
                key='Total number of reads'
            elif path.split('/')[-1].split('_')[-1]=='filter.stat':
                key='Number mapped of reads'

            with open(path) as infile:
                for line in infile:
                    dic[key]=re.search(r'\d+',line.strip()).group()
                    break
    elif tech=='ATAC':
        for path in glob.glob(os.path.join(Alignment_3_folder_abs,'*.stat')):
            temp=[]
            with open(path) as infile:
                for line in infile:
                    temp.append(re.search(r'\d+',line.strip()).group())
            dic['Total number of reads']=temp[0]
            dic['Number mapped of reads']=temp[4]

    #4-MACS
    MACS_4_folder_abs=folder+'/4-MACS/'
    number_of_peaks=0
    for path in glob.glob(os.path.join(MACS_4_folder_abs,'*_peaks.narrowPeak')):
        with open(path) as infile:
            for line in infile:
                number_of_peaks+=1

    #/6-Annotation/
    Annotation_6_folder_abs=folder+'/6-Annotation/'
    png_paths_dic=pdf2png(Annotation_6_folder_abs)

    #7-Motif
    Motif_7_folder_abs=folder+'/7-Motif/'
    homer_result_df=get_homer_result(glob.glob(Motif_7_folder_abs+'/*/homerResults.html')[0],Motif_7_folder_abs)

    known_result_df=get_known_result(glob.glob(Motif_7_folder_abs+'/*/knownResults.html')[0],Motif_7_folder_abs)

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
            logo_wrapper.append(StandAloneGraphic(image_options="width=400px",filename=logo_file))

    #Add document title
    with first_page.create(Head("R")) as right_header:
        with right_header.create(MiniPage(width=NoEscape(r"0.49\textwidth"),pos='c',align='r')) as title_wrapper:
            title_wrapper.append(LargeText(bold("Analysis Report")))


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
    #Quality control
    with doc.create(Section('Quality control')):
        doc.append('The figures below show the quality control of bases before and after trimming.')

        with doc.create(LongTabu("X[2l] X[2l]",
                                 row_height=1.5)) as data_table:
            [data_table.add_empty_row() for i in range(2)]
            data_table.add_row(['Before trimming',''])
            data_table.add_row([StandAloneGraphic(path,image_options="width=180px") for path in QC_1_image_paths])
            data_table.add_row(['After trimming',''])
            data_table.add_row([StandAloneGraphic(path,image_options="width=180px") for path in Trim_2_image_paths])
        
    doc.append(NewPage())

    #Mapping  and peaks calling statistic
    with doc.create(Section('Mapping statistic')):
        doc.append('The table below shows the total reads number and mapped reads number for reads mapping and number of peaks for peaks calling.')

        with doc.create(LongTabu("X[1.5l] X[2l] X[2l] X[2l]",
                                 row_height=1.5)) as data_table:
            [data_table.add_empty_row() for i in range(2)]
            data_table.add_hline()
            data_table.add_row(['Sample_ID',
                                'Total number of reads',
                                'Number mapped of reads',
                                'Total number of peaks'],
                                mapper=bold,color="lightgray")
            data_table.add_hline()

            data_table.add_row([name,dic['Total number of reads'],dic['Number mapped of reads'],number_of_peaks])
            data_table.add_hline()  
    
    doc.append(NewPage())
    #Peaks annotation    
    with doc.create(Section('Peaks annotations')):
        doc.append('The figures below shows genomic location annotation and functional annotation of the identified peaks.')
        with doc.create(LongTabu("X[3l] X[3l] X[3l]",
                                 row_height=2)) as data_table:
            [data_table.add_empty_row() for i in range(1)]
            data_table.add_row(['Peak distribution over chromosomes','',''])
            data_table.add_row([StandAloneGraphic(png_paths_dic['Peak_distribution_over_chromosomes'][0],image_options="width=150px"),'',''])
            data_table.add_row(['Peak distribution over genome location','',''])
            data_table.add_row([StandAloneGraphic(png_paths_dic['Peak_distribution_over_genome_location'][0],image_options="width=150px"),StandAloneGraphic(png_paths_dic['Peak_distribution_over_genome_location'][1],image_options="width=150px"),StandAloneGraphic(png_paths_dic['Peak_distribution_over_genome_location'][2],image_options="width=150px")])
            data_table.add_row(['Peak function annotation','',''])
            data_table.add_row([StandAloneGraphic(png_paths_dic['Peak_function_annotation'][0],image_options="width=150px"),'',''])

    doc.append(NewPage())  


    #Motif enrichment
    with doc.create(Section('Known motif enrichment')):
        doc.append('The table below shows the known motif enrichment reuslt using Homer.')

        with doc.create(LongTabu("X[3l] X[0.6l] X[3l] X[1l] X[1.5l] X[0.8l] X[0.8l] X[0.8l] X[0.8l] X[0.8l]",
                                 row_height=1.5)) as data_table:
            [data_table.add_empty_row() for i in range(2)]
            data_table.add_hline()
            data_table.add_row(['Motif']+known_result_df.columns.tolist()[:-1],
                               mapper=bold,color="lightgray")
            data_table.add_hline() 
            for i in known_result_df.index:
                cairosvg.svg2png(url=known_result_df.loc[i,'SVG'],write_to=known_result_df.loc[i,'SVG'].split('.')[0]+'.png')
                if i%2==0:                    
                    data_table.add_row([StandAloneGraphic(known_result_df.loc[i,'SVG'].split('.')[0]+'.png',image_options="width=80px")]+['%s'%j for j in known_result_df.iloc[i,:-1]],
                                           color="white")
                else:
                    data_table.add_row([StandAloneGraphic(known_result_df.loc[i,'SVG'].split('.')[0]+'.png',image_options="width=80px")]+['%s'%j for j in known_result_df.iloc[i,:-1]],
                                           color="lightgray")
         
  
    doc.append(NewPage())    
    #Motif enrichment        
    with doc.create(Section('de novo motif enrichment')):
        doc.append('The table below shows the de novo motif enrichment reuslt using homer.(* denotes possible false positive)')  
        with doc.create(LongTabu("X[3l] X[1l] X[1l] X[1l] X[1l] X[1l] X[1l] X[3l]",
                                 row_height=1.5)) as data_table:
            [data_table.add_empty_row() for i in range(2)]
            data_table.add_hline()
            data_table.add_row(['Motif']+homer_result_df.columns.tolist()[:-1],
                               mapper=[bold,SmallText,SmallText,SmallText,SmallText],color="lightgray")
            data_table.add_hline() 
            for i in homer_result_df.index:
                cairosvg.svg2png(url=homer_result_df.loc[i,'Motif File'],write_to=homer_result_df.loc[i,'Motif File'].split('.')[0]+'.png')
                if i%2==0:                    
                    data_table.add_row([StandAloneGraphic(homer_result_df.loc[i,'Motif File'].split('.')[0]+'.png',image_options="width=80px")]+[SmallText(SmallText('%s'%j)) for j in homer_result_df.iloc[i,:-1]],
                                           color="white")
                else:
                    data_table.add_row([StandAloneGraphic(homer_result_df.loc[i,'Motif File'].split('.')[0]+'.png',image_options="width=80px")]+['%s'%j for j in homer_result_df.iloc[i,:-1]],
                                           color="lightgray")          
        
    doc.append(NewPage()) 
    #References  
    with doc.create(Section('References')):
        doc.append('[1]Fast gapped-read alignment with Bowtie 2.\n[2]Improved peak-calling with MACS2.\n[3]ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization.\n[4]Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities.\n')
    doc.generate_pdf(os.path.join(out_dir,output_file_name), clean_tex=False)

def main():
    name='Test'
    out_dir=False
    output_file_name=False
    folder=''
    tech=''
    cover_image='/home/chenhaojie/Pictures/cover.png'
    logo_image='/home/chenhaojie/Pictures/1.png'
    try:
        opts,args=getopt(argv[1:],'h',['name=','input=','output_file_name=','tech=','outdir=','cover_image=','logo_image=','help'])
        for i,j in opts:   
            if i=="-h" or i=="--help":
                stdout.write(__doc__)
                exit(0)
            elif i=='--name':
                name=j                             
            elif i=='--input':
                folder=j
            elif i=='--tech':
                if j in ['ChIP','ATAC']:
                    tech=j
                else:
                    raise Exception("--tech parameter only takes ChIP or ATAC.")
            elif i=='--output_file_name':
                output_file_name=j
            elif i=='--outdir':
                out_dir=j
            elif i=='--cover_image':
                cover_image=j
            elif i=='--logo_image':
                logo_image=j
            else:
                raise Exception("Internal errors occur when parsing command line arguments.")
    except Exception as e:
        stderr.write("%s\n" % e)
        stderr.write("Type 'python generate_report.py --help' for more information.\n")
        exit(1)
    print(cover_image)
    print(logo_image)
    if not out_dir:
        out_dir=os.getcwd()

    if not output_file_name:
        output_file_name='temp'        
    generate(name=name,
             folder=folder,
             out_dir=out_dir,
             output_file_name=output_file_name,
             cover_image_path=cover_image,
             logo_image_path=logo_image,
             tech=tech)


if __name__=='__main__':
    main()
