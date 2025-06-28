# generate_report.py
# 2021-10-08
# Haojie Chen

"""
generate_report.py --name=name --reads_mapping=path --peaks_calling=path --reads_counting=path
                   [--output_file_name=<str>] [--outdir=output_directory]

--name=<str>                   Name for the analysis.
--reads_mapping=<str>          Reads mapping log file.
--peaks_calling=<str>          Peaks calling log file.
--reads_counting=<str>         Reads counting log file.
[--output_file_name=<str>]     Output file name.
                               Default: temp
[--outdir=<str>]               Output directory for the processed result.
                               Default: current directory

--help/-h              print this page.                       

"""
from sys import argv, stderr, stdin, stdout
from getopt import getopt
import os
from pylatex import Document,PageStyle,Head,Foot,MiniPage,StandAloneGraphic,MultiColumn,Tabu,LongTabu,LargeText,MediumText,LineBreak,NewPage, Tabularx,TextColor,simple_page_number,Section,Subsection
from pylatex.utils import bold,NoEscape
from datetime import date
import pandas as pd

def main():
    reads_mapping=''
    peaks_calling=''
    reads_counting=''
    name='Lung cancer ATAC-seq data set'
    out_dir=False
    output_file_name=False
    try:
        opts,args=getopt(argv[1:],'h',['name=','reads_mapping=','peaks_calling=','reads_counting=',
                                       'output_file_name=','outdir=','help'])
        for i,j in opts:   
            if i=="-h" or i=="--help":
                stdout.write(__doc__)
                exit(0)
            elif i=='--name':
                name=j                
            elif i=='--reads_mapping':
                reads_mapping=j
            elif i=='--peaks_calling':
                peaks_calling=j
            elif i=='--reads_counting':
                reads_counting=j
            elif i=='--output_file_name':
                output_file_name=j
            elif i=='--outdir':
                out_dir=j
            else:
                raise Exception("Internal errors occur when parsing command line arguments.")
    except Exception as e:
        stderr.write("%s\n" % e)
        stderr.write("Type 'python generate_report.py --help' for more information.\n")
        exit(1)

    if not out_dir:
        out_dir=os.getcwd()

    if not output_file_name:
        output_file_name='temp'        
    print(name,reads_mapping,peaks_calling,reads_counting,output_file_name,out_dir)
    generate_unique(name,reads_mapping,peaks_calling,reads_counting,output_file_name,out_dir)

def generate_unique(name,reads_mapping,peaks_calling,reads_counting,output_file_name,out_dir):
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
            logo_file='/home/chenhaojie/Pictures/logo.png'
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
        cheque_file='/home/chenhaojie/Pictures/cover.png'
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
    #Quality control for reads mapping
    with doc.create(Section('Quality control for reads mapping')):
        doc.append('The table below shows the number of reads keeped after reads mapping and removing PCR duplicates.')

    with doc.create(LongTabu("X[2.5l] X[2l] X[2l] X[2l]",
                             row_height=1.5)) as data_table:
        data_table.add_hline()
        data_table.add_row(["Sample_ID",
                            "Total number of reads",
                            "Unique mapped reads",
                            "Deduplicated reads"],
                           mapper=bold,color="lightgray")
        data_table.add_hline()

        df=pd.read_csv(reads_mapping,sep='\t')
        print(df)
        for i in df.index:
            if i%2==0:
                data_table.add_row(['%s'%j for j in df.loc[i,:]],
                                   color="white")
            else:
                data_table.add_row(['%s'%j for j in df.loc[i,:]],
                                   color="lightgray")
        data_table.add_hline()       
        
        
    doc.append(NewPage())
    #Quality control for peak calling
    with doc.create(Section('Quality control for peak calling')):
        with doc.create(Subsection('The table below shows the number of peaks and TSS enrichment score in each sample.')):
            with doc.create(LongTabu("X[3l] X[2l] X[2l]",
                                     row_height=1.5)) as data_table:
                data_table.add_hline()
                data_table.add_row(["Sample_ID",
                                    "Total number of peak",
                                    "TSS enrichment score"],
                                   mapper=bold,color="lightgray"
                                   )
                data_table.add_hline()
                df=pd.read_csv(peaks_calling,sep='\t')
                print(df)
                for i in df.index:
                    if i%2==0:
                        data_table.add_row(['%s'%j for j in df.loc[i,:]]+['NA'],
                                           color="white")
                    else:
                        data_table.add_row(['%s'%j for j in df.loc[i,:]]+['NA'],
                                           color="lightgray")
                data_table.add_hline()  
                data_table.add_empty_row()
        with doc.create(Subsection('The figures below show TSS enrichment profiles in each sample.')):
            doc.append('TSS enrichment score is a raio between aggregate distribution of reads centered on TSSs and that flanking the corresponding TSSs.')
            with doc.create(LongTabu("X[2c] X[2c] X[2c]",
                                     row_height=2.5)) as data_table:
                cheque_file = '/home/chenhaojie/Pictures/TSS_score.png'
                cheque = StandAloneGraphic(cheque_file, image_options="width=180px")
                [data_table.add_empty_row() for i in range(2)]
                data_table.add_row([cheque,cheque,cheque])
    
    doc.append(NewPage())              
    #Quality control for reads counting
    with doc.create(Section('Quality control for reads counting')):
        doc.append('The table below shows the number of reads within peaks for each sample.')
        with doc.create(LongTabu("X[2l] X[2l] X[2l]",
                                    row_height=1.5)) as data_table:
            data_table.add_hline()
            data_table.add_row(["Sample_ID",
                                    "Reads within bins",
                                    "Within ratio"],
                                   mapper=bold,color="lightgray")
            data_table.add_hline()
            df=pd.read_csv(reads_counting,sep='\t')
            print(df)
            for i in df.index:
                if i%2==0:
                    data_table.add_row(['%s'%j for j in df.loc[i,:]],color="white")
                else:
                    data_table.add_row(['%s'%j for j in df.loc[i,:]],color="lightgray")
            data_table.add_hline()  
            data_table.add_empty_row()    
            
    doc.append(NewPage())            
    #Downstream analysis
    with doc.create(Section('Downstream analysis')):
        with doc.create(Subsection('The volcano plot below shows the signal changes between two conditions.')):
            with doc.create(LongTabu("X[1c] X[3c] X[1c]",
                                     row_height=2.5)) as data_table:
                cheque_file = '/home/chenhaojie/Pictures/DMPs.png'
                cheque = StandAloneGraphic(cheque_file, image_options="width=200px")
                data_table.add_row(['',cheque,''])            
        with doc.create(Subsection('The scatter plot below shows dimensionality reduction representation of each sample.')):            
            with doc.create(LongTabu("X[1c] X[3c] X[1c]",
                                     row_height=2.5)) as data_table:
                cheque_file = '/home/chenhaojie/Pictures/TSNE.png'
                cheque = StandAloneGraphic(cheque_file, image_options="width=200px")
                data_table.add_row(['',cheque,''])        
        
    doc.append(NewPage()) 
    #References  
    with doc.create(Section('References')):
        doc.append('[1]Fast gapped-read alignment with Bowtie 2.\n[2]Improved peak-calling with MACS2.\n[3]MAnorm2 for quantitatively comparing groups of ChIP-seq samples.\n[4]HyperChIP for identifying hypervariable signals across ChIP/ATAC-seq samples.\n')

    doc.generate_pdf(os.path.join(out_dir,output_file_name), clean_tex=False)


if __name__ == '__main__':
    main()
