import copepodTCR as cpp
import pandas as pd
import numpy as np
import seaborn as sn
from itertools import combinations
from math import comb
import re
from matplotlib import pyplot as plt
import ast
from shiny import App, render, ui, reactive, Inputs, Outputs, Session
import shiny.experimental as x
from shiny.types import ImgData
import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from io import StringIO
import zipfile
from io import BytesIO

question_circle_fill = ui.HTML(
    '<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-question-circle-fill mb-1" viewBox="0 0 16 16"><path d="M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zM5.496 6.033h.825c.138 0 .248-.113.266-.25.09-.656.54-1.134 1.342-1.134.686 0 1.314.343 1.314 1.168 0 .635-.374.927-.965 1.371-.673.489-1.206 1.06-1.168 1.987l.003.217a.25.25 0 0 0 .25.246h.811a.25.25 0 0 0 .25-.25v-.105c0-.718.273-.927 1.01-1.486.609-.463 1.244-.977 1.244-2.056 0-1.511-1.276-2.241-2.673-2.241-1.267 0-2.655.59-2.75 2.286a.237.237 0 0 0 .241.247zm2.325 6.443c.61 0 1.029-.394 1.029-.927 0-.552-.42-.94-1.029-.94-.584 0-1.009.388-1.009.94 0 .533.425.927 1.01.927z"/></svg>'
)

def read_markdown_file(path):
    with open(path, 'r') as file:
        return file.read()

app_ui = ui.page_navbar(
    ui.nav(
            ui.output_image('logo', inline=True),
           ui.output_image('logo2', inline=True),
           ui.markdown(read_markdown_file('/Users/hmeyer/myapp/www/copepodTCR.md')),
           ui.output_image('scheme1', inline=True),
           ui.markdown(read_markdown_file('/Users/hmeyer/myapp/www/instructions.md')),
           ),
    ui.nav(ui.output_image('Instructions', inline=True),
           ui.markdown(read_markdown_file('/Users/hmeyer/myapp/www/instructions_2.md')),
           ui.output_table('simulation_example'),
           ui.markdown(read_markdown_file('/Users/hmeyer/myapp/www/instructions_3.md')),
           ui.output_table('simulation_do_example'),
           ui.markdown(read_markdown_file('/Users/hmeyer/myapp/www/instructions_4.md')),
           ui.output_image('scheme2', inline=True),
           ui.markdown(read_markdown_file('/Users/hmeyer/myapp/www/instructions_5.md')),
           ui.output_image('scheme3', inline=True),
           ui.markdown(read_markdown_file('/Users/hmeyer/myapp/www/instructions_6.md')),
    ),
    ui.nav(ui.output_image('Balance_check', inline=True),
        ui.panel_title("Selection of parameters", "Balance check"),
        x.ui.tooltip(ui.span("Number of pools ", question_circle_fill),
                "Number of pools in the experiment",
                placement="right",
                id="card_tooltip",),
        ui.input_numeric("n_pools1", '', value=9),
        x.ui.tooltip(ui.span("Number of peptides ", question_circle_fill),
                "Number of tested peptides",
                placement="right",
                id="card_tooltip",),
        ui.input_numeric("len_lst1", "", value=100),
        x.ui.tooltip(ui.span("Possible peptide occurence across pools ", question_circle_fill),
                "Number of pools to which a single peptide would be added",
                placement="right",
                id="card_tooltip",),
        ui.output_ui('iters_slider'),
        x.ui.tooltip(ui.span("Resulting balance ", question_circle_fill),
                "Resulting number of unique peptides per pool",
                placement="right",
                id="card_tooltip",), 
        ui.output_text_verbatim('balance'),
        ui.output_plot('balance_plot')
    ),
    ui.nav(ui.output_image('Overlap_check', inline=True),
           ui.panel_title("Check for overlap length of the peptides", "Overlap check"),
           x.ui.tooltip(ui.span("Choose a file with peptides to upload: ", question_circle_fill),
                "TSV or CSV table with peptides in a column without header",
                placement="right",
                id="card_tooltip",),
           ui.input_file("file1", "", multiple=False),
           x.ui.tooltip(ui.span("Or enter them here: ", question_circle_fill),
                "List of peptides separated by comma",
                placement="right",
                id="card_tooltip",),
            ui.input_text_area('text_peptides1', "", 'Peptides', resize = 'both'),
            ui.input_action_button('sample_data1', 'Load sample peptides',
                                   class_="btn", style="color: #015396; background-color: #fff; border-color: #015396"),
            ui.input_action_button('get_overlap', 'Check',
                                   class_="btn-primary", style="color: #fff; background-color: #015396; border-color: #015396"),
            ui.output_text_verbatim('file_content'),
            ui.output_text_verbatim('strange_overlaps')
           ),
    ui.nav(ui.output_image('Pooling_scheme', inline=True),
           ui.panel_title("Get a peptide pooling scheme ", "Pooling scheme"),
           x.ui.tooltip(ui.span("Number of pools ", question_circle_fill),
                "Number of pools in the experiment",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric("n_pools", '', value = 20),
           x.ui.tooltip(ui.span("Peptide occurrence across pools ", question_circle_fill),
                "Number of pools to which a single peptide would be added",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric('iters', '', value = 4),
           x.ui.tooltip(ui.span("Expected epitope length ", question_circle_fill),
                "The length in amino acids of the epitope for simulation step",
                placement="right",
                id="card_tooltip",),
           ui.input_slider('ep_length', '', min = 8, max = 15, value=8),
           x.ui.tooltip(ui.span("Select regime ", question_circle_fill),
                "Account for experimental mistakes (dropped pools) in simulation step or not",
                placement="right",
                id="card_tooltip",),
           ui.input_select('drops', '',
                           {'without dropouts': 'Without drop-outs',
                            'with dropouts': 'With drop-outs'}),
            x.ui.tooltip(ui.span("Choose a file with peptides to upload: ", question_circle_fill),
                "TSV or CSV table with peptides in a column without header",
                placement="right",
                id="card_tooltip",),
           ui.input_file("file2", "",
                         multiple=False),
            x.ui.tooltip(ui.span("Or enter them here: ", question_circle_fill),
                "List of peptides separated by comma",
                placement="right",
                id="card_tooltip",),
            ui.input_text_area('text_peptides2', "", 'Peptides',
                               resize = 'both'),
            ui.input_action_button('sample_data', 'Load sample peptides',
                                   class_="btn", style="color: #015396; background-color: #fff; border-color: #015396"),
            ui.input_text("letter", "Enter your email", placeholder="example@mail.ru"),
            ui.input_action_button('get_results', 'Get scheme and simulation results',
                                   class_="btn-primary", style="color: #fff; background-color: #015396; border-color: #015396"),
            ui.input_action_button("send", "Send results", class_="btn-success", style="color: #015396; background-color: #fff; border-color: #015396"),
            ui.output_text('notification'),
            ui.download_button('download1', 'Download results', class_="btn-success", style="color: #015396; background-color: #fff; border-color: #015396"),
            ui.output_ui('sim_text'),
            ui.output_plot('dropouts'),
            ui.output_data_frame('check_results'),
    ),
    ui.nav(ui.output_image('STL_files', inline=True),
           ui.panel_title("Download zip file with STL ", "STL files"),
           x.ui.tooltip(ui.span("Choose a file with pools: ", question_circle_fill),
                "Produced by copepodTCR in previous steps",
                placement="right",
                id="card_tooltip",),
           ui.input_file('file4', '', multiple=False),
           x.ui.tooltip(ui.span("Choose a file with peptide arrangement in the plate: ", question_circle_fill),
                "TSV table with one peptide per cell, rows and columns are corresponding rows and columns of the plate",
                placement="right",
                id="card_tooltip",),
           ui.input_file('file5', '', multiple=False),
           ui.output_image('scheme4', inline=False),
           x.ui.tooltip(ui.span("Number of rows ", question_circle_fill),
                "Number of rows in the plate",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric("rows", '', value = 16),
           x.ui.tooltip(ui.span("Number of columns ", question_circle_fill),
                "Number of columns in the plate",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric("cols", '', value = 24),
           x.ui.tooltip(ui.span("Length (mm) ", question_circle_fill),
                "Length of the plate in mm",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric("length", '', value = 122.10),
           x.ui.tooltip(ui.span("Width (mm) ", question_circle_fill),
                "Width of the plate in mm",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric("width", '', value = 79.97),
           x.ui.tooltip(ui.span("Thickness (mm) ", question_circle_fill),
                "Thickness of the plate in mm",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric("thickness", '', value = 1.5),
           x.ui.tooltip(ui.span("Hole radius (mm) ", question_circle_fill),
                "Diameter of the hole divided by two",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric("hole_radius", '', value = 2),
           x.ui.tooltip(ui.span("X offset (mm) ", question_circle_fill),
                "Margin along the X axis for the A1 hole, in mm",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric("x_offset", '', value = 9.05),
           x.ui.tooltip(ui.span("Y offset (mm) ", question_circle_fill),
                "Margin along the Y axis for the A1 hole, in mm",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric("y_offset", '', value = 6.20),
           x.ui.tooltip(ui.span("Well spacing (mm) ", question_circle_fill),
                "Distance between wells",
                placement="right",
                id="card_tooltip",),
           ui.input_numeric("well_spacing", '', value = 4.5),
           x.ui.tooltip(ui.span("Add well with coordinates m-k ", question_circle_fill),
                "To better orient tip pattern, you can add the hole with these coordinates, but it should be done only if there is no peptide in the plate in respective well ",
                placement="right",
                id="card_tooltip",),
           ui.input_switch('well24', '', value=True),
           ui.input_action_button('get_stl', 'Generate STL files', class_="btn-primary",
                                  style="color: #fff; background-color: #015396; border-color: #015396"),
           ui.output_text('stl_generation'),
           ui.input_text("stl_letter", "Enter your email", placeholder="example@mail.ru"),
           ui.input_action_button("send_stl", "Send results", class_="btn-success",
                                  style="color: #015396; background-color: #fff; border-color: #015396"),
           ui.output_text('stl_notification'),
           ui.download_button('download_stl', 'Download zip file with STL files',
                              style="color: #015396; background-color: #fff; border-color: #015396")
    ),

    ui.nav(ui.output_image('Interpretation', inline=True),
           ui.panel_title("Analyze results ", "Activated pools"),
           x.ui.tooltip(ui.span("Choose a file with simulation: ", question_circle_fill),
                "Produced during pooling scheme step with or without drop-outs",
                placement="right",
                id="card_tooltip",),
           ui.input_file('file6', '', multiple=False),
            x.ui.tooltip(ui.span("Experiment read out ", question_circle_fill),
                "Numeric results of the experiment, in percentage of activated T cells, a table with columns 'Pool', 'Percentage'",
                placement="right",
                id="card_tooltip",),
            ui.input_file('file7', '', multiple=False),
            x.ui.tooltip(ui.span("Number of pools ", question_circle_fill),
                "Number of pools in the experiment",
                placement="right",
                id="card_tooltip",),
            ui.input_numeric('n_pools3', '', value=0),
            ui.input_action_button('get_model', 'Analyze', class_="btn-primary",
                                   style="color: #fff; background-color: #015396; border-color: #015396"),
            ui.output_plot('fitting'),
            ui.output_text_verbatim('activated_pools'),
            ui.output_text_verbatim('resulting_peptides2'),
            ui.download_button('download_peptide_probs', 'Download results',  style="color: #015396; background-color: #fff; border-color: #015396"),
            ui.input_text("letter2", "Enter your email:", placeholder="example@mail.ru"),
            ui.input_action_button("send_probs", "Send results", class_="btn-success", style="color: #015396; background-color: #fff; border-color: #015396"),
            ui.output_text('send_probs_note'),
            
    )
)

def server(input: Inputs, output: Outputs, session: Session):

    @output
    @render.image
    def logo():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/logo.svg", "width": "130px"}
        return img
    @output
    @render.image
    def logo2():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/logo.svg", "width": "400px"}
        return img
    @output
    @render.image
    def Instructions():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/Instructions.svg", "width": "130px"}
        return img
    @output
    @render.image
    def scheme1():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/scheme1.svg"}
        return img
    @output
    @render.table
    def simulation_example():
        table = pd.read_csv('/Users/hmeyer/myapp/www/simulation_example.tsv', sep = "\t")
        return table
    @output
    @render.table
    def simulation_do_example():
        table = pd.read_csv('/Users/hmeyer/myapp/www/simulation_do_example.tsv', sep = "\t")
        return table
    @output
    @render.image
    def scheme2():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/scheme2.svg"}
        return img
    @output
    @render.image
    def scheme3():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/scheme3.svg"}
        return img
    @output
    @render.image
    def Balance_check():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/Balance_check.svg", "width": "130px"}
        return img
    @output
    @render.image
    def Overlap_check():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/Overlap_check.svg", "width": "130px"}
        return img
    @output
    @render.image
    def Pooling_scheme():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/Pooling_scheme.svg", "width": "130px"}
        return img
    @output
    @render.image
    def STL_files():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/STL_files.svg", "width": "130px"}
        return img
    @output
    @render.image
    def scheme4():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/scheme2.svg"}
        return img
    @output
    @render.image
    def Interpretation():
        img: ImgData = {"src": "/Users/hmeyer/myapp/www/Interpretation.svg", "width": "130px"}
        return img
    @output
    @render.ui
    def iters_slider():
        x1 = input.n_pools1()
        x2 = input.len_lst1()
        if type(x1) != int or type(x2) != int:
            return "Please enter a valid integer"
        range_iters = []
        for i in range(x1):
            if comb(x1, i) >= x2:
                range_iters.append(i)

        if range_iters != []:
            res = ui.input_slider("iters1",
                              "",
                              min=range_iters[0],
                              max=range_iters[-1],
                              value=range_iters[0],
                              step=1)
            return res
        else:
            return "Warning: Increase number of pools"
    @reactive.Calc
    def balance_calc():
        if type(input.n_pools1()) == int and type(input.len_lst1()) == int:
            if comb(input.n_pools1(), input.iters1()) >= input.len_lst1():
                b, lines = cpp.address_rearrangement_AU(n_pools=input.n_pools1(),
                                                iters=input.iters1(),
                                                len_lst=input.len_lst1())
                return b
    @output
    @render.text
    def balance():
        if balance_calc() is not None:
            return balance_calc()
    @output
    @render.plot
    def balance_plot():
        b = balance_calc()
        if b is not None:
            ideal = [int(- (-(input.len_lst1()*input.iters1()/input.n_pools1()) // 1))]*input.n_pools1()
            fig = sn.scatterplot(y=b, x = list(range(1, len(b)+1)), label = 'current balance',
                           color = '#999999ff')
            fig = sn.lineplot(y = ideal, x = list(range(1, len(b)+1)),
                         color = '#0c5394', label = 'perfect balance')
            plt.ylim(0, max(b)+10)
            plt.xlim(0, len(b)+1)
            plt.yticks([int(tick) for tick in plt.yticks()[0]])
            plt.xticks([int(tick) for tick in plt.xticks()[0]][:-1])
            plt.ylabel('Number of peptides per pool')
            plt.xlabel('Pool index')
            plt.title('Peptides per pool')
            plt.legend()
            return fig
    
    @reactive.Effect
    @reactive.event(input.sample_data1)
    def _():
        ui.update_text_area('text_peptides1', value = 'MFVFLVLLPLVSSQCVN, VLLPLVSSQCVNLTTRT, VSSQCVNLTTRTQLPPA, VNLTTRTQLPPAYTNSF, RTQLPPAYTNSFTRGVY, PAYTNSFTRGVYYPDKV, SFTRGVYYPDKVFRSSV, VYYPDKVFRSSVLHSTQ, KVFRSSVLHSTQDLFLP, SVLHSTQDLFLPFFSNV, TQDLFLPFFSNVTWFHA, LPFFSNVTWFHAIHVSG, NVTWFHAIHVSGTNGTK, HAIHVSGTNGTKRFDNP, SGTNGTKRFDNPVLPFN, TKRFDNPVLPFNDGVYF, NPVLPFNDGVYFASTEK, FNDGVYFASTEKSNIIR, YFASTEKSNIIRGWIFG, EKSNIIRGWIFGTTLDS, IRGWIFGTTLDSKTQSL, FGTTLDSKTQSLLIVNN, DSKTQSLLIVNNATNVV, SLLIVNNATNVVIKVCE, NNATNVVIKVCEFQFCN, VVIKVCEFQFCNDPFLG, CEFQFCNDPFLGVYYHK, CNDPFLGVYYHKNNKSW, LGVYYHKNNKSWMESEF, HKNNKSWMESEFRVYSS, SWMESEFRVYSSANNCT, EFRVYSSANNCTFEYVS, SSANNCTFEYVSQPFLM, CTFEYVSQPFLMDLEGK, VSQPFLMDLEGKQGNFK, LMDLEGKQGNFKNLREF, GKQGNFKNLREFVFKNI, FKNLREFVFKNIDGYFK, EFVFKNIDGYFKIYSKH, NIDGYFKIYSKHTPINL, FKIYSKHTPINLVRDLP, KHTPINLVRDLPQGFSA, NLVRDLPQGFSALEPLV, LPQGFSALEPLVDLPIG, SALEPLVDLPIGINITR, LVDLPIGINITRFQTLL, IGINITRFQTLLALHRS, TRFQTLLALHRSYLTPG, LLALHRSYLTPGDSSSG, RSYLTPGDSSSGWTAGA, PGDSSSGWTAGAAAYYV, SGWTAGAAAYYVGYLQP, GAAAYYVGYLQPRTFLL, YVGYLQPRTFLLKYNEN, QPRTFLLKYNENGTITD, LLKYNENGTITDAVDCA, ENGTITDAVDCALDPLS, TDAVDCALDPLSETKCT, CALDPLSETKCTLKSFT, LSETKCTLKSFTVEKGI, CTLKSFTVEKGIYQTSN, FTVEKGIYQTSNFRVQP, GIYQTSNFRVQPTESIV, SNFRVQPTESIVRFPNI, QPTESIVRFPNITNLCP, IVRFPNITNLCPFGEVF, NITNLCPFGEVFNATRF, CPFGEVFNATRFASVYA, VFNATRFASVYAWNRKR, RFASVYAWNRKRISNCV, YAWNRKRISNCVADYSV, KRISNCVADYSVLYNSA, CVADYSVLYNSASFSTF, SVLYNSASFSTFKCYGV, SASFSTFKCYGVSPTKL, TFKCYGVSPTKLNDLCF, GVSPTKLNDLCFTNVYA, KLNDLCFTNVYADSFVI, CFTNVYADSFVIRGDEV, YADSFVIRGDEVRQIAP, VIRGDEVRQIAPGQTGK, EVRQIAPGQTGKIADYN, APGQTGKIADYNYKLPD, GKIADYNYKLPDDFTGC, YNYKLPDDFTGCVIAWN, PDDFTGCVIAWNSNNLD, GCVIAWNSNNLDSKVGG, WNSNNLDSKVGGNYNYL, LDSKVGGNYNYLYRLFR, GGNYNYLYRLFRKSNLK, YLYRLFRKSNLKPFERD, FRKSNLKPFERDISTEI, LKPFERDISTEIYQAGS, RDISTEIYQAGSTPCNG, EIYQAGSTPCNGVEGFN, GSTPCNGVEGFNCYFPL, NGVEGFNCYFPLQSYGF, FNCYFPLQSYGFQPTNG, PLQSYGFQPTNGVGYQP, GFQPTNGVGYQPYRVVV, NGVGYQPYRVVVLSFEL, QPYRVVVLSFELLHAPA, VVLSFELLHAPATVCGP, ELLHAPATVCGPKKSTN, PATVCGPKKSTNLVKNK, GPKKSTNLVKNKCVNFN, TNLVKNKCVNFNFNGLT, NKCVNFNFNGLTGTGVL, FNFNGLTGTGVLTESNK, LTGTGVLTESNKKFLPF, VLTESNKKFLPFQQFGR, NKKFLPFQQFGRDIADT, PFQQFGRDIADTTDAVR, GRDIADTTDAVRDPQTL, DTTDAVRDPQTLEILDI, VRDPQTLEILDITPCSF, TLEILDITPCSFGGVSV, DITPCSFGGVSVITPGT, SFGGVSVITPGTNTSNQ, SVITPGTNTSNQVAVLY, GTNTSNQVAVLYQDVNC, NQVAVLYQDVNCTEVPV, LYQDVNCTEVPVAIHAD, NCTEVPVAIHADQLTPT, PVAIHADQLTPTWRVYS, ADQLTPTWRVYSTGSNV, PTWRVYSTGSNVFQTRA, YSTGSNVFQTRAGCLIG, NVFQTRAGCLIGAEHVN, RAGCLIGAEHVNNSYEC, IGAEHVNNSYECDIPIG, VNNSYECDIPIGAGICA, ECDIPIGAGICASYQTQ, IGAGICASYQTQTNSPR, CASYQTQTNSPRRARSV, TQTNSPRRARSVASQSI, PRRARSVASQSIIAYTM, SVASQSIIAYTMSLGAE, SIIAYTMSLGAENSVAY, TMSLGAENSVAYSNNSI, AENSVAYSNNSIAIPTN, AYSNNSIAIPTNFTISV, SIAIPTNFTISVTTEIL, TNFTISVTTEILPVSMT, SVTTEILPVSMTKTSVD, ILPVSMTKTSVDCTMYI, MTKTSVDCTMYICGDST, VDCTMYICGDSTECSNL, YICGDSTECSNLLLQYG, STECSNLLLQYGSFCTQ, NLLLQYGSFCTQLNRAL, YGSFCTQLNRALTGIAV, TQLNRALTGIAVEQDKN, ALTGIAVEQDKNTQEVF, AVEQDKNTQEVFAQVKQ, KNTQEVFAQVKQIYKTP, VFAQVKQIYKTPPIKDF, KQIYKTPPIKDFGGFNF, TPPIKDFGGFNFSQILP, DFGGFNFSQILPDPSKP, NFSQILPDPSKPSKRSF, LPDPSKPSKRSFIEDLL, KPSKRSFIEDLLFNKVT, SFIEDLLFNKVTLADAG, LLFNKVTLADAGFIKQY, VTLADAGFIKQYGDCLG, AGFIKQYGDCLGDIAAR, QYGDCLGDIAARDLICA, LGDIAARDLICAQKFNG, ARDLICAQKFNGLTVLP, CAQKFNGLTVLPPLLTD, NGLTVLPPLLTDEMIAQ, LPPLLTDEMIAQYTSAL, TDEMIAQYTSALLAGTI, AQYTSALLAGTITSGWT, ALLAGTITSGWTFGAGA, TITSGWTFGAGAALQIP, WTFGAGAALQIPFAMQM, GAALQIPFAMQMAYRFN, IPFAMQMAYRFNGIGVT, QMAYRFNGIGVTQNVLY, FNGIGVTQNVLYENQKL, VTQNVLYENQKLIANQF, LYENQKLIANQFNSAIG, KLIANQFNSAIGKIQDS, QFNSAIGKIQDSLSSTA, IGKIQDSLSSTASALGK, DSLSSTASALGKLQDVV, TASALGKLQDVVNQNAQ, GKLQDVVNQNAQALNTL, VVNQNAQALNTLVKQLS, AQALNTLVKQLSSNFGA, TLVKQLSSNFGAISSVL, LSSNFGAISSVLNDILS, GAISSVLNDILSRLDKV, VLNDILSRLDKVEAEVQ, LSRLDKVEAEVQIDRLI, KVEAEVQIDRLITGRLQ, VQIDRLITGRLQSLQTY, LITGRLQSLQTYVTQQL, LQSLQTYVTQQLIRAAE, TYVTQQLIRAAEIRASA, QLIRAAEIRASANLAAT, AEIRASANLAATKMSEC, SANLAATKMSECVLGQS, ATKMSECVLGQSKRVDF, ECVLGQSKRVDFCGKGY, QSKRVDFCGKGYHLMSF, DFCGKGYHLMSFPQSAP, GYHLMSFPQSAPHGVVF, SFPQSAPHGVVFLHVTY, APHGVVFLHVTYVPAQE, VFLHVTYVPAQEKNFTT, TYVPAQEKNFTTAPAIC, QEKNFTTAPAICHDGKA, TTAPAICHDGKAHFPRE, ICHDGKAHFPREGVFVS, KAHFPREGVFVSNGTHW, REGVFVSNGTHWFVTQR, VSNGTHWFVTQRNFYEP, HWFVTQRNFYEPQIITT, QRNFYEPQIITTDNTFV, EPQIITTDNTFVSGNCD, TTDNTFVSGNCDVVIGI, FVSGNCDVVIGIVNNTV, CDVVIGIVNNTVYDPLQ, GIVNNTVYDPLQPELDS, TVYDPLQPELDSFKEEL, LQPELDSFKEELDKYFK, DSFKEELDKYFKNHTSP, ELDKYFKNHTSPDVDLG, FKNHTSPDVDLGDISGI, SPDVDLGDISGINASVV, LGDISGINASVVNIQKE, GINASVVNIQKEIDRLN, VVNIQKEIDRLNEVAKN, KEIDRLNEVAKNLNESL, LNEVAKNLNESLIDLQE, KNLNESLIDLQELGKYE, SLIDLQELGKYEQYIKW, QELGKYEQYIKWPWYIW, YEQYIKWPWYIWLGFIA, KWPWYIWLGFIAGLIAI, IWLGFIAGLIAIVMVTI, IAGLIAIVMVTIMLCCM, AIVMVTIMLCCMTSCCS, TIMLCCMTSCCSCLKGC, CMTSCCSCLKGCCSCGS, CSCLKGCCSCGSCCKFD, GCCSCGSCCKFDEDDSE, GSCCKFDEDDSEPVLKG, FDEDDSEPVLKGVKLHY')

        
    @output
    @render.text
    @reactive.event(input.get_overlap, ignore_none=False)
    def file_content():
        file_infos = input.file1()
        if not file_infos:
            lst = re.split(r'\W+', input.text_peptides1())
        else:
            for file_info in file_infos:
                input_df=pd.read_csv(file_info["datapath"], header = None)
                lst = list(input_df[0])
        overlaps = cpp.all_overlaps(lst)
        overlaps_text = ''
        for key in overlaps:
            t = str(overlaps[key]+1) + ' peptides with overlap ' + str(key) +  ' amino acids' + '\n'
            overlaps_text = overlaps_text + t

        if len(overlaps) > 1:
            m = ui.modal(
            "Peptide set with different overlaps hampens results interpretation.",
            title="Warning",
            easy_close=True,
            footer=None)
            ui.modal_show(m)

        return overlaps_text
    @output
    @render.text
    @reactive.event(input.get_overlap, ignore_none=False)
    def strange_overlaps():
        file_infos = input.file1()
        if not file_infos:
            lst = re.split(r'\W+', input.text_peptides1())
        else:
            for file_info in file_infos:
                input_df=pd.read_csv(file_info["datapath"], header = None)
                lst = list(input_df[0])
        overlaps = cpp.all_overlaps(lst)
        strange_overlaps = [i for i in overlaps if overlaps[i]!=max(overlaps.values())]
        pair_text = ''
        for overlap in strange_overlaps:
            pair = cpp.find_pair_with_overlap(lst, overlap)
            pair_text = pair_text + 'Overlap of ' + str(overlap) + ' amino acids:' + '\n'
            for p in pair:
                for peptide in p:
                    pair_text = pair_text + peptide + ' '
                pair_text = pair_text + '\n'
        return pair_text
    
    @reactive.Effect
    @reactive.event(input.sample_data)
    def _():
        ui.update_text_area('text_peptides2', value = 'MFVFLVLLPLVSSQCVN, VLLPLVSSQCVNLTTRT, VSSQCVNLTTRTQLPPA, VNLTTRTQLPPAYTNSF, RTQLPPAYTNSFTRGVY, PAYTNSFTRGVYYPDKV, SFTRGVYYPDKVFRSSV, VYYPDKVFRSSVLHSTQ, KVFRSSVLHSTQDLFLP, SVLHSTQDLFLPFFSNV, TQDLFLPFFSNVTWFHA, LPFFSNVTWFHAIHVSG, NVTWFHAIHVSGTNGTK, HAIHVSGTNGTKRFDNP, SGTNGTKRFDNPVLPFN, TKRFDNPVLPFNDGVYF, NPVLPFNDGVYFASTEK, FNDGVYFASTEKSNIIR, YFASTEKSNIIRGWIFG, EKSNIIRGWIFGTTLDS, IRGWIFGTTLDSKTQSL, FGTTLDSKTQSLLIVNN, DSKTQSLLIVNNATNVV, SLLIVNNATNVVIKVCE, NNATNVVIKVCEFQFCN, VVIKVCEFQFCNDPFLG, CEFQFCNDPFLGVYYHK, CNDPFLGVYYHKNNKSW, LGVYYHKNNKSWMESEF, HKNNKSWMESEFRVYSS, SWMESEFRVYSSANNCT, EFRVYSSANNCTFEYVS, SSANNCTFEYVSQPFLM, CTFEYVSQPFLMDLEGK, VSQPFLMDLEGKQGNFK, LMDLEGKQGNFKNLREF, GKQGNFKNLREFVFKNI, FKNLREFVFKNIDGYFK, EFVFKNIDGYFKIYSKH, NIDGYFKIYSKHTPINL, FKIYSKHTPINLVRDLP, KHTPINLVRDLPQGFSA, NLVRDLPQGFSALEPLV, LPQGFSALEPLVDLPIG, SALEPLVDLPIGINITR, LVDLPIGINITRFQTLL, IGINITRFQTLLALHRS, TRFQTLLALHRSYLTPG, LLALHRSYLTPGDSSSG, RSYLTPGDSSSGWTAGA, PGDSSSGWTAGAAAYYV, SGWTAGAAAYYVGYLQP, GAAAYYVGYLQPRTFLL, YVGYLQPRTFLLKYNEN, QPRTFLLKYNENGTITD, LLKYNENGTITDAVDCA, ENGTITDAVDCALDPLS, TDAVDCALDPLSETKCT, CALDPLSETKCTLKSFT, LSETKCTLKSFTVEKGI, CTLKSFTVEKGIYQTSN, FTVEKGIYQTSNFRVQP, GIYQTSNFRVQPTESIV, SNFRVQPTESIVRFPNI, QPTESIVRFPNITNLCP, IVRFPNITNLCPFGEVF, NITNLCPFGEVFNATRF, CPFGEVFNATRFASVYA, VFNATRFASVYAWNRKR, RFASVYAWNRKRISNCV, YAWNRKRISNCVADYSV, KRISNCVADYSVLYNSA, CVADYSVLYNSASFSTF, SVLYNSASFSTFKCYGV, SASFSTFKCYGVSPTKL, TFKCYGVSPTKLNDLCF, GVSPTKLNDLCFTNVYA, KLNDLCFTNVYADSFVI, CFTNVYADSFVIRGDEV, YADSFVIRGDEVRQIAP, VIRGDEVRQIAPGQTGK, EVRQIAPGQTGKIADYN, APGQTGKIADYNYKLPD, GKIADYNYKLPDDFTGC, YNYKLPDDFTGCVIAWN, PDDFTGCVIAWNSNNLD, GCVIAWNSNNLDSKVGG, WNSNNLDSKVGGNYNYL, LDSKVGGNYNYLYRLFR, GGNYNYLYRLFRKSNLK, YLYRLFRKSNLKPFERD, FRKSNLKPFERDISTEI, LKPFERDISTEIYQAGS, RDISTEIYQAGSTPCNG, EIYQAGSTPCNGVEGFN, GSTPCNGVEGFNCYFPL, NGVEGFNCYFPLQSYGF, FNCYFPLQSYGFQPTNG, PLQSYGFQPTNGVGYQP, GFQPTNGVGYQPYRVVV, NGVGYQPYRVVVLSFEL, QPYRVVVLSFELLHAPA, VVLSFELLHAPATVCGP, ELLHAPATVCGPKKSTN, PATVCGPKKSTNLVKNK, GPKKSTNLVKNKCVNFN, TNLVKNKCVNFNFNGLT, NKCVNFNFNGLTGTGVL, FNFNGLTGTGVLTESNK, LTGTGVLTESNKKFLPF, VLTESNKKFLPFQQFGR, NKKFLPFQQFGRDIADT, PFQQFGRDIADTTDAVR, GRDIADTTDAVRDPQTL, DTTDAVRDPQTLEILDI, VRDPQTLEILDITPCSF, TLEILDITPCSFGGVSV, DITPCSFGGVSVITPGT, SFGGVSVITPGTNTSNQ, SVITPGTNTSNQVAVLY, GTNTSNQVAVLYQDVNC, NQVAVLYQDVNCTEVPV, LYQDVNCTEVPVAIHAD, NCTEVPVAIHADQLTPT, PVAIHADQLTPTWRVYS, ADQLTPTWRVYSTGSNV, PTWRVYSTGSNVFQTRA, YSTGSNVFQTRAGCLIG, NVFQTRAGCLIGAEHVN, RAGCLIGAEHVNNSYEC, IGAEHVNNSYECDIPIG, VNNSYECDIPIGAGICA, ECDIPIGAGICASYQTQ, IGAGICASYQTQTNSPR, CASYQTQTNSPRRARSV, TQTNSPRRARSVASQSI, PRRARSVASQSIIAYTM, SVASQSIIAYTMSLGAE, SIIAYTMSLGAENSVAY, TMSLGAENSVAYSNNSI, AENSVAYSNNSIAIPTN, AYSNNSIAIPTNFTISV, SIAIPTNFTISVTTEIL, TNFTISVTTEILPVSMT, SVTTEILPVSMTKTSVD, ILPVSMTKTSVDCTMYI, MTKTSVDCTMYICGDST, VDCTMYICGDSTECSNL, YICGDSTECSNLLLQYG, STECSNLLLQYGSFCTQ, NLLLQYGSFCTQLNRAL, YGSFCTQLNRALTGIAV, TQLNRALTGIAVEQDKN, ALTGIAVEQDKNTQEVF, AVEQDKNTQEVFAQVKQ, KNTQEVFAQVKQIYKTP, VFAQVKQIYKTPPIKDF, KQIYKTPPIKDFGGFNF, TPPIKDFGGFNFSQILP, DFGGFNFSQILPDPSKP, NFSQILPDPSKPSKRSF, LPDPSKPSKRSFIEDLL, KPSKRSFIEDLLFNKVT, SFIEDLLFNKVTLADAG, LLFNKVTLADAGFIKQY, VTLADAGFIKQYGDCLG, AGFIKQYGDCLGDIAAR, QYGDCLGDIAARDLICA, LGDIAARDLICAQKFNG, ARDLICAQKFNGLTVLP, CAQKFNGLTVLPPLLTD, NGLTVLPPLLTDEMIAQ, LPPLLTDEMIAQYTSAL, TDEMIAQYTSALLAGTI, AQYTSALLAGTITSGWT, ALLAGTITSGWTFGAGA, TITSGWTFGAGAALQIP, WTFGAGAALQIPFAMQM, GAALQIPFAMQMAYRFN, IPFAMQMAYRFNGIGVT, QMAYRFNGIGVTQNVLY, FNGIGVTQNVLYENQKL, VTQNVLYENQKLIANQF, LYENQKLIANQFNSAIG, KLIANQFNSAIGKIQDS, QFNSAIGKIQDSLSSTA, IGKIQDSLSSTASALGK, DSLSSTASALGKLQDVV, TASALGKLQDVVNQNAQ, GKLQDVVNQNAQALNTL, VVNQNAQALNTLVKQLS, AQALNTLVKQLSSNFGA, TLVKQLSSNFGAISSVL, LSSNFGAISSVLNDILS, GAISSVLNDILSRLDKV, VLNDILSRLDKVEAEVQ, LSRLDKVEAEVQIDRLI, KVEAEVQIDRLITGRLQ, VQIDRLITGRLQSLQTY, LITGRLQSLQTYVTQQL, LQSLQTYVTQQLIRAAE, TYVTQQLIRAAEIRASA, QLIRAAEIRASANLAAT, AEIRASANLAATKMSEC, SANLAATKMSECVLGQS, ATKMSECVLGQSKRVDF, ECVLGQSKRVDFCGKGY, QSKRVDFCGKGYHLMSF, DFCGKGYHLMSFPQSAP, GYHLMSFPQSAPHGVVF, SFPQSAPHGVVFLHVTY, APHGVVFLHVTYVPAQE, VFLHVTYVPAQEKNFTT, TYVPAQEKNFTTAPAIC, QEKNFTTAPAICHDGKA, TTAPAICHDGKAHFPRE, ICHDGKAHFPREGVFVS, KAHFPREGVFVSNGTHW, REGVFVSNGTHWFVTQR, VSNGTHWFVTQRNFYEP, HWFVTQRNFYEPQIITT, QRNFYEPQIITTDNTFV, EPQIITTDNTFVSGNCD, TTDNTFVSGNCDVVIGI, FVSGNCDVVIGIVNNTV, CDVVIGIVNNTVYDPLQ, GIVNNTVYDPLQPELDS, TVYDPLQPELDSFKEEL, LQPELDSFKEELDKYFK, DSFKEELDKYFKNHTSP, ELDKYFKNHTSPDVDLG, FKNHTSPDVDLGDISGI, SPDVDLGDISGINASVV, LGDISGINASVVNIQKE, GINASVVNIQKEIDRLN, VVNIQKEIDRLNEVAKN, KEIDRLNEVAKNLNESL, LNEVAKNLNESLIDLQE, KNLNESLIDLQELGKYE, SLIDLQELGKYEQYIKW, QELGKYEQYIKWPWYIW, YEQYIKWPWYIWLGFIA, KWPWYIWLGFIAGLIAI, IWLGFIAGLIAIVMVTI, IAGLIAIVMVTIMLCCM, AIVMVTIMLCCMTSCCS, TIMLCCMTSCCSCLKGC, CMTSCCSCLKGCCSCGS, CSCLKGCCSCGSCCKFD, GCCSCGSCCKFDEDDSE, GSCCKFDEDDSEPVLKG, FDEDDSEPVLKGVKLHY')
    
    @reactive.Calc
    @reactive.event(input.get_results)
    def pooling_scheme_react():
        file_infos = input.file2()
        if not file_infos:
            lst = re.split(r'\W+', input.text_peptides2())
        else:
            for file_info in file_infos:
                input_df=pd.read_csv(file_info["datapath"], header = None)
            lst = list(input_df[0])
        b, lines = cpp.address_rearrangement_AU(n_pools=input.n_pools(),
                                                iters=input.iters(),
                                                len_lst=len(lst))
        pools, peptide_address = cpp.pooling(lst=lst, addresses=lines, n_pools=input.n_pools())

        act_profile = cpp.epitope_pools_activation(peptide_address, lst, input.ep_length())
    
        check_results = pd.DataFrame(columns = ['Peptide', 'Address', 'Epitope', 'Act Pools',
                                        '# of pools', '# of epitopes', '# of peptides', 'Remained', '# of lost',
                                           'Right peptide', 'Right epitope'])

        with ui.Progress(min=1, max=len(lst)) as p:
            p.set(message="Simulating...", detail="This may take a while...")

            for peptide in lst:
                p.set(lst.index(peptide), message="Computing")
                for i in range(len(peptide)):
                    ep = peptide[i:i+input.ep_length()]
                    if len(ep) == input.ep_length():
                        act = cpp.pools_activation(pools, ep)
                        if input.drops() == 'without dropouts':
                            peps, eps = cpp.peptide_search(lst=lst, act_profile=act_profile,
                                           act_pools = act,
                                           iters = input.iters(), n_pools = input.n_pools(),
                                           regime = 'without dropouts')
                            right_pep = str(peptide in peps)
                            right_ep = str(ep in eps)
                            row = {'Peptide':peptide, 'Address':str(peptide_address[peptide]), 'Epitope':ep,
                           'Act Pools':str(sorted(list(act))), '# of pools':len(act),
                           '# of epitopes':len(eps), '# of peptides':len(peps), 'Remained':'-', '# of lost':0,
                           'Right peptide':right_pep, 'Right epitope':right_ep}
                            check_results = pd.concat([check_results, pd.DataFrame(row, index = [0])])
                        elif input.drops() == 'with dropouts':
                            l = len(act)
                            for i in range(1, l+1):
                                lost = len(act) - i
                                lost_combs = list(combinations(act, i))
                                for lost_comb in lost_combs:
                                    peps, eps = cpp.peptide_search(lst=lst, act_profile=act_profile,
                                           act_pools = list(lost_comb),
                                           iters = input.iters(), n_pools = input.n_pools(),
                                           regime = 'with dropouts')
                                    right_pep = str(peptide in peps)
                                    right_ep = str(ep in eps)
                
                                    row = {'Peptide':peptide, 'Address':str(peptide_address[peptide]), 'Epitope':ep,
                                   'Act Pools':str(sorted(list(act))), '# of pools':len(act),
                                   '# of epitopes':len(eps), '# of peptides':len(peps),
                                   'Remained':str(list(lost_comb)), '# of lost':lost,
                                   'Right peptide':right_pep, 'Right epitope':right_ep}
                                    check_results = pd.concat([check_results, pd.DataFrame(row, index = [0])])

            
        #check_results = cpp.run_experiment(lst=lst, peptide_address=peptide_address,
        #                                   ep_length=input.ep_length(),
        #                                   pools=pools, iters=input.iters(), n_pools=input.n_pools(),
        #                                   regime=input.drops())
        
        if len(check_results[check_results['# of peptides'] > 2]) > 0 and input.drops() == 'without dropouts':
            m = ui.modal(
            "Your peptide set has different overlap lengths, it might hampen result interpretation.",
            title="Warning",
            easy_close=True,
            footer=None)
            ui.modal_show(m)
        
        return [b, pools, peptide_address, check_results, lst]
    
    @output
    @render.text()
    @reactive.event(input.send)
    def notification():
        pools = pooling_scheme_react()[1]
        pools_df = dict()
        for key in pools:
            pools_df[key] = ';'.join(pools[key])
        pools_df = pd.DataFrame.from_dict(pools_df, orient = 'index')
        pools_df = pools_df.reset_index()
        pools_df.columns = ['Pool', 'Peptides']

        peptide_address = pooling_scheme_react()[2]
        peptide_address_df = pd.DataFrame(list(peptide_address.items()), columns=['Peptide', 'Address'])

        check_results_df = pooling_scheme_react()[3]

        send_dataframe_via_email([input.letter()], [pools_df, peptide_address_df, check_results_df,
                                                    [input.n_pools(), input.iters(), len(pooling_scheme_react()[4])]])
        return 'Sent!'
    
    @output
    @render.ui
    @reactive.event(input.get_results)
    def sim_text():
        res = x.ui.tooltip(ui.span("Simulation results ", question_circle_fill),
                "Results of the experiment for any possible epitope of provided length",
                placement="right",
                id="card_tooltip",)
        return res
    
    @output
    @render.plot
    @reactive.event(input.get_results)
    def dropouts():
        check_results = pooling_scheme_react()[3]
        
        if input.drops() == 'with dropouts':
            sn.lineplot(data = check_results, x = '# of lost', y = '# of peptides', hue = 'Right peptide',
                        legend = True)
            fig1 = plt.gca()
            fig1.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        else:
            sn.histplot(check_results, x="# of peptides", hue="Right peptide", legend = True, binwidth=0.1)
            fig1 = plt.gca()
            fig1.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        return fig1

    @output
    @render.data_frame
    @reactive.event(input.get_results)
    def check_results():
        table = pooling_scheme_react()[3]
        fr_table = int(np.floor(len(table)*0.03))
        return table.sample(frac = 1).head(fr_table)
    @session.download(filename='Results.zip')
    def download1():
        pools = pooling_scheme_react()[1]
        pools_df = dict()
        for key in pools:
            pools_df[key] = ';'.join(pools[key])
        pools_df = pd.DataFrame.from_dict(pools_df, orient = 'index')
        pools_df = pools_df.reset_index()
        pools_df.columns = ['Pool', 'Peptides']

        peptide_address = pooling_scheme_react()[2]
        peptide_address_df = pd.DataFrame(list(peptide_address.items()), columns=['Peptide', 'Address'])

        table = pooling_scheme_react()[3]

        zip_buffer = dataframes_to_zip([pools_df, peptide_address_df, table],
                          ['pools.tsv', 'peptide_addresses.tsv', 'simulation.tsv'])
        yield zip_buffer.getvalue()
    
    @reactive.Calc
    @reactive.event(input.get_stl)
    def stl_files_react():
        file_infos = input.file4()
        if not file_infos:
            return
        for file_info in file_infos:
            pools_stl=pd.read_csv(file_info["datapath"], sep = "\t", index_col = 'Pool')
        file_infos = input.file5()
        if not file_infos:
            return
        for file_info in file_infos:
            peptides_table_stl=pd.read_csv(file_info["datapath"], sep = "\t", index_col = 0)
        
        with ui.Progress(min=1, max=len(pools_stl.index)+1) as p:
            p.set(message="Generating...", detail="This may take a while...")
            count_stl = 1
            meshes_list = dict()

            for pool_N in set(pools_stl.index):
                p.set(count_stl, message="Rendering")
                coordinates = []
                for peptide in pools_stl['Peptides'].iloc[pool_N].split(';'):
                    row_value = int([(x, peptides_table_stl.columns[y]) for x, y in zip(*np.where(peptides_table_stl.values == peptide))][0][0]+1)
                    column_value = int([(x, peptides_table_stl.columns[y]) for x, y in zip(*np.where(peptides_table_stl.values == peptide))][0][1])
                    coordinates.append([row_value, column_value])

                if input.well24():
                    coordinates = coordinates + [[16, 24]]
        
                name = 'pool' + str(pool_N)
        
                m = cpp.stl_generator(input.rows(), input.cols(), input.length(), input.width(),
                                      input.thickness(), input.hole_radius(), input.x_offset(),
                                      input.y_offset(), input.well_spacing(), coordinates)
                
                meshes_list[name] = m
                count_stl = count_stl +1

        return meshes_list
    
    @output
    @render.text()
    @reactive.event(input.get_stl)
    def stl_generation():
        meshes_list = stl_files_react()
        return 'Generation is complete'
    
    @output
    @render.text()
    @reactive.event(input.send_stl)
    def stl_notification():
        zip_buffer = cpp.zip_meshes(stl_files_react())
        send_zip_via_email(receivers_list=[input.stl_letter()], zip_buffer = zip_buffer)
        return 'Sent!'

    @output
    @session.download(filename='STL_files.zip')
    def download_stl():
        zip_buffer = cpp.zip_meshes(stl_files_react())
        yield zip_buffer.getvalue()
    
    @reactive.Calc
    @reactive.event(input.get_model)
    def model_fitting():
        file_infos = input.file7()
        if not file_infos:
            return
        for file_info in file_infos:
            exp_results=pd.read_csv(file_info["datapath"])

        with ui.Progress(min=1, max=2) as p:
            p.set(message="Model fitting started.", detail="It will take some time...")
            cells = list(exp_results['Percentage'])
            #cells = cells/np.max(cells)
            inds = list(exp_results['Pool'])
            fig, probs = cpp.activation_model(cells, input.n_pools3(), inds)
            progress = 2

        return fig, probs
    
    @output
    @render.plot
    def fitting():
        fig, _ = model_fitting()
        return fig
    
    @reactive.Calc
    def peptide_probs_react():
        file_infos = input.file6()
        if not file_infos:
            return
        for file_info in file_infos:
            sim=pd.read_csv(file_info["datapath"], sep = "\t")
        _, probs = model_fitting()
        peptide_probs = cpp.peptide_probabilities(sim, probs)
        return peptide_probs

    @output
    @render.text()
    def activated_pools():
        act_pools = []
        _, probs = model_fitting()
        for i in range(len(probs)):
            if probs['assign'].iloc[i] < 0.5:
                act_pools.append(i)
        text = 'Activated pools: '
        for item in act_pools:
            text = text + str(item) + ', '
        return text[:-2]

    @output
    @render.text()
    def resulting_peptides2():
        file_infos = input.file6()
        if not file_infos:
            return
        for file_info in file_infos:
            sim=pd.read_csv(file_info["datapath"], sep = "\t")
        _, probs = model_fitting()
        peptide_probs = peptide_probs_react()
        note, most, possible = cpp.results_analysis(peptide_probs, probs, sim)
        if len(most) != 0:
            text = note + '\nResulting peptides: '
            for p in most:
                text = text + p + ', '
            text = text[:-2]
            text = text + '\n\nOther possible peptides: '
            for p in possible:
                text = text + p + ', '
        elif len(possible) != 0 and len(most) == 0:
            text = note
            text = text + '\n\nPossible peptides: '
            for p in possible:
                text = text + p + ', '
        else:
            text = 'Nothing was found. '
        return text[:-2]
    
    @session.download(filename='Peptide_probs.tsv')
    def download_peptide_probs():
        peptide_probs = peptide_probs_react()
        peptide_probs = peptide_probs.sort_values(by = 'Probability', ascending=False)
        yield peptide_probs.to_csv(sep = "\t", index=False)
    
    @output
    @render.text()
    @reactive.event(input.send_probs)
    def send_probs_note():
        file_infos = input.file6()
        if not file_infos:
            return
        for file_info in file_infos:
            sim=pd.read_csv(file_info["datapath"], sep = "\t")
        _, probs = model_fitting()
        peptide_probs = peptide_probs_react()
        note, most, possible = cpp.results_analysis(peptide_probs, probs, sim)
        if len(most) != 0:
            text = note + '\nResulting peptides: '
            for p in most:
                text = text + p + ', '
            text = text[:-2]
            text = text + '\n\nOther possible peptides: '
            for p in possible:
                text = text + p + ', '
        elif len(possible) != 0:
            text = note
            text = text + '\n\nPossible peptides: '
            for p in possible:
                text = text + p + ', '
        else:
            text = 'Nothing was found. '
        peptide_probs = peptide_probs.sort_values(by = 'Probability', ascending=False)

        send_probs_via_email([input.letter2()], peptide_probs, text[:-2])
        return 'Sent!'
        

app = App(app_ui, server, debug=True)

def send_dataframe_via_email(receivers_list, dataframe):
    lines = []
    with open('/Users/hmeyer/myapp/pass.txt') as f:
        for line in f:
            lines.append(line.strip())
    password = lines[0]
    login = lines[1]
    host_address = lines[2]
    channel = int(lines[3])

    with smtplib.SMTP(host=host_address, port=channel) as server:
        server.login(login, password)
        
        msg = MIMEMultipart()
        msg['From'] = login
        msg['To'] = ', '.join(receivers_list)
        msg['Subject'] = 'copepodTCR: combinatorial peptide pooling scheme'

        # Pools
        csv_buffer = StringIO()
        dataframe[0].to_csv(csv_buffer, index=False, sep = "\t")
        csv_str = csv_buffer.getvalue()
        part = MIMEBase('text', 'tsv')
        part.set_payload(csv_str)
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="pools.tsv"')
        msg.attach(part)

        # Peptide_address
        csv_buffer = StringIO()
        dataframe[1].to_csv(csv_buffer, index=False, sep = "\t")
        csv_str = csv_buffer.getvalue()
        part = MIMEBase('text', 'tsv')
        part.set_payload(csv_str)
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="peptide_address.tsv"')
        msg.attach(part)
        
        # Simulation
        csv_buffer = StringIO()
        dataframe[2].to_csv(csv_buffer, index=False, sep = "\t")
        csv_str = csv_buffer.getvalue()
        part = MIMEBase('text', 'tsv')
        part.set_payload(csv_str)
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="simulation_results.tsv"')
        msg.attach(part)

        # README.txt
        n_pools = dataframe[3][0]
        iters = dataframe[3][1]
        l = dataframe[3][2]
        readme_content = f"n_pools = {n_pools}\npeptide occurence across pools (to how many pools one peptide is added) = {iters}\nnumber of peptides = {l}\nThank you for using copepodTCR! Made by CSHL."
        part = MIMEBase('text', 'plain')
        part.set_payload(readme_content)
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="README.txt"')
        msg.attach(part)
        
        server.sendmail(
            from_addr=login,
            to_addrs=receivers_list,
            msg=msg.as_string()
        )

def send_probs_via_email(receivers_list, dataframe, message):
    lines = []
    with open('/Users/hmeyer/myapp/pass.txt') as f:
        for line in f:
            lines.append(line.strip())
    password = lines[0]
    login = lines[1]
    host_address = lines[2]
    channel = int(lines[3])

    with smtplib.SMTP(host=host_address, port=channel) as server:
        server.login(login, password)
        
        msg = MIMEMultipart()
        msg['From'] = login
        msg['To'] = ', '.join(receivers_list)
        msg['Subject'] = 'copepodTCR: results interpretation'

        # Peptide probabilities
        csv_buffer = StringIO()
        dataframe.to_csv(csv_buffer, index=False, sep = "\t")
        csv_str = csv_buffer.getvalue()
        part = MIMEBase('text', 'tsv')
        part.set_payload(csv_str)
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="Peptide_probabilities.tsv"')
        msg.attach(part)

        # README.txt
        readme_content = f"{message}\nThank you for using copepodTCR! Made by CSHL."
        part = MIMEBase('text', 'plain')
        part.set_payload(readme_content)
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="README.txt"')
        msg.attach(part)
        
        server.sendmail(
            from_addr=login,
            to_addrs=receivers_list,
            msg=msg.as_string()
        )

def dataframes_to_zip(dataframes, filenames):
    zip_buffer = BytesIO()

    with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED) as zip_file:
        for dataframe, filename in zip(dataframes, filenames):
            df_buffer = BytesIO()
            dataframe.to_csv(df_buffer, index=False, sep='\t')
            zip_file.writestr(filename, df_buffer.getvalue())
    zip_buffer.seek(0)

    return zip_buffer

def send_zip_via_email(receivers_list, zip_buffer):
    
    lines = []
    with open('/Users/hmeyer/myapp/pass.txt') as f:
        for line in f:
            lines.append(line.strip())
    password = lines[0]
    login = lines[1]
    host_address = lines[2]
    channel = int(lines[3])

    with smtplib.SMTP(host=host_address, port=channel) as server:
        #server.starttls()
        server.login(login, password)
        
        msg = MIMEMultipart()
        msg['From'] = login
        msg['To'] = ', '.join(receivers_list)
        msg['Subject'] = 'copepodTCR: STL files'

        # Zips
        attachment_part = MIMEBase("application", "zip")
        attachment_part.set_payload(zip_buffer.read())
        encoders.encode_base64(attachment_part)
        attachment_part.add_header("Content-Disposition", 'attachment; filename="STL_files.zip"')
        
        msg.attach(attachment_part)
        
        server.sendmail(
            from_addr=login,
            to_addrs=receivers_list,
            msg=msg.as_string()
        )