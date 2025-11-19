import copepodTCR as cpp
import codepub as cdp
from itertools import combinations
from collections import Counter
from math import comb
from math import ceil

import pandas as pd
import numpy as np
from Bio import SeqIO

#import pymc as pm
#import arviz as az

import seaborn as sn
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import plotly.graph_objects as go

import re

from shiny import App, render, ui, reactive, Inputs, Outputs, Session
#import shiny.experimental as x
from shiny.types import ImgData
from shinywidgets import output_widget
from shinywidgets import render_plotly

import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from io import StringIO
import zipfile
from io import BytesIO
import ssl

path_user = ''

question_circle_fill = ui.HTML(
    '<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-question-circle-fill mb-1" viewBox="0 0 16 16"><path d="M16 8A8 8 0 1 1 0 8a8 8 0 0 1 16 0zM5.496 6.033h.825c.138 0 .248-.113.266-.25.09-.656.54-1.134 1.342-1.134.686 0 1.314.343 1.314 1.168 0 .635-.374.927-.965 1.371-.673.489-1.206 1.06-1.168 1.987l.003.217a.25.25 0 0 0 .25.246h.811a.25.25 0 0 0 .25-.25v-.105c0-.718.273-.927 1.01-1.486.609-.463 1.244-.977 1.244-2.056 0-1.511-1.276-2.241-2.673-2.241-1.267 0-2.655.59-2.75 2.286a.237.237 0 0 0 .241.247zm2.325 6.443c.61 0 1.029-.394 1.029-.927 0-.552-.42-.94-1.029-.94-.584 0-1.009.388-1.009.94 0 .533.425.927 1.01.927z"/></svg>'
)

def read_markdown_file(path):
    with open(path, 'r') as file:
        return file.read()

app_ui = ui.page_navbar(

# Main page
    ui.nav_panel(
        ##  page title and subtitle
        ui.output_image('logo', inline=True),
        ## markdown file with text
        ui.markdown(read_markdown_file(f'{path_user}www/copepodTCR.md')),

    ),

# Instructions
    ui.nav_panel(
        ## page title and subtitle
        ui.output_image('Instructions', inline=True),
        ui.markdown(read_markdown_file(f'{path_user}www/instructions1.md')),
        ui.output_table('simulation_example'),
        ui.markdown(read_markdown_file(f'{path_user}www/instructions2.md')),
        ui.output_image('peptide_plate_scheme_instr', inline=True),
        ui.markdown(read_markdown_file(f'{path_user}www/instructions3.md')),
        ui.output_image('activation_model', inline=True),
        ui.markdown(read_markdown_file(f'{path_user}www/instructions4.md'),),
        ui.output_table('cell_f1'),
        ui.markdown(read_markdown_file(f'{path_user}www/instructions5.md'),),
    ),

# Balance check

    ui.nav_panel(
        ## page title and subtitle
        ui.output_image('Balance_check', inline=True),
        ui.panel_title('Pooling scheme settings \n', 'Balance check'),

        ## enter n_pools (num input)
        ui.tooltip(ui.span('Number of pools ', question_circle_fill),
                     'Number of pools in the experiment', 
                     placement = 'right',
                     id = 'card_tooltip1_bc'),
        ui.input_numeric('npools_bc', '', value = 10),

        ## enter n_peptides (num input)
        ui.tooltip(ui.span('Number of peptides ', question_circle_fill),
                     'Number of tested peptides',
                     placement = 'right',
                     id = 'card_tooltip2_bc'),
        ui.input_numeric('lenlst_bc', '', value = 100),

        ## peptide occurrence (slider)
        ui.tooltip(ui.span('Possible peptide occurrence ', question_circle_fill),
                     'Number of pools to which one peptide is added',
                     placement = 'right',
                     id = 'card_tooltip3_bc'),
        ui.output_ui('iters_slider_bc'),

        ## resulting negative share (text)
        ui.tooltip(ui.span('Negative share ', question_circle_fill),
                     'Expected number of non-activated pools',
                     placement = 'right',
                     id = 'card_tooltip4_bc'),
        ui.output_text_verbatim('negshare_bc'),

        ## resulting balance (text)
        ui.tooltip(ui.span('Resulting balance ', question_circle_fill),
                      'Distribution of peptides per pool',
                      placement = 'right',
                      id = 'card_tooltip5_bc'),
        ui.output_text_verbatim('balance_bc'),

        ## resulting balance (plot)
        ui.output_plot('balance_plot_bc', width="700px", height="400px")
    ),

# Peptide check
    ui.nav_panel(
        ## page title and subtitle
        ui.output_image('Peptide_check', inline=True),
        ui.panel_title('Generate peptides \n', 'Peptide check'),

        ## upload proteins (file)
        ui.tooltip(ui.span('Choose a fasta file with your protein or proteins ', question_circle_fill),
                     'FASTA file with amino acid sequences',
                     placement = 'right',
                     id = 'card_tooltip1_pc'
                     ),
        ui.input_file('proteins_file_pc', '', multiple = False, accept = '.fasta'),

        ## upload proteins (enter)
        ui.tooltip(ui.span('Or enter a protein sequence here ', question_circle_fill),
                     'Amino acid sequences of proteins separated by comma',
                     placement = 'right',
                     id = 'card_tooltip2_pc'
                     ),
        ui.input_text_area('proteins_text_pc', '', placeholder = 'protein1,', resize = 'both'),

        ## generation parameters
        ### peptide length
        ui.tooltip(ui.span('Peptide length ', question_circle_fill),
                     'Length of the generated peptides',
                     placement = 'right',
                     id = 'card_tooltip3_pc'
                     ),
        ui.input_numeric('pep_length_pc', '', value = 17, min = 10, max = 50),
        ### shift
        ui.tooltip(ui.span('Shift between peptides ', question_circle_fill),
                     'Number of amino acids to shift when generating peptides by sliding window',
                     placement = 'right',
                     id = 'card_tooltip4_pc'
                     ),
        ui.input_numeric('shift_pc', '', value = 5, min = 1, max = 50),
        ### protein end
        ui.tooltip(ui.span('Padd protein end? ', question_circle_fill),
                     'If True, last peptide will be padded by upstream sequence if shorter than expected',
                     placement = 'right',
                     id = 'card_tooltip5_pc'
                     ),
        ui.input_switch('protein_end_pc', '', value = True),

        ## button to generate peptides
        ui.input_action_button('generate_peptides_pc', 'Generate peptides',
                               class_="btn-primary",
                               style="color: #fff; background-color: #015396; border-color: #015396"),

        ## generated peptides text box
        ui.br(),
        ui.br(),
        ui.output_text_verbatim('generated_peptides_pc'),

        ## button to download resulting peptides
        ui.download_button('download_peptides_pc', 'Download generated peptides',
                           class_="btn-success",
                           style="color: #015396; background-color: #fff; border-color: #015396"),

        ## peptide check title
        ui.br(),
        ui.br(),
        ui.br(),
        ui.panel_title('Check peptides \n', 'Peptide check'),

        ## upload peptides
        ui.tooltip(ui.span('Choose a file with peptides to upload ', question_circle_fill),
                     'TSV or CSV table with peptides in a column without a header',
                     placement = 'right',
                     id = 'card_tooltip6_pc'
                     ),
        ui.input_file('file_peptides_pc', '', multiple = False, accept = '.tsv,.csv'),

        ## enter peptides
        ui.tooltip(ui.span('Or enter peptides here ', question_circle_fill),
                     'List of peptides separated by comma',
                     placement = 'right',
                     id = 'card_tooltip7_pc'
                     ),
        ui.input_text_area('text_peptides_pc', '', placeholder = 'peptide1, peptide2', resize = 'both'),

        ### button to get sample peptides
        ui.input_action_button('sample_peptides_pc', 'Load sample peptides',
                               class_="btn",
                               style="color: #015396; background-color: #fff; border-color: #015396"
                               ),

        ## button to check for overlap length
        ui.input_action_button('get_overlap_pc', 'Check',
                               class_="btn-primary",
                               style="color: #fff; background-color: #015396; border-color: #015396"),
        
        ## output (text)
        ui.br(),
        ui.br(),
        ui.output_text_verbatim('peptides_info_pc'),
        ui.output_text_verbatim('strange_overlaps_pc'),
    ),

# Pooling
    ui.nav_panel(
        ## page title and subtitle
        ui.output_image('Pooling_scheme', inline=True),
        ui.panel_title('Pooling scheme settings \n', 'Pooling'),

        ## number of pools (num input)
        ui.tooltip(ui.span("Number of pools ", question_circle_fill),
                     "Number of pools in the experiment",
                     placement="right",
                     id="card_tooltip1_po"),
        ui.input_numeric("npools_po", '', value = 20),

        ## upload peptides
        ui.tooltip(ui.span('Choose a file with peptides to upload ', question_circle_fill),
                     'TSV or CSV table with peptides in a column without a header',
                     placement = 'right',
                     id = 'card_tooltip2_po'
                     ),
        ui.input_file('file_peptides_po', '', multiple = False, accept = '.tsv,.csv'),

        ## enter peptides
        ui.tooltip(ui.span('Or enter peptides here ', question_circle_fill),
                     'List of peptides separated by comma',
                     placement = 'right',
                     id = 'card_tooltip3_po'
                     ),
        ui.input_text_area('text_peptides_po', '', placeholder = 'peptide1, peptide2', resize = 'both'),

        ### button to get sample peptides
        ui.input_action_button('sample_peptides_po', 'Load sample peptides',
                               class_="btn",
                               style="color: #015396; background-color: #fff; border-color: #015396"
                               ),

        ## peptide occurrence (slider)
        ui.br(),
        ui.br(),
        ui.tooltip(ui.span('Possible peptide occurrence ', question_circle_fill),
                     'Number of pools to which one peptide is added',
                     placement = 'right',
                     id = 'card_tooltip4_po'),
        ui.output_ui('iters_slider_po'),

        ## epitope length (slider)
        ui.tooltip(ui.span("Expected epitope length ", question_circle_fill),
                     "The length in amino acids of the epitope for simulation step",
                     placement="right",
                     id="card_tooltip5_po"),
        ui.input_slider('ep_length_po', '', min = 8, max = 15, value=8),

        ## resulting negative share (text)
        ui.tooltip(ui.span('Negative share ', question_circle_fill),
                     'Expected number of non-activated pools, the value will be displayed after you enter peptides',
                     placement = 'right',
                     id = 'card_tooltip6_po'),
        ui.output_text_verbatim('negshare_po', placeholder=True),

        ## email input
        ui.input_text("letter_po", "Enter your email", placeholder="example@mail.ru"),

        ## button to generate pooling scheme
        ui.input_action_button('generate_scheme_po', 'Generate a pooling scheme',
                               class_="btn-primary",
                               style="color: #fff; background-color: #015396; border-color: #015396"),

        ## button to send generated scheme
        ui.span("   "),
        ui.input_action_button("send_po", "Send the scheme",
                               class_="btn-success",
                               style="color: #015396; background-color: #fff; border-color: #015396"),
        
        ## generated scheme was sent (notification)
        ui.output_text('notification_po'),

        ## button to download generated scheme
        ui.br(),
        ui.download_button('download_po', 'Download the scheme',
                           class_="btn-success",
                           style="color: #015396; background-color: #fff; border-color: #015396"),

        ## notification about resulting scheme (text)
        ui.output_ui('scheme_result_po'),

        ## resulting balance (plot)
        ui.output_plot('balance_plot_po', width="700px", height="400px"),

        ## resulting scheme (dataframe)
        ui.output_data_frame('check_results_po'),
    ),

# STL files
    ui.nav_panel(
        ## page title and subtitle
        ui.output_image('STL_files', inline=True),
        ui.panel_title('Download zip with STL files \n', 'STL files'),

        ## upload a file with pools (pools_df)
        ui.tooltip(ui.span("Choose a file with pools ", question_circle_fill),
                "Produced by copepodTCR in pooling step, pools",
                placement="right",
                id="card_tooltip1_stl"),
        ui.input_file('pools_df_stl', '', multiple=False, accept = '.tsv'),

        ## upload a files with peptide arrangement in the plate
        ui.tooltip(ui.span("Choose a file with peptide arrangement in the plate ", question_circle_fill),
                "TSV table with one peptide per cell, rows and columns are corresponding rows and columns of the plate",
                placement="right",
                id="card_tooltip2_stl"),
        ui.input_file('peptide_plate_stl', '', multiple=False, accept = '.TSV'),

        ## scheme with parameters
        ui.output_image('peptide_plate_scheme_stl', inline=False),

        ## parameters
        ### rows
        ui.tooltip(ui.span("Number of rows ", question_circle_fill),
                "Number of rows in the plate",
                placement="right",
                id="card_tooltip3_stl",),
        ui.input_numeric("rows_stl", '', value = 16),
        ### columns
        ui.tooltip(ui.span("Number of columns ", question_circle_fill),
                "Number of columns in the plate",
                placement="right",
                id="card_tooltip4_stl",),
        ui.input_numeric("cols_stl", '', value = 24),
        ### length
        ui.tooltip(ui.span("Length (mm) ", question_circle_fill),
                "Length of the plate in mm",
                placement="right",
                id="card_tooltip5_stl",),
        ui.input_numeric("length_stl", '', value = 122.10),
        ### width
        ui.tooltip(ui.span("Width (mm) ", question_circle_fill),
                "Width of the plate in mm",
                placement="right",
                id="card_tooltip6_stl",),
        ui.input_numeric("width_stl", '', value = 79.97),
        ### thickness
        ui.tooltip(ui.span("Thickness (mm) ", question_circle_fill),
                "Thickness of the plate in mm",
                placement="right",
                id="card_tooltip7_stl",),
        ui.input_numeric("thickness_stl", '', value = 1.5),
        ### hole radius
        ui.tooltip(ui.span("Hole radius (mm) ", question_circle_fill),
                "Diameter of the hole divided by two",
                placement="right",
                id="card_tooltip8_stl",),
        ui.input_numeric("hole_radius_stl", '', value = 2),
        ### X offset
        ui.tooltip(ui.span("X offset (mm) ", question_circle_fill),
                "Margin along the X axis for the A1 hole, in mm",
                placement="right",
                id="card_tooltip9_stl",),
        ui.input_numeric("x_offset_stl", '', value = 9.05),
        ### Y offset
        ui.tooltip(ui.span("Y offset (mm) ", question_circle_fill),
                "Margin along the Y axis for the A1 hole, in mm",
                placement="right",
                id="card_tooltip10_stl",),
        ui.input_numeric("y_offset_stl", '', value = 6.20),
        ### well spacing
        ui.tooltip(ui.span("Well spacing (mm) ", question_circle_fill),
                "Distance between wells",
                placement="right",
                id="card_tooltip11_stl",),
        ui.input_numeric("well_spacing_stl", '', value = 4.5),
        ### well with coordinates 16-24
        ui.tooltip(ui.span("Add well with coordinates m-k? ", question_circle_fill),
                "To better orient tip pattern, you can add the hole with these coordinates, but it should be done only if there is no peptide in the plate in respective well. ",
                placement="right",
                id="card_tooltip12_stl",),
        ui.input_switch('well24_stl', '', value=True),
        ### marks
        ui.tooltip(ui.span("Add pool indices to the masks? ", question_circle_fill),
                "To distinguish between masks for different pools, you can mark them.",
                placement="right",
                id="card_tooltip13_stl",),
        ui.input_switch('marks_stl', '', value=True),

        ## button to generate stl files
        ui.input_action_button('generate_stl', 'Generate STL files',
                               class_="btn-primary",
                               style="color: #fff; background-color: #015396; border-color: #015396"),
        
        ## notification about generated stl files
        ui.br(),
        ui.output_text('stl_generation'),
        ui.br(),

        ## enter email to send stl files
        ui.input_text("letter_stl", "Enter your email", placeholder="example@mail.ru"),

        ## button to send stl files
        ui.input_action_button("send_stl", "Send results",
                               class_="btn-success",
                               style="color: #015396; background-color: #fff; border-color: #015396"),

        ## notification that files were sent
        ui.output_text('notification_stl'),

        ## button to download generated STL files
        ui.br(),
        ui.download_button('download_stl', 'Download zip file with STL files',
                           style="color: #015396; background-color: #fff; border-color: #015396")
    ),

# Results interpretation
    ui.nav_panel(
        ## page title and subtitle
        ui.output_image('Interpretation', inline=True),
        ui.panel_title('Analyze results \n', 'Interpretation'),

        ## upload check_results file
        ui.tooltip(ui.span("Choose a file with simulation ", question_circle_fill),
                "Produced during pooling step, simulation.tsv",
                placement="right",
                id="card_tooltip1_in",),
        ui.input_file('check_results_in', '', multiple=False, accept=".tsv"),

        ## upload experimental read out
        ui.tooltip(ui.span("Experiment read out ", question_circle_fill),
                "Numeric results of the experiment, in percentage of activated T cells, a table with columns 'Pool', 'Percentage'",
                placement="right",
                id="card_tooltip2_in",),
        ui.input_file('exp_in', '', multiple=False, accept=".csv"),

        ### button to show example of the experimental read-out
        ui.input_action_button('exp_example_in', 'Show example read out',
                               class_="btn",
                               style="color: #015396; background-color: #fff; border-color: #015396"
                               ),
        ui.br(),
        ui.br(),

        ## example dataframe
        ui.output_data_frame('exp_example_out'),

        ## pooling settings
        ### number of pools
        ui.tooltip(ui.span("Number of pools ", question_circle_fill),
                "Number of pools in the experiment",
                placement="right",
                id="card_tooltip3_in",),
        ui.input_numeric('npools_in', '', value=0),
        ### peptide occurrence
        ui.tooltip(ui.span("Peptide occurrence ", question_circle_fill),
                "Number of pools to which a single peptide was added",
                placement="right",
                id="card_tooltip4_in",),
        ui.input_numeric('iters_in', '', value=0),
        ### epitope length
        ui.tooltip(ui.span("Epitope length ", question_circle_fill),
                "Length of the peptide epitope",
                placement="right",
                id="card_tooltip5_in",),
        ui.input_slider('ep_length_in', '', value=8, min = 8, max = 15),

        ## negative controls
        ui.tooltip(ui.span("Negative controls ", question_circle_fill),
                "Negative control values, in percentage of activated T cells, different replicates separated by comma.",
                placement="right",
                id="card_tooltip6_in",),
        ui.input_text('negcontrol_in', '', placeholder = 'replicate1,replicate2', value=''),

        ## button to start the model
        ui.input_action_button('start_model_in', 'Analyze',
                               class_="btn-primary",
                               style="color: #fff; background-color: #015396; border-color: #015396"),
        
        ## enter email
        ui.br(),
        ui.br(),
        ui.input_text("letter_in", "Enter your email ", placeholder="example@mail.ru"),

        ## button to send peptide probabilities
        ui.input_action_button("send_in", "Send results",
                              class_="btn-success",
                              style="color: #015396; background-color: #fff; border-color: #015396"),

        ## notification that they were sent
        ui.output_text('notification_in'),

        ## button to download peptide probabilities
        ui.br(),
        ui.download_button('download_in', 'Download results', 
                           style="color: #015396; background-color: #fff; border-color: #015396"),

        ## model fit
        ui.output_plot('fitting_in', width="700px", height="400px"),

        ## activated pools (text)
        ui.output_text_verbatim('activated_pools_in'),

        ## activated pools (poolplot)
        ui.output_plot('poolplot_in', width="700px", height="400px"),

        ## cognate peptides (text)
        ui.output_text_verbatim('cognate_peptides_in'),

        ## cognate peptides (bubbleplot)
        output_widget("bubbleplot_in"),
    )

)

def server(input: Inputs, output: Outputs, session: Session):

# Header
    @output
    @render.image
    def logo():
        img: ImgData = {"src": f"{path_user}www/intro.svg", "width": "130px"}
        return img
    
    @output
    @render.image
    def Instructions():
        img: ImgData = {"src": f"{path_user}www/Instructions.svg", "width": "130px"}
        return img

    @output
    @render.image
    def Balance_check():
        img: ImgData = {"src": f'{path_user}www/Balance_check.svg', "width": "130px"}
        return img
    @output
    @render.image
    def Peptide_check():
        img: ImgData = {"src": f"{path_user}www/Peptide_check.svg", "width": "130px"}
        return img
    @output
    @render.image
    def Pooling_scheme():
        img: ImgData = {"src": f"{path_user}www/Pooling_scheme.svg", "width": "130px"}
        return img
    @output
    @render.image
    def STL_files():
        img: ImgData = {"src": f"{path_user}www/STL_files.svg", "width": "130px"}
        return img
    @output
    @render.image
    def Interpretation():
        img: ImgData = {"src": f"{path_user}www/Interpretation.svg", "width": "130px"}
        return img
    
# Instructions
    @output
    @render.table
    def simulation_example():
        table = pd.read_csv(f'{path_user}www/simulation_example.tsv', sep = "\t")
        return table
    @output
    @render.image
    def peptide_plate_scheme_instr():
        img: ImgData = {"src": f"{path_user}www/peptide_plate_scheme.svg"}
        return img
    @output
    @render.image
    def activation_model():
        img: ImgData = {"src": f"{path_user}www/activation_model.svg"}
        return img
    @output
    @render.table
    def cell_f1():
        table = pd.read_csv(f'{path_user}www/cell_f1.csv')
        return table

# Balance check
    
    ## peptide occurrence (slider)
    @output
    @render.ui
    def iters_slider_bc():
        m = input.npools_bc()
        n = input.lenlst_bc()

        if type(m) != int or type(n) != int:
            return 'Please enter an integer'
        
        range_iters = []
        for r in range(m):
            if comb(m, r)*0.8 >= n and comb(m, r+1)*0.8 >= n:
                range_iters.append(r)

        if len(range_iters) != 0:
            result_iters_bc = ui.input_slider('iters_bc',
                                              '',
                                              min = range_iters[0],
                                              max = range_iters[-1],
                                              value = range_iters[0],
                                              step = 1)
            return result_iters_bc
        else:
            return 'Warning: InvalidParameters. Increase number of pools or decrease number of peptides.'
        
    ## negative share (text)
    @reactive.Calc
    def negshare_calculation():
        m = input.npools_bc()
        r = input.iters_bc()
        n = input.lenlst_bc()
        if type(m) == int and r < m and type(r) == int:
            if comb(m, r)*0.8 >= n and comb(m, r+1)*0.8 >= n:
                negshare_result_bc = (m - r - 1)/m
                return negshare_result_bc
    
    @output
    @render.text
    def negshare_bc():
        if negshare_calculation() is not None:
            return negshare_calculation()

    ## balance (text)
    @reactive.Calc
    def balance_calculation():
        m = input.npools_bc()
        r = input.iters_bc()
        n = input.lenlst_bc()
        if type(m) == int and type(n) == int:
            if comb(m, r)*0.8 >= n and comb(m, r+1)*0.8 >= n:
                if n < 200:
                    b, lines = cdp.bba(m = m, r = r, n = n)
                elif n >= 200:
                    b, lines = cdp.rcbba(m = m, r = r, n = n)
                return b
            else:
                return 'Warning: InvalidParameters. Increase number of pools.'
            
    @output
    @render.text
    def balance_bc():
        if balance_calculation() is not None:
            return balance_calculation()
        
    ## balance (plot)
    @output
    @render.plot
    def balance_plot_bc():
        m = input.npools_bc()
        r = input.iters_bc()
        n = input.lenlst_bc()

        b = balance_calculation()
        if b is not None and type(b) != str:
            perfect = [ceil(n*r/m)]*m

            fig, axes = plt.subplots(figsize = (6, 4))
            sn.scatterplot(y = b, x = list(range(1, m+1)), label = 'balance',
                           facecolor = 'none', s = 50, edgecolor = '#015396', linewidth = 2)
            sn.lineplot(y = perfect, x = list(range(1, m+1)),
                        color = '#E2E2ED', linestyle = '--', label = 'perfect balance', linewidth = 2)

            plt.ylim(0, max(b) + 10)
            plt.xlim(0, m + 1)
            plt.xticks(range(1, m + 1))
            plt.yticks([int(tick) for tick in plt.yticks()[0]])

            plt.ylabel('Number of peptides per pool')
            plt.xlabel('Pool index')
            sn.despine()

            return fig
        
# Peptide check
    
    ## generate peptides
    @reactive.Calc
    @reactive.event(input.generate_peptides_pc)
    def peptide_generation_react():

        file_infos = input.proteins_file_pc()

        if not file_infos:
            protein_seqs = re.split(r'\W+', input.proteins_text_pc())
        else:
            for file_info in file_infos:
                file_proteins_path = file_info['datapath']
            protein_ids, protein_seqs = parse_fasta(file_proteins_path)
        
        peptides = cpp.peptide_generation(protein_seqs,
                                          input.pep_length_pc(),
                                          input.shift_pc(), input.protein_end_pc())
        return peptides
    
    ## return generated peptides (text)
    @output
    @render.text
    @reactive.event(input.generate_peptides_pc)
    def generated_peptides_pc():
        peptides = peptide_generation_react()
        text = ''
        for p in peptides:
            text = text + p + ','
        return text
    
    ## return generated peptides (download)
    @output
    @session.download(filename = 'peptides.zip')
    def download_peptides_pc():
        peptides = peptide_generation_react()
        if len(peptides) > 0:
            peptides_df = pd.DataFrame(peptides)

            zip_buffer = dataframes_to_zip([peptides_df],
                                           [f'peptides_{input.pep_length_pc()}aa_{input.shift_pc()}shift.tsv'])
            yield zip_buffer.getvalue()

# Overlap check

    ## load sample peptides
    @reactive.Effect
    @reactive.event(input.sample_peptides_pc)
    def _():
        ui.update_text_area('text_peptides_pc', value = 'MFVFLVLLPLVSSQCVN, VLLPLVSSQCVNLTTRT, VSSQCVNLTTRTQLPPA, VNLTTRTQLPPAYTNSF, RTQLPPAYTNSFTRGVY, PAYTNSFTRGVYYPDKV, SFTRGVYYPDKVFRSSV, VYYPDKVFRSSVLHSTQ, KVFRSSVLHSTQDLFLP, SVLHSTQDLFLPFFSNV, TQDLFLPFFSNVTWFHA, LPFFSNVTWFHAIHVSG, NVTWFHAIHVSGTNGTK, HAIHVSGTNGTKRFDNP, SGTNGTKRFDNPVLPFN, TKRFDNPVLPFNDGVYF, NPVLPFNDGVYFASTEK, FNDGVYFASTEKSNIIR, YFASTEKSNIIRGWIFG, EKSNIIRGWIFGTTLDS, IRGWIFGTTLDSKTQSL, FGTTLDSKTQSLLIVNN, DSKTQSLLIVNNATNVV, SLLIVNNATNVVIKVCE, NNATNVVIKVCEFQFCN, VVIKVCEFQFCNDPFLG, CEFQFCNDPFLGVYYHK, CNDPFLGVYYHKNNKSW, LGVYYHKNNKSWMESEF, HKNNKSWMESEFRVYSS, SWMESEFRVYSSANNCT, EFRVYSSANNCTFEYVS, SSANNCTFEYVSQPFLM, CTFEYVSQPFLMDLEGK, VSQPFLMDLEGKQGNFK, LMDLEGKQGNFKNLREF, GKQGNFKNLREFVFKNI, FKNLREFVFKNIDGYFK, EFVFKNIDGYFKIYSKH, NIDGYFKIYSKHTPINL, FKIYSKHTPINLVRDLP, KHTPINLVRDLPQGFSA, NLVRDLPQGFSALEPLV, LPQGFSALEPLVDLPIG, SALEPLVDLPIGINITR, LVDLPIGINITRFQTLL, IGINITRFQTLLALHRS, TRFQTLLALHRSYLTPG, LLALHRSYLTPGDSSSG, RSYLTPGDSSSGWTAGA, PGDSSSGWTAGAAAYYV, SGWTAGAAAYYVGYLQP, GAAAYYVGYLQPRTFLL, YVGYLQPRTFLLKYNEN, QPRTFLLKYNENGTITD, LLKYNENGTITDAVDCA, ENGTITDAVDCALDPLS, TDAVDCALDPLSETKCT, CALDPLSETKCTLKSFT, LSETKCTLKSFTVEKGI, CTLKSFTVEKGIYQTSN, FTVEKGIYQTSNFRVQP, GIYQTSNFRVQPTESIV, SNFRVQPTESIVRFPNI, QPTESIVRFPNITNLCP, IVRFPNITNLCPFGEVF, NITNLCPFGEVFNATRF, CPFGEVFNATRFASVYA, VFNATRFASVYAWNRKR, RFASVYAWNRKRISNCV, YAWNRKRISNCVADYSV, KRISNCVADYSVLYNSA, CVADYSVLYNSASFSTF, SVLYNSASFSTFKCYGV, SASFSTFKCYGVSPTKL, TFKCYGVSPTKLNDLCF, GVSPTKLNDLCFTNVYA, KLNDLCFTNVYADSFVI, CFTNVYADSFVIRGDEV, YADSFVIRGDEVRQIAP, VIRGDEVRQIAPGQTGK, EVRQIAPGQTGKIADYN, APGQTGKIADYNYKLPD, GKIADYNYKLPDDFTGC, YNYKLPDDFTGCVIAWN, PDDFTGCVIAWNSNNLD, GCVIAWNSNNLDSKVGG, WNSNNLDSKVGGNYNYL, LDSKVGGNYNYLYRLFR, GGNYNYLYRLFRKSNLK, YLYRLFRKSNLKPFERD, FRKSNLKPFERDISTEI, LKPFERDISTEIYQAGS, RDISTEIYQAGSTPCNG, EIYQAGSTPCNGVEGFN, GSTPCNGVEGFNCYFPL, NGVEGFNCYFPLQSYGF, FNCYFPLQSYGFQPTNG, PLQSYGFQPTNGVGYQP, GFQPTNGVGYQPYRVVV, NGVGYQPYRVVVLSFEL, QPYRVVVLSFELLHAPA, VVLSFELLHAPATVCGP, ELLHAPATVCGPKKSTN, PATVCGPKKSTNLVKNK, GPKKSTNLVKNKCVNFN, TNLVKNKCVNFNFNGLT, NKCVNFNFNGLTGTGVL, FNFNGLTGTGVLTESNK, LTGTGVLTESNKKFLPF, VLTESNKKFLPFQQFGR, NKKFLPFQQFGRDIADT, PFQQFGRDIADTTDAVR, GRDIADTTDAVRDPQTL, DTTDAVRDPQTLEILDI, VRDPQTLEILDITPCSF, TLEILDITPCSFGGVSV, DITPCSFGGVSVITPGT, SFGGVSVITPGTNTSNQ, SVITPGTNTSNQVAVLY, GTNTSNQVAVLYQDVNC, NQVAVLYQDVNCTEVPV, LYQDVNCTEVPVAIHAD, NCTEVPVAIHADQLTPT, PVAIHADQLTPTWRVYS, ADQLTPTWRVYSTGSNV, PTWRVYSTGSNVFQTRA, YSTGSNVFQTRAGCLIG, NVFQTRAGCLIGAEHVN, RAGCLIGAEHVNNSYEC, IGAEHVNNSYECDIPIG, VNNSYECDIPIGAGICA, ECDIPIGAGICASYQTQ, IGAGICASYQTQTNSPR, CASYQTQTNSPRRARSV, TQTNSPRRARSVASQSI, PRRARSVASQSIIAYTM, SVASQSIIAYTMSLGAE, SIIAYTMSLGAENSVAY, TMSLGAENSVAYSNNSI, AENSVAYSNNSIAIPTN, AYSNNSIAIPTNFTISV, SIAIPTNFTISVTTEIL, TNFTISVTTEILPVSMT, SVTTEILPVSMTKTSVD, ILPVSMTKTSVDCTMYI, MTKTSVDCTMYICGDST, VDCTMYICGDSTECSNL, YICGDSTECSNLLLQYG, STECSNLLLQYGSFCTQ, NLLLQYGSFCTQLNRAL, YGSFCTQLNRALTGIAV, TQLNRALTGIAVEQDKN, ALTGIAVEQDKNTQEVF, AVEQDKNTQEVFAQVKQ, KNTQEVFAQVKQIYKTP, VFAQVKQIYKTPPIKDF, KQIYKTPPIKDFGGFNF, TPPIKDFGGFNFSQILP, DFGGFNFSQILPDPSKP, NFSQILPDPSKPSKRSF, LPDPSKPSKRSFIEDLL, KPSKRSFIEDLLFNKVT, SFIEDLLFNKVTLADAG, LLFNKVTLADAGFIKQY, VTLADAGFIKQYGDCLG, AGFIKQYGDCLGDIAAR, QYGDCLGDIAARDLICA, LGDIAARDLICAQKFNG, ARDLICAQKFNGLTVLP, CAQKFNGLTVLPPLLTD, NGLTVLPPLLTDEMIAQ, LPPLLTDEMIAQYTSAL, TDEMIAQYTSALLAGTI, AQYTSALLAGTITSGWT, ALLAGTITSGWTFGAGA, TITSGWTFGAGAALQIP, WTFGAGAALQIPFAMQM, GAALQIPFAMQMAYRFN, IPFAMQMAYRFNGIGVT, QMAYRFNGIGVTQNVLY, FNGIGVTQNVLYENQKL, VTQNVLYENQKLIANQF, LYENQKLIANQFNSAIG, KLIANQFNSAIGKIQDS, QFNSAIGKIQDSLSSTA, IGKIQDSLSSTASALGK, DSLSSTASALGKLQDVV, TASALGKLQDVVNQNAQ, GKLQDVVNQNAQALNTL, VVNQNAQALNTLVKQLS, AQALNTLVKQLSSNFGA, TLVKQLSSNFGAISSVL, LSSNFGAISSVLNDILS, GAISSVLNDILSRLDKV, VLNDILSRLDKVEAEVQ, LSRLDKVEAEVQIDRLI, KVEAEVQIDRLITGRLQ, VQIDRLITGRLQSLQTY, LITGRLQSLQTYVTQQL, LQSLQTYVTQQLIRAAE, TYVTQQLIRAAEIRASA, QLIRAAEIRASANLAAT, AEIRASANLAATKMSEC, SANLAATKMSECVLGQS, ATKMSECVLGQSKRVDF, ECVLGQSKRVDFCGKGY, QSKRVDFCGKGYHLMSF, DFCGKGYHLMSFPQSAP, GYHLMSFPQSAPHGVVF, SFPQSAPHGVVFLHVTY, APHGVVFLHVTYVPAQE, VFLHVTYVPAQEKNFTT, TYVPAQEKNFTTAPAIC, QEKNFTTAPAICHDGKA, TTAPAICHDGKAHFPRE, ICHDGKAHFPREGVFVS, KAHFPREGVFVSNGTHW, REGVFVSNGTHWFVTQR, VSNGTHWFVTQRNFYEP, HWFVTQRNFYEPQIITT, QRNFYEPQIITTDNTFV, EPQIITTDNTFVSGNCD, TTDNTFVSGNCDVVIGI, FVSGNCDVVIGIVNNTV, CDVVIGIVNNTVYDPLQ, GIVNNTVYDPLQPELDS, TVYDPLQPELDSFKEEL, LQPELDSFKEELDKYFK, DSFKEELDKYFKNHTSP, ELDKYFKNHTSPDVDLG, FKNHTSPDVDLGDISGI, SPDVDLGDISGINASVV, LGDISGINASVVNIQKE, GINASVVNIQKEIDRLN, VVNIQKEIDRLNEVAKN, KEIDRLNEVAKNLNESL, LNEVAKNLNESLIDLQE, KNLNESLIDLQELGKYE, SLIDLQELGKYEQYIKW, QELGKYEQYIKWPWYIW, YEQYIKWPWYIWLGFIA, KWPWYIWLGFIAGLIAI, IWLGFIAGLIAIVMVTI, IAGLIAIVMVTIMLCCM, AIVMVTIMLCCMTSCCS, TIMLCCMTSCCSCLKGC, CMTSCCSCLKGCCSCGS, CSCLKGCCSCGSCCKFD, GCCSCGSCCKFDEDDSE, GSCCKFDEDDSEPVLKG, FDEDDSEPVLKGVKLHY')

    ## information about loaded peptides (text)
    @output
    @render.text
    @reactive.event(input.get_overlap_pc, ignore_none = True)
    def peptides_info_pc():
        file_infos = input.file_peptides_pc()
        if not file_infos:
            peptides = re.split(r'\W+', input.text_peptides_pc())
        else:
            for file_info in file_infos:
                input_peptides = pd.read_csv(file_info['datapath'], header = None)
                peptides = list(input_peptides[0])

        p_lengths = [len(p) for p in peptides]
        p_length_counts = Counter(p_lengths)

        overlaps = cpp.all_overlaps(peptides)
        overlaps_text = ''

        for key in overlaps:
            t = f'{str(overlaps[key]+1)} peptides with overlap {str(key)} amino acids\n'
            overlaps_text = overlaps_text +  t
        
        overlaps_text = overlaps_text + '\n'

        for key in p_length_counts:
            tpl = f'{str(p_length_counts[key])} peptides with length {str(key)} amino acids\n'
            overlaps_text = overlaps_text + tpl

        if len(overlaps) > 1:
            m = ui.modal(
            "Variable overlap lengths in the peptide set hinder result interpretation.",
            title="Warning",
            easy_close=True,
            footer=None)
            ui.modal_show(m)

        return overlaps_text
    
    ## overlap information (text)
    @output
    @render.text
    @reactive.event(input.get_overlap_pc, ignore_none = False)
    def strange_overlaps_pc():
        file_infos = input.file_peptides_pc()
        if not file_infos:
            peptides = re.split(r'\W+', input.text_peptides_pc())
        else:
            for file_info in file_infos:
                input_peptides = pd.read_csv(file_info['datapath'], header = None)
                peptides = list(input_peptides[0])
        
        overlaps = cpp.all_overlaps(peptides)
        strange_overlaps = [i for i in overlaps if overlaps[i]!=max(overlaps.values())]
        pair_text = ''
        for overlap in strange_overlaps:
            pair = cpp.find_pair_with_overlap(peptides, overlap)
            pair_text = pair_text + 'Overlap of ' + str(overlap) + ' amino acids:' + '\n'
            for p in pair:
                for peptide in p:
                    pair_text = pair_text + peptide + ' '
                pair_text = pair_text + '\n'
        return pair_text
    
# Pooling scheme
    
    ## load sample peptides
    @reactive.Effect
    @reactive.event(input.sample_peptides_po)
    def _():
        ui.update_text_area('text_peptides_po', value = 'MFVFLVLLPLVSSQCVN, VLLPLVSSQCVNLTTRT, VSSQCVNLTTRTQLPPA, VNLTTRTQLPPAYTNSF, RTQLPPAYTNSFTRGVY, PAYTNSFTRGVYYPDKV, SFTRGVYYPDKVFRSSV, VYYPDKVFRSSVLHSTQ, KVFRSSVLHSTQDLFLP, SVLHSTQDLFLPFFSNV, TQDLFLPFFSNVTWFHA, LPFFSNVTWFHAIHVSG, NVTWFHAIHVSGTNGTK, HAIHVSGTNGTKRFDNP, SGTNGTKRFDNPVLPFN, TKRFDNPVLPFNDGVYF, NPVLPFNDGVYFASTEK, FNDGVYFASTEKSNIIR, YFASTEKSNIIRGWIFG, EKSNIIRGWIFGTTLDS, IRGWIFGTTLDSKTQSL, FGTTLDSKTQSLLIVNN, DSKTQSLLIVNNATNVV, SLLIVNNATNVVIKVCE, NNATNVVIKVCEFQFCN, VVIKVCEFQFCNDPFLG, CEFQFCNDPFLGVYYHK, CNDPFLGVYYHKNNKSW, LGVYYHKNNKSWMESEF, HKNNKSWMESEFRVYSS, SWMESEFRVYSSANNCT, EFRVYSSANNCTFEYVS, SSANNCTFEYVSQPFLM, CTFEYVSQPFLMDLEGK, VSQPFLMDLEGKQGNFK, LMDLEGKQGNFKNLREF, GKQGNFKNLREFVFKNI, FKNLREFVFKNIDGYFK, EFVFKNIDGYFKIYSKH, NIDGYFKIYSKHTPINL, FKIYSKHTPINLVRDLP, KHTPINLVRDLPQGFSA, NLVRDLPQGFSALEPLV, LPQGFSALEPLVDLPIG, SALEPLVDLPIGINITR, LVDLPIGINITRFQTLL, IGINITRFQTLLALHRS, TRFQTLLALHRSYLTPG, LLALHRSYLTPGDSSSG, RSYLTPGDSSSGWTAGA, PGDSSSGWTAGAAAYYV, SGWTAGAAAYYVGYLQP, GAAAYYVGYLQPRTFLL, YVGYLQPRTFLLKYNEN, QPRTFLLKYNENGTITD, LLKYNENGTITDAVDCA, ENGTITDAVDCALDPLS, TDAVDCALDPLSETKCT, CALDPLSETKCTLKSFT, LSETKCTLKSFTVEKGI, CTLKSFTVEKGIYQTSN, FTVEKGIYQTSNFRVQP, GIYQTSNFRVQPTESIV, SNFRVQPTESIVRFPNI, QPTESIVRFPNITNLCP, IVRFPNITNLCPFGEVF, NITNLCPFGEVFNATRF, CPFGEVFNATRFASVYA, VFNATRFASVYAWNRKR, RFASVYAWNRKRISNCV, YAWNRKRISNCVADYSV, KRISNCVADYSVLYNSA, CVADYSVLYNSASFSTF, SVLYNSASFSTFKCYGV, SASFSTFKCYGVSPTKL, TFKCYGVSPTKLNDLCF, GVSPTKLNDLCFTNVYA, KLNDLCFTNVYADSFVI, CFTNVYADSFVIRGDEV, YADSFVIRGDEVRQIAP, VIRGDEVRQIAPGQTGK, EVRQIAPGQTGKIADYN, APGQTGKIADYNYKLPD, GKIADYNYKLPDDFTGC, YNYKLPDDFTGCVIAWN, PDDFTGCVIAWNSNNLD, GCVIAWNSNNLDSKVGG, WNSNNLDSKVGGNYNYL, LDSKVGGNYNYLYRLFR, GGNYNYLYRLFRKSNLK, YLYRLFRKSNLKPFERD, FRKSNLKPFERDISTEI, LKPFERDISTEIYQAGS, RDISTEIYQAGSTPCNG, EIYQAGSTPCNGVEGFN, GSTPCNGVEGFNCYFPL, NGVEGFNCYFPLQSYGF, FNCYFPLQSYGFQPTNG, PLQSYGFQPTNGVGYQP, GFQPTNGVGYQPYRVVV, NGVGYQPYRVVVLSFEL, QPYRVVVLSFELLHAPA, VVLSFELLHAPATVCGP, ELLHAPATVCGPKKSTN, PATVCGPKKSTNLVKNK, GPKKSTNLVKNKCVNFN, TNLVKNKCVNFNFNGLT, NKCVNFNFNGLTGTGVL, FNFNGLTGTGVLTESNK, LTGTGVLTESNKKFLPF, VLTESNKKFLPFQQFGR, NKKFLPFQQFGRDIADT, PFQQFGRDIADTTDAVR, GRDIADTTDAVRDPQTL, DTTDAVRDPQTLEILDI, VRDPQTLEILDITPCSF, TLEILDITPCSFGGVSV, DITPCSFGGVSVITPGT, SFGGVSVITPGTNTSNQ, SVITPGTNTSNQVAVLY, GTNTSNQVAVLYQDVNC, NQVAVLYQDVNCTEVPV, LYQDVNCTEVPVAIHAD, NCTEVPVAIHADQLTPT, PVAIHADQLTPTWRVYS, ADQLTPTWRVYSTGSNV, PTWRVYSTGSNVFQTRA, YSTGSNVFQTRAGCLIG, NVFQTRAGCLIGAEHVN, RAGCLIGAEHVNNSYEC, IGAEHVNNSYECDIPIG, VNNSYECDIPIGAGICA, ECDIPIGAGICASYQTQ, IGAGICASYQTQTNSPR, CASYQTQTNSPRRARSV, TQTNSPRRARSVASQSI, PRRARSVASQSIIAYTM, SVASQSIIAYTMSLGAE, SIIAYTMSLGAENSVAY, TMSLGAENSVAYSNNSI, AENSVAYSNNSIAIPTN, AYSNNSIAIPTNFTISV, SIAIPTNFTISVTTEIL, TNFTISVTTEILPVSMT, SVTTEILPVSMTKTSVD, ILPVSMTKTSVDCTMYI, MTKTSVDCTMYICGDST, VDCTMYICGDSTECSNL, YICGDSTECSNLLLQYG, STECSNLLLQYGSFCTQ, NLLLQYGSFCTQLNRAL, YGSFCTQLNRALTGIAV, TQLNRALTGIAVEQDKN, ALTGIAVEQDKNTQEVF, AVEQDKNTQEVFAQVKQ, KNTQEVFAQVKQIYKTP, VFAQVKQIYKTPPIKDF, KQIYKTPPIKDFGGFNF, TPPIKDFGGFNFSQILP, DFGGFNFSQILPDPSKP, NFSQILPDPSKPSKRSF, LPDPSKPSKRSFIEDLL, KPSKRSFIEDLLFNKVT, SFIEDLLFNKVTLADAG, LLFNKVTLADAGFIKQY, VTLADAGFIKQYGDCLG, AGFIKQYGDCLGDIAAR, QYGDCLGDIAARDLICA, LGDIAARDLICAQKFNG, ARDLICAQKFNGLTVLP, CAQKFNGLTVLPPLLTD, NGLTVLPPLLTDEMIAQ, LPPLLTDEMIAQYTSAL, TDEMIAQYTSALLAGTI, AQYTSALLAGTITSGWT, ALLAGTITSGWTFGAGA, TITSGWTFGAGAALQIP, WTFGAGAALQIPFAMQM, GAALQIPFAMQMAYRFN, IPFAMQMAYRFNGIGVT, QMAYRFNGIGVTQNVLY, FNGIGVTQNVLYENQKL, VTQNVLYENQKLIANQF, LYENQKLIANQFNSAIG, KLIANQFNSAIGKIQDS, QFNSAIGKIQDSLSSTA, IGKIQDSLSSTASALGK, DSLSSTASALGKLQDVV, TASALGKLQDVVNQNAQ, GKLQDVVNQNAQALNTL, VVNQNAQALNTLVKQLS, AQALNTLVKQLSSNFGA, TLVKQLSSNFGAISSVL, LSSNFGAISSVLNDILS, GAISSVLNDILSRLDKV, VLNDILSRLDKVEAEVQ, LSRLDKVEAEVQIDRLI, KVEAEVQIDRLITGRLQ, VQIDRLITGRLQSLQTY, LITGRLQSLQTYVTQQL, LQSLQTYVTQQLIRAAE, TYVTQQLIRAAEIRASA, QLIRAAEIRASANLAAT, AEIRASANLAATKMSEC, SANLAATKMSECVLGQS, ATKMSECVLGQSKRVDF, ECVLGQSKRVDFCGKGY, QSKRVDFCGKGYHLMSF, DFCGKGYHLMSFPQSAP, GYHLMSFPQSAPHGVVF, SFPQSAPHGVVFLHVTY, APHGVVFLHVTYVPAQE, VFLHVTYVPAQEKNFTT, TYVPAQEKNFTTAPAIC, QEKNFTTAPAICHDGKA, TTAPAICHDGKAHFPRE, ICHDGKAHFPREGVFVS, KAHFPREGVFVSNGTHW, REGVFVSNGTHWFVTQR, VSNGTHWFVTQRNFYEP, HWFVTQRNFYEPQIITT, QRNFYEPQIITTDNTFV, EPQIITTDNTFVSGNCD, TTDNTFVSGNCDVVIGI, FVSGNCDVVIGIVNNTV, CDVVIGIVNNTVYDPLQ, GIVNNTVYDPLQPELDS, TVYDPLQPELDSFKEEL, LQPELDSFKEELDKYFK, DSFKEELDKYFKNHTSP, ELDKYFKNHTSPDVDLG, FKNHTSPDVDLGDISGI, SPDVDLGDISGINASVV, LGDISGINASVVNIQKE, GINASVVNIQKEIDRLN, VVNIQKEIDRLNEVAKN, KEIDRLNEVAKNLNESL, LNEVAKNLNESLIDLQE, KNLNESLIDLQELGKYE, SLIDLQELGKYEQYIKW, QELGKYEQYIKWPWYIW, YEQYIKWPWYIWLGFIA, KWPWYIWLGFIAGLIAI, IWLGFIAGLIAIVMVTI, IAGLIAIVMVTIMLCCM, AIVMVTIMLCCMTSCCS, TIMLCCMTSCCSCLKGC, CMTSCCSCLKGCCSCGS, CSCLKGCCSCGSCCKFD, GCCSCGSCCKFDEDDSE, GSCCKFDEDDSEPVLKG, FDEDDSEPVLKGVKLHY')

    ## peptide occurrence (slider)
    @output
    @render.ui
    def iters_slider_po():
        m = input.npools_po()

        file_infos = input.file_peptides_po()
        if not file_infos:
            peptides = re.split(r'\W+', input.text_peptides_po())
        else:
            for file_info in file_infos:
                input_peptides = pd.read_csv(file_info['datapath'], header = None)
                peptides = list(input_peptides[0])

        n = len(peptides)

        if type(m) != int or type(n) != int:
            return 'Please enter an integer'
        
        range_iters = []
        for r in range(m):
            if comb(m, r)*0.8 >= n and comb(m, r+1)*0.8 >= n:
                range_iters.append(r)

        if len(range_iters) != 0:
            result_iters_po = ui.input_slider('iters_po',
                                              '',
                                              min = range_iters[0],
                                              max = range_iters[-1],
                                              value = range_iters[0],
                                              step = 1)
            return result_iters_po
        else:
            return 'Warning: InvalidParameters. Increase number of pools or decrease number of peptides.'
    
    ## negative share (calculation)
    @reactive.Calc
    def negshare_calculation_po():

        file_infos = input.file_peptides_po()
        if not file_infos:
            peptides = re.split(r'\W+', input.text_peptides_po())
        else:
            for file_info in file_infos:
                input_peptides = pd.read_csv(file_info['datapath'], header = None)
                peptides = list(input_peptides[0])

        n = len(peptides)

        m = input.npools_po()
        r = input.iters_po()
        ep_length = input.ep_length_po()

        if len(peptides) > 1:
            t, _ = cpp.how_many_peptides(peptides, ep_length)
            e = max(t, key=t.get)

            if type(m) == int and r < m:
                negshare_result_po = (m - r - e + 1)/m
                return negshare_result_po
            
    ## negative share (text)
    @output
    @render.text
    def negshare_po():
        if negshare_calculation_po() is not None:
            return negshare_calculation_po()
    
    ## peptide pooling (reactive)
    @reactive.Calc
    @reactive.event(input.generate_scheme_po)
    def pooling_scheme_react():

        m = input.npools_po()
        r = input.iters_po()
        ep_length = input.ep_length_po()

        file_infos = input.file_peptides_po()
        if not file_infos:
            lst = re.split(r'\W+', input.text_peptides_po())
        else:
            for file_info in file_infos:
                input_df=pd.read_csv(file_info['datapath'], header = None)
            lst = list(input_df[0])

        n = len(lst)

        if type(m) == int:
            if comb(m, r) >= n:
                if n < 200:
                    b, lines = cdp.bba(m = m, r = r, n = n)
                else:
                    b, lines = cdp.rcbba(m = m, r = r, n = n)

                pools, peptide_address = cpp.pooling(lst=lst, addresses=lines, n_pools=m)

                act_profile = cpp.epitope_pools_activation(peptide_address, lst, ep_length)

                check_results = pd.DataFrame(columns = ['Peptide', 'Address', 'Epitope', 'Act Pools',
                                        '# of act pools', '# of epitopes with these act pools', '# of peptides with these act pools'])

                with ui.Progress(min=1, max=len(lst)) as p:
                    p.set(message="Analyzing generated scheme...", detail="This may take a while...")

                    for peptide in lst:
                        p.set(lst.index(peptide), message="Computing")
                        for i in range(len(peptide)):
                            ep = peptide[i:i+ep_length]
                            if len(ep) == ep_length:
                                act = cpp.pools_activation(pools, ep)
                                peps, eps = cpp.peptide_search(lst=lst,
                                                               act_profile=act_profile,
                                                               act_pools = act,
                                                               iters = r,
                                                               n_pools = m,
                                                               regime = 'without dropouts')
                                right_pep = str(peptide in peps)
                                right_ep = str(ep in eps)
                                row = {'Peptide':peptide,
                                       'Address':str(peptide_address[peptide]),
                                       'Epitope':ep,
                                       'Act Pools':str(sorted(list(act))),
                                       '# of act pools':len(act),
                                       '# of epitopes with these act pools':len(eps),
                                       '# of peptides with these act pools':len(peps)}
                                check_results = pd.concat([check_results, pd.DataFrame(row, index = [0])])
        
                if len(check_results[check_results['# of peptides with these act pools'] > 2]) > 0 :
                    m = ui.modal(
                    "Peptides in your peptide set have variable overlap lengths. Potentially, this may hinder result interpretation.",
                    title="Warning",
                    easy_close=True,
                    footer=None)
                    ui.modal_show(m)
        
                return [b, pools, peptide_address, check_results, lst]
            else:
                return 'Impossible set of parameters'
        else:
            return 'Impossible set of parameters'
        
    ## notification about resulting scheme (text)
    @output
    @render.ui
    @reactive.event(input.generate_scheme_po)
    def scheme_result_po():
        if pooling_scheme_react() == 'Impossible set of parameters':
            return 'Impossible set of parameters'
        else:
            res = ui.tooltip(ui.span("Simulation results ", question_circle_fill),
                "Results of the experiment for any possible epitope of provided length",
                placement="right",
                id="card_tooltip1_r_po",)
            return res
    
    ## send generated scheme, notification that it was sent
    @output
    @render.text()
    @reactive.event(input.send_po)
    def notification_po():
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

        n_pools = len(pools_df)
        iters = input.iters_po()
        n_peptides = len(pooling_scheme_react()[4])

        send_dataframe_via_email([input.letter_po()], [pools_df, peptide_address_df, check_results_df, [n_pools, iters, n_peptides]])
        return 'Sent!'

    ## resulting balance (plot)
    @output
    @render.plot
    @reactive.event(input.generate_scheme_po)
    def balance_plot_po():
        result = pooling_scheme_react()
        if result != 'Impossible set of parameters':
            b = result[0]

            if b is not None:
                m = input.npools_po()
                r = input.iters_po()
                n = len(result[4])

                perfect = [ceil(n*r/m)]*m

                fig, axes = plt.subplots(figsize = (6, 4))
                sn.scatterplot(y = b, x = list(range(1, m+1)), label = 'balance',
                            facecolor = 'none', s = 50, edgecolor = '#015396', linewidth = 2)
                sn.lineplot(y = perfect, x = list(range(1, m+1)),
                            color = '#E2E2ED', linestyle = '--', label = 'perfect balance', linewidth = 2)

                plt.ylim(0, max(b) + 10)
                plt.xlim(0, m + 1)
                plt.xticks(range(1, m + 1))
                plt.yticks([int(tick) for tick in plt.yticks()[0]])

                plt.ylabel('Number of peptides per pool')
                plt.xlabel('Pool index')
                sn.despine()

                return fig

    ## resulting scheme (dataframe)
    @output
    @render.data_frame
    @reactive.event(input.generate_scheme_po)
    def check_results_po():
        if pooling_scheme_react() != 'Impossible set of parameters':
            table = pooling_scheme_react()[3]
            fr_table = int(np.floor(len(table)*0.03))
            return table.sample(frac = 1).head(fr_table)
    ## download resulting scheme
    @session.download(filename='Results.zip')
    def download_po():
        if pooling_scheme_react() != 'Impossible set of parameters':
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

# STL files

    @output
    @render.image
    def peptide_plate_scheme_stl():
        img: ImgData = {"src": f"{path_user}www/peptide_plate_scheme.svg"}
        return img

    ## generation of STL files (reactive)
    @reactive.Calc
    @reactive.event(input.generate_stl)
    def stl_files_react():

        file_infos = input.pools_df_stl()
        if not file_infos:
            return
        for file_info in file_infos:
            pools_stl=pd.read_csv(file_info["datapath"], sep = "\t", index_col = 'Pool')

        file_infos = input.peptide_plate_stl()
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

                if input.well24_stl():
                    coordinates = coordinates + [[input.rows_stl(), input.cols_stl()]]
        
                name = 'pool' + str(pool_N)
        
                m = cpp.stl_generator(input.rows_stl(), input.cols_stl(), input.length_stl(), input.width_stl(),
                                      input.thickness_stl(), input.hole_radius_stl(), input.x_offset_stl(),
                                      input.y_offset_stl(), input.well_spacing_stl(), coordinates, engine = 'manifold', marks = input.marks_stl())
                
                meshes_list[name] = m
                count_stl = count_stl +1

        return meshes_list
    
    ## notification that generation was complete
    @output
    @render.text()
    @reactive.event(input.generate_stl)
    def stl_generation():
        meshes_list = stl_files_react()
        return 'Generation is complete!'
    
    @output
    @render.text()
    @reactive.event(input.send_stl)
    def notification_stl():
        zip_buffer = cpp.zip_meshes(stl_files_react())
        send_zip_via_email(receivers_list=[input.letter_stl()], zip_buffer = zip_buffer)
        return 'Sent!'
    
    @output
    @session.download(filename='STL_files.zip')
    def download_stl():
        zip_buffer = cpp.zip_meshes(stl_files_react())
        yield zip_buffer.getvalue()

# Interpretation

    message_visible = False

    @render.data_frame
    def exp_example_out():
        if input.exp_example_in():  # Detect button clicks
            nonlocal message_visible
            message_visible = not message_visible # Toggle state

        if message_visible:
            table = pd.read_csv(f'{path_user}www/cell_f1.csv')
            return table
        else:
            return 

    ## activation model (calculation)
    @reactive.Calc
    @reactive.event(input.start_model_in, ignore_none = True)
    def model_fitting_react():
        file_infos = input.exp_in()
        if not file_infos:
            return
        for file_info in file_infos:
            exp_results=pd.read_csv(file_info["datapath"])

        file_infos = input.check_results_in()
        if not file_infos:
            return
        for file_info in file_infos:
            check_results=pd.read_csv(file_info["datapath"], sep = "\t")

        if input.negcontrol_in():
            neg_control = input.negcontrol_in()
            neg_control = neg_control.split(',')
            neg_control = np.array(neg_control).astype(float)
        else:
            neg_control = None

        lst = check_results[['Peptide']].drop_duplicates().values.flatten().tolist()
        n_pools = input.npools_in()
        iters  = input.iters_in()
        ep_length = input.ep_length_in()

        with ui.Progress(min=1, max=2) as p:
            p.set(message="Model fitting started.", detail="It will take some time...")
            cells = list(exp_results['Percentage'])
            inds = list(exp_results['Pool'])

            t, _ = cpp.how_many_peptides(lst, ep_length)
            e = max(t, key=t.get)
            neg_share = (n_pools - iters - e + 1)/n_pools

            model, fig, probs, n_c, pp, pn = cpp.activation_model(cells, n_pools, inds, neg_control, neg_share = neg_share)
            plt.xlabel('Activation signal')
            plt.title('Model')
            plt.legend()

        return fig, probs
    
    ## fit plot
    @output
    @render.plot
    def fitting_in():
        res = model_fitting_react()
        if not res:
            return
        fig_like, _ = res

        if isinstance(fig_like, Figure):
            fig = fig_like
        elif hasattr(fig_like, "figure"):
            fig = fig_like.figure
        elif isinstance(fig_like, np.ndarray) and fig_like.size:
            fig = fig_like.flat[0].figure
        else:
            return

        return fig
    
    ## peptide probabilities (calculation)
    @reactive.Calc
    def peptide_probs_react():

        res = model_fitting_react()
        if not res:
            return

        file_infos = input.check_results_in()
        if not file_infos:
            return
        for file_info in file_infos:
            sim=pd.read_csv(file_info["datapath"], sep = "\t")

        _, probs = model_fitting_react()
        peptide_probs = cpp.peptide_probabilities(sim, probs)
        return peptide_probs
    
    ## resulting peptides (calculation)
    @reactive.Calc
    def resulting_peptides_react():

        res = model_fitting_react()
        if not res:
            return

        file_infos = input.check_results_in()
        if not file_infos:
            return
        for file_info in file_infos:
            sim=pd.read_csv(file_info["datapath"], sep = "\t")

        _, probs = model_fitting_react()
        peptide_probs = peptide_probs_react()

        n_act_pools, message, most, possible = cpp.results_analysis(peptide_probs, probs, sim)
        
        return message, most, possible
    
    ## activated pools
    @output
    @render.text()
    def activated_pools_in():
        res = model_fitting_react()
        if not res:
            return

        act_pools = []
        _, probs = model_fitting_react()
        for i in range(len(probs)):
            if probs['assign'].iloc[i] < 0.5:
                act_pools.append(i)
        text = 'Activated pools: '
        for item in act_pools:
            text = text + str(item) + ', '
        return text[:-2]
    
    @output
    @render.text()
    def cognate_peptides_in():

        res = resulting_peptides_react()
        if not res:
            return

        message, most, possible = resulting_peptides_react()
        if len(most) != 0:
            text = message + '\nResulting peptides: '
            for p in most:
                text = text + p + ', '
            text = text[:-2]
            text = text + '\n\nOther possible peptides: '
            for p in possible:
                text = text + p + ', '
        elif len(possible) != 0 and len(most) == 0:
            text = message
            text = text + '\n\nPossible peptides: '
            for p in possible:
                text = text + p + ', '
        else:
            text = 'Nothing was found.   '
        return text
    
    ## poolplot
    @output
    @render.plot
    def poolplot_in():

        res = resulting_peptides_react()
        if not res:
            return

        file_infos = input.exp_in()
        if not file_infos:
            return
        for file_info in file_infos:
            exp_results=pd.read_csv(file_info["datapath"])
        
        cells = list(exp_results['Percentage'])
        inds = list(exp_results['Pool'])
        most = resulting_peptides_react()[1]

        _, probs = model_fitting_react()
        fig = cpp.poolplot(probs, cells, inds, most)
        return fig
    
    ## bubbleplot (widget)
    @output
    @render_plotly
    def bubbleplot_in():
        res = peptide_probs_react()
        if res is None or len(res) == 0:
            return go.Figure()
        
        # Defensive copy and cleaning
        df = res.copy()
        num_cols = df.select_dtypes(include=["number"]).columns

        # Replace infinities with NaN, then drop any rows that have NaNs in numeric columns
        df[num_cols] = df[num_cols].replace([np.inf, -np.inf], np.nan)
        df = df.dropna(subset=num_cols, how="any")

        # Optional: if everything was dropped, return empty figure instead of exploding
        if df.empty:
            return go.Figure()

        return cpp.hover_bubbleplot(df)
    
    ## download peptide probabilities
    @session.download(filename='Peptide_probs.tsv')
    def download_in():
        peptide_probs = peptide_probs_react()
        peptide_probs = peptide_probs.sort_values(by = 'Probability', ascending=False)
        yield peptide_probs.to_csv(sep = "\t", index=False)

    @output
    @render.text()
    @reactive.event(input.send_in)
    def notification_in():
        file_infos = input.check_results_in()
        if not file_infos:
            return
        for file_info in file_infos:
            sim=pd.read_csv(file_info["datapath"], sep = "\t")
            
        _, probs = model_fitting_react()
        peptide_probs = peptide_probs_react()
        note, most, possible = resulting_peptides_react()
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

        send_probs_via_email([input.letter_in()], peptide_probs, text[:-2])
        return 'Sent!'
    


app = App(app_ui, server, debug=True)

# Supporting functions

## to parse FASTA files
def parse_fasta(filepath):
    ids = []
    seqs = []
    with open(filepath, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ids.append(record.id)
            seqs.append(str(record.seq))
    return ids, seqs

## to download dataframes with zips
def dataframes_to_zip(dataframes, filenames):
    zip_buffer = BytesIO()

    with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED) as zip_file:
        for dataframe, filename in zip(dataframes, filenames):
            df_buffer = BytesIO()
            dataframe.to_csv(df_buffer, index=False, sep='\t')
            zip_file.writestr(filename, df_buffer.getvalue())
    zip_buffer.seek(0)

    return zip_buffer

## to send zip by email
def send_zip_via_email(receivers_list, zip_buffer):
    lines = []
    with open(f"{path_user}pass.txt", "r", encoding="utf-8") as f:
        for line in f:
            lines.append(line.strip())

    login = lines[0]
    host_address = lines[1]
    channel = int(lines[2])

    with smtplib.SMTP(host=host_address, port=channel, timeout=30) as server:
        server.ehlo()
        server.starttls()
        server.ehlo()
        #server.login(login, "")

        msg = MIMEMultipart()
        msg["From"] = login
        msg["To"] = ", ".join(receivers_list)
        msg["Subject"] = "copepodTCR: STL files"

        # ZIP attachment
        zip_buffer.seek(0)
        attachment_part = MIMEBase("application", "zip")
        attachment_part.set_payload(zip_buffer.read())
        encoders.encode_base64(attachment_part)
        attachment_part.add_header(
            "Content-Disposition",
            'attachment; filename="STL_files.zip"'
        )
        msg.attach(attachment_part)

        server.sendmail(
            from_addr=login,
            to_addrs=receivers_list,
            msg=msg.as_string()
        )


## to send dataframes by email
def send_dataframe_via_email(receivers_list, dataframe):
    lines = []
    with open(f"{path_user}pass.txt", "r", encoding="utf-8") as f:
        for line in f:
            lines.append(line.strip())

    login = lines[0]
    host_address = lines[1]
    channel = int(lines[2])

    with smtplib.SMTP(host=host_address, port=channel, timeout=30) as server:
        server.ehlo()
        server.starttls()
        server.ehlo()
        #server.login(login, "")

        msg = MIMEMultipart()
        msg["From"] = login
        msg["To"] = ", ".join(receivers_list)
        msg["Subject"] = "copepodTCR: combinatorial peptide pooling scheme"

        # Pools (TSV)
        csv_buffer = StringIO()
        dataframe[0].to_csv(csv_buffer, index=False, sep="\t")
        csv_str = csv_buffer.getvalue()
        part = MIMEBase("text", "tsv")
        part.set_payload(csv_str.encode("utf-8"))
        encoders.encode_base64(part)
        part.add_header("Content-Disposition", 'attachment; filename="pools.tsv"')
        msg.attach(part)

        # Peptide_address (TSV)
        csv_buffer = StringIO()
        dataframe[1].to_csv(csv_buffer, index=False, sep="\t")
        csv_str = csv_buffer.getvalue()
        part = MIMEBase("text", "tsv")
        part.set_payload(csv_str.encode("utf-8"))
        encoders.encode_base64(part)
        part.add_header("Content-Disposition", 'attachment; filename="peptide_address.tsv"')
        msg.attach(part)

        # Simulation (TSV)
        csv_buffer = StringIO()
        dataframe[2].to_csv(csv_buffer, index=False, sep="\t")
        csv_str = csv_buffer.getvalue()
        part = MIMEBase("text", "tsv")
        part.set_payload(csv_str.encode("utf-8"))
        encoders.encode_base64(part)
        part.add_header("Content-Disposition", 'attachment; filename="simulation_results.tsv"')
        msg.attach(part)

        # README.txt
        n_pools = dataframe[3][0]
        iters = dataframe[3][1]
        l = dataframe[3][2]
        readme_content = (
            f"n_pools = {n_pools}\n"
            f"peptide occurence across pools (to how many pools one peptide is added) = {iters}\n"
            f"number of peptides = {l}\n"
            "Thank you for using copepodTCR! Made at CSHL."
        )
        part = MIMEBase("text", "plain")
        part.set_payload(readme_content.encode("utf-8"))
        encoders.encode_base64(part)
        part.add_header("Content-Disposition", 'attachment; filename="README.txt"')
        msg.attach(part)

        server.sendmail(
            from_addr=login,
            to_addrs=receivers_list,
            msg=msg.as_string(),
        )

## to send peptide probabilities with email
def send_probs_via_email(receivers_list, dataframe, message):
    lines = []
    with open(f"{path_user}pass.txt", "r", encoding="utf-8") as f:
        for line in f:
            lines.append(line.strip())

    login = lines[0]
    host_address = lines[1]
    channel = int(lines[2])

    with smtplib.SMTP(host=host_address, port=channel, timeout=30) as server:
        server.ehlo()
        server.starttls()
        server.ehlo()
        #server.login(login, "")

        msg = MIMEMultipart()
        msg["From"] = login
        msg["To"] = ", ".join(receivers_list)
        msg["Subject"] = "copepodTCR: results interpretation"

        # Peptide probabilities (TSV)
        csv_buffer = StringIO()
        dataframe.to_csv(csv_buffer, index=False, sep="\t")
        csv_str = csv_buffer.getvalue()
        part = MIMEBase("text", "tsv")
        part.set_payload(csv_str.encode("utf-8"))
        encoders.encode_base64(part)
        part.add_header("Content-Disposition", 'attachment; filename="Peptide_probabilities.tsv"')
        msg.attach(part)

        # README.txt
        readme_content = f"{message}\nThank you for using copepodTCR! Made at CSHL."
        part = MIMEBase("text", "plain")
        part.set_payload(readme_content.encode("utf-8"))
        encoders.encode_base64(part)
        part.add_header("Content-Disposition", 'attachment; filename="README.txt"')
        msg.attach(part)

        server.sendmail(
            from_addr=login,
            to_addrs=receivers_list,
            msg=msg.as_string()
        )