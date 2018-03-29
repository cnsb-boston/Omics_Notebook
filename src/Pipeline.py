#!/usr/bin/env python3

# Python script to take variables and write Variables.R and pipeline.sh scripts to run Omics_Notebook

import tkinter as tk
from tkinter import ttk, filedialog
import os
import subprocess

app = tk.Tk()
app.title('Omics Notebook Setup')
app.geometry('900x600')

lbl = tk.Label(app, text="Enter variables below to configure Omics Notebook")
lbl.grid(column=0, row=0, columnspan=3, pady=(20,20))

# DIRECTORY
fileDir_lbl = tk.Label(app, text="Files Directory:")
fileDir_lbl.grid(column=0, row=2)
fileDir_field = tk.Entry(app, width=60)
fileDir_field.grid(column=1, row=2, columnspan=2)
def fileDir_clicked():
    global fileDir
    fileDir = filedialog.askdirectory(initialdir="/projectnb/cnsbomic", title='Choose directory')
    fileDir_field.configure(textvariable=tk.StringVar(app, value=fileDir))
fileDir_btn= tk.Button(app, text="Choose Directory", command=lambda : fileDir_clicked())
fileDir_btn.grid(column=3, row=2)

# DATA TYPES
isprot_state = tk.BooleanVar()
isprot_state.set(True)
isprot = ttk.Checkbutton(app, text="Proteomics", var=isprot_state)
isprot.grid(column=0, row=3, pady=(10,10))

isphos_state = tk.BooleanVar()
isphos_state.set(True)
isphos = ttk.Checkbutton(app, text="Phospho", var=isphos_state)
isphos.grid(column=1, row=3, pady=(10,10))

isPTM2_state = tk.BooleanVar()
isPTM2_state.set(False)
isPTM2 = ttk.Checkbutton(app, text="PTM 2", var=isPTM2_state, state="disabled")
isPTM2.grid(column=2, row=3, pady=(10,10))

isPTM3_state = tk.BooleanVar()
isPTM3_state.set(False)
isPTM3 = ttk.Checkbutton(app, text="PTM 3", var=isPTM3_state, state="disabled")
isPTM3.grid(column=3, row=3, pady=(10,10))

# FILE NAMES
nameAnn_lbl = tk.Label(app, text="Annotation file:")
nameAnn_lbl.grid(column=0, row=5)
nameAnn_var = tk.StringVar(app, value="Annotation.xlsx")
nameAnn = tk.Entry(app, width=20, textvariable=nameAnn_var)
nameAnn.grid(column=1, row=5)

fileProt_lbl = tk.Label(app, text="Protein Groups file:")
fileProt_lbl.grid(column=0, row=6)
fileProt_var = tk.StringVar(app, value="proteinGroups.txt")
fileProt = tk.Entry(app, width=20, textvariable=fileProt_var)
fileProt.grid(column=1, row=6)

filePhos_lbl = tk.Label(app, text="Phospho Sites file:")
filePhos_lbl.grid(column=0, row=7)
filePhos_var = tk.StringVar(app, value="Phospho (STY)Sites.txt")
filePhos = tk.Entry(app, width=20, textvariable=filePhos_var)
filePhos.grid(column=1, row=7, pady=(0,20))

# RUN PARAMETERS
runEnrichr_lbl = tk.Label(app, text="Run EnrichR:")
runEnrichr_lbl.grid(column=0, row=10)
runEnrichr = ttk.Combobox(app)
runEnrichr['values']=('TRUE', 'FALSE')
runEnrichr.current(0)
runEnrichr.grid(column=1, row=10)

runGSEA_lbl = tk.Label(app, text="Run GSEA/EnrichmentMap:")
runGSEA_lbl.grid(column=0, row=11, padx=(10,0))
runGSEA = ttk.Combobox(app)
runGSEA['values']=('TRUE', 'FALSE')
runGSEA.current(0)
runGSEA.grid(column=1, row=11)

txtFolder_lbl = tk.Label(app, text="txt Folder report:")
txtFolder_lbl.grid(column=0, row=16)
txtFolder_lbl2 = tk.Label(app, text="Copy mqpar file into txt folder for complete results.")
txtFolder_lbl2.grid(column=2, row=16, columnspan=1)
txtFolder = ttk.Combobox(app)
txtFolder['values']=('TRUE', 'FALSE')
txtFolder.current(1)
txtFolder.grid(column=1, row=16)

normMethod_lbl = tk.Label(app, text="Normalization Method:")
normMethod_lbl.grid(column=0, row=12, pady=(20,0))
normMethod = ttk.Combobox(app)
normMethod['values']=('loess', 'quantile')
normMethod.current(0)
normMethod.grid(column=1, row=12, pady=(20,0))

newcontrast_lbl = tk.Label(app, text="New contrast only?")
newcontrast_lbl.grid(column=0, row=15, pady=(20,0))
newcontrast = ttk.Combobox(app)
newcontrast['values']=('TRUE', 'FALSE')
newcontrast.current(1)
newcontrast.grid(column=1, row=15, pady=(20,0))

heatcolors_lbl = tk.Label(app, text="Heatmap color scale:")
heatcolors_lbl.grid(column=0, row=14)
heatcolors = ttk.Combobox(app)
heatcolors['values']=('viridis', 'RdYlBu', 'RdBu')
heatcolors.current(0)
heatcolors.grid(column=1, row=14)

isInteractive_state = tk.BooleanVar()
isInteractive_state.set(True)
isInteractive = ttk.Checkbutton(app, text="Interactive", var=isInteractive_state, state='disabled')
isInteractive.grid(column=0, row=18, pady=(10,10))

# SAVE VARIABLES FILE AND CLOSE
def clicked():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    variable_file = dir_path + "/Variables.R"
    outputfile = open(variable_file, "w")
    outputfile.write("working.dir <- '" + str(tk.StringVar(app, value=fileDir).get()) + "'; # Directory for analysis and where data is \n")
    outputfile.write("isthereProteo <- " + str(isprot_state.get()).upper() + ";\n")
    outputfile.write("istherePhospho <- " + str(isphos_state.get()).upper() + ";\n")
    
    outputfile.write("maxq.global.filename <- '" + fileProt_var.get() +"';\n")
    outputfile.write("maxq.phospho.filename <- '" + filePhos_var.get() +"';\n")
    outputfile.write("annotation.filename <- '" + nameAnn_var.get() +"';\n")

    outputfile.write("enrichr_section <- " + str(runEnrichr.get()).upper() + ";\n")
    outputfile.write("gsea_section <- " + str(runGSEA.get()).upper() + ";\n")
    outputfile.write("enrichment_map <- " + str(isInteractive_state.get()).upper() + ";\n")

    outputfile.write("map.color <- '" + str(heatcolors.get()) + "';\n")
    outputfile.write("norm.method <- '" + str(normMethod.get()) + "';\n")
    outputfile.write("newcontrastonly <- " + str(newcontrast.get()).upper() + ";\n")
    outputfile.write("txtFolder <- " + str(txtFolder.get()).upper() + ";\n")
    outputfile.write("isInteractive <- " + str(isInteractive_state.get()).upper() + ";\n")
    outputfile.write("notebook_dir <- '" + dir_path + "';\n")

    outputfile.close()
    lbl.configure(text="Variables file created.")
    
    if isInteractive_state.get()==True:
        pipeline_file = dir_path + "/pipeline.sh"
        outputfile = open(pipeline_file, "w")
        outputfile.write("#!/bin/bash\n")
        outputfile.write("# run Interactive pipeline\n")
        outputfile.write("qrsh -P cnsbomic -pe omp 10 -l gpus=0.2 -b y \"" + dir_path + "/pipelineInteractive.sh\"")
        outputfile.close()
    elif isInteractive_state.get()==False:
        pipeline_file = dir_path + "/pipeline.sh"
        outputfile = open(pipeline_file, "w")
        outputfile.write("#!/bin/bash\n")
        outputfile.write("# run pipeline as batch job\n")
        outputfile.write("qsub -P cnsbomic -pe omp 16 -l gpus=0.125 -b y \"" + dir_path + "/pipelineBatch.sh\"")
        outputfile.close()

    quit()

btn= tk.Button(app, text="Enter", command=clicked)
btn.grid(column=0, row=20, columnspan=3, pady=(20,0))

app.mainloop()
