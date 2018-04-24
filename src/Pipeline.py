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

nameProj_lbl = tk.Label(app, text="Project Name:")
nameProj_lbl.grid(column=0, row=1)
nameProj_var = tk.StringVar(app, value="Analysis")
nameProj = tk.Entry(app, width=20, textvariable=nameProj_var)
nameProj.grid(column=1, row=1)

# DIRECTORY
fileDir_lbl = tk.Label(app, text="Files Directory:")
fileDir_lbl.grid(column=0, row=2)
fileDir_field = tk.Entry(app, width=60)
fileDir_field.grid(column=1, row=2, columnspan=2)
def fileDir_clicked():
    global fileDir
    fileDir = filedialog.askdirectory(initialdir="/projectnb/cnsbomic", title='Choose directory') # change default based on install
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
incDifferential_state = tk.BooleanVar()
incDifferential_state.set(True)
incDifferential = ttk.Checkbutton(app, text="Differential Analysis", var=incDifferential_state)
incDifferential.grid(column=2, row=5)

newcontrast_state = tk.BooleanVar()
newcontrast_state.set(False)
newcontrast = ttk.Checkbutton(app, text="New Contrast Only", var=newcontrast_state)
newcontrast.grid(column=2, row=6)

runXlsx_state = tk.BooleanVar()
runXlsx_state.set(True)
runXlsx = ttk.Checkbutton(app, text="Excel Report", var=runXlsx_state)
runXlsx.grid(column=0, row=11, padx=(2,0))

runEnrichr_state = tk.BooleanVar()
runEnrichr_state.set(False)
runEnrichr = ttk.Checkbutton(app, text="Run EnrichR", var=runEnrichr_state)
runEnrichr.grid(column=0, row=10)

runGSEA_state = tk.BooleanVar()
runGSEA_state.set(True)
isInteractive_state = tk.BooleanVar()
isInteractive_state.set(True)

runGSEA = ttk.Checkbutton(app, text="Run GSEA", var=runGSEA_state, command=lambda:isInteractive_state.set(isInteractive_state.get() and runGSEA_state.get()))
runGSEA.grid(column=1, row=10)
isInteractive = ttk.Checkbutton(app, text="Interactive (include EnrichmentMap)", var=isInteractive_state, command=lambda:isInteractive_state.set(isInteractive_state.get() and runGSEA_state.get()))
isInteractive.grid(column=2, row=10, columnspan=1)

normMethod_lbl = tk.Label(app, text="Normalization Method:")
normMethod_lbl.grid(column=0, row=12, padx=(10,0),pady=(20,0))
normMethod = ttk.Combobox(app)
normMethod['values']=('loess', 'quantile')
normMethod.current(0)
normMethod.grid(column=1, row=12, pady=(20,0))

heatcolors_lbl = tk.Label(app, text="Heatmap color scale:")
heatcolors_lbl.grid(column=0, row=14)
heatcolors = ttk.Combobox(app)
heatcolors['values']=('viridis', 'RdYlBu', 'RdBu')
heatcolors.current(0)
heatcolors.grid(column=1, row=14)

txtFolder_lbl = tk.Label(app, text="txt Folder report:")
txtFolder_lbl.grid(column=0, row=16)
txtFolder_lbl2 = tk.Label(app, text="Copy mqpar file into txt folder for complete results.")
txtFolder_lbl2.grid(column=2, row=16, columnspan=1)
txtFolder = ttk.Combobox(app)
txtFolder['values']=('TRUE', 'FALSE')
txtFolder.current(1)
txtFolder.grid(column=1, row=16)

# SAVE VARIABLES FILE AND CLOSE
def clicked():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    variable_file = dir_path + "/Variables.R"
    outputfile = open(variable_file, "w")
    outputfile.write("project_name <- '" + nameProj_var.get() +"';\n")
    outputfile.write("working_dir <- '" + str(tk.StringVar(app, value=fileDir).get()) + "'; # Directory for analysis and where data is \n")
    outputfile.write("isthereProteo <- " + str(isprot_state.get()).upper() + ";\n")
    outputfile.write("istherePhospho <- " + str(isphos_state.get()).upper() + ";\n")
    
    outputfile.write("maxq_global_filename <- '" + fileProt_var.get() +"';\n")
    outputfile.write("maxq_phospho_filename <- '" + filePhos_var.get() +"';\n")
    outputfile.write("annotation_filename <- '" + nameAnn_var.get() +"';\n")

    outputfile.write("saveXlsx <- " + str(runXlsx_state.get()).upper() + ";\n")
    outputfile.write("enrichr_section <- " + str(runEnrichr_state.get()).upper() + ";\n")
    outputfile.write("gsea_section <- " + str(runGSEA_state.get()).upper() + ";\n")
    outputfile.write("enrichment_map <- " + str(isInteractive_state.get()).upper() + ";\n")

    outputfile.write("map_color <- '" + str(heatcolors.get()) + "';\n")
    outputfile.write("norm_method <- '" + str(normMethod.get()) + "';\n")
    outputfile.write("newcontrastonly <- " + str(newcontrast_state.get()).upper() + ";\n")
    outputfile.write("txtFolder <- " + str(txtFolder.get()).upper() + ";\n")
    outputfile.write("isInteractive <- " + str(isInteractive_state.get()).upper() + ";\n")
    outputfile.write("notebook_dir <- '" + dir_path + "';\n")
    outputfile.write("runDifferential <- " + str(incDifferential_state.get()).upper() + ";\n")

    outputfile.close()
    lbl.configure(text="Variables file created.")
    
    if isInteractive_state.get()==True:
        pipeline_file = dir_path + "/pipeline.sh"
        outputfile = open(pipeline_file, "w")
        outputfile.write("#!/bin/bash\n")
        outputfile.write("# run Interactive pipeline\n")
        outputfile.write("source \"" + dir_path + "/pipelineInteractive.sh\"")
        outputfile.close()
    elif isInteractive_state.get()==False:
        pipeline_file = dir_path + "/pipeline.sh"
        outputfile = open(pipeline_file, "w")
        outputfile.write("#!/bin/bash\n")
        outputfile.write("# run pipeline as batch job\n")
        outputfile.write("source \"" + dir_path + "/pipelineBatch.sh\"")
        outputfile.close()

    quit()

btn= tk.Button(app, text="Enter", command=clicked)
btn.grid(column=0, row=20, columnspan=3, pady=(20,0))

app.mainloop()
