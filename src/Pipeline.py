#!/usr/bin/env python3

# Python script to take variables and write Variables.R and run Omics_Notebook

import tkinter as tk
from tkinter import ttk, filedialog
import os
import subprocess
import sys

class GUI(tk.Frame):
  def __init__(self, parent):
    tk.Frame.__init__(self,parent)
    lbl = tk.Label(self, text="Enter variables below to configure Omics Notebook")
    lbl.grid(column=0, row=0, columnspan=3, pady=(20,20))
    # Project Name:
    nameProj_lbl = tk.Label(self, text="Project Name:") 
    nameProj_lbl.grid(column=0, row=1)
    nameProj_var = tk.StringVar(self, value="Analysis")
    nameProj = tk.Entry(self, width=20, textvariable=nameProj_var)
    nameProj.grid(column=1, row=1)
    # DIRECTORY
    fileDir_lbl = tk.Label(self, text="Files Directory:")
    fileDir_lbl.grid(column=0, row=2)
    fileDir_field = tk.Entry(self, width=60)
    fileDir_field.grid(column=1, row=2, columnspan=2)
    startDir = "/projectnb/cnsbomic" # change default based on install
    def fileDir_clicked():
      global fileDir
      fileDir = filedialog.askdirectory(initialdir=startDir, title='Choose directory') 
      fileDir_field.configure(textvariable=tk.StringVar(self, value=fileDir))
    fileDir_btn= tk.Button(self, text="Choose Directory", command= lambda : fileDir_clicked())
    fileDir_btn.grid(column=3, row=2)
    # FILE NAME
    nameAnn_lbl = tk.Label(self, text="Annotation file:")
    nameAnn_lbl.grid(column=0, row=5)
    nameAnn_field = tk.Entry(self, width=60)
    nameAnn_field.grid(column=1, row=5, columnspan=2)
    def nameAnn_clicked():
       global nameAnn
       if 'fileDir' in globals():
         start_dir = fileDir
       else:
         start_dir=startDir
       nameAnn = filedialog.askopenfilename(initialdir=start_dir, title='Choose file')
       nameAnn_field.configure(textvariable=tk.StringVar(self, value=nameAnn))
    nameAnn_btn= tk.Button(self, text="Choose Annotation File", command=lambda : nameAnn_clicked())
    nameAnn_btn.grid(column=3, row=5)
    # RUN PARAMETERS
    incDifferential_state = tk.BooleanVar()
    incDifferential_state.set(True)
    incDifferential = ttk.Checkbutton(self, text="Differential Analysis", var=incDifferential_state)
    incDifferential.grid(column=1, row=6, pady=(10,10))

    newcontrast_state = tk.BooleanVar()
    newcontrast_state.set(False)
    newcontrast = ttk.Checkbutton(self, text="New Contrast Only", var=newcontrast_state)
    newcontrast.grid(column=2, row=6, pady=(10,10))

    runXlsx_state = tk.BooleanVar()
    runXlsx_state.set(True)
    runXlsx = ttk.Checkbutton(self, text="Excel Report", var=runXlsx_state)
    runXlsx.grid(column=0, row=11, padx=(2,0))

    runEnrichr_state = tk.BooleanVar()
    runEnrichr_state.set(False)
    runEnrichr = ttk.Checkbutton(self, text="Run EnrichR", var=runEnrichr_state)
    runEnrichr.grid(column=0, row=10)

    runGSEA_state = tk.BooleanVar()
    runGSEA_state.set(True)
    isInteractive_state = tk.BooleanVar()
    isInteractive_state.set(True)

    runGSEA = ttk.Checkbutton(self, text="Run GSEA", var=runGSEA_state, command=lambda:isInteractive_state.set(isInteractive_state.get() and runGSEA_state.get()))
    runGSEA.grid(column=1, row=10)
    isInteractive = ttk.Checkbutton(self, text="Interactive (include EnrichmentMap)", var=isInteractive_state, command=lambda:isInteractive_state.set(isInteractive_state.get() and runGSEA_state.get()))
    isInteractive.grid(column=2, row=10, columnspan=1)

    normMethod_lbl = tk.Label(self, text="Normalization Method:")
    normMethod_lbl.grid(column=0, row=12, padx=(10,0),pady=(20,0))
    normMethod = ttk.Combobox(self)
    normMethod['values']=('loess', 'quantile')
    normMethod.current(0)
    normMethod.grid(column=1, row=12, pady=(20,0))

    heatcolors_lbl = tk.Label(self, text="Heatmap color scale:")
    heatcolors_lbl.grid(column=0, row=14)
    heatcolors = ttk.Combobox(self)
    heatcolors['values']=('viridis', 'RdYlBu', 'RdBu')
    heatcolors.current(0)
    heatcolors.grid(column=1, row=14)
    
    species_lbl = tk.Label(self, text="Species:")
    species_lbl.grid(column=0, row=15)
    species = ttk.Combobox(self)
    species['values']=('Human', 'Mouse', 'Other')
    species.current(0)
    species.grid(column=1, row=15)

    txtFolder_lbl = tk.Label(self, text="txt Folder report:")
    txtFolder_lbl.grid(column=0, row=16)
    txtFolder_lbl2 = tk.Label(self, text="Copy mqpar file into txt folder for complete results.")
    txtFolder_lbl2.grid(column=2, row=16, columnspan=1)
    txtFolder = ttk.Combobox(self)
    txtFolder['values']=('TRUE', 'FALSE')
    txtFolder.current(1)
    txtFolder.grid(column=1, row=16)

    # NUMERIC ENTRIES
    vcmd = (self.register(self.onValidate), '%d','%i','%P','%s','%S','%v','%V','%W')

    emapVar_lbl = tk.Label(self, text="EnrichmentMap Parameters:")
    emapVar_lbl.grid(column=2, columnspan=1, row=17, pady=(20,0))
    emapVar_p_lbl = tk.Label(self, text="P Value:")
    emapVar_p_lbl.grid(column=2, row=18)
    emapVar_p_start = tk.StringVar(self)
    emapVar_p_start.set("0.05")
    enrichmentmap_p_val = tk.Spinbox(self, from_=0.01, to=0.50, increment=0.01, textvariable=emapVar_p_start, validate='all', validatecommand=vcmd)
    enrichmentmap_p_val.grid(column=3, row=18)

    emapVar_q_lbl = tk.Label(self, text="Q Value:")
    emapVar_q_lbl.grid(column=2, row=19)
    emapVar_q_start = tk.StringVar(self)
    emapVar_q_start.set("0.30")
    enrichmentmap_q_val = tk.Spinbox(self, from_=0.01, to=0.50, increment=0.01, textvariable=emapVar_q_start, validate='all', validatecommand=vcmd)
    enrichmentmap_q_val.grid(column=3, row=19)

    cutoff_lbl = tk.Label(self, text="Cut Offs:")
    cutoff_lbl.grid(column=0, columnspan=1, row=17, pady=(20,0))
    cutoff_zero_lbl = tk.Label(self, text="Zero Percent:")
    cutoff_zero_lbl.grid(column=0, row=18)
    cutoff_zero_start = tk.StringVar(self)
    cutoff_zero_start.set("0.70")
    cutoff_zero_val = tk.Spinbox(self, from_=0.01, to=1.00, increment=0.01, textvariable=cutoff_zero_start, validate='all', validatecommand=vcmd)
    cutoff_zero_val.grid(column=1, row=18)

    cutoff_fdr_lbl = tk.Label(self, text="FDR:")
    cutoff_fdr_lbl.grid(column=0, row=19)
    cutoff_fdr_start = tk.StringVar(self)
    cutoff_fdr_start.set("0.05")
    cutoff_fdr_val = tk.Spinbox(self, from_=0.01, to=0.50, increment=0.01, textvariable=cutoff_fdr_start, validate='all', validatecommand=vcmd)
    cutoff_fdr_val.grid(column=1, row=19)

    # SAVE VARIABLES FILE AND CLOSE COMMAND
    def clicked_button():
      dir_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__)))
      analysis_dir = str(tk.StringVar(self, value=fileDir).get()) 
      
      variable_file = analysis_dir + "/Variables.R"
      outputfile = open(variable_file, "w")
      
      outputfile.write("project_name <- '" + nameProj_var.get() +"';\n")
      outputfile.write("working_dir <- '" + analysis_dir + "'; # Directory for analysis and where data is \n")
      outputfile.write("annotation_filename <- '" + str(tk.StringVar(self, value=nameAnn).get()) + "';\n")

      outputfile.write("saveXlsx <- " + str(runXlsx_state.get()).upper() + ";\n")
      outputfile.write("enrichr_section <- " + str(runEnrichr_state.get()).upper() + ";\n")
      outputfile.write("gsea_section <- " + str(runGSEA_state.get()).upper() + ";\n")
      outputfile.write("enrichment_map <- " + str(isInteractive_state.get()).upper() + ";\n")

      outputfile.write("map_color <- '" + str(heatcolors.get()) + "';\n")
      outputfile.write("species <- '" + str(species.get()) + "';\n")
      outputfile.write("norm_method <- '" + str(normMethod.get()) + "';\n")
      outputfile.write("newcontrastonly <- " + str(newcontrast_state.get()).upper() + ";\n")
      outputfile.write("txtFolder <- " + str(txtFolder.get()).upper() + ";\n")
      outputfile.write("isInteractive <- " + str(isInteractive_state.get()).upper() + ";\n")
      outputfile.write("notebook_dir <- '" + dir_path + "';\n")
      outputfile.write("runDifferential <- " + str(incDifferential_state.get()).upper() + ";\n")

      outputfile.write("enrichmentmap_p_val <- " + str(enrichmentmap_p_val.get()) + ";\n")
      outputfile.write("enrichmentmap_q_val <- " + str(enrichmentmap_q_val.get()) + ";\n")
    
      outputfile.write("zero_percent <- " + str(cutoff_zero_val.get()) + ";\n")
      outputfile.write("adjpcutoff <- " + str(cutoff_fdr_val.get()) + ";\n")
      outputfile.write("BUSCC <- TRUE;\n")

      outputfile.close()
      lbl.configure(text="Variables file created.")
      sys.exit(analysis_dir)
    
    btn= tk.Button(self, text="Enter", command=clicked_button)
    btn.grid(column=0, row=20, columnspan=3, pady=(20,0))

  # Validate numbers command
  def onValidate(self, d, i, P, s, S, v, V, W):
    minval=int(self.nametowidget(W).config('from')[4])
    maxval=int(self.nametowidget(W).config('to')[4])
    if P.isdigit() and P in range(minval, maxval):
      return True
    else:
      return False

root=tk.Tk()
root.title('Omics Notebook Setup')
GUI(root).pack(fill="both", expand=True)
root.mainloop()


