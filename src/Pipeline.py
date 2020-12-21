#!/usr/bin/env python3

# Python script to take variables and write Paramters.R and run Omics_Notebook

# Load modules
import tkinter as tk
from tkinter import ttk, filedialog
import os
import subprocess
import sys

############################################

# Define GUI instance
class GUI(tk.Frame):
  def __init__(self, parent):
    tk.Frame.__init__(self,parent)
    
    ############################################
    ##
    ##  SET FILE PATHS FOR LOCAL INSTALL
    ##
    ############################################

    # R libraries path
    libraries_path = '/project/cnsbomic/Tools/R3.6' # change default based on install

    # Pandoc path
    #pandoc_path = SET IF NEEDED

    # Inherit paths
    inherit_paths = "TRUE"

    startDir = "/projectnb/cnsbomic" # change default based on install

    ############################################
    
    # WINDOW TITLE AND DIMENSIONS
    lbl = tk.Label(self, text="Enter variables below to configure Omics Notebook")
    lbl.grid(column=0, row=0, columnspan=3, pady=(20,20))
    
    # PROJECT NAME
    nameProj_lbl = tk.Label(self, text="Project Name:") 
    nameProj_lbl.grid(column=0, row=1, sticky=tk.E)
    nameProj_var = tk.StringVar(self, value="Analysis")
    nameProj = tk.Entry(self, width=20, textvariable=nameProj_var)
    nameProj.grid(column=1, row=1, sticky=tk.W)
    
    # DIRECTORY
    fileDir_lbl = tk.Label(self, text="Files Directory:")
    fileDir_lbl.grid(column=0, row=2, sticky=tk.E)
    fileDir_field = tk.Entry(self, width=60, state="disabled")
    fileDir_field.grid(column=1, row=2, columnspan=3, sticky=tk.W)
    if os.path.isdir(os.path.normpath(startDir) ):
      startDir = os.path.normpath(startDir)
    else :
      startDir = os.path.normpath("./")
    def fileDir_clicked():
      global fileDir
      fileDir = filedialog.askdirectory(initialdir=startDir, title='Choose directory') 
      fileDir_field.configure(textvariable=tk.StringVar(self, value=fileDir))
    fileDir_btn= tk.Button(self, text="Choose Directory", command= lambda : fileDir_clicked())
    fileDir_btn.grid(column=4, row=2, sticky=tk.W)
    
    # ANNOTATION FILE NAME
    nameAnn_lbl = tk.Label(self, text="Annotation file:")
    nameAnn_lbl.grid(column=0, row=5, sticky=tk.E)
    nameAnn_field = tk.Entry(self, width=60, state="disabled")
    nameAnn_field.grid(column=1, row=5, columnspan=3, sticky=tk.W)
    def nameAnn_clicked():
       global nameAnn
       if 'fileDir' in globals():
         start_dir = fileDir
       else:
         start_dir=startDir
       nameAnn = filedialog.askopenfilename(initialdir=start_dir, title='Choose file')
       nameAnn_field.configure(textvariable=tk.StringVar(self, value=nameAnn))
    nameAnn_btn= tk.Button(self, text="Choose Annotation File", command=lambda : nameAnn_clicked())
    nameAnn_btn.grid(column=4, row=5, sticky=tk.W)
    
    # RUN PARAMETERS
    incDifferential_state = tk.BooleanVar()
    incDifferential_state.set(True)
    incDifferential = ttk.Checkbutton(self, text="Differential Analysis", var=incDifferential_state)
    incDifferential.grid(column=1, row=6, pady=(10,10), sticky=tk.W)

    queryWeb_state = tk.BooleanVar()
    queryWeb_state.set(True)
    queryWeb = ttk.Checkbutton(self, text="Query Uniprot", var=queryWeb_state)
    queryWeb.grid(column=2, row=6, pady=(10,10), sticky=tk.W)

    runXlsx_state = tk.BooleanVar()
    runXlsx_state.set(True)
    runXlsx = ttk.Checkbutton(self, text="Excel Report", var=runXlsx_state)
    runXlsx.grid(column=3, row=6, pady=(10,10), sticky=tk.W)
    
    isIntHM_state = tk.BooleanVar()
    isIntHM_state.set(False)
    isIntHM = ttk.Checkbutton(self, text="Interactive Heatmap", var=isIntHM_state)
    isIntHM.grid(column=4, row=6, pady=(10,10), sticky=tk.W)
    
    #showSampleName_state = tk.BooleanVar()
    #showSampleName_state.set(True)
    #showSampleName = ttk.Checkbutton(self, text="Show Names", var=showSampleName_state)
    #showSampleName.grid(column=2, row=7, pady=(0,10), sticky=tk.W)
    
    isIntVl_state = tk.BooleanVar()
    isIntVl_state.set(True)
    isIntVl = ttk.Checkbutton(self, text="Interactive Volcano", var=isIntVl_state)
    isIntVl.grid(column=4, row=7, pady=(0,10), sticky=tk.W)

    runEnrichr_state = tk.BooleanVar()
    runEnrichr_state.set(True)
    runEnrichr = ttk.Checkbutton(self, text="Run EnrichR", var=runEnrichr_state)
    runEnrichr.grid(column=1, row=10, pady=(10,10), sticky=tk.W)
    
    runGSEA_state = tk.BooleanVar()
    runGSEA_state.set(True)
    runGSEA = ttk.Checkbutton(self, text="Run GSEA", var=runGSEA_state)
    runGSEA.grid(column=2, row=10, pady=(10,10), sticky=tk.W)
    
    runMOMENTA_state = tk.BooleanVar()
    runMOMENTA_state.set(True)
    runMOMENTA = ttk.Checkbutton(self, text="Run MOMENTA", var=runMOMENTA_state)
    runMOMENTA.grid(column=3, row=10, pady=(10,10), sticky=tk.W)
    
    runActPath_state = tk.BooleanVar()
    runActPath_state.set(False)
    runActPath = ttk.Checkbutton(self, text="Run ActivePathways", var=runActPath_state)
    runActPath.grid(column=4, row=10, pady=(10,10), sticky=tk.W)
    
    allcomparisons_state = tk.BooleanVar()
    allcomparisons_state.set(True)
    allcomparisons = ttk.Checkbutton(self, text="All Comparisons (if False, All to First)", var=allcomparisons_state)
    allcomparisons.grid(column=1, row=7, columnspan=1, pady=(10,10), sticky=tk.W)
    
    useSiteNorm_state = tk.BooleanVar()
    useSiteNorm_state.set(False)
    useSiteNorm = ttk.Checkbutton(self, text="Use Normalized to 1st Data", var=useSiteNorm_state)
    useSiteNorm.grid(column=3, row=11,columnspan=2, pady=(20,0), sticky=tk.W)
    
    incShinyNorm_state = tk.BooleanVar()
    incShinyNorm_state.set(False) # adjust when working
    incShinyNorm = ttk.Checkbutton(self, text="Interactive Normalization", var=incShinyNorm_state, state=tk.DISABLED)
    incShinyNorm.grid(column=1, row=14, sticky=tk.W)

    normMethod_lbl2 = tk.Label(self, text="When not using interactive normalization:")
    normMethod_lbl2.grid(column=2, row=12, columnspan=2 ,pady=(40,0), sticky=tk.E)
    normMethod_lbl = tk.Label(self, text="Normalization Method:")
    normMethod_lbl.grid(column=2, row=14, sticky=tk.E)
    normMethod = ttk.Combobox(self)
    normMethod['values']=('loess', 'quantile', 'median', 'z transform', 'none')
    normMethod.current(0)
    normMethod.grid(column=3, row=14)
    
    logtransform_state = tk.BooleanVar()
    logtransform_state.set(True)
    logtransform = ttk.Checkbutton(self, text="Log Transform", var=logtransform_state)
    logtransform.grid(column=4, row=14)

    species_lbl = tk.Label(self, text="Species:")
    species_lbl.grid(column=0, row=16, pady=(30,0), sticky=tk.E)
    species = ttk.Combobox(self)
    species['values']=('Human (9606)', 'Mouse (10090)','Other', 'Yeast (559292)', 'E. coli (511145)','Zebrafish (7955)', 'C. elegans (6239)', 'Fruit Fly (7227)','Rat (10116)')
    species.current(0)
    species.grid(column=1, row=16, pady=(30,0))

    heatcolors_lbl = tk.Label(self, text="Heatmap color scale:")
    heatcolors_lbl.grid(column=0, row=17, pady=(0,10))
    heatcolors = ttk.Combobox(self)
    heatcolors['values']=('viridis', 'RdYlBu', 'RdBu')
    heatcolors.current(0)
    heatcolors.grid(column=1, row=17, pady=(0,10))

    # NUMERIC ENTRIES
    vcmd = (self.register(self.onValidate), '%d','%i','%P','%s','%S','%v','%V','%W')

    #emapVar_lbl = tk.Label(self, text="EnrichmentMap Parameters:")
    #emapVar_lbl.grid(column=3, columnspan=1, sticky=tk.W, row=18, pady=(20,0))
    
    #emapVar_p_lbl = tk.Label(self, text="P Value:")
    #emapVar_p_lbl.grid(column=2, row=19, sticky=tk.E)
    #emapVar_p_start = tk.StringVar(self)
    #emapVar_p_start.set("0.05")
    #enrichmentmap_p_val = tk.Spinbox(self, from_=0.01, to=0.50, increment=0.01, textvariable=emapVar_p_start, validate='all', validatecommand=vcmd)
    #enrichmentmap_p_val.grid(column=3, row=19)

    #emapVar_q_lbl = tk.Label(self, text="Q Value:")
    #emapVar_q_lbl.grid(column=2, row=20, sticky=tk.E)
    #emapVar_q_start = tk.StringVar(self)
    #emapVar_q_start.set("0.30")
    #enrichmentmap_q_val = tk.Spinbox(self, from_=0.01, to=0.50, increment=0.01, textvariable=emapVar_q_start, validate='all', validatecommand=vcmd)
    #enrichmentmap_q_val.grid(column=3, row=20)

    cutoff_lbl = tk.Label(self, text="Cut Offs:")
    cutoff_lbl.grid(column=0, columnspan=1, row=18, pady=(20,0))
    
    cutoff_zero_lbl = tk.Label(self, text="Zero Value Filter:")
    cutoff_zero_lbl.grid(column=2, row=15, sticky=tk.E)
    cutoff_zero_start = tk.StringVar(self)
    cutoff_zero_start.set("0.70")
    cutoff_zero_val = tk.Spinbox(self, from_=0.00, to=1.00, increment=0.01, textvariable=cutoff_zero_start, validate='all', validatecommand=vcmd)
    cutoff_zero_val.grid(column=3, row=15)

    cutoff_fdr_lbl = tk.Label(self, text="FDR:")
    cutoff_fdr_lbl.grid(column=0, row=19, sticky=tk.E)
    cutoff_fdr_start = tk.StringVar(self)
    cutoff_fdr_start.set("0.05")
    cutoff_fdr_val = tk.Spinbox(self, from_=0.01, to=0.50, increment=0.01, textvariable=cutoff_fdr_start, validate='all', validatecommand=vcmd)
    cutoff_fdr_val.grid(column=1, row=19)
    
    cutoff_percent_lbl = tk.Label(self, text="Top %:")
    cutoff_percent_lbl.grid(column=0, row=20, sticky=tk.E)
    cutoff_percent_start = tk.StringVar(self)
    cutoff_percent_start.set("0.010")
    cutoff_percent_val = tk.Spinbox(self, from_=0.001, to=0.50, increment=0.001, textvariable=cutoff_percent_start, validate='all', validatecommand=vcmd)
    cutoff_percent_val.grid(column=1, row=20)
    
    txtFolder_lbl = tk.Label(self, text="txt Folder report:")
    txtFolder_lbl.grid(column=0, row=24, pady=(20,0), sticky=tk.E)
    txtFolder_lbl2 = tk.Label(self, text="Copy mqpar file into txt folder for complete results.")
    txtFolder_lbl2.grid(column=2, row=24, columnspan=2, pady=(20,0), sticky=tk.W)
    txtFolder = ttk.Combobox(self)
    txtFolder['values']=('TRUE', 'FALSE')
    txtFolder.current(1)
    txtFolder.grid(column=1, row=24, pady=(20,0))
    
    # Check if inherit paths directories are there
    
    if os.path.isdir(os.path.normpath(libraries_path) ):
      libraries_path = os.path.normpath(libraries_path)
    else :
      libraries_path = os.path.normpath("./")
      inherit_paths = "FALSE"

    ############################################
    # SAVE VARIABLES FILE AND CLOSE COMMAND
    def clicked_button():
      analysis_dir = str(tk.StringVar(self, value=fileDir).get()) 
      
      variable_file = analysis_dir + "/Parameters.R"
      outputfile = open(variable_file, "w")
      
      outputfile.write("project_name <- '" + nameProj_var.get() +"';\n")
      outputfile.write("working_dir <- '" + analysis_dir + "'; # Directory for analysis and where data is \n")
      outputfile.write("annotation_filename <- '" + str(tk.StringVar(self, value=nameAnn).get()) + "';\n")

      outputfile.write("query_web <- " + str(queryWeb_state.get()).upper() + ";\n")
      #outputfile.write("show_names <- " + str(showSampleName_state.get()).upper() + ";\n")
      outputfile.write("map_color <- '" + str(heatcolors.get()) + "';\n")
      
      outputfile.write("zero_percent <- " + str(cutoff_zero_val.get()) + ";\n")
      outputfile.write("norm_method <- '" + str(normMethod.get()) + "';\n")
      outputfile.write("log_transform <- "+ str(logtransform_state.get()).upper() + ";\n")
      outputfile.write("shinyNorm <- " + str(incShinyNorm_state.get()).upper() +  ";\n")
      
      outputfile.write("saveXlsx <- " + str(runXlsx_state.get()).upper() + ";\n")
      outputfile.write("runDifferential <- " + str(incDifferential_state.get()).upper() + ";\n")
      outputfile.write("adjpcutoff <- " + str(cutoff_fdr_val.get()) + ";\n")
      outputfile.write("sig_percent <- " + str(cutoff_percent_val.get()) + ";\n")
      
      outputfile.write("enrichr_section <- " + str(runEnrichr_state.get()).upper() + ";\n")
      outputfile.write("gsea_section <- " + str(runGSEA_state.get()).upper() + ";\n")
      outputfile.write("gsea_MOMENTA_section <- " + str(runMOMENTA_state.get()).upper() + ";\n")
      outputfile.write("activePathways_section <- " + str(runActPath_state.get()).upper() + ";\n")
      outputfile.write("all_comparisons <- "+ str(allcomparisons_state.get()).upper() + ";\n")
      
      outputfile.write("species <- '" + str(species.get()) + "';\n")
      #outputfile.write("enrichmentmap_p_val <- " + str(enrichmentmap_p_val.get()) + ";\n")
      #outputfile.write("enrichmentmap_q_val <- " + str(enrichmentmap_q_val.get()) + ";\n")
      outputfile.write("use_site_norm <- " + str(useSiteNorm_state.get()).upper() + ";\n")
      
      outputfile.write("int_heatmap_section <- " + str(isIntHM_state.get()).upper() + ";\n")
      outputfile.write("int_volcano_section <- " + str(isIntVl_state.get()).upper() + ";\n")
      outputfile.write("txtFolder <- " + str(txtFolder.get()).upper() + ";\n")
      
      outputfile.write("inherit_paths <- " + inherit_paths +  ";\n")
      outputfile.write("libraries_path <- '" + libraries_path +  "';\n")
      #outputfile.write("pandoc_path <- '" + pandoc_path +  "';\n")

      outputfile.close()
      lbl.configure(text="Parameters file created.")
      sys.exit(analysis_dir)
      
    ############################################
    
    # ENTER BUTTON
    btn= tk.Button(self, text="Enter", command=clicked_button)
    btn.grid(column=0, row=30, columnspan=3, pady=(20,0))

  ###############################################
  
  # Validate numbers function
  def onValidate(self, d, i, P, s, S, v, V, W):
    minval=int(self.nametowidget(W).config('from')[4])
    maxval=int(self.nametowidget(W).config('to')[4])
    if P.isdigit() and P in range(minval, maxval):
      return True
    else:
      return False

###############################################

# Run App

root=tk.Tk()
root.title('Omics Notebook Setup')
GUI(root).pack(fill="both", expand=True)
root.mainloop()

###############################################

