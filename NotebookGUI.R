library(tcltk)

tclbool = function(tclx) as.logical(as.integer(tclvalue(tclx)))
tclnum = function(tclx) as.numeric(tclvalue(tclx))
tclstring = function(tclx) paste0("'",tclvalue(tclx),"'")
tcltype=function(default,tfun){ x=tclVar(default); attr(x,"type")=tfun; x }
tclBool=function(default=T) tcltype(default,tclbool)
tclStr=function(default="") tcltype(default,tclstring)
tclNum=function(default=0) tcltype(default,tclnum)

make.gui = function(startdir="/projectnb/cnsb2022", extra_v=list(), extra_widgets=NULL, param_file=NULL){
  write.params = function(vars){
    fname=paste0(tclvalue(vars$working_dir), "/Parameters.R")
    write(paste("#OmicsNotebook"), fname, append=F)
    for(i in 1:length(vars)){
      f=attr(vars[[i]],"type")
      write(paste(names(vars)[i],"=",f(vars[[i]])), fname, append=T)
    }
  }

  tt <- tktoplevel()

  v=list(project_name=tclStr("Analysis"),
         working_dir=tclStr(""),
         annotation_filename=tclStr(""),
         query_web=tclBool(TRUE),
         map_color=tclStr('viridis'),
         zero_percent=tclNum(0.70),
         norm_method = tclStr('loess'),
         log_transform = tclBool(TRUE),
         shinyNorm = tclBool(FALSE),
         saveXlsx = tclBool(TRUE),
         runDifferential = tclBool(TRUE),
         adjpcutoff = tclNum(0.05),
         sig_percent = tclNum(0.010),
         enrichr_section = tclBool(FALSE),
         gsea_section = tclBool(TRUE),
         gsea_MOMENTA_section = tclBool(TRUE),
         activePathways_section = tclBool(TRUE),
         all_comparisons = tclBool(TRUE),
         combined_met = tclBool(TRUE),
         species = tclStr('Human (9606)'),
         use_site_norm = tclBool(FALSE),
         int_heatmap_section = tclBool(FALSE),
         int_volcano_section = tclBool(TRUE),
         txtFolder = tclBool(FALSE),
         inherit_paths = tclBool(FALSE),
         libraries_path = tclStr('.'))
  v=c(v,extra_v)

  # set default values from previous Parameters.R
  if(!is.null(param_file) && file.exists(param_file)){
    p.env = new.env()
    sys.source(param_file, envir=p.env, toplevel.env=p.env)
    for(p in names(p.env)){
      if(!is.null(v[[p]]))
        tclvalue(v[[p]]) = p.env[[p]]
    }
  }

  tkgrid(tklabel(tt, text="Project Name:"), column=0, row=1, sticky="E") 
  tkgrid(tkentry(tt, textvariable=v$project_name), column=1, row=1, sticky="W")

  annotation_validate=tclStr("")
  tkgrid(tklabel(tt, textvariable=annotation_validate), column=1, row=2, columnspan=3, sticky="E")

  fileDir_clicked=function()tclvalue(v$working_dir)=tkchooseDirectory(initialdir=startdir, title="Choose directory")
  tkgrid(tklabel(tt, text="Files Directory:"), column=0, row=3, sticky="E")
  tkgrid(tkentry(tt, textvariable=v$working_dir, width=60, state="disabled"), column=1, row=3, columnspan=3, sticky="W")
  tkgrid(tkbutton(tt, text="Choose Directory", command=fileDir_clicked), column=4, row=3, sticky="W")

  anno_clicked=function(){
    fname=tkgetOpenFile(initialdir=tclvalue(v$working_dir), title="Choose Annotation file", filetypes="{ {Excel Files} {.xlsx} } { {All Files} * }")
    fi=file.path(tclvalue(v$working_dir),basename(tclvalue(fname)))
    if(file.exists(fi)){
      tclvalue(v$annotation_filename)=fi
      tclvalue(annotation_validate)=""
    } else {
      tclvalue(annotation_validate)="Annotation file must be in the project directory."
    }
  }
  tkgrid(tklabel(tt, text="Annotation File:"), column=0, row=5, sticky="E")
  tkgrid(tkentry(tt, textvariable=v$annotation_filename, width=60, state="disabled"), column=1, row=5, columnspan=3, sticky="W")
  tkgrid(tkbutton(tt, text="Choose Annotation File", command=anno_clicked), column=4, row=5, sticky="W")

  gencheck=function(col,row,text,var, colspan=1){
    tkgrid(tkcheckbutton(tt, text=text, variable=v[[var]]), column=col, row=row, columnspan=colspan, sticky="W")
  }

  gencheck(1, 6, "Differential Analysis", "runDifferential")
  gencheck(2, 6, "Query Uniprot", "query_web")
  gencheck(3, 6, "Excel Report", "saveXlsx")
  gencheck(4, 6, "Interactive Heatmap", "int_heatmap_section")

  gencheck(1, 7, "All Comparisons", "all_comparisons", colspan=1)
  gencheck(4, 7, "Interactive Volcano", "int_volcano_section")

  gencheck(1, 10, "Run EnrichR", "enrichr_section")
  gencheck(2, 10, "Run GSEA", "gsea_section")
  gencheck(3, 10, "Run MOMENTA", "gsea_MOMENTA_section")
  gencheck(4, 10, "Run ActivePathways", "activePathways_section")

  gencheck(1, 11, "Combined Metabolomics", "combined_met", colspan=1)
  gencheck(3, 11, "Use Normalized to 1st Data", "use_site_norm", colspan=2)

  #gencheck(1, 14, "Interactive Normalization", "shinyNorm")
  gencheck(4, 14, "Log Transform", "log_transform")

  norm_methods = c('loess', 'quantile', 'median', 'z transform', 'none')
  #tkgrid(tklabel(tt, text="When not using interactive normalization:"), column=2, row=12)
  tkgrid(tklabel(tt, text="Normalization Method"), column=2, row=14)
  tkgrid(ttkcombobox(tt, textvariable=v$norm_method, values=norm_methods), column=3, row=14)

  tkgrid(tklabel(tt, text="Zero Value Filter:"), column=2, row=15)
  tkgrid(ttkspinbox(tt, textvariable=v$zero_percent, from=0, to=1, increment=.01), column=3, row=15)

  species = c('Human (9606)', 'Mouse (10090)','Other', 'Yeast (559292)', 'E. coli (511145)','Zebrafish (7955)', 'C. elegans (6239)', 'Fruit Fly (7227)','Rat (10116)')
  tkgrid(tklabel(tt, text="Species"), column=0, row=16)
  tkgrid(ttkcombobox(tt, textvariable=v$species, values=species), column=1, row=16)

  map_color = c('viridis', 'RdYlBu', 'RdBu')
  tkgrid(tklabel(tt, text="Heatmap color scale:"), column=0, row=17)
  tkgrid(ttkcombobox(tt, textvariable=v$map_color, values=map_color), column=1, row=17)

  tkgrid(tklabel(tt, text="Cut Offs:"), column=0, row=18, columnspan=1)

  tkgrid(tklabel(tt, text="FDR:"), column=0, row=19)
  tkgrid(ttkspinbox(tt, textvariable=v$adjpcutoff, from=0.01, to=0.5, increment=.01), column=1, row=19)
  tkgrid(tklabel(tt, text="Top %:"), column=0, row=20)
  tkgrid(ttkspinbox(tt, textvariable=v$sig_percent, from=0.001, to=0.5, increment=.001), column=1, row=20)

  gencheck(0, 25, "txt Folder report", "txtFolder")
  tkgrid(tklabel(tt, text="Copy mqpar file into txt folder for complete results."), column=0, row=24, columnspan=2)

  tkgrid(tklabel(tt, text="Cut Offs:"), column=0, row=18, columnspan=1)

  if(is.function(extra_widgets)) extra_widgets(tt,v)

  tkgrid(tkbutton(tt, text="Enter", command=function(){write.params(v); tkdestroy(tt)}), column=1, row=30)

  tkwait.window(tt)
  Sys.sleep(3) # sometimes the window hangs if we don't wait?
  tclvalue(v$working_dir)
}
