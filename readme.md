## Software requirements

* Install R 3.5 and Python 3.6.

* Install Rstudio.

* Install all required R packages. (See install.R files)


#### For Enrichment Map

* Install cytoscape 3.5.1 and java 8. 

* Install enrichmentmap pipeline collection.

Open cytoscape. Default port for REST API is 1234. If your computer already has this port in use, open cytoscape from the command line and specity another port. 

Note: it may be easier to run the scripts and import enrichment results into Cytoscape manually.

---

## Structure and Instructions

"Notebook.py" is a python script that will call two scripts. The first is "src/Pipeline.py", which will generate a Parameters.R file in the directory with the data. This can be helpful to assure correct variable formatting. The second is "src/Pipeline.R", which will knit the R markdown files.

It may be easier to use "Pipeline.py" to generate a "Parameters.R". Place the "Parameters.R" with the directory with the data and annotation files, and then use R studio to walk through the .Rmd files. This can be done with the umbrella "Notebook.Rmd", exploratory "Notebook_1_Expl.Rmd", differential "Notebook_2_Diff.Rmd", enrichment "Notebook_3_Enrch.Rmd", integrative "Notebook_4_Integ.Rmd", and "notebookFunctions.R".

Data files should be either .txt or .csv, with data column names specified in the annotation file. All other columns will be used as annotation. Columns with HGNC gene symbols should be called "Gene", columns with uniprot ID's should be called "Protein." The script will parse the MaxQuant output to create Gene and Protein columns if the fasta file was correctly configured.

---
