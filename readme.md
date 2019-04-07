## Software requirements

* R 3.5,  Python 3.6, Java (for GSEA .jar and Cytoscape)

* Install Rstudio.

* Install all required R packages. (See install.R)

## For Enrichment Analysis

#### GSEA
* Download the GSEA .jar file from the Broad Institute. Update directory to appropriate script headers. 


#### Enrichment Map

* Install cytoscape>3.5.1.

* Install enrichmentmap pipeline collection.

Open cytoscape. Default port for REST API is 1234. If your computer already has this port in use, open cytoscape from the command line and specity another port. 

Note: it may be easier to run the scripts and import enrichment results into Cytoscape manually.

---

## Structure and Instructions

"Notebook.py" is a python script that will call two scripts. The first is "src/Pipeline.py", which will generate a Parameters.R file in the directory with the data. This can be helpful to assure correct variable formatting. The second is "src/Pipeline.R", which will knit the R markdown files.

It may be easier to use "Pipeline.py" to generate a "Parameters.R". Place the "Parameters.R" with the directory with the data and annotation files, and then use R studio to walk through the .Rmd files. This can be done with the umbrella "Notebook.Rmd", exploratory "Notebook_1_Expl.Rmd", differential "Notebook_2_Diff.Rmd", enrichment "Notebook_3_Enrch.Rmd", integrative "Notebook_4_Integ.Rmd", and functions in the "R" directory.

Data files should be either .txt or .csv, with data column names specified in the annotation file. All other columns will be used as annotation. Columns with HGNC gene symbols should be called "Gene", columns with uniprot ID's should be called "Protein." The script will parse the MaxQuant output to create Gene and Protein columns if the fasta file was correctly configured.

The makeEset function should be modified to parse additional input data formats. 

---

## Annotation file

The annotation file provides a standard way to read in the experimental data. The top section has information about what differential analyses to perform, what Group each sample belongs to, and the sample names for each sample. Datasets names, "type" and filename should be specified in the bottom left, with no empty rows. Each row in the lower section is a different omics dataset. For a given row, the columns under the samplenames should correspond to the samples. Other columns will be kepts as annotation information.

In the second sheet of the annotation file, the sample names are repeated. Additional rows can specify additional annotation information. A couple will control additional functionality. 

* Batch : if "Batch" is found, the row will be used to try to perform batch correction. The default is to use ComBat. Care should be take not to overfit data.
* ColorsHex : if "ColorsHex" is used, and the row has functional hex colors values corresponding to "Group", these hex colors will be used in figures where possible.
* Group2: if "Group2" is used, this will specify shapes in the PCA plot. More functionality may be added later.
* TimeSeries: if a numeric value, will run limma as for a time series analysis. Will only work with 1 or 2 Groups. Code may be edited for other analysis.



---
