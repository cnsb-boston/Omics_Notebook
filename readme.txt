Install cytoscape and java 8. 

Install enrichmentmap pipeline collection.

Install R 3.5 and Python 3.6.

Install Rstudio.

Install all required R packages.

Open cytoscape. Default port for REST API is 1234. If your computer already has this port in use, open cytoscape from the command line and specity another port. 

Run Notebook.py

Data files should be either .txt or .csv, with data column names specified in the annotation file. All other columns will be used as annotation. Columns with HGNC gene symbols should be called "Gene", columns with uniprot ID's should be called "Protein." The script will parse the MaxQuant output to create Gene and Protein columns if the fasta file was correctly configured.
