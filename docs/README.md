## Instructions

#### 1. Install git. Clone github repository:
`git clone https://github.com/cnsb-boston/Omics_Notebook.git`


#### 2 A. Run Natively.

##### i. Install specified versions of R and Python and all packages in install.R file (see below).

##### ii. Run Notebook.py script, which will automate entire pipeline.
`python3 Notebook.py`


#### 2 B. Run with Docker.
Given complexity with R package dependencies, it may be easiest to run with the assistance of docker on local computers or cloud resources. 

##### i. Install Python (3.6), tkinter, and Docker. 

##### ii. From Omics Notebook Directory, run:
`docker pull bblum/omics_notebook:latest`
Or build from Dockerfile.

##### iii. Run Notebook.py with "Docker" argument.
`python3 Notebook.py Docker`


Or run each component on its own, which may be easier for integrating into other workflows.

##### i. Install Python (3.6), tkinter, and Docker. 

##### ii. From Omics Notebook Directory, run:
`docker pull bblum/omics_notebook:latest`
Or build from Dockerfile.

##### iii. Generate Parameters.R file.
The GUI component automates the creation of the Parameters.R file and should be run natively with Python3 and tkinter. While there are solutions for GUI in docker (e.g. VNC or configuring X11 socket), they may be difficult to configure on all systems.
`python3 Omics_Notebook/src/Pipeline.py` Adjust path for file location.

##### iv. Run R analysis using docker:
```
docker run -it --rm \
  -u docker \
  -v /PATH/TO/OMICS NOTEBOOK:/home:rw \
  -v /PATH/TO/DATA ANALYSIS DIR:/data:rw \
  bblum/omics_notebook Rscript /home/src/Pipeline.R "/home" "/data"
```
The Analysis Directory is the directory where the Parameters.R, Annotation, and Data files are and where the output will be saved. Remember, files are relative to the docker container. The third and fourth -v lines mount the local Omics_Notebook directory to the home directory in the container, and mount the Analysis Directory to the data directory.

The modified docker run command may look like:
```
docker run -it --rm \
  -u docker \
  -v ~/Omics_Notebook:/home:rw \
  -v ~/Omics_Notebook/example:/data:rw \
  bblum/omics_notebook Rscript "/home/src/Pipeline.R" "/home" "/data"
```


#### 2 C. Run with Singularity.
For use with shared HPC linux environments.

##### i. Install Python (3.6), tkinter, and Singularity. 

##### ii. Create singularity image from Docker.
`singularity build "ON.simg" "docker://bblum/omics_notebook:latest"`

##### iii. Run Notebook.py with "Docker" argument.
`python3 Notebook.py Singularity '/Path/to/SingularityImage.simg'`


Or run each component on its own, which may be easier for integrating into other workflows.

##### i. Install Python (3.6), tkinter, and Singularity. 

##### ii. Create singularity image from Docker.
`singularity build "ON.simg" "docker://bblum/omics_notebook:latest"`

##### iii. Generate Parameters.R file.
The GUI component automates the creation of the Parameters.R file and should be run natively with Python3 and tkinter. While there are solutions for GUI in docker (e.g. VNC or configuring X11 socket), they may be difficult to configure on all systems.
`python3 Omics_Notebook/src/Pipeline.py` Adjust path for file location.

##### iv. Run R analysis using docker:
```
singularity run \
  --bind /PATH/TO/OMICS NOTEBOOK:/home:rw \
  --bind /PATH/TO/DATA ANALYSIS DIR:/data:rw \
  /PATH/TO/ON.simg Rscript /home/src/Pipeline.R "/home" "/data"
```


---

#### Annotation file

The annotation file provides a standard way to read in the experimental data. See supplementaryinformation.pdf for extended explanations with screenshots.

The top section has information about what differential analyses to perform, what Group each sample belongs to, and the sample names for each sample. Data set names, "type" and filename should be specified in the bottom left, with no empty rows. Each row in the lower section is a different omics dataset. For a given row, the columns next to SampleName should correspond to the specific samples in the data. Other columns will be kept as annotation information.

In the second sheet of the annotation file, the sample names are repeated. Additional rows can specify additional annotation information. A couple will control additional functionality. 

* Batch : if "Batch" is found, the row will be used to try to perform batch correction. The default is to use ComBat. Care should be take not to overfit data.
* ColorsHex : if "ColorsHex" is used, and the row has functional hex colors values corresponding to "Group", these hex colors will be used in figures where possible.
* Group2: if "Group2" is used, this will specify shapes in the PCA plot. More functionality may be added later.
* TimeSeries: if a numeric value, will run limma as for a time series analysis. Will only work with 2 Groups. Code may be edited for other analyses. If there is 1 group, enter the time point as "Group."
* Pairs: if "Pairs" is found, will try to run differential analysis accounting for paired samples.

This list can be extended for custom analysis.

---

#### Input Data

Data files should be either .txt or .csv, with data column names specified in the annotation file. All other columns will be used as annotation. Columns with HGNC gene symbols should be called "Gene", columns with uniprot ID's should be called "Protein." The script will parse the MaxQuant output to create Gene and Protein columns if the fasta file was correctly configured. Example data, annotation file, and output are provided in "example/".

The makeEset function should be modified to parse additional input data formats. 

---

#### Custom Input Parameters

Omics Notebook will automatically search for the following in the Analysis Directory to perform additional analysis or offer customization.

* "Paramaters.R" file: permits the over-writing of variables specified in Parameters.R on in Notebook.RMD. Note proper formatting is required but provides many options. For example, if different normalization methods are desired for different datasets, the variable "norm_method" may be defined as a vector with the appropriate definitions. See the Options.R chunk in Notebook.RMD for additional variables that may be modified.
* "Subset_Lists" directory: If a directory called "Subset_Lists" exists, with .txt files specifying lists of feature identifiers or Gene symbols provided on separate lines, they will be incorporated in several of the outputs, including heatmaps of only those features, or labels on volcano plots.
* "Gene_Sets" directory: If a directory called "Gene_Sets" with gene sets formatted as .gmt files inside, additional GSEA analysis will be performed with the provided .gmt files and output in the "GSEA_Custom" directory.

---

#### Output

The output of the Omics Notebook R markdown is an html file and analysis directory with results saved as text or image files. The output adapts depending on the input data but includes exploratory, differential, and enrichment analysis.

A formatted excel file facilitates the sharing of normalized feature-level data with annotation from the analysis.

The complete pipeline should run on a normal desktop computer in a few hours, possibly more or much less depending on options and size of the input data set.

---

## Structure

"Notebook.py" is a python script that will call two scripts. The first is "src/Pipeline.py", which will generate a Parameters.R file in the directory with the data (Analysis Directory). This can be helpful to assure correct variable formatting. The second is "src/Pipeline.R", which will knit the R markdown files.

For cutsome analysis, it may be easier to use "Pipeline.py" to generate a "Parameters.R". Place the "Parameters.R" in the directory with the data and annotation files, and then use R studio to walk through the .Rmd files. This can be done with the umbrella "Notebook.Rmd", exploratory "Notebook_1_Expl.Rmd", differential "Notebook_2_Diff.Rmd", enrichment "Notebook_3_Enrch.Rmd", integrative "Notebook_4_Integ.Rmd", and functions in the "R" directory.

In this way, the r markdown reports can be generated automatically with wrapper scripts or manually with R studio. 

---

## GUI Instructions

Read /docs/SupplementaryInformation.pdf for screen shots and instructions.

---

## Software requirements for running natively.

In general, this software is designed to make use of containerization for managing the many R package installation requirements and dependencies. It also permits use of this software on local computers or shared (e.g., HPC or cloud) computing resources with Docker or Singularity. However, the following will be updated as needed.

* R 3.6,  Python 3.6, Tkinter.

* Install Rstudio.

* Install all required R packages. (See install.R and respective R packages for compatibility requirements).

This software has been tested for use on Linux (Cent OS 6 and 7, Ubuntu 18.04, 20.04) and MacOS (10.14 and 10.15)

---

