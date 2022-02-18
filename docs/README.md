## Instructions

For additional information, see:
https://github.com/cnsb-boston/Omics_Notebook_Docs

#### 1. Install git. Clone github repository:
`git clone https://github.com/cnsb-boston/Omics_Notebook.git`

Edit Config.R configuration variables to match your installation paths.

#### 2 A. Run Natively.

##### i. Install specified versions of R and the OmicsNotebook package and its dependencies
It's recommended to use something like [remotes](https://cran.r-project.org/web/packages/remotes/) to install this package.
If this code was downloaded into the "Omics_Notebook" directory:

`R -e "remotes::install_local('Omics_Notebook')"`

##### ii. Run Notebook.R script, which will automate entire pipeline.
`Rscript Notebook.R`


#### 2 B. Run in a container.
Given complexity with R package dependencies, it may be easiest to run with the assistance of a container on local computers or cloud resources. 

##### i. Install R and the container system (Docker or Singularity). 

##### ii. From Omics Notebook Directory, run:
*Docker:*
`docker pull cnsbboston/omicsnotebook:latest`

Or build from Dockerfile.

`docker build -t omicsnotebook-base -f Dockerfile-base .` base OS with required system dependencies (libraries, R, etc)

`docker build -t cnsbboston/omicsnotebook -f Dockerfile .` full image with R packages, built on top of the base image

*Singularity:*
`singularity build "ON.simg" "docker://cnsbboston/omicsnotebook:latest"`

The base image can be used if it's preferable to keep R packages outside the container. In this case the R library directory can be mapped onto the container with -v (Docker) or -B (Singularity), or using the libdir variable in Config.R.

##### iii. Run Notebook.R with a container argument.
*Docker:*
`Rscript Notebook.R Docker`

`-c` can also be passed here (like Singularity below) to specify a different docker image tag to run

*Singularity:*
`Rscript Notebook.R Singularity -c '/Path/to/SingularityImage.simg'`

Or run each component on its own, which may be easier for integrating into other workflows.

##### iii. Generate Parameters.R file.
The GUI component automates the creation of the Parameters.R file and should be run natively with R. While there are solutions for GUI in docker (e.g. VNC or configuring X11 socket), they may be difficult to configure on all systems.
`Rscript Notebook.R GUI` Adjust path for file location.

##### iv. Run R analysis using docker:
```
docker run -it --rm \
  -u docker \
  -v /PATH/TO/OMICS NOTEBOOK:/home:rw \
  -v /PATH/TO/DATA ANALYSIS DIR:/data:rw \
  cnsbboston/omicsnotebook Rscript /home/src/Pipeline.R "/home" "/data"
```
The Analysis Directory is the directory where the Parameters.R, Annotation, and Data files are and where the output will be saved. Remember, files are relative to the docker container. The third and fourth -v lines mount the local Omics_Notebook directory to the home directory in the container, and mount the Analysis Directory to the data directory.

The modified docker run command may look like:
```
docker run -it --rm \
  -u docker \
  -v ~/Omics_Notebook:/home:rw \
  -v ~/Omics_Notebook/example:/data:rw \
  cnsbboston/omicsnotebook Rscript "/home/src/Pipeline.R" "/home" "/data"
```

##### iv. Run R analysis using Singularity:
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

This list can be extended for custom analysis. See /docs/SupplementaryInformation.pdf for screen shots and additional instructions.

---

#### Input Data

Data files should be either .txt or .csv, with data column names, corresponding to samples, specified in the annotation file. All other columns in the data files will be used as annotation. Standardized annotation currently configured for analysis is as follows:

* Columns with HGNC gene symbols should be called "Gene". 
* Columns with uniprot ID's should be called "Protein." 

The script will parse the MaxQuant output to create Gene and Protein columns if the fasta file was correctly configured. Example data, annotation file, and output are provided in "example/".

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

"Notebook.R" is an R script that will call two scripts. The first is "NotebookGUI.R", which will generate a Parameters.R file in the directory with the data (Analysis Directory). This can be helpful to assure correct variable formatting. The second is "src/Pipeline.R", which will knit the R markdown files.

For custom analysis, it may be easier to call `Rscript Notebook.R GUI` to generate a "Parameters.R". Place the "Parameters.R" in the directory with the data and annotation files, and then use R studio to walk through the .Rmd files. This can be done with the umbrella "Notebook.Rmd", exploratory "Notebook_1_Expl.Rmd", differential "Notebook_2_Diff.Rmd", enrichment "Notebook_3_Enrch.Rmd", integrative "Notebook_4_Integ.Rmd", and functions in the "R" directory.

In this way, the r markdown reports can be generated automatically with wrapper scripts or manually with R studio. 

---

## GUI Instructions

Read /docs/SupplementaryInformation.pdf for screen shots and instructions.

---

## Software requirements for running natively.

In general, this software is designed to make use of containerization for managing the many R package installation requirements and dependencies. It also permits use of this software on local computers or shared (e.g., HPC or cloud) computing resources with Docker or Singularity. However, the following will be updated as needed.

* R 4.1

* Install Rstudio.

* Install all required R packages. (See install.R and respective R packages for compatibility requirements).

This software has been tested for use on Linux (Cent OS 6 and 7, Ubuntu 18.04, 20.04) and MacOS (10.14 and 10.15)

---

