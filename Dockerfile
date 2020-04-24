# hash:sha256:8d6d14519232bd92b52eb66fd7715e34728802f891f279fc2492dc30d9fd8c68
FROM r-base:3.6.1

ARG DEBIAN_FRONTEND=interactive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libnetcdf-dev=1:4.6.0-2build1 \
        libxt-dev=1:1.1.5-1 \
        pandoc=1.19.2.4~dfsg-1build4 \
        python3=3.6.7-1~18.04 \
        python3-tk=3.6.9-1~18.04 \
        r-base-dev=3.6.2-1bionic \
    && rm -rf /var/lib/apt/lists/*

RUN export R_REMOTES_NO_ERRORS_FROM_WARNINGS=true \
    && Rscript -e 'remotes::install_version("BiocManager", "1.30.10")' \
    && Rscript -e 'remotes::install_version("Cairo", "1.5-10")' \
    && Rscript -e 'remotes::install_version("PTXQC", "1.0.0")' \
    && Rscript -e 'remotes::install_version("R.utils", "2.9.2")' \
    && Rscript -e 'remotes::install_version("RColorBrewer", "1.1-2")' \
    && Rscript -e 'remotes::install_version("RCurl", "1.95-4.13")' \
    && Rscript -e 'remotes::install_version("RJSONIO", "1.3-1.4")' \
    && Rscript -e 'remotes::install_version("Rcpp", "1.0.3")' \
    && Rscript -e 'remotes::install_version("UpSetR", "1.4.0")' \
    && Rscript -e 'remotes::install_version("VennDiagram", "1.6.20")' \
    && Rscript -e 'remotes::install_version("calibrate", "1.7.5")' \
    && Rscript -e 'remotes::install_version("corrplot", "0.84")' \
    && Rscript -e 'remotes::install_version("devtools", "2.2.1")' \
    && Rscript -e 'remotes::install_version("dplyr", "0.8.3")' \
    && Rscript -e 'remotes::install_version("enrichR", "2.1")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.2.1")' \
    && Rscript -e 'remotes::install_version("ggrepel", "0.8.1")' \
    && Rscript -e 'remotes::install_version("gridExtra", "2.3")' \
    && Rscript -e 'remotes::install_version("haven", "2.2.0")' \
    && Rscript -e 'remotes::install_version("heatmaply", "1.0.0")' \
    && Rscript -e 'remotes::install_version("httr", "1.4.1")' \
    && Rscript -e 'remotes::install_version("knitr", "1.27")' \
    && Rscript -e 'remotes::install_version("lattice", "0.20-38")' \
    && Rscript -e 'remotes::install_version("ncdf4", "1.17")' \
    && Rscript -e 'remotes::install_version("openxlsx", "4.1.4")' \
    && Rscript -e 'remotes::install_version("pls", "2.7-2")' \
    && Rscript -e 'remotes::install_version("reshape", "0.8.8")' \
    && Rscript -e 'remotes::install_version("reshape2", "1.4.3")' \
    && Rscript -e 'remotes::install_version("rlang", "0.4.2")' \
    && Rscript -e 'remotes::install_version("rmarkdown", "2.0")' \
    && Rscript -e 'remotes::install_version("scatterplot3d", "0.3-41")' \
    && Rscript -e 'remotes::install_version("som", "0.3-5.1")' \
    && Rscript -e 'remotes::install_version("spls", "2.2-3")' \
    && Rscript -e 'remotes::install_version("stringr", "1.4.0")' \
    && Rscript -e 'remotes::install_version("tidyr", "1.0.0")' \
    && Rscript -e 'remotes::install_version("viridis", "0.5.1")'

RUN Rscript -e 'remotes::install_github( \
        "cnsb-boston/KSEAapp", \
        ref = "805fb0c78bbc3e01636142d6738fe29f5930064a")'

RUN Rscript -e 'install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install(c( \
        "Biobase", \
        "BiocParallel", \
        "CAMERA", \
        "ComplexHeatmap", \
        "Glimma", \
        "GlobalAncova", \
        "Hmisc", \
        "KEGGgraph", \
        "MSnbase", \
        "ROCR", \
        "Rgraphviz", \
        "Rserve", \
        "SSPA", \
        "caTools", \
        "car", \
        "caret", \
        "data.table", \
        "e1071", \
        "edgeR", \
        "ellipse", \
        "fgsea", \
        "fitdistrplus", \
        "genefilter", \
        "globaltest", \
        "gplots", \
        "igraph", \
        "impute", \
        "lars", \
        "limma", \
        "magrittr", \
        "metap", \
        "mzR", \
        "pROC", \
        "pcaMethods", \
        "pheatmap", \
        "plotly", \
        "preprocessCore", \
        "randomForest", \
        "siggenes", \
        "sva", \
        "xcms" \
    ))' # Original versions: 2.46.0 1.20.1 1.42.0 2.2.0 1.14.0 4.4.0 4.3-0 1.46.0 2.12.0 1.0-7 2.30.0 1.7-3.1 2.26.0 1.18.0 3.0-6 6.0-85 1.12.8 1.7-3 3.28.0 0.4.1 1.12.0 1.0-14 1.68.0 5.40.0 3.0.1.2 1.2.4.2 1.60.0 1.2 3.42.0 1.5 1.2 2.20.0 1.16.1 1.78.0 1.0.12 4.9.1 1.48.0 4.6-14 1.60.0 3.34.0 3.8.1

RUN Rscript -e 'devtools::install_github("cnsb-boston/MetaboAnalystR", build_vignettes=FALSE);'


