FROM rocker/r-ver:4.1.1

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libnetcdf-dev \
        libxml2-dev \
        libssl-dev \
        libxt-dev \
        libcairo2-dev \
        pandoc \
        python3 \
        python3-tk \
        r-base-dev \
        libglpk-dev \
        libzmq5 \
        && rm -rf /var/lib/apt/lists/*

RUN R -e "options(repos = c(CRAN = 'https://cran.rstudio.com')); install.packages('remotes')"
RUN R -e "remotes::install_github('cnsb-boston/Omics_Notebook',repos='https://cran.rstudio.com')"

# Create user
 RUN useradd -ms /bin/bash docker
 USER docker

# Start at prompt
CMD ["/bin/bash"]
