FROM rocker/r-ver:3.6.1

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
        && rm -rf /var/lib/apt/lists/*


# Install R libraries
COPY install.R /home/install.R

RUN Rscript home/install.R


# Start at prompt
CMD ["/bin/bash"]