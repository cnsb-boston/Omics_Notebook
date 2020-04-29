FROM debian:sid-slim

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libnetcdf-dev \
        libxml2-dev \
        libssl-dev \
        libxt-dev \
        libcairo2-dev \
        pandoc \
        python3=3.6.7 \
        python3-tk \
        r-base-dev \
        && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get -y install --no-install-recommends --no-install-suggests \
        ca-certificates software-properties-common gnupg2 gnupg1 \
      && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
      && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
      && apt-get install r-base=3.6.3

# Install R libraries
COPY install.R /home/install.R

RUN Rscript home/install.R



# Create user
RUN useradd -ms /bin/bash docker
USER docker


# Start at prompt
CMD ["/bin/bash"]