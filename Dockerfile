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
        r-base-dev


# Install R libraries
COPY install.R /home/install.R

RUN Rscript home/install.R


# Install vnc, xvfb in order to create a 'fake' display and firefox
RUN     apt-get install -y x11vnc xvfb
RUN     mkdir ~/.vnc
# Setup a password
RUN     x11vnc -storepasswd 1234 ~/.vnc/passwd


# Create user
RUN useradd -ms /bin/bash docker
USER newuser


# Start at prompt
CMD ["/bin/bash"]