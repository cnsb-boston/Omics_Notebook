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

RUN apt-get update && apt-get install -y \
	ca-certificates \
	curl \
	dirmngr \
	gnupg \
	libasound2 \
	libdbus-glib-1-2 \
	libgtk-3-0 \
	libxrender1 \
	libx11-xcb-dev \
	libx11-xcb1 \
	libxt6 \
	xz-utils \
	file \
	--no-install-recommends \
	&& rm -rf /var/lib/apt/lists/*


# Create user
RUN useradd -ms /bin/bash docker
USER docker


# Start at prompt
CMD ["/bin/bash"]