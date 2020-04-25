FROM r-base:3.6.1

COPY install.R /home/install.R

RUN Rscript home/install.R


