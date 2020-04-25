FROM r-base:3.6.1

CMD apt-get install libxml2-dev



COPY install.R /home/install.R

RUN Rscript home/install.R


