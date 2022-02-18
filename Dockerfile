FROM omicsnotebook-base

COPY DESCRIPTION /tmp/ON/
RUN R -e "remotes::install_deps('/tmp/ON',repos='https://cran.rstudio.com')"
COPY ./ /home
RUN R -e "remotes::install_local('/home', dependencies=F, force=T)"
RUN rm -rf /tmp/ON /home/Notebook.py

# Create user
RUN useradd -ms /bin/bash docker
USER docker

# Start at prompt
CMD ["/bin/bash"]
