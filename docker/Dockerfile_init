# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:3.6.3

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    gfortran \ 
    libgdal-dev \
    gdal-bin\
    libgsl0-dev\
    libudunits2-dev

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean


## renv.lock file
#COPY ../renv.lock ./renv.lock

# install renv & restore packages
#RUN Rscript -e 'install.packages("renv")'
#RUN Rscript -e 'renv::restore()'

#  TBD: remove 'tools' from list; 'XML' requires R-4.0
RUN R -e "utils::install.packages(c('devtools','Hmisc','shinyBS', 'shinycssloaders', 'shinyFiles','DT','tools', 'igraph', 'dHSIC', 'fda.usc', 'rlist','xtable', 'inlmisc','CHNOSZ','networkD3','promises','future','future.apply','shinythemes','XML'), dependencies = TRUE)"
# 'HiveR' # Not used at this stage; would require openGL (see CRAN)....

# Add necessary files

## data
ADD Reactor_Runs/ChemDB.tar.gz /

# Mounting point for local Projects folder
RUN mkdir /Projects

## reactor
COPY Reactor /Reactor/
RUN  cd Reactor \
  && ./compil \
  && mv reactor /ReactorCodes \
  && cd /

## app folder
COPY ReactorUI/global.R /app/
COPY ReactorUI/ui.R /app/
COPY ReactorUI/server.R /app/
COPY ReactorUI/ui_files /app/ui_files/
COPY ReactorUI/server_files /app/server_files/
COPY ReactorUI/R ./app/R/

# expose port
EXPOSE 3838

# run app on container start
WORKDIR /app
CMD ["R", "-e", "shiny::runApp(host = '0.0.0.0', port = 3838)"]
