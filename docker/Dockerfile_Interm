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

# Probably useless, to be removed...
RUN apt-get -y install libgfortran4

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean

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

## renv
ENV RENV_VERSION 0.11.0-17
RUN Rscript -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN Rscript -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY ReactorUI/.Rprofile /ReactorUI/
COPY ReactorUI/renv/activate.R /ReactorUI/renv/activate.R
COPY ReactorUI/renv.lock /ReactorUI/renv.lock
WORKDIR /ReactorUI
RUN R -e 'renv::consent(TRUE)' \
 && R -e 'renv::restore()'

## app files
COPY ReactorUI/global.R /ReactorUI/
COPY ReactorUI/ui.R /ReactorUI/
COPY ReactorUI/server.R /ReactorUI/
COPY ReactorUI/ui_files /ReactorUI/ui_files/
COPY ReactorUI/server_files /ReactorUI/server_files/
COPY ReactorUI/R /ReactorUI/R/

# expose port
EXPOSE 3838

# run app on container start
WORKDIR /ReactorUI
CMD ["R", "-e", "shiny::runApp(host = '0.0.0.0', port = 3838)"]
