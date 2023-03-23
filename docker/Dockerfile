FROM  ppernot1/shiny_base

WORKDIR /

# Mounting point for local Projects folder
RUN mkdir /Projects

# Unpack misc scripts
ADD Reactor_Runs/ReacCodes.tar.gz /

# Install code
COPY Reactor /Reactor/
WORKDIR /Reactor
RUN  ./compil \
  && mv reactor /ReactorCodes
WORKDIR /

# Install packages

## renv
ENV RENV_VERSION 0.11.0-17
RUN Rscript -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN Rscript -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

## restore packages
COPY ReactorUI/.Rprofile /ReactorUI/
COPY ReactorUI/renv/activate.R /ReactorUI/renv/activate.R
COPY ReactorUI/renv.lock /ReactorUI/renv.lock
WORKDIR /ReactorUI
RUN R -e 'renv::consent(TRUE)' \
 && R -e 'renv::restore()'

# Install App
COPY ReactorUI/global.R /ReactorUI/
COPY ReactorUI/ui.R /ReactorUI/
COPY ReactorUI/server.R /ReactorUI/
COPY ReactorUI/ui_files /ReactorUI/ui_files/
COPY ReactorUI/server_files /ReactorUI/server_files/
COPY ReactorUI/R /ReactorUI/R/

# Chemical database
# ADD Reactor_Runs/ChemDB.tar.gz /

# Add mounting point for local ChemDB
RUN mkdir /ChemDBLocal

# Expose port & run app on container start
EXPOSE 3838
WORKDIR /ReactorUI
CMD ["R", "-e", "shiny::runApp(host = '0.0.0.0', port = 3838)"]






