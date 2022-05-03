# CRAN Libraries ####
libs <- c(
  # "devtools",
  "Hmisc",
  "shiny",
  "shinyBS",
  "shinycssloaders",
  "shinyFiles",
  "DT",
  "igraph",
  "dHSIC",
  "fda.usc",
  "xtable",
  "inlmisc",
  "CHNOSZ",
  "networkD3",
  "promises",
  "future",
  "future.apply",
  "rlist",
  "data.table",
  "visNetwork"
)
for  (lib in libs) {
  # Installation is managed in the container
  # if (!require(lib, character.only = TRUE, quietly = TRUE)) {
  #   install.packages(
  #     lib,
  #     dependencies = TRUE,
  #     repos = "https://cran.rstudio.com"
  #   )
  #  }
  library(lib, character.only = TRUE, quietly = TRUE)
}

# Github
# if (!require("nml", character.only = TRUE, quietly = TRUE)) {
#   devtools::install_github("jsta/nml")
#   library("nml", quietly = TRUE)
# }

# Parallelism by future/promises
plan(multiprocess)
