# Options ####
Sys.setlocale(category = "LC_NUMERIC", locale = "C")

options(
  shiny.maxRequestSize = 20 * 1024^2,
  width = 60,
  warn = 0
)

# options(shiny.json.digits=32)

# CRAN Libraries ####
libs <- c(
  "devtools",
  "shiny", "shinyBS", "shinycssloaders", "shinyFiles",
  "DT", "tools", "igraph", "dHSIC", "fda.usc", "rlist",
  "xtable", "inlmisc","CHNOSZ","networkD3",
  "promises","future","future.apply"
)
for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(
      lib,
      dependencies = TRUE,
      repos = "https://cran.univ-paris1.fr"
    )
    library(lib, quietly = TRUE)
  }
}

# Github
# if (!require("nml", character.only = TRUE, quietly = TRUE)) {
#   devtools::install_github("jsta/nml")
#   library("nml", quietly = TRUE)
# }

# Parallelism by future/promises
plan(multiprocess)
