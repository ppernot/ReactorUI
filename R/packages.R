# CRAN Libraries ####
libs <- c(
  "devtools","Hmisc",
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
    library(lib, character.only = TRUE, quietly = TRUE)
  }
}

# Github
# if (!require("nml", character.only = TRUE, quietly = TRUE)) {
#   devtools::install_github("jsta/nml")
#   library("nml", quietly = TRUE)
# }

# Parallelism by future/promises
plan(multiprocess)
