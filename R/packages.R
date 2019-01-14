# Options ####
Sys.setlocale(category = "LC_NUMERIC", locale = "C")

options(
  shiny.maxRequestSize = 20 * 1024^2,
  width = 60,
  warn = 0
)

# options(shiny.json.digits=32)

# Libraries ####
libs <- c(
  "shiny", "shinyBS", "shinycssloaders", "shinyFiles",
  "DT", "tools", "igraph","dHSIC","fda.usc",
  "xtable", "inlmisc","CHNOSZ",
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

# Parallelism by future/promises
plan(multiprocess)
