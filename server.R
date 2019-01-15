function(input, output, session) {

  # Initialize ####
  if (!dir.exists("outputDir")) {
    dir.create("outputDir", showWarnings = FALSE)
  }

  # Load Server files ####
  files <- c(
    "project.R" ,
    "model.R",
    "analysis.R",
    "fluxes.R"
  )

  for (f in files)
    source(
      file.path("server_files", f),
      local = TRUE
    )
}
