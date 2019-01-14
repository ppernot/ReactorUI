function(request) {
  source_ui <- function(...) {
    source(
      file.path("ui_files", ...),
      local = TRUE
    )$value
  }

  navbarPage(
    "Reactor",
    theme = shinythemes::shinytheme(
      c("cosmo", "cerulean", "spacelab", "yeti")[3]
    ),
    tabPanel(
      title = "Project",
      source_ui("project.R")
    ),
    tabPanel(
      title = "Analysis",
      source_ui("analysis.R")
    ),
    tabPanel(
      title = "Fluxes",
      source_ui("fluxes.R")
    ),
    # tabPanel(
    #   title = "Downloads",
    #   source_ui("downloads.R")
    # ),
    tabPanel(
      title = "About",
      source_ui("about.R")
    )
  )
}
