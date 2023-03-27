function(request) {
  source_ui <- function(...) {
    source(
      file.path("ui_files", ...),
      local = TRUE
    )$value
  }

  navbarPage(strong(paste0("Reactor ",version)),
    theme = shinythemes::shinytheme(
      c("cosmo", "cerulean", "spacelab", "yeti")[3]
    ),
    selected = "Project", # The tags section takes the first entry !?!?
    tags$head(
      tags$style(
        HTML("hr { height: 1px; margin-top: 0.0em; background: #666;}")
      )
    ),
    tabPanel(
      title = "Project",
      source_ui("project.R")
    ),
    tabPanel(
      title = "Model",
      source_ui("model.R")
    ),
    tabPanel(
      title = "Run",
      source_ui("run.R")
    ),
    tabPanel(
      title = "Analysis",
      source_ui("analysis.R")
    ),
    tabPanel(
      title = "Fluxes",
      source_ui("fluxes.R")
    )
  )
}
