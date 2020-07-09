sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    h4("Select Project Directory"),
    hr( style="border-color: #666;"),
    shinyDirButton(
      id = "projectDir",
      label = "Choose",
      title = "Choose Project Directory",
      buttonType = "default",
      class = NULL,
      icon = shiny::icon('folder'),
      style = NULL)
  ),
  mainPanel(
    width = mainWidth,
    verbatimTextOutput("selectMsg"),
    verbatimTextOutput("contentsMsg")
  )
)
