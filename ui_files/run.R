sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    h4("Run Reactor"),
    hr( style="border-color: #666;"),
    uiOutput("nMCRunSelect"),
    actionButton(
      'reactorRun',
      'Start',
      icon = icon('gear')
    )
  ),
  mainPanel(
    width = mainWidth
  )
)
