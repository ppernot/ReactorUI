sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    h4("Run Reactor"),
    hr( style="border-color: #666;"),
    uiOutput("nMCRunSelect"),
    fluidRow(
      column(
        width=12,
        actionButton(
          'reactorRun',
          'Start',
          icon = icon('gear')
        )#,
        # actionButton(
        #   'updateOut',
        #   'Update',
        #   icon = icon('file')
        # )
      )
    )
  ),
  mainPanel(
    width = mainWidth,
    wellPanel(
      tagAppendAttributes(
        verbatimTextOutput(
          "reactorOutput"
        ),
        style="white-space:pre-wrap; text-align: left;
               overflow-y:scroll; max-height: 400px;"
      ),
      tagAppendAttributes(
        verbatimTextOutput(
          "reactorErrors"
        ),
        style="white-space:pre-wrap; text-align: left;
               overflow-y:scroll; max-height: 100px;"
      )
    )
  )
)
