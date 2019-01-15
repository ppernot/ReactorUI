sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    h4("Model controle"),
    hr( style="border-color: #666;"),
    verbatimTextOutput("contentsNmlMsg")
  ),
  mainPanel(
    width = mainWidth,
    wellPanel(
      tabsetPanel(
        tabPanel(
          title=h4("Reactor"),
          br(),
          verbatimTextOutput("reactorParams")
        ),
        tabPanel(
          title=h4("Irradiation"),
          br(),
          fixedRow(
            column(
              width = 4,
              verbatimTextOutput("irradParams"),
              uiOutput("irradUI")
            ),
            column(
              width = 8,
              plotOutput("irradSpectrum",
                         height = plotHeight)
            )
          )
        ),
        tabPanel(
          title=h4("Chemistry"),
          br(),
          verbatimTextOutput("chemistryParams")
        )
      )
    )
  )
)
