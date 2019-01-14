sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    fluidRow(
      column(
        width = 6,
        actionButton(
          'calcFlux',
          'Compute fluxes',
          icon = icon('gear')
        )
      )
    ),
    br(),
    fluidRow(
      column(
        width = 12,
        textInput(
          "flSpec",
          "Choose species",
          width = "50%",
          value = "CH4",
          placeholder = "Type species"
        ),
        sliderInput(
          "topShow",
          "Threshold",
          min   = 0.001,
          max   = 0.5,
          step  = 0.01,
          value = 0.05
        )
      )
    )
  ),
  mainPanel(
    width = mainWidth,
    wellPanel(
      tabsetPanel(
        tabPanel(
          title=h4("ViewFlow"),
          br(),
          withSpinner(
            plotOutput(
              "viewFlow",
              height = 700
            ),
            type=4
          )
        ),
        tabPanel(
          title=h4("Budget/Target"),
          br(),
          withSpinner(
            verbatimTextOutput(
              "viewBudget"
            ),
            type=4
          ),
          withSpinner(
            verbatimTextOutput(
              "viewTarget"
            ),
            type=4
          )
        ),
        id="fluxesTabset"
      )
    )
  )
)
