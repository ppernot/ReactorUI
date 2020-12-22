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
          value = "H",
          placeholder = "Type species"
        ),
        sliderInput(
          "topShow",
          "Threshold",
          min   = 0.001,
          max   = 1.00,
          step  = 0.01,
          value = 0.05
        ),
        checkboxInput(
          "curvedArrow",
          "Curved arrows",
          value = FALSE
        ),
        checkboxInput(
          "level",
          "Level = 2 (WIP)",
          value = FALSE
        ),
        selectInput(
          "fluxGraphAlgo",
          "Graph Algorithm",
          choices = list(
            "GEM"        = "GEM",
            "FR"         = "FR",
            # "DrL"        = "DrL",
            "LGL"        = "LGL",
            "Bipartite"  = "Bipartite"
            ),
          selected = "GEM",
          width = "50%"
        ),
        sliderInput(
          'forceNetChargeFlux',
          'Nodes attraction',
          min   = -160,
          max   =    0,
          step  =   20,
          value = -100
        ),
        shiny::checkboxInput(
          'scaleDistD3',
          label = 'Weigted distances',
          value = FALSE
        )
      )
    )
  ),
  mainPanel(
    width = mainWidth,
    wellPanel(
      tabsetPanel(
        tabPanel(
          title = h4("Budget/Target"),
          br(),
          h5("Productions & Losses"),
          verbatimTextOutput("viewBudget"),
          br(),
          h5("Main production path"),
          verbatimTextOutput("viewTarget")
        ),
        tabPanel(
          title=h4("ViewFlow-D3"),
          br(),
          # withSpinner(
          forceNetworkOutput(
            "viewFlowD3",
            height = 700
            # ),
            # type=4
          )
        ),
        tabPanel(
          title=h4("ViewFlow (WIP)"),
          br(),
          # withSpinner(
            plotOutput(
              "viewFlow",
              height = 700
            # ),
            # type=4
          )
        ),
        id="fluxesTabset"
      )
    )
  )
)
