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
    textInput(
      "flSpec",
      "Choose species",
      width = "50%",
      value = "",
      placeholder = "Type species"
    ),
    h4('Budget / Target'),
    hr(),
    sliderInput(
      "topShow",
      "Log-Threshold",
      min   = -15,
      max   =   0,
      step  =   1,
      value =  -1
    ),
    fluidRow(
      column(
        width = 6,
        checkboxInput(
          "loopAvoid",
          "Avoid loops",
          value = TRUE
        )
      ),
      column(
        width = 5,
        checkboxInput(
          "maxVolpert",
          "Shortest path",
          value = TRUE
        )
      )
    ),
    h4('ViewFlow'),
    hr(),
    fluidRow(
      column(
        width = 6,
        checkboxInput(
          "graphSimplify",
          "Simplify",
          value = FALSE
        )
      ),
      column(
        width = 6,
        shiny::checkboxInput(
          'scaleDistD3',
          label = 'Weigted distances',
          value = FALSE
        )
      )
    ),
    sliderInput(
      'forceNetChargeFlux',
      'Nodes attraction',
      min   = -160,
      max   =    0,
      step  =   20,
      value = -100
    ),
    sliderInput(
      'flFontSize',
      'Font Size',
      min   =  1,
      max   = 30,
      step  =  1,
      value = 10
    )
  ),
  mainPanel(
    width = mainWidth,
    wellPanel(
      tabsetPanel(
        tabPanel(
          title = h4("Budget/Target"),
          br(),
          fluidRow(
            column(
              width = 6,
              h5("Productions & Losses"),
              verbatimTextOutput("viewBudget")
            ),
            column(
              width = 6,
              h5("Main production path"),
              verbatimTextOutput("viewTarget")
            )
          )
        ),
        tabPanel(
          title = h4("ViewFlow"),
          style = " background-color: white;",
          br(),
          shinycssloaders::withSpinner(
            forceNetworkOutput("viewFlowD3", height = 700),
            type = 7
          ),
        ),
        # tabPanel(
        #   title=h4("ViewFlow (WIP)"),
        #   br(),
        #   shinycssloaders::withSpinner(
        #     plotOutput(
        #       "viewFlow",
        #       height = 700
        #     ),
        #     type = 7
        #   )
        # ),
        id="fluxesTabset"
      )
    )
  )
)
