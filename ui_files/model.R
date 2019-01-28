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
              uiOutput("irradUI"),
              verbatimTextOutput("irradParams")
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
          fixedRow(
            column(
              width = 3,
              uiOutput("chemistryParams"),
              actionButton(
                'generateNetwork',
                'Generate Reactions',
                icon = icon('gear')
              ),
              hr(),
              sliderInput(
                'forceNetCharge',
                'Nodes attraction',
                min   = -150,
                max   =    0,
                step  =   10,
                value = -100
              ),
              sliderInput(
                'vlpMax',
                'Max Volpert Index',
                min   =  0,
                max   = 10,
                step  =  1,
                value = 10
              )
            ),
            column(
              width = 9,
              # plotOutput(
              forceNetworkOutput(
                "plotScheme",
                height = plotHeight
              )
            )
          )
        )
      )
    )
  )
)
