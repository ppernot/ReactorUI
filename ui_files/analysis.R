sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    fluidRow(
      column(
        width = 12,
        actionButton(
          'loadMC',
          'Load MC results',
          icon = icon('clone')
        ),
        br(),
        h4("Summary"),
        withSpinner(
          verbatimTextOutput(
            "loadMsg"
          ),
          type=4
        )
      )
    )
  ),
  mainPanel(
    width = mainWidth,
    wellPanel(
      tabsetPanel(
        tabPanel(
          title=h4("Kinetics"),
          br(),
          fluidRow(
            column(
              width = 3,
              wellPanel(
                checkboxInput(
                  "mcPlot",
                  "Show error bands",
                  value = TRUE
                ),
                checkboxInput(
                  "colSel",
                  "Fixed colors",
                  value = TRUE
                ),
                checkboxInput(
                  "ppscale",
                  "Draw PPM scale",
                  value = FALSE
                ),
                checkboxGroupInput(
                  "categsPlot",
                  "Categories",
                  choices = c(
                    "Neutrals"     = "neutrals",
                    "Ions"         = "ions",
                    "Radicals only"= "radicals",
                    "Hydrocarbons" = "hydrocarbons",
                    "N-bearing"    = "N-bearing",
                    "O-bearing"    = "O-bearing",
                    "C0"           = "C0",
                    "C1"           = "C1",
                    "C2"           = "C2",
                    "C3"           = "C3",
                    "C4"           = "C4",
                    "C5"           = "C5",
                    "C6"           = "C6",
                    "C>6"          = "Cmore"
                  ),
                  # selected = c("neutrals","ions","hydrocarbons","N-bearing",
                  #              "C0","C1","C2","C3","C4","C5","C6","Cmore")
                  selected = c("neutrals","hydrocarbons","C0","C1","C2")
                )
              )
            ),
            column(
              width = 9,
              withSpinner(
                plotOutput(
                  "kinetics",
                  height = plotHeight,
                  dblclick = "kinetics_dblclick",
                  brush = brushOpts(
                    id = "kinetics_brush",
                    resetOnNew = TRUE
                  )
                ),
                type=4
              )
            )
          )
        ),
        tabPanel(
          title=h4("Sensitivity"),
          br(),
          fluidRow(
            column(
              width = 3,
              wellPanel(
                radioButtons(
                  "anaType",
                  "Sensitivity indices",
                  choices = c(
                    "Rank Correl." = "spearman",
                    "dCorr"        = "dcorr",
                    "dHSIC"        = "hsic"
                  )
                ),
                textInput(
                  "SASpecies",
                  "Choose species",
                  width = "50%",
                  value = NULL,
                  placeholder = "Type species"
                ),
                actionButton(
                  'doSA',
                  'Run SA',
                  icon = icon('gear')
                ),
                hr(),
                radioButtons(
                  "SAPlotType",
                  "Plot Type",
                  choices = c(
                    "Bar Plot"    = "barplot",
                    "Scatterplot" = "scatterplot"
                  )
                )
              )
            ),
            column(
              width = 9,
              withSpinner(
                plotOutput(
                  "sensitivity",
                  height = plotHeight
                ),
                type=4
              )
            )
          )
        ),
        id="analysisTabset"
      )
    )
  )
)
