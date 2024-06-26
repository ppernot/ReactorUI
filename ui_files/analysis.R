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
        shinycssloaders::withSpinner(
          verbatimTextOutput(
            "loadMsg"
          ),
          type = 7
        )
      )
    )
  ),
  mainPanel(
    width = mainWidth,
    wellPanel(
      tabsetPanel(
        tabPanel(
          title=h4("Sanity checks"),
          wellPanel(
            tabsetPanel(
              tabPanel(
                "Outputs",
                shinycssloaders::withSpinner(
                  verbatimTextOutput(
                    "sanityOutputs"
                  ),
                  type = 7
                )
              ),
              tabPanel(
                "Integration",
                plotOutput(
                  "sanityInteg",
                  dblclick = "sanity_dblclick",
                  brush = brushOpts(
                    id = "sanity_brush",
                    resetOnNew = TRUE
                  ),
                  height = plotHeight
                )
              ),
              tabPanel(
                "Spectral radius",
                plotOutput(
                  "sanitySR",
                  dblclick = "sanitySR_dblclick",
                  brush = brushOpts(
                    id = "sanitySR_brush",
                    resetOnNew = TRUE
                  ),
                  height = plotHeight
                )
              )
            )
          )
        ),
        tabPanel(
          title=h4("Kinetics"),
          br(),
          fluidRow(
            column(
              width = 3,
              wellPanel(
                uiOutput(
                  "categsPlot"
                ),
                strong("Plot"),
                hr(),
                checkboxInput(
                  "mcPlot",
                  "Show error bands",
                  value = FALSE
                ),
                checkboxInput(
                  "colSel",
                  "Fixed colors",
                  value = TRUE
                ),
                checkboxInput(
                  "colByCompo",
                  "Color by composition",
                  value = TRUE
                ),
                checkboxInput(
                  "ppscale",
                  "Draw PPM scale",
                  value = FALSE
                )
              )
            ),
            column(
              width = 9,
              shinycssloaders::withSpinner(
                plotOutput(
                  "kinetics",
                  height = plotHeight,
                  dblclick = "kinetics_dblclick",
                  brush = brushOpts(
                    id = "kinetics_brush",
                    resetOnNew = TRUE
                  )
                ),
                type = 7
              )
            )
          )
        ),
        tabPanel(
          title=h4("Pseudo-MS"),
          br(),
          fluidRow(
            column(
              width = 3,
              wellPanel(
                uiOutput(
                  "categsPlotMS"
                ),
                strong("Plot"),
                hr(),
                checkboxInput(
                  "mcPlotMS",
                  "Show error bars",
                  value = TRUE
                ),
                checkboxInput(
                  "ppScaleMS",
                  "Draw PPM scale",
                  value = FALSE
                ),
                sliderInput(
                  "timeMS",
                  "Log sampling time",
                  min   = -10,
                  max   =  10,
                  step  =   1,
                  value =  10
                ),
                sliderInput(
                  "threshMS",
                  "MS Amplitude Range",
                  min = -30,
                  max = 0,
                  step = 1,
                  value = c(-15, 0)
                ),
                sliderInput(
                  "widthMS",
                  "Adjust bar width",
                  min = 1,
                  max = 10,
                  step = 1,
                  value = 3
                ),
                sliderInput(
                  "MStext.cex",
                  "Size of names",
                  min = 0,
                  max = 9,
                  step = 1,
                  value = 5
                )
              )
            ),
            column(
              width = 9,
              shinycssloaders::withSpinner(
                plotOutput(
                  "pseudoMS",
                  height = plotHeight
                ),
                type = 7
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
                fluidRow(
                  column(6,
                         textInput(
                           "SASpecies",
                           "Target",
                           width = "100%",
                           value = NULL,
                           placeholder = "Species"
                         )
                  ),
                  column(6,
                         actionButton(
                           'doSA',
                           'Run SA',
                           icon = icon('gear')
                         ),
                         tags$style(
                           type = 'text/css',
                           "#doSA { width:100%; margin-top: 27px;}"
                         )
                  )
                ),
                radioButtons(
                  "anaType",
                  "Sensitivity indices",
                  choices = c(
                    "Rank Correl." = "spearman",
                    "dCorr"        = "dcorr",
                    "dHSIC"        = "hsic"
                  )
                ),
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
              shinycssloaders::withSpinner(
                plotOutput(
                  "sensitivity",
                  height = plotHeight
                ),
                type = 7
              )
            )
          )
        ),
        id="analysisTabset"
      )
    )
  )
)
