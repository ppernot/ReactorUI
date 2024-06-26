sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    h4("Model control"),
    # hr( style="border-color: #666;"),
    verbatimTextOutput("contentsNmlMsg")
  ),
  mainPanel(
    width = mainWidth,
    wellPanel(
      tabsetPanel(
        tabPanel(
          title=h4("Chemistry"),
          wellPanel(
            tabsetPanel(
              tabPanel(
                title=h4("Generate"),
                br(),
                fixedRow(
                  column(
                    width = 3,
                    uiOutput("chemistryParams"),
                    actionButton(
                      'generateNetwork',
                      'Generate Reactions',
                      icon = icon('gear')
                    )
                  ),
                  column(
                    width = 9,
                    shinycssloaders::withSpinner(
                      tagAppendAttributes(
                        verbatimTextOutput(
                          "summaryScheme"
                        ),
                        style="white-space:pre-wrap; text-align: left;
                               overflow-y:scroll; max-height: 450px;"
                      ),
                      type = 7
                    )
                  )
                )
              ),
              tabPanel(
                title=h4("Simplify"),
                br(),
                uiOutput("chemistrySimplify")
              ),
              tabPanel(
                title=h4("Visualize"),
                br(),
                fixedRow(
                  column(
                    width = 3,
                    checkboxGroupInput(
                      'netCtrl',
                      '',
                      choiceNames = list(
                        'Digraph','Legend'
                      ),
                      choiceValues = list(
                        'digraph','legend'
                      ),
                      inline = TRUE,
                      selected = 'legend'
                    ),
                    sliderInput(
                      'forceNetCharge',
                      'Nodes attraction',
                      min   = -200,
                      max   =  -10,
                      step  =   10,
                      value =  -30
                    ),
                    sliderInput(
                      'vlpMax',
                      'Max Volpert Index',
                      min   =  0,
                      max   = 10,
                      step  =  1,
                      value = 10
                    ),
                    shiny::radioButtons(
                      'netColoring',
                      'Color scheme',
                      choiceNames = list(
                        'Volpert','Charge','Radicals','Composition','Heavy atoms'
                      ),
                      choiceValues = list(
                        'volpert','charge','radicals','compo','mass'
                      )
                    ),
                    # sliderInput(
                    #   'fontSizeNet',
                    #   'Labels size',
                    #   min   = 6,
                    #   max   = 24,
                    #   step  = 1,
                    #   value = 12
                    # ),
                    sliderInput(
                      'linkDensNet',
                      'Links opacity',
                      min   = 0,
                      max   = 1,
                      step  = 0.1,
                      value = 0.2
                    )
                  ),
                  column(
                    width = 9,
                    visNetworkOutput(
                      "netScheme",
                      height = plotHeight
                    )
                  )
                )
              ),
              tabPanel(
                title=h4("Reactions"),
                br(),
                fixedRow(
                  column(
                    width = 3,
                    textInput(
                      "targetSpecies",
                      label = NULL,
                      value = NA,
                      placeholder = "Target sp."
                    )
                  ),
                  column(
                    width = 9,
                    shinycssloaders::withSpinner(
                      dataTableOutput(
                        "tabScheme",
                        height = "auto"
                      ),
                      type = 7
                    )
                  )
                )
              ),
              tabPanel(
                title=h4("Sinks"),
                br(),
                fixedRow(
                  column(
                    width = 12,
                    tagAppendAttributes(
                      verbatimTextOutput(
                        "quality"
                      ),
                      style="white-space:pre-wrap; text-align: left;
                               overflow-y:scroll; max-height: 450px;"
                    )
                  )
                )
              ),
              tabPanel(
                title=h4("Sample"),
                br(),
                fixedRow(
                  column(
                    width = 3,
                    uiOutput("nMCButton"),
                    actionButton(
                      'sampleChem',
                      'Generate Samples',
                      icon = icon('gear')
                    )
                  )
                )
              )
            )
          )
        ),
        tabPanel(
          title=h4("Irradiation"),
          br(),
          fixedRow(
            column(
              width = 4,
              uiOutput("irradUI"),
              uiOutput("irradUI2"),
              verbatimTextOutput("irradParams")
            ),
            column(
              width = 8,
              plotOutput("irradSpectrum",
                         height = 600)
            )
          )
        ),
        tabPanel(
          title=h4("Reactor"),
          br(),
          uiOutput("reactorUI"),
          verbatimTextOutput("reactorParams")
        ),
        tabPanel(
          title=h4("ChemDB Versions"),
          br(),
          fixedRow(
            column(
              width = 4,
              selectInput(
                'phoVers',
                label    = 'PhotoProcs',
                choices  = list(0,1),
                selected = 0,
                width    = '350px'
              ),
              selectInput(
                'neuVers',
                label    = 'Neutrals',
                choices  = list(0,1),
                selected = 0,
                width    = '350px'
              ),
              selectInput(
                'ionVers',
                label    = 'Ions',
                choices  = list(0,1),
                selected = 0,
                width    = '350px'
              )
            ),
            column(
              width = 4,
              selectInput(
                'speReso',
                label    = 'Spectral Resolution',
                choices  = list(0,1),
                selected = 0,
                width    = '350px'
              )
            )
          ),
          verbatimTextOutput("checkChanges")
        )
      )
    )
  )
)
