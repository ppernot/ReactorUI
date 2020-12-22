sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    h4("Model control"),
    hr( style="border-color: #666;"),
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
                    # withSpinner(
                      tagAppendAttributes(
                        verbatimTextOutput(
                          "summaryScheme"
                        ),
                        style="white-space:pre-wrap; text-align: left;
                               overflow-y:scroll; max-height: 450px;"
                      # ),
                      # type=4
                    )
                  )
                )
              ),
              tabPanel(
                title=h4("Network"),
                br(),
                fixedRow(
                  column(
                    width = 3,
                    checkboxGroupInput(
                      'netCtrl',
                      '',
                      choiceNames = list(
                        'Digraph','Full names','Legend'
                      ),
                      choiceValues = list(
                        'digraph','showNames','legend'
                      ),
                      inline = TRUE,
                      selected = c('showNames','legend')
                    ),
                    sliderInput(
                      'forceNetCharge',
                      'Nodes attraction',
                      min   = -500,
                      max   =    0,
                      step  =   25,
                      value = -100
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
                      'Nodes colors',
                      choiceNames = list(
                        'Volpert','Charge','Radicals','C/N/O','Mass'
                      ),
                      choiceValues = list(
                        'volpert','charge','radicals','C/N/O','Mass'
                      )
                    ),
                    sliderInput(
                      'fontSizeNet',
                      'Labels size',
                      min   = 6,
                      max   = 24,
                      step  = 1,
                      value = 12
                    ),
                    sliderInput(
                      'linkDensNet',
                      'Links opacity',
                      min   = 1,
                      max   = 9,
                      step  = 1,
                      value = 3
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
                    # withSpinner(
                      dataTableOutput(
                        "tabScheme",
                        height = "auto"
                      ) #,
                    #   type=4
                    # )
                  )
                )
              ),
              tabPanel(
                title=h4("Checks"),
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
          uiOutput("chemDBVersions"),
          verbatimTextOutput("checkChanges")
        )
      )
    )
  )
)
