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
          wellPanel(
            tabsetPanel(
              tabPanel(
                title=h4("Initial mix"),
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
                    withSpinner(
                      tagAppendAttributes(
                        verbatimTextOutput(
                          "summaryScheme"
                        ),
                        style="white-space:pre-wrap; text-align: left;
                               overflow-y:scroll; max-height: 600px;"
                      ),
                      type=4
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
              ),
              tabPanel(
                title=h4("Reac. List"),
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
                    withSpinner(
                      dataTableOutput(
                        "tabScheme",
                        height = "auto"
                      ),
                      # tagAppendAttributes(
                      #   verbatimTextOutput(
                      #     "listScheme"
                      #   ),
                      #   style="white-space:pre-wrap; text-align: left;
                      #   overflow-y:scroll; max-height: 600px;"
                      # ),
                      type=4
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)
