sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    h4("Run Reactor"),
    hr(),
    uiOutput("nMCRunSelect"),
    tabsetPanel(
      tabPanel(
        title = h5('Local'),
        h3(' '),
        fluidRow(
          column(
            width=6,
            checkboxInput(
              'logTail',
              'Log tail only',
              value = TRUE
            )
          ),
          column(
            width=6,
            actionButton(
              'reactorRun',
              'Start',
              icon = icon('gear')
            )
          )
        )
      ),
      tabPanel(
        title = h5('Cloud (TBD)'),
        # h3(' '),
        # uiOutput("cloudSelect"),
        # fluidRow(
        #   column(
        #     width=6,
        #     actionButton(
        #       'createCluster',
        #       'Create Cluster',
        #       icon = icon('gear')
        #     ),
        #     actionButton(
        #       'cleanCluster',
        #       'Clean Cluster',
        #       icon = icon('gear')
        #     ),
        #     actionButton(
        #       'deleteCluster',
        #       'Delete Cluster',
        #       icon = icon('gear')
        #     )
        #   ),
        #   column(
        #     width=6,
        #     actionButton(
        #       'runCluster',
        #       'Start Runs',
        #       icon = icon('gear')
        #     ),
        #     actionButton(
        #       'checkCluster',
        #       'Update',
        #       icon = icon('gear')
        #     )
        #   )
        # )
      )
    )
  ),
  mainPanel(
    width = mainWidth,
    wellPanel(
      tagAppendAttributes(
        verbatimTextOutput(
          "reactorOutput"
        ),
        style="white-space:pre-line;
               text-align: left;
               height: 450px;
               overflow-y:scroll;
               overflow-x:scroll;"
      ),
      tagAppendAttributes(
        verbatimTextOutput(
          "reactorErrors"
        ),
        style="white-space:pre-wrap;
               text-align: left;
               overflow-y:scroll;
               height: 100px;
               overflow-x:scroll;"
      )
    )
  )
)
