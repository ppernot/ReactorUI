sidebarLayout(
  sidebarPanel(
    width = sideWidth,
    tabsetPanel(
      tabPanel(
        title = h4('Open/New'),
        h3(' '),
        checkboxInput(
          'newProj',
          label = 'New Project'
        ),
        conditionalPanel(
          condition = "input.newProj != 0",
          selectInput(
            'newProjTemplate',
            label = 'Template',
            choices = c('APSIS','Titan')
          )
        ),
        shinyDirButton(
          id = "projectDir",
          label = "Select a Project",
          title = "Select a Project Directory",
          buttonType = "default",
          class = NULL,
          icon = shiny::icon('folder'),
          style = NULL)
      ),
      tabPanel(
        title = h4('Save'),
        h3(' '),
        actionButton(
          'saveProj',
          'Save Project',
          icon = icon('gear')
        )

      )
    )
  ),
  mainPanel(
    width = mainWidth,
    verbatimTextOutput("selectMsg"),
    verbatimTextOutput("contentsMsg")
  )
)
