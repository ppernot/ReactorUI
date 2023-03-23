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
      ),
      br(),br(),br(),
      wellPanel(
        h4("About"),
        hr(),
        h5("Author      : P. Pernot"),
        h5("Affiliation : ",a(href="https://www.cnrs.fr/","CNRS")),
        h5(paste0("Version     : ",version)),
        # h5("Date        : 2022/09/19"),
        br(),
        a(href="https://github.com/ppernot/ReactorUI","How to cite..."),
        br(),
        a(href="https://github.com/ppernot/ReactorUI","code@github"),
        br(),
        a(href="https://github.com/ppernot/ReactorUI/issues",
          "Bugs report, Features request"),
        br(),
        a(href="https://ppernot.github.io/ReactorUI",
          "User's manual")
      )
    )
  ),
  mainPanel(
    width = mainWidth,
    verbatimTextOutput("selectMsg"),
    verbatimTextOutput("contentsMsg")
  )
)
