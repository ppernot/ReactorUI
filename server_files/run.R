output$nMCRunSelect <- renderUI({

  # Max runs to max MC data samples
  projectDir  = ctrlPars$projectDir
  InputMCDir = paste0(projectDir,'/MC_Input/Reactions')
  maxMC = length(
    list.files(path = InputMCDir, pattern = '.csv'))
  maxMC = maxMC - 1 # Do not count nominal run

  choices = c(0, maxMC)
  if (maxMC > 10 & maxMC <= 100)
    choices = c(0, seq(10, maxMC, by = 10))
  if (maxMC > 100)
    choices = c(0, 10, seq(100, maxMC, by = 100))

  ui <- list(
    selectInput(
      'nMCRun',
      label    = '# MC Runs (0: nominal)',
      choices  = choices,
      width    = '200px'
    )
  )
  ui
})

observeEvent(
  input$reactorRun, {

    # Nb MC runs
    nMC = as.numeric(input$nMCRun)
    print(nMC)
  }
)
