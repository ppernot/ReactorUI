output$nMCRunSelect <- renderUI({

  # Max runs to max MC data samples
  InputMCDir = paste0(projectDir(),'/MC_Input/Reactions')
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
    ),
    numericInput(
      'nbSnap',
      label = '# snapshots per run ',
      value = reacData()$nbSnapshots
    )
  )
  ui
})
observeEvent(
  # Update nb snapshots
  input$nbSnap, ({
    ll = reacData()
    ll$nbSnapshots = input$nbSnap
    reacData(ll)
  })
)

running = reactiveVal(NULL)
observeEvent(
  input$reactorRun, {

    projectDir  = ctrlPars$projectDir

    # Clean standard outputs
    stdout = paste0(projectDir,'/Run/runOut.txt')
    if(file.exists(stdout))
      file.remove(stdout)
    stderr = paste0(projectDir,'/Run/runErr.txt')
    if(file.exists(stderr))
      file.remove(stderr)

    # Nb MC runs
    nMC = as.numeric(input$nMCRun)

    if (nMC == 0) # Nominal run
      running(
        try(
          system2(
            command = paste0(projectDir,'/Scripts/OneRun_Loc.sh'),
            args    = c('0',projectDir),
            stdout  = stdout,
            stderr  = stderr,
            wait    = FALSE ),
          silent = TRUE
        )
      )
    else
      running(
        try(
          system2(
            command = paste0(projectDir,'/Scripts/MCRun_Loc.sh'),
            args    = c(nMC,projectDir),
            stdout  = stdout,
            stderr  = stderr,
            wait    = FALSE ),
          silent = TRUE
        )
      )

  }
)
stdErr = reactiveFileReader(
  500, session,
  paste0(ctrlPars$projectDir,'/Run/runErr.txt'),
  readLines
)
output$reactorErrors <- renderPrint({
  if (is.null(running())) {
    cat('Please run code ...')
    return()
  }
  cat('Error messages:\n')
  runMsg =''
  if(running() != 0)
    runMsg = cat(running(),'\n')
  cat(
    runMsg,
    paste(stdErr(),collapse ='\n')
  )
})
stdOut = reactiveFileReader(
  500, session,
  paste0(ctrlPars$projectDir,'/Run/runOut.txt'),
  readLines
)
output$reactorOutput <- renderPrint({
  if (is.null(running())) {
    cat('Please run code ...')
    return()
  }
  stdOut()
})

