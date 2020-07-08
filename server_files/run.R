form = function(x) {

  if(is.logical(x)){
    val = '.false.'
    if(x)
      val = '.true.'
    return(val)
  }

  if(is.vector(x) & length(x) > 1) {
    if(is.numeric(x[1])) {
      sel = !is.na(x)
      count = 10 - sum(sel)
      val = paste(
        paste(signif(x[sel],5),collapse=','),
        ',',
        paste0(count,'*0')
      )
    } else {
       val = paste(x,collapse=',')
    }
    return(val)
  }

  if(is.numeric(x) | is.character(x))
    return(x)



}
generateUpdatedControl = function(ll) {
  namelist = '&REAC_DATA\n'
  for(n in names(reacData()) )
    namelist = paste(
      namelist,
      n, '=',  form( reacData()[[n]] ), ',\n'
    )
  namelist = paste(namelist,'/')
print(namelist)
}



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

    # Generate updated control.dat
    generateUpdatedControl(reacData()); stop()

    # Clean standard outputs
    stdout = paste0(projectDir(),'/Run/runOut.txt')
    if(file.exists(stdout))
      file.remove(stdout)
    stderr = paste0(projectDir(),'/Run/runErr.txt')
    if(file.exists(stderr))
      file.remove(stderr)

    # Nb MC runs
    nMC = as.numeric(input$nMCRun)

    if (nMC == 0) # Nominal run
      running(
        try(
          system2(
            command = paste0(projectDir(),'/Scripts/OneRun_Loc.sh'),
            args    = c('0',projectDir()),
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
            command = paste0(projectDir(),'/Scripts/MCRun_Loc.sh'),
            args    = c(nMC,projectDir()),
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
  paste0(ctrlPars$ProjectDir,'/Run/runErr.txt'),
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
  paste0(ctrlPars$ProjectDir,'/Run/runOut.txt'),
  readLines
)
output$reactorOutput <- renderPrint({
  if (is.null(running())) {
    cat('Please run code ...')
    return()
  }
  stdOut()
})

