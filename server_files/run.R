form = function(x) {
  if (is.logical(x)) {
    val = ifelse(x, '.true.', '.false.')
    return(val)
  }

  if (is.vector(x) & length(x) > 1) {
    val = paste(x, collapse = ',')
    return(val)
  }

  if (is.character(x)) {
    val = paste0("'", x, "'")
    return(val)
  }

  if (is.numeric(x))
    return(x)

}
formSp = function(x, nsp) {
  count = 10 - nsp

  if (is.numeric(x[1])) {
    val = paste0(paste0(signif(x[1:nsp], 5), collapse = ','),
                 ',', paste0(count, '*0'))
  } else {
    val = strsplit(x, ',')[[1]][1:nsp]
    val = paste0(paste0("'", val, "'", collapse = ','),
                 ',',
                 paste0(count, "*''"))
  }

  return(val)

}
generateUpdatedControl = function() {
  ll = reacData()
  namelist = '&REAC_DATA\n'
  for (n in names(ll)) {
    if (n == 'reactantsSpecies' |
        n == 'reactantsComposition') {
      nsp = sum(is.finite(ll$reactantsComposition))
      val = formSp(ll[[n]], nsp)
    } else {
      val = form(ll[[n]])
    }
    namelist = paste0(namelist,
                      ' ', n, '=',  val, ',\n')
  }
  namelist = paste0(namelist, '/')

  ll = chemDBData()
  namelist = paste0(namelist,'\n','&DB_DATA\n')
  for (n in names(ll)) {
    val = form(ll[[n]])
    namelist = paste0(namelist,
                      ' ', n, '=',  val, ',\n')
  }
  namelist = paste0(namelist, '/')

  # Save old control.dat, if any
  oldCtrl = paste0(projectDir(), '/Run/control.dat')
  if (file.exists(oldCtrl)) {
    savCtrl = paste0(projectDir(), '/Run/control.dat_sav')
    file.copy(from = oldCtrl,
              to = savCtrl,
              overwrite = TRUE)
  }
  writeLines(namelist, oldCtrl)
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
    conditionalPanel(
      condition = "input.nMCRun != 0",
      checkboxInput(
        'appendMC',
        label = 'Append to existing MC runs',
        value = FALSE
      )
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
    generateUpdatedControl()

    # Clean standard outputs
    stdout = paste0(projectDir(),'/Run/runOut.txt')
    if(file.exists(stdout))
      file.remove(stdout)
    stderr = paste0(projectDir(),'/Run/runErr.txt')
    if(file.exists(stderr))
      file.remove(stderr)

    # Nb MC runs
    nMC = as.numeric(input$nMCRun)

    if (nMC == 0) {
      # Nominal run
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
    } else {
      if(input$appendMC) {
        mcf = list.files(
          path = paste0(projectDir(),'/MC_Output'),
          pattern = 'fracmol_'
        )
        nmc = length(mcf)
        first = nmc # Run 0 always present
        nrun  = nMC - 1 # if append, do not do nominal
        nSamples = chemDBData()$nMC
        if(first+nrun > nSamples) {
          nrunMax = nSamples - first
          showModal(modalDialog(
            title = ">>>> Too many runs <<<< ",
            paste0('The cumulated number of runs (',first + nrun,
                   ') exceeds the number of chemDB samples.\n',
                   'This session will be truncated to ',nrunMax,' runs.\n',
                   'If you need more, generate more chemDB samples.'),
            easyClose = TRUE,
            footer = modalButton("OK"),
            size = 's'
          ))
          nrun = nrunMax
        }
      } else {
        first = 0
        nrun  = nMC
      }
      running(
        try(
          system2(
            command = paste0(projectDir(),'/Scripts/MCRun_Loc1.sh'),
            args    = c(first,nrun,projectDir()),
            stdout  = stdout,
            stderr  = stderr,
            wait    = FALSE ),
          silent = TRUE
        )
      )
    }
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

