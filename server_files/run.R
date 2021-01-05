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
  oldCtrl = file.path(projectDir(), 'Run','control.dat')
  if (file.exists(oldCtrl)) {
    savCtrl = file.path(projectDir(), 'Run','control.dat_sav')
    file.copy(from = oldCtrl,
              to = savCtrl,
              overwrite = TRUE)
  }
  writeLines(namelist, oldCtrl)
}
output$nMCRunSelect <- renderUI({

  # Max runs to max MC data samples
  InputMCDir = file.path(projectDir(),'MC_Input','Reactions')
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
    # actionButton(
    #   'mcClean',
    #   'Clean Runs',
    #   icon = icon('eraser')
    # ),
    h3(''),
    numericInput(
      'nbSnap',
      label = '# snapshots per run ',
      value = reacData()$nbSnapshots,
      width = '200px'
    ),
    h3(''),
    sliderInput(
      "logRelErr",
      "Log relative error",
      min = -12,
      max =  -2,
      step =  1,
      value = log10(
        ifelse(
          is.null(reacData()$relativeError),
          REAC_DATA_default$relativeError,
          reacData()$relativeError
        )
      )
    ),
    h3(''),
    sliderInput(
      "logAbsErr",
      "Log absolute error",
      min = -12,
      max =  -2,
      step =  1,
      value = log10(
        ifelse(
          is.null(reacData()$absoluteError),
          REAC_DATA_default$absoluteError,
          reacData()$absoluteError
        )
      )
    ),
    h3('')
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
observeEvent(
  # Update integrator relative error
  input$logRelErr, ({
    ll = reacData()
    ll$relativeError = 10^input$logRelErr
    reacData(ll)
  })
)
observeEvent(
  # Update integrator absolute error
  input$logAbsErr, ({
    ll = reacData()
    ll$absoluteError = 10^input$logAbsErr
    reacData(ll)
  })
)
# observeEvent(
#   # Empty MC_Output dir
#   input$mcClean, ({
#     showModal(modalDialog(
#       title = ">>>> WARNING <<<< ",
#       paste0('You are about to remove all results in MC_Output'),
#       easyClose = FALSE,
#       footer = tagList(
#         modalButton("Cancel"),
#         actionButton("mcCleanOK", "OK")
#       ),
#       size = 's'
#     ))
#   })
# )
# observeEvent(input$mcCleanOK, {
#   removeModal()
#   allFiles = list.files(
#     path = paste0(projectDir(),'/MC_Output'),
#     pattern = '.dat',
#     full.names = TRUE
#   )
#   file.remove(allFiles)
# })

# Local ####
running = reactiveVal(NULL)
observeEvent(
  input$reactorRun, {

    # Generate updated control.dat
    generateUpdatedControl()

    # Clean standard outputs
    stdout = file.path(projectDir(),'Run','runOut.txt')
    if(file.exists(stdout))
      file.remove(stdout)
    stderr = file.path(projectDir(),'Run','runErr.txt')
    if(file.exists(stderr))
      file.remove(stderr)

    # Nb MC runs
    nMC = as.numeric(input$nMCRun)

    if (nMC == 0) {
      # Nominal run
      running(
        try(
          system2(
            command = file.path(projectDir(),
                                'Scripts',
                                'OneRun_Loc.sh'),
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
          path = file.path(projectDir(),'MC_Output'),
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
        # Clean existing runs to avoid mixing with older runs
        allFiles = list.files(
          path = file.path(projectDir(),'MC_Output'),
          pattern = '.dat',
          full.names = TRUE
        )
        file.remove(allFiles)
      }

      running(
        try(
          system2(
            command = file.path(projectDir(),'Scripts','MCRun_Loc1.sh'),
            args    = c(first,nrun,projectDir()),
            stdout  = stdout,
            stderr  = stderr,
            wait    = FALSE ),
          silent = TRUE
        )
      )

      # # Tried to use progressbar by managing the loop internally,
      # # PB: have to put "wait = TRUE" to avoid process mingling
      # #     and the output of the runs is not displayed in RT
      # withProgress(
      #   message = 'Reactor',
      #   {
      #     for (iMC in first:(first+nrun)) {
      #       incProgress(
      #         1/nrun,
      #         detail = paste('Run #', iMC, '/', nrun))
      #       running(
      #         try(
      #           system2(
      #             command = paste0(projectDir(),'/Scripts/OneRun_Loc.sh'),
      #             args    = c(iMC,projectDir()),
      #             stdout  = stdout,
      #             stderr  = stderr,
      #             wait    = TRUE ), # PB: does not display code output
      #           silent = TRUE
      #         )
      #       )
      #     }
      #   }
      # )
    }
  }
)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# PB: reactiveFileReader works hapazardly in container
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
stdErr = reactiveFileReader(
  intervalMillis = 500,
  session  = session,
  filePath = file.path(ctrlPars$projectDir,'Run','runErr.txt'),
  readFunc = readLines
)
output$reactorErrors <- renderPrint({
  if (is.null(running())) {
    cat('Please run code ...')
    return(base::invisible())
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
  intervalMillis = 500,
  session  = session,
  filePath = file.path(ctrlPars$projectDir,'Run','runOut.txt'),
  readFunc = readLines
)
output$reactorOutput <- renderPrint({
  if (is.null(running())) {
    cat('Please run code ...')
    return(base::invisible())
  }

  if(input$logTail)
    cat(tail(stdOut(), n = 23), sep = '\n')
  else
    cat(stdOut(), sep = '\n'
    )

})

# Cloud ####
#
# RQ: requires a solid process and errors management
# (errors in the nova scripts might generate infinite
# loops of error messages)
#
# output$cloudSelect <- renderUI({
#
#   # Max runs to max MC data samples
#   InputMCDir = paste0(projectDir(),'/MC_Input/Reactions')
#   maxMC = length(
#     list.files(path = InputMCDir, pattern = '.csv'))
#   maxMC = maxMC - 1 # Do not count nominal run
#
#   # Rudimentary dispatch strategy
#   nMC = as.numeric(input$nMCRun)
#   if(nMC <= 10) {
#     CL_SIZE = max(nMC,1) # Size of cluster
#     NB_CORE = 1          # Nb cores on VM
#     RUN_BY_CORE = 1
#   } else if(nMC > 10 & nMC <= 100) {
#     CL_SIZE = 10
#     NB_CORE = nMC / CL_SIZE
#     RUN_BY_CORE = 1
#   } else{
#     CL_SIZE = 10    # Size of cluster
#     NB_CORE = 10    # Nb cores on VM
#     RUN_BY_CORE = nMC/(CL_SIZE*NB_CORE)
#   }
#
#   ui <- tagList(
#     fluidRow(
#       column(
#         width = 6,
#         numericInput(
#           'clSize',
#           label = 'Cluster size',
#           value = CL_SIZE,
#           min   =  1,
#           max   = 10,
#           width = '200px'
#         ),
#         numericInput(
#           'nbCore',
#           label = 'Cores per VM',
#           value = NB_CORE,
#           min   = 1,
#           max   = 20,
#           width = '200px'
#         )
#       ),
#       column(
#         width = 6,
#         numericInput(
#           'runByCore',
#           label = '# runs by core',
#           value = RUN_BY_CORE,
#           min   =  1,
#           max   = 10,
#           width = '200px'
#         )
#       )
#     ),
#     h3('')
#   )
#   ui
# })
# observeEvent(
#   input$createCluster, {
#
#     # Generate updated control.dat
#     generateUpdatedControl()
#
#     # Clean standard outputs
#     stdout = paste0(projectDir(), '/Run/runOut.txt')
#     if (file.exists(stdout))
#       file.remove(stdout)
#     stderr = paste0(projectDir(), '/Run/runErr.txt')
#     if (file.exists(stderr))
#       file.remove(stderr)
#
#     # Generate cluster config file used by all c-scripts
#     running('')
#     config = paste0(
#       'export CL_SIZE=',input$clSize,'\n',
#       'export NB_CORE=',input$nbCore,'\n',
#       'export MC_RUNS=',max(input$nMCRun,1),'\n',
#       'export VM_KEY=PPkey\n',
#       'export VM_FLAVOR=os.',input$nbCore,'\n',
#       'export VM_IMAGE=titan_2018-07-04'
#     )
#     writeLines(paste0('>>> Generate cl_config.rc <<<\n',config), stdout)
#     writeLines('',stderr)
#     oldCtrl = paste0(projectDir(), '/Run/cl_config.rc')
#     if (file.exists(oldCtrl)) {
#       savCtrl = paste0(projectDir(), '/Run/cl_config.rc_sav')
#       file.copy(from = oldCtrl,
#                 to = savCtrl,
#                 overwrite = TRUE)
#     }
#     writeLines(config, oldCtrl)
#
#     # Create cluster
#     writeLines('Creating cluster', stdout)
#     running(try(system2(
#       command = paste0(projectDir(), '/../../Scripts/claunch.sh'),
#       args    = projectDir(),
#       stdout  = stdout,
#       stderr  = stderr,
#       wait    = FALSE
#     ),
#     silent = TRUE))
#
#   }
# )
