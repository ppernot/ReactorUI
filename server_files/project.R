# Global reactive values
projectDir   = reactiveVal(NULL)
chemDBDir    = reactiveVal(NULL)
reacData     = reactiveVal(NULL) # Content of namelist REAC_DATA
chemDBData   = reactiveVal(NULL) # Content of namelist DB_DATA
spectrumData = reactiveVal(NULL)
reacScheme   = reactiveVal(NULL)
concList     = reactiveVal(NULL)
ratesList    = reactiveVal(NULL)
fluxesList   = reactiveVal(NULL)
graphsList   = reactiveVal(NULL)

# Open ####
shinyFiles::shinyDirChoose(
  input,
  'projectDir',
  roots = roots
)
observeEvent(
  ignoreInit = TRUE,
  input$projectDir, {

  dir = shinyFiles::parseDirPath(
    roots = roots,
    input$projectDir)

  # dir = path.expand(dir)
  req(dir)
  projectDir(dir)
  # Prevent msg from irradUI when switching projects
  reacData(NULL)

  if(input$newProj) {

    # Check that new project dir is in 'Projects' folder
    if( !dir.exists(file.path(dir,'..','..','Projects')) ) {
      showModal(modalDialog(
        title = ">>>> Wrong projects folder <<<< ",
        paste0('The selected project (',dir,
               ') is not in the "Projects" folder.\n',
               'Please select another one.'),
        easyClose = FALSE,
        footer = modalButton("OK"),
        size = 's'
      ))
      projectDir(NULL)
      REAC_DATA = NULL
      DB_DATA   = NULL

    } else {
      # Populate project
      for (dname in c(
        'Run',
        file.path('Run', 'Photo'),
        'MC_Output',
        'MC_Input',
        file.path('MC_Input', 'Photoprocs'),
        file.path('MC_Input', 'Reactions'),
        'Scripts',
        'Outputs'
      ))
        dir.create(file.path(dir, dname), showWarnings = FALSE)

      from = file.path(dir, '..', '..', 'ReactorCodes', 'reactor')
      to   = file.path(dir, 'Run', 'reactor')
      file.copy(from, to, overwrite = TRUE)

      for (script in c(
        'OneRun_Loc.sh'  ,
        'MCRun_Loc1.sh'  ,
        'OneRun_Cloud.sh',
        'MCRun_Cloud.sh'
      )) {
        from = file.path(dir, '..', '..', 'ReactorCodes', script)
        to   = file.path(dir, 'Scripts', script)
        file.copy(from, to, overwrite = TRUE)
      }

      DB_DATA   = DB_DATA_default
      REAC_DATA = REAC_DATA_default
      if(input$newProjTemplate != 'APSIS')
        REAC_DATA = REAC_DATA_Titan

      # Copy spectrum file
      from =  file.path(
        dir,'..','..','ChemDBPublic','BeamSpectrumFiles','1nm',
        REAC_DATA$beamSpectrumFile
      )
      to   = file.path(
        dir,'Run',
        REAC_DATA$beamSpectrumFile)
      file.copy(from, to)
    }

  } else {

    # Attempt at reading "control.dat"
    ctrlList = readCtrl(dir)

    if(is.null(ctrlList$REAC_DATA))
      REAC_DATA = REAC_DATA_default
    else
      REAC_DATA = ctrlList$REAC_DATA

    if(is.null(ctrlList$DB_DATA))
      DB_DATA = DB_DATA_default
    else
      DB_DATA = ctrlList$DB_DATA
  }

  # Populate/Update reactive values
  reacData(REAC_DATA)
  chemDBData(DB_DATA)
  spectrumData(NULL) # Reinit

  # Save status to file
  if(!is.null(projectDir())) {
    ctrlPars[['projectDir']] <<- dir
    if (!is.null(ctrlPars$projectDir) &
        is.character(ctrlPars$projectDir))
      rlist::list.save(ctrlPars, 'ctrlParams.yaml')
  }

  # Set ChemDBDir
  dir = file.path(projectDir(),'..','..','ChemDBPublic')
  if(!dir.exists(dir)) # When run in container, dirs are at the top...
    dir = '/ChemDBPublic'
  chemDBDir(dir)

})
# Save ####

# Outputs ####
output$selectMsg <- renderPrint({
  if (is.null(projectDir())) {
    cat('Select/Create a project directory\n')

  } else {
    cat('Selected: ', projectDir())

  }

})

output$contentsMsg <- renderPrint({
  if (is.null(projectDir())) {
    cat('')

  } else {
    # Count MC outputs
    mcf = list.files(
      path = paste0(projectDir(),'/MC_Output'),
      pattern = 'fracmol_'
    )
    nmc = length(mcf)

    if (nmc == 0) {
      cat('The project has not been run yet !')

    } else {
      cat('Contents: ', nmc,' MC runs\n')

      # Consistency checks (data MC runs
      # posterior to control.dat ...)
      fmtime = sapply(
        file.path(projectDir(),'MC_Output',mcf),
        file.mtime
      )
      ctrl = file.path(projectDir(),'Run','control.dat')
      if( min(fmtime) < file.mtime(ctrl) )
        cat('WARNING : some MC runs are older than ',
            'control.dat. \n',
            '          They might be outdated...')
    }

  }

})
observeEvent(
  input$saveProj, {
    generateUpdatedControl()
})
