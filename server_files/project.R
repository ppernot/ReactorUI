shinyFiles::shinyDirChoose(
  input,
  'projectDir',
  roots = roots
)

# projectDir = reactiveVal(ctrlPars$projectDir) # Restart ???
projectDir   = reactiveVal(NULL)
reacData     = reactiveVal(NULL)
chemDBData   = reactiveVal(NULL)
spectrumData = reactiveVal(NULL)

observeEvent(
  input$projectDir, {

  dir = shinyFiles::parseDirPath(
    roots = roots,
    input$projectDir)

  projectDir(dir)

  # Save status to file
  ctrlPars[['projectDir']] <<- dir
  if (!is.null(ctrlPars$projectDir) &
      is.character(ctrlPars$projectDir))
    rlist::list.save(ctrlPars, 'ctrlParams.yaml')

  # Attempt reading "control.dat"
  ctrlList = readCtrl(dir)

  if(is.null(ctrlList)) {
    # Use default data
    REAC_DATA = REAC_DATA_default
    DN_DATA   = DB_DATA_default
    # writeCtrl(REAC_DATA,DB_DATA) # Generate control.dat

  } else {

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
  reacData(ctrlList$REAC_DATA)
  chemDBData(ctrlList$DB_DATA)
  spectrumData(NULL) # Reinit

})
output$selectMsg <- renderPrint({
  if (is.null(projectDir())) {
    cat('Select a project directory\n')

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
      samp = paste0(projectDir(),'/MC_Output/',mcf[1])
      ctrl = paste0(projectDir(),'/Run/control.dat')
      if( file.mtime(samp) < file.mtime(ctrl) )
        cat('WARNING : MC runs older than control.dat. \n',
            '          Might be outdated...')
    }

  }

})
