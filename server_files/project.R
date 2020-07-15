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

# Open ####
observeEvent(
  input$projectDir, {

  dir = shinyFiles::parseDirPath(
    roots = roots,
    input$projectDir)
  projectDir(dir)
  req(projectDir())

  if(input$newProj) {
    for (dname in c('/Run',
                    '/Run/Photo',
                    '/MC_Output',
                    '/MC_Input',
                    '/MC_Input/Photoprocs',
                    '/MC_Input/Reactions',
                    '/Scripts') )
      dir.create(paste0(dir,dname),showWarnings = FALSE)

    file.copy(
      from = paste0(dir,'/../../ReactorCodes/reactor'),
      to   = paste0(dir,'/Run/reactor'),
      overwrite = TRUE)
    file.copy(
      from = paste0(dir,'/../../ReactorCodes/OneRun_Loc.sh'),
      to   = paste0(dir,'/Scripts/OneRun_Loc.sh'),
      overwrite = TRUE)
    file.copy(
      from = paste0(dir,'/../../ReactorCodes/MCRun_Loc1.sh'),
      to   = paste0(dir,'/Scripts/MCRun_Loc1.sh'),
      overwrite = TRUE)
    file.copy(
      from = paste0(dir,'/../../ReactorCodes/OneRun_Cloud.sh'),
      to   = paste0(dir,'/Scripts/OneRun_Loc.sh'),
      overwrite = TRUE)
    file.copy(
      from = paste0(dir,'/../../ReactorCodes/MCRun_Cloud.sh'),
      to   = paste0(dir,'/Scripts/MCRun_Loc1.sh'),
      overwrite = TRUE)

    DB_DATA   = DB_DATA_default
    REAC_DATA = REAC_DATA_default
    if(input$newProjTemplate != 'APSIS')
      REAC_DATA = REAC_DATA_Titan

    # Copy spectrum file
    fromFile =  paste0(
      dir,
      '/../../ChemDBPublic/BeamSpectrumFiles/1nm/',
      REAC_DATA$beamSpectrumFile
    )
    toFile   = paste0(
      dir,
      '/Run/',
      REAC_DATA$beamSpectrumFile)
    file.copy(from = fromFile, to = toFile)

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
  ctrlPars[['projectDir']] <<- dir
  if (!is.null(ctrlPars$projectDir) &
      is.character(ctrlPars$projectDir))
    rlist::list.save(ctrlPars, 'ctrlParams.yaml')

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
        paste0(projectDir(),'/MC_Output/',mcf),
        file.mtime
      )
      ctrl = paste0(projectDir(),'/Run/control.dat')
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
