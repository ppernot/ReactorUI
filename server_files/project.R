shinyFiles::shinyDirChoose(
  input,
  'projectDir',
  roots = roots
)

output$selectMsg <- renderPrint({
  if (!"path" %in% names(input$projectDir)) {
    cat('Select a project directory\n')

  } else {

    # Get selection
    ctrlPars[['projectDir']] <<- shinyFiles::parseDirPath(
      roots = roots,
      input$projectDir
    )

    # Save selection to file
    if(!is.null(ctrlPars$projectDir) &
       is.character(ctrlPars$projectDir))
      rlist::list.save(ctrlPars,'ctrlParams.yaml')

    cat('Selected: ', ctrlPars$projectDir)

  }

})

output$contentsMsg <- renderPrint({
  if (!"path" %in% names(input$projectDir)) {
    cat('')

  } else {
    # Count MC outputs
    mcf = list.files(
      path = paste0(ctrlPars$projectDir,'/MC_Output'),
      pattern = 'fracmol_'
    )
    nmc = length(mcf)

    if (nmc==0)
      cat('The project has not been run yet !')
    else
      cat('Contents: ', nmc,' MC runs')

  }

})
