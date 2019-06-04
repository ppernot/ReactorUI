# Functions ####
calcFluxes = function(C,R) {

  # Extract conc et al.
  for (n in names(C))
    assign(n,rlist::list.extract(C,n))
  # Extract rates et al.
  for (n in names(R))
    assign(n,rlist::list.extract(R,n))

  # Concentrations at stationnary state (last snapshot)
  conc = conc[,nt,]

  # MC fluxes
  allRates   = cbind(photoRates,rates)
  reacNames  = colnames(allRates)
  nbReacs    = ncol(allRates)
  nf         = nrow(allRates)
  flux       = matrix(NA,nrow=nf,ncol=nbReacs)
  nbSpecies  = nsp

  ## get D and L matrices
  D = L = matrix(0,ncol=nbSpecies,nrow=nbReacs)
  DL  = readLines(
    con = paste0(ctrlPars$projectDir,'/Run/photo_DL.dat'),
    n   = 4
  )
  nnz = as.numeric(DL[1])
  Dsp = as.numeric(unlist(strsplit(DL[2],split=' ')))
  Dsp = matrix(Dsp,nrow=nnz, ncol=3)
  for (i in 1:nrow(Dsp))
    D[Dsp[i,1],Dsp[i,2]]=Dsp[i,3]
  nnz = as.numeric(DL[3])
  Lsp = as.numeric(unlist(strsplit(DL[4],split=' ')))
  Lsp = matrix(Lsp,nrow=nnz, ncol=3)
  for (i in 1:nrow(Lsp))
    L[Lsp[i,1],Lsp[i,2]]=Lsp[i,3]
  nPh = ncol(photoRates)

  DL  = readLines(
    con = paste0(ctrlPars$projectDir,'/Run/reac_DL.dat'),
    n   = 4
  )
  nnz = as.numeric(DL[1])
  Dsp = as.numeric(unlist(strsplit(DL[2],split=' ')))
  Dsp = matrix(Dsp,nrow=nnz, ncol=3)
  for (i in 1:nrow(Dsp))
    D[nPh+Dsp[i,1],Dsp[i,2]]=Dsp[i,3]
  nnz = as.numeric(DL[3])
  Lsp = as.numeric(unlist(strsplit(DL[4],split=' ')))
  Lsp = matrix(Lsp,nrow=nnz, ncol=3)
  for (i in 1:nrow(Lsp))
    L[nPh+Lsp[i,1],Lsp[i,2]]=Lsp[i,3]

  R = D + L

  for (i in 1:nf)
    flux[i,] = allRates[i,] * apply(L,1,function(x) prod(conc[i,]^x))

  logFlux = ifelse(flux==0, NA, log10(flux))

  # Summaries
  flMean  = apply(logFlux,2,function(x) 10^mean(x,na.rm=TRUE))
  flMed   = apply(logFlux,2,function(x) 10^median(x,na.rm=TRUE))
  flF     = apply(logFlux,2,function(x) 10^sd(x,na.rm=TRUE))
  flLow   = apply(logFlux,2,
                function(x) 10^quantile(x,probs = 0.025,na.rm=TRUE))
  flSup   = apply(logFlux,2,
                function(x) 10^quantile(x,probs = 0.975,na.rm=TRUE))

  names(flMean) = reacNames
  names(flMed)  = reacNames
  names(flF)    = reacNames
  names(flLow)  = reacNames
  names(flSup)  = reacNames

  return(
    list(
      flMean = flMean,
      flMed  = flMed,
      flF    = flF,
      flLow  = flLow,
      flSup  = flSup,
      L      = L,
      R      = R
    )
  )
}

# Interactive ####

fluxesList <- reactiveVal()
observeEvent(
  input$calcFlux,
  {
    if(is.null(concList())) {
      C = getConc()
      concList(C)
    }
    C = concList()

    if(is.null(ratesList())){
      R = getRates()
      ratesList(R)
    }
    R = ratesList()

    future({calcFluxes(C,R)}) %...>% fluxesList()
  }
)

output$viewFlow <- renderPlot({
  if(is.null(concList()) |
     is.null(ratesList()) |
     is.null(fluxesList()) |
     is.null(input$flSpec)  )
    return(NULL)

  # Extract data from lists
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))
  for (n in names(ratesList()))
    assign(n,rlist::list.extract(ratesList(),n))
  for (n in names(fluxesList()))
    assign(n,rlist::list.extract(fluxesList(),n))

  # Remove dummy species
  sel  = ! species %in% spDummy
  spec = species[sel]
  LR   = L[,sel]
  RR   = R[,sel]

  reacTypeNames=c("Ph","R")
  g = viewFlow(
    input$flSpec,
    LR,
    RR,
    spec,
    reacs    = names(flMean),
    reacType = c(rep(1, ncol(photoRates)),
                 rep(2, ncol(rates))),
    reacTypeNames = reacTypeNames,
    flMean = flMean,
    spInit = spInit,
    topShow = input$topShow,
    level = ifelse(input$level, 2, 1),
    showLegend = TRUE,
    curved = input$curvedArrow
  )
  layout = switch(
    input$fluxGraphAlgo,
    FR  = layout_with_fr(g),
    LGL = layout_with_lgl(g,
                          root = input$flSpec,
                          repulserad = vcount(g)^3 * input$topShow),
    Bipartite = layout_as_bipartite(g,vgap = 100 * input$topShow),
    layout_with_gem(g) # Default
  )
  plot(g, layout=layout)
  # legend(
  #   'topright',cex=1,bty='n',
  #   legend=reacTypeNames,
  #   lty=1,lwd=2,
  #   col=[1:length(reacTypeNames)]
  # )

})

output$viewBudget <- renderPrint({
  if(is.null(concList()) |
     is.null(ratesList()) |
     is.null(fluxesList()) |
     is.null(input$flSpec)  )
    return(NULL)


  # Extract data from lists
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))
  for (n in names(ratesList()))
    assign(n,rlist::list.extract(ratesList(),n))
  for (n in names(fluxesList()))
    assign(n,rlist::list.extract(fluxesList(),n))

  flSpec = input$flSpec
  if(flSpec != "") # Check validity
    validate(
      need(flSpec %in% species, 'Invalid species !')
    )

  # Remove dummy species
  sel  = ! species %in% spDummy
  spec = species[sel]
  LR   = L[,sel]
  RR   = R[,sel]

  budget(
    flSpec,
    LR,RR,spec,
    names(flMean),
    flMean,
    weightThresh=input$topShow
  )
})

output$viewTarget <- renderPrint({
  if(is.null(concList()) |
     is.null(ratesList()) |
     is.null(fluxesList()) |
     is.null(input$flSpec)  )
    return(NULL)

  # Extract data from lists
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))
  for (n in names(ratesList()))
    assign(n,rlist::list.extract(ratesList(),n))
  for (n in names(fluxesList()))
    assign(n,rlist::list.extract(fluxesList(),n))

  flSpec = input$flSpec
  if(flSpec != "") # Check validity
    validate(
      need(flSpec %in% species, 'Invalid species !')
    )

  # Remove dummy species
  sel  = ! species %in% spDummy
  spec = species[sel]
  LR   = L[,sel]
  RR   = R[,sel]

  traceBack(
    flSpec,
    LR, RR, spec, spInit,
    names(flMean),
    flMean
  )
})
