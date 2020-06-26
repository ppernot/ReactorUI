# Function  ####
## Chemistry
generateNetwork <- function(spInit) {

  # Locations of reactions database
  reacSourceDirs = c('/../../ChemDBPublic/Neutrals/',
                     '/../../ChemDBPublic/Ions/')
  photoSourceDir = '/../../ChemDBPublic/PhotoProcs/'

  # Kinetic Parser
  nbReac=0
  reactants = products = params = type = orig = locnum = list()

  ## Photo processes
  photoData = c(
    paste0(ctrlPars$projectDir,photoSourceDir,'PhotoScheme.dat'),
    paste0(ctrlPars$projectDir,photoSourceDir,'PhotoIonScheme.dat')
  )
  for (filename in photoData) {
    if(file.exists(filename)){ # PhotoIonScheme.dat is optional
      scheme  = read.fwf(file=filename, widths= rep(11,12))
      scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
      for (i in 1:nrow(scheme)) {
        nbReac = nbReac + 1
        terms=scheme[i,1:2]
        reactants[[nbReac]] = terms[!is.na(terms) & terms!="" & terms!="HV"]
        terms=scheme[i,3:6]
        products[[nbReac]]  = terms[!is.na(terms) & terms!="" & terms!="HV"]
        terms=scheme[i,7:12]
        params[[nbReac]]    = terms[!is.na(terms) & terms!=""]
        type[[nbReac]]      = 'photo'
        locnum[[nbReac]]    = i
        orig[[nbReac]]      = filename
      }
    }
  }

  ## Reactions
  for (sourceDir in reacSourceDirs) {
    filename=paste0(ctrlPars$projectDir,sourceDir,'run_0000.csv')
    scheme  = read.csv(file=filename,header=FALSE,sep=';')
    scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
    for (i in 1:nrow(scheme)) {
      nbReac = nbReac + 1
      terms=scheme[i,1:3]
      reactants[[nbReac]] = terms[!is.na(terms) & terms!=""]
      terms=scheme[i,4:7]
      products[[nbReac]]  = terms[!is.na(terms) & terms!="" & terms!="HV"]
      terms=scheme[i,8:ncol(scheme)]
      params[[nbReac]]    = terms[!is.na(terms) & terms!=""]
      type[[nbReac]]      = terms[length(terms)]
      locnum[[nbReac]]    = i
      orig[[nbReac]]      = sourceDir
    }
  }

  # Build species list from reactants and products
  species   = levels(as.factor(unlist(c(reactants,products))))
  nbSpecies = length(species)

  # Elements in initial species
  elInit = names(CHNOSZ::count.elements(paste0(spInit,collapse='')))

  # Stoechiometry
  compo = t(apply(as.matrix(species,ncol=1),1,get.atoms))
  colnames(compo)=elements
  mass  = apply(compo,1,massFormula)
  names(mass) = species
  ## Attribute mass to dummy species
  dummySpecies = spDummy # From global variables
  dummyMass   = round(max(mass,na.rm=TRUE)+2)
  mass[dummySpecies] = dummyMass
  # orderByMass=order(mass)
  # species1=species[orderByMass]
  # mass1=mass[orderByMass]

  # Full R, L, D matrices
  L = R = D = matrix(0,ncol=nbSpecies,nrow=nbReac)
  for (m in 1:nbReac) {
    reac = unlist(reactants[m])
    prod = unlist(products[m] )
    for (n in 1:nbSpecies) {
      search=species[n]
      L[m,n] = length(which( search == reac )) # Loss
      R[m,n] = length(which( search == prod )) # Prod
    }
  }
  D = R - L # Step change matrix
  colnames(L)=species
  colnames(R)=species
  colnames(D)=species

  # Volpert analysis
  # Build species list and reactions list
  # from initial species list by
  # iterations over the reaction network
  reacs    = 1:nbReac
  spInit0  = spInit
  lReacs   =
    rowSums(
      as.matrix(
        L[,spInit] == 1,
        ncol = length(spInit)
        )
      ) & type=='photo'
  reacList = reacs[lReacs]
  spProds  = species[colSums(R[lReacs,]) != 0 ]
  spProds  = spProds[!(spProds %in% spInit)]

  vlpInd  = rep(NA,length(species))
  names(vlpInd)   = species
  vlpInd[spInit]  = 0
  vlpInd[spProds] = 1
  ncount = 1
  while ( length(spProds)!=0 ) {
    spInit   = c(spInit,spProds)
    lReacsU  = rowSums(L[,spInit] == 1) & type=='photo'
    lReacsB  = rowSums(L[,! colnames(L) %in% spInit]) == 0 & type!='photo'
    lReacs   = lReacsU | lReacsB
    reacList = c(reacList,reacs[lReacs])
    spProds  = species[colSums(R[lReacs,]) != 0]
    spProds  = spProds[!(spProds %in% spInit)]
    ncount   = ncount + 1
    vlpInd[spProds] = ncount
  }
  # Orphan species: species without parents in scheme
  ## For info only. At the moment, they are not rejected from scheme
  # reject = apply(compo, 1, function(x)
  #   sum(x[elInit],na.rm=TRUE) != sum(x,na.rm=TRUE) )
  # orphans = species[!species %in% spInit &
  #                     !reject &
  #                     !species %in% dummySpecies ]
  # Reduce reactions and species list to accessible species
  species   = spInit
  nbSpecies = length(species)
  mass      = mass[species]
  vlpInd    = vlpInd[species]

  reacList  = sort(unique(reacList))
  nbReac    = length(reacList)
  L         = L[reacList,species]
  R         = R[reacList,species]
  params    = params[reacList]
  type      = type[reacList]
  locnum    = locnum[reacList]
  orig      = orig[reacList]

  # Build links matrix for network plots
  linksR=matrix(0,ncol = nbSpecies, nrow = nbSpecies)
  colnames(linksR) = species
  rownames(linksR) = species
  for (i in 1:nbReac) {
    reacts = which(L[i,] != 0)
    prodts = which(R[i,] != 0)
    linksR[reacts, prodts] = 1
    # linksR[prodts, reacts] = 1
  }

  return(
    list(
      spInit  = spInit0,
      species = species,
      reacs   = reacList,
      L       = L,
      R       = R,
      linksR  = linksR,
      params  = params,
      type    = type,
      locnum  = locnum,
      orig    = orig,
      mass    = mass,
      vlpInd  = vlpInd
    )
  )
}

# Interactive ####

reacData <- reactiveVal()
reacProc <- reactiveValues()

output$contentsNmlMsg <- renderPrint({
  if (is.null(input$projectDir)) {
    return(NULL)
  }
  ctrlList = readCtrl(ctrlPars$projectDir)
  # Update reactive value
  reacData(ctrlList)
})
output$reactorParams <- renderPrint({
  if (is.null(reacData())) {
    return(NULL)
  }

  listPars <- c(
    "reactorLength",
    "reactorSection",
    "gasTemperature",
    "electronsTemperature",
    "totalPressure",
    "reactantsPressure",
    "reactantsFlux",
    "reactionTime"
  )

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))

  for (n in listPars)
    cat(n, " = ", unlist(mget(n, ifnotfound = NA)), "\n")
})

# Load spectrum file
spectrumData <- reactive({
  if (is.null(reacData())) {
    return(NULL)
  }

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))

  # Get spectrum data
  sp <- read.csv(
    file = paste0(ctrlPars$projectDir, "/Run/", beamSpectrumFile),
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE
  )

  list(
    wavelength = sp[, 1],
    photonFlux = sp[, 2]
  )
})
output$irradUI <- renderUI({
  if (is.null(spectrumData())) {
    return(NULL)
  }

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))
  for (n in names(spectrumData()))
    assign(n, rlist::list.extract(spectrumData(), n))

  tagList(
    fileInput(
      "beamSpectrumFileName",
      label = "Beam spectrum file",
      placeholder = paste0(ctrlPars$projectDir, "/Run/", beamSpectrumFile)
    ),
    sliderInput(
      "spectrumRange",
      "Spectrum Range (nm)",
      min = min(wavelength), max = max(wavelength),
      value = range(wavelength),
      step = 10
    ),
    numericInput(
      "beamIntensity",
      "Beam Intensity (ph.cm^-2.s^-1)",
      value = beamIntensity,
      width = "75%",
      step  = beamIntensity/10
    ),
    numericInput(
      "beamSection",
      "Beam Section (cm^2)",
      value = beamSection,
      min   = 0,
      step  = 0.01,
      width = "75%"
    )
  )
})
output$irradParams <- renderPrint({
  if (is.null(spectrumData())) {
    return(NULL)
  }

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))
  for (n in names(spectrumData()))
    assign(n, rlist::list.extract(spectrumData(), n))

  # Select in range
  selw <- wavelength >= input$spectrumRange[1] &
          wavelength <= input$spectrumRange[2]

  # Rescale if adequate
  if (input$beamIntensity >= 0) {
    photonFlux <- photonFlux / sum(photonFlux[selw]) *
                  input$beamIntensity
  }
  # Spread beam over full reactor cross-section
  photonFlux <- photonFlux * (input$beamSection / reactorSection)

  cat("\n")
  cat("Statistics:\n")
  cat("Photons **************************\n")
  cat(
    "flux  / ph.cm^-2.s^-1 =",
    signif(sum(photonFlux), 2), "\n"
  )
  cat(
    "flux  / ph.s^-1       =",
    signif(sum(photonFlux) * reactorSection, 2), "\n"
  )
  cat(
    "integ / ph            =",
    signif(sum(photonFlux) * reactorSection * reactionTime, 2), "\n"
  )
  cat("**********************************")

})
output$irradSpectrum <- renderPlot({
  if (is.null(spectrumData())) {
    return(NULL)
  }

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))
  for (n in names(spectrumData()))
    assign(n, rlist::list.extract(spectrumData(), n))
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  selw <- wavelength >= input$spectrumRange[1] &
          wavelength <= input$spectrumRange[2]

  # Rescale if adequate
  if (input$beamIntensity >= 0) {
    photonFlux <- photonFlux /
      sum(photonFlux[selw]) * input$beamIntensity
  }
  # Spread beam over full reactor cross-section
  photonFlux <- photonFlux * (input$beamSection / reactorSection)

  par(
    mfrow = c(1, 1),
    cex = cex, cex.main = cex, mar = mar,
    mgp = mgp, tcl = tcl, pty = pty, lwd = lwd
  )
  plot(
    wavelength[selw], photonFlux[selw],
    type = "l", col = cols[5],
    xlab = "Wavelength / nm",
    ylab = "Intensity / ph.cm^-2.s^-1.nm^-1",
    main = "Photon flux"
  )
  grid()
  box()
})

output$chemistryParams <- renderUI({
  if (is.null(reacData())) {
    return(NULL)
  }

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))

  csp = reactantsComposition
  sel = which(is.finite(csp))
  csp = csp[sel]/sum(csp[sel])
  rsp = reactantsSpecies
  rsp = strsplit(rsp,',')[[1]][sel]

  width = 6

  ui <- list(
    fluidRow(
      column(
        width,
        h5("Species")
      ),
      column(
        width,
        h5("Compo.")
      )
    ),
    hr(style = "border-color: #666;")
  )
  for (i in 1:length(rsp)) {
    ui[[2+i]] <-
      fixedRow(
        column(
          width,
          textInput(
            paste0("sp_", i),
            label = NULL,
            value = rsp[i]
          )
        ),
        column(
          width,
          numericInput(
            paste0("c0_", i),
            label = NULL,
            value = csp[i],
            min   = 0,
            max   = 1,
            step  = 0.1
          )
        )
      )
  }
  # Additional empty lines
  for (i in 1:2) {
    ui[[2+length(rsp)+i]] <-
      fixedRow(
        column(
          width,
          textInput(
            paste0("sp_", length(rsp)+i),
            label = NULL,
            value = NA,
            placeholder = "New sp."
          )
        ),
        column(
          width,
          numericInput(
            paste0("c0_", length(rsp)+i),
            label = NULL,
            value = 0,
            min   = 0,
            max   = 1,
            step  = 0.1
          )
        )
      )
  }
  ui
})

reacScheme <- reactiveVal()
observeEvent(
  input$generateNetwork, {
    spInit = c()
    for(sp in 1:5)
      spInit = c(spInit, input[[paste0('sp_',sp)]])
    spInit = spInit[!is.na(spInit) & spInit != ""]
    future({generateNetwork(spInit)}) %...>% reacScheme()
  }
)
output$summaryScheme <- renderPrint({
  if (is.null(reacScheme())) {
    return(cat('Please Generate Reactions...'))
  }

  cat('Summary\n')
  cat('-------\n\n')
  cat('Nb species         = ', length(reacScheme()$species),'\n')
  cat('Nb reactions       = ', length(reacScheme()$reacs)  ,'\n\n')

  vlpInd =reacScheme()$vlpInd
  maxVlpInd = max(vlpInd)
  cat('Max. Volpert Index = ', maxVlpInd,'\n\n')
  for (i in 0:maxVlpInd)
    cat('VlpI = ',i,' / Species : ',names(vlpInd[vlpInd == i]),'\n\n')

})
output$listScheme <- renderPrint({
  if (is.null(reacScheme())) {
    return(cat('Please Generate Reactions...'))
  }

  reacs   = reacScheme()$reacs
  species = reacScheme()$species
  params  = reacScheme()$params
  L       = reacScheme()$L
  R       = reacScheme()$R

  formatReac <- function(i) {
    # Arrange reactants list
    sel = L[i,]!=0
    react = colnames(L)[sel]
    stoec = L[i,][sel]
    for (j in 1:length(react))
      if(stoec[j] != 1)
        react[j] = paste(rep(react[j],stoec[j]),collapse = ' + ')
    if(reacScheme()$type[[i]] == 'photo')
      react = c(react,'hv')
    react = paste(react, collapse = ' + ')
    # Arrange products list
    sel = R[i,]!=0
    prods = colnames(R)[sel]
    stoec = R[i,][sel]
    for (j in 1:length(prods))
      prods[j] = paste(rep(prods[j],stoec[j]),collapse = ' + ')
    prods = paste(prods, collapse = ' + ')
    # Print reaction
    cat(react,' --> ',prods,unlist(params[[i]]),'\n')
  }

  if(is.na(input$targetSpecies) | input$targetSpecies == "") {
    # Full reaction list
    cat('Summary\n')
    cat('-------\n\n')
    cat('Nb species         = ', length(species),'\n')
    cat('Nb reactions       = ', length(reacs)  ,'\n\n')

    for (i in 1: length(reacs))
      formatReac(i)

  } else {
    # Species-specific reaction list
    if(! input$targetSpecies %in% species )
      return(cat('Species not in list : ',input$targetSpecies))

    indx = which(species == input$targetSpecies)
    # Reactions
    cat('Loss\n')
    cat('----\n')
    sel = which(L[,indx]!=0)
    for (i in sel)
      formatReac(i)
    cat('\n')
    # Productions
    cat('Production\n')
    cat('----------\n')
    sel = which(R[,indx]!=0)
    for (i in sel)
      formatReac(i)
    cat('\n')
  }

})

output$tabScheme <- renderDataTable({
  if (is.null(reacScheme())) {
    return(cat('Please Generate Reactions...'))
  }

  reacs   = reacScheme()$reacs
  species = reacScheme()$species
  params  = reacScheme()$params
  type    = reacScheme()$type
  L       = reacScheme()$L
  R       = reacScheme()$R

  formatReac <- function(i) {
    # Arrange reactants list
    sel = L[i,]!=0
    react = colnames(L)[sel]
    stoec = L[i,][sel]
    for (j in 1:length(react))
      if(stoec[j] != 1)
        react[j] = paste(rep(react[j],stoec[j]),collapse = ' + ')
    if(reacScheme()$type[[i]] == 'photo')
      react = c(react,'hv')
    react = paste(react, collapse = ' + ')
    # Arrange products list
    sel = R[i,]!=0
    prods = colnames(R)[sel]
    stoec = R[i,][sel]
    for (j in 1:length(prods))
      prods[j] = paste(rep(prods[j],stoec[j]),collapse = ' + ')
    prods = paste(prods, collapse = ' + ')
    # Print reaction
    typ  = unlist(type[[i]])
    pars = unlist(params[[i]])
    if(typ != 'photo')
      pars = pars[-length(pars)]

    return(
      data.frame(
        Id        = i,
        Reactants = react,
        Products  = prods,
        Params    = paste0(pars, collapse =', '),
        Type      = typ
      )
    )
  }

  dat = data.frame(
    Id = NA,
    Reactants = NA,
    Products = NA,
    Params = NA,
    Type = NA
  )
  if(is.na(input$targetSpecies) | input$targetSpecies == "") {
    for (i in 1: length(reacs))
      dat = rbind(dat,formatReac(i))

  } else {
    # Species-specific reaction list
    if(! input$targetSpecies %in% species )
      return(cat('Species not in list : ',input$targetSpecies))
    indx = which(species == input$targetSpecies)
    # Losses
    sel = which(L[,indx]!=0)
    for (i in sel)
      dat = rbind(dat,formatReac(i))
    # Productions
    sel = which(R[,indx]!=0)
    for (i in sel)
      dat = rbind(dat,formatReac(i))
  }

  return(dat[-1,])
},
rownames = FALSE,
extensions = c('Buttons','Scroller'),
options = list(
  dom = 'Btip',
  buttons =
    list('copy', list(
      extend = 'collection',
      buttons = c('csv', 'excel', 'pdf'),
      text = 'Download')
    ),
  deferRender = TRUE,
  scrollY = 550,
  scroller = TRUE
)
)

output$plotScheme <- renderForceNetwork({
  if (is.null(reacScheme())) {
    return(NULL)
  }

  for (n in names(reacScheme()))
    assign(n, rlist::list.extract(reacScheme(), n))
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  # Select max Volpert index to display
  sel = vlpInd <= input$vlpMax
  linksR = linksR[sel,sel]

  g = simplify(
    graph_from_adjacency_matrix(
      linksR,
      mode = "undirected",
      weighted = TRUE
    )
  )

  graph_d3 <- igraph_to_networkD3(g, group = vlpInd[sel])

  forceNetwork(
    Links   = graph_d3$links,
    Nodes   = graph_d3$nodes,
    Source  = 'source',
    Target  = 'target',
    NodeID  = 'name',
    Group   = 'group',
    bounded = TRUE,
    zoom    = TRUE,
    arrows  = FALSE,
    opacity = 0.8,
    colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"),
    # clickAction = 'alert(d.name );'
    linkColour = "#DDDDDD",
    legend = TRUE,
    opacityNoHover = 0.5,
    charge = input$forceNetCharge
  )

})

# Bottom ####
