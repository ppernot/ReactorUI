# Function  ####
## Chemistry
generateNetwork <- function(
  spInit,
  photoSourceDir    = NULL,
  neutralsSourceDir = NULL,
  ionsSourceDir     = NULL
  ) {

  # Kinetic Parser
  nbReac=0
  reactants = products = params =
    type = orig = locnum = reacTagFull = list()

  ## Photo processes
  photoData = file.path(photoSourceDir,'..','PhotoScheme.dat')
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
        rr = paste(reactants[[nbReac]],collapse = ' + ')
        pp = paste(products[[nbReac]],collapse = ' + ')
        reacTagFull[[nbReac]]   = paste0(rr,' --> ',pp)
      }
    }
  }

  ## Reactions
  for (sourceDir in c(neutralsSourceDir, ionsSourceDir)) {
    filename = file.path(sourceDir,'run_0000.csv')
    # scheme   = read.csv(file = filename, header = FALSE, sep = ';')
    scheme   = as.data.frame(
      data.table::fread(file = filename, header = FALSE, sep = ';')
    )
    scheme   = t(apply(scheme, 1, function(x) gsub(" ", "", x)))
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
      rr = paste(reactants[[nbReac]],collapse = ' + ')
      pp = paste(products[[nbReac]],collapse = ' + ')
      reacTagFull[[nbReac]]   = paste0(rr,' --> ',pp)
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
      ) & type == 'photo'
  reacList = reacs[lReacs]

  spProds  = species[colSums(R[lReacs,]) != 0 ]
  spProds  = spProds[!(spProds %in% spInit)]

  vlpInd  = rep(NA,length(species))
  names(vlpInd)   = species
  vlpInd[spInit]  = 0
  vlpInd[spProds] = 1
  ncount = 1
  while ( length(spProds) != 0 ) {
    spInit   = c(spInit,spProds)
    lReacsU  = rowSums(L[,spInit] == 1) & type == 'photo'
    lReacsB  = rowSums(L[,! colnames(L) %in% spInit]) == 0 &
      type != 'photo'
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
  L         = L[reacList, species]
  R         = R[reacList, species]
  D         = D[reacList, species]
  params    = params[reacList]
  type      = type[reacList]
  locnum    = locnum[reacList]
  orig      = orig[reacList]
  reacTagFull = reacTagFull[reacList]

  # Short tags
  nbPhoto = sum(type == 'photo')
  reacTags    = c(
    paste0('Ph', 1:nbPhoto),
    paste0('R',  1:(nbReac-nbPhoto))
  )
  rownames(L) = reacTags
  rownames(R) = reacTags
  rownames(D) = reacTags

  return(
    list(
      spInit  = spInit0,
      species = species,
      reacs   = reacList,
      nbReac  = nbReac,
      nbSpecies = nbSpecies,
      L       = L,
      R       = R,
      D       = D,
      params  = params,
      type    = type,
      locnum  = locnum,
      orig    = orig,
      mass    = mass,
      vlpInd  = vlpInd,
      reacTagFull = reacTagFull,
      reacTags = reacTags
    )
  )
}
generateLinks = function (LR, weights = NULL) {

  # Chemical network params
  L = LR$L
  R = LR$R

  nbSpecies = ncol(R)
  species   = colnames(R)
  nbReacs   = nrow(R)
  reacTags  = rownames(R)

  # Build species connectivity  matrix for network plots
  linksR=matrix(0,ncol = nbSpecies, nrow = nbSpecies)
  colnames(linksR) = species
  rownames(linksR) = species
  for (i in 1:nbReacs) {
    reacts = which(L[i,] != 0)
    prodts = which(R[i,] != 0)
    linksR[reacts, prodts] = linksR[reacts, prodts] + 1
  }

  # Build bi-partite/digraph connectivity matrix
  nbNodes = nbReacs + nbSpecies
  lSpecies = (nbReacs + 1):(nbReacs + nbSpecies)
  nodeNames = c(reacTags,species)
  linksR2 = matrix(0, ncol = nbNodes, nrow = nbNodes)
  colnames(linksR2) = nodeNames
  rownames(linksR2) = nodeNames
  for (i in 1:nbReacs) {
    linksR2[lSpecies, i] = linksR2[lSpecies, i] + L[i, 1:nbSpecies]
    linksR2[i, lSpecies] = linksR2[i, lSpecies] + R[i, 1:nbSpecies]
  }

  return(
    list(
      linksR  = linksR,
      linksR2 = linksR2
    )
  )
}
writeDL = function(outputDir,RS) {

  # Chemical network params
  for (n in names(RS))
    assign(n, rlist::list.extract(RS, n))

  # Save species list, mass [& charge (TBD)]
  sp_aux =paste(
    nbSpecies,'\n',
    paste(species,collapse = ' '), '\n',
    paste(mass,   collapse = ' '), '\n',
    paste(rep(0,nbSpecies), collapse = ' ')
  )
  writeLines(
    sp_aux,
    file.path(outputDir,'Run','species_aux.dat')
  )

  # Split photodissociations and reactions
  photo = type == 'photo'
  Dphoto = matrix(D[ photo,],nrow=sum(photo), ncol=nbSpecies)
  D      = matrix(D[!photo,],nrow=sum(!photo),ncol=nbSpecies)
  Lphoto = matrix(L[ photo,],nrow=sum(photo), ncol=nbSpecies)
  L      = matrix(L[!photo,],nrow=sum(!photo),ncol=nbSpecies)

  # Treat Reactions ###
  # 1/ Convert D & L to triplets (sparse matrices)
  Dsp = cbind(which(D!=0,arr.ind=TRUE),D[D!=0])
  Lsp = cbind(which(L!=0,arr.ind=TRUE),L[L!=0])

  reac_DL = paste(nrow(Dsp),'\n',
                  trimws(paste(Dsp, collapse = " ")),'\n',
                  nrow(Lsp),'\n',
                  trimws(paste(Lsp, collapse = " ")))
  writeLines(
    reac_DL,
    file.path(outputDir,'Run','reac_DL.dat')
  )

  # Reactions list
  reac_list = ''
  sel = (1:nbReac)[!photo]
  for (i in sel) {
    reac_list = paste(reac_list, reacTagFull[[i]])
    if( i != sel[length(sel)])
      reac_list = paste(reac_list, '\n')
  }
  writeLines(
    reac_list,
    file.path(outputDir,'Run','reac_list.dat')
  )

  # Single output file from nominal data
  reac_params = ''
  for (i in sel) {
    reac_params = paste(
      reac_params,
      paste(params[[i]], collapse = " "))
    if( i != sel[length(sel)])
      reac_params = paste(reac_params,'\n')
  }
  writeLines(
    reac_params,
    file.path(outputDir,'Run','reac_params.dat')
  )

  # Treat Nominal Photo-processes ###

  # Reactions list
  reac_list = ''
  sel = (1:nbReac)[photo]
  for (i in sel) {
    reac_list = paste(reac_list, reacTagFull[[i]])
    if( i != sel[length(sel)])
      reac_list = paste(reac_list, '\n')
  }
  writeLines(
    reac_list,
    file.path(outputDir,'Run','photo_list.dat')
  )

  if (dim(Dphoto)[1]!=0) {
    # Convert D & L to triplets (sparse matrices)
    Dsp = cbind(which(Dphoto!=0,arr.ind=TRUE),Dphoto[Dphoto!=0])
    Lsp = cbind(which(Lphoto!=0,arr.ind=TRUE),Lphoto[Lphoto!=0])

    reac_DL = paste(nrow(Dsp),'\n',
                    trimws(paste(Dsp, collapse = " ")),'\n',
                    nrow(Lsp),'\n',
                    trimws(paste(Lsp, collapse = " ")))
    writeLines(
      reac_DL,
      file.path(outputDir,'Run','photo_DL.dat')
    )

    reac_params = ''
    for (i in sel) {
      sp=species[which(Lphoto[i,]!=0)]
      if( params[[i]][1] == 0 )
        reac_params = paste0(reac_params,'se',sp,'.dat')
      else
        reac_params = paste(
          reac_params,
          paste0('se',sp,'.dat'),
          paste0('qy',sp,'_',params[[i]][1],'.dat'))
      if( i != sel[length(sel)])
        reac_params = paste(reac_params, '\n')
    }
    writeLines(
      reac_params,
      file.path(outputDir,'Run','photo_params.dat')
    )
  }
}
# Interactive ####
output$contentsNmlMsg <- renderPrint({
  if (is.null(reacData()) |
      is.null(chemDBData())) {
    return(NULL)
  }

  # Printout
  cat('_ REAC_DATA:\n')
  for (n in names(reacData())) {
    assign(n, rlist::list.extract(reacData(), n))
    cat(n, " = ", unlist(mget(n, ifnotfound = NA)), "\n")
  }
  cat('_ DB_DATA:\n')
  for (n in names(chemDBData())) {
    assign(n, rlist::list.extract(chemDBData(), n))
    cat(n, " = ", unlist(mget(n, ifnotfound = NA)), "\n")
  }

})

# ChemDB ####
output$chemDBVersions <- renderUI({
  # Generate UI to display and select database versions
  if (is.null(chemDBData()) |
      is.null(reacData())) {
    return(NULL)
  }

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))

  for (n in names(chemDBData()))
    assign(n, rlist::list.extract(chemDBData(), n))

  allAvailable = list.dirs(path = chemDBDir(),recursive = FALSE)
  allVersions  = basename(allAvailable)

  photoDir = file.path(chemDBDir(), 'PhotoProcs_latest')
  allReso  = list.dirs(path = photoDir,recursive = FALSE)
  allReso  = basename(allReso)

  phoVers = photoVersion
  if(is.null(phoVers))
    phoVers = DB_DATA_default$photoVersion
  phoVers = ifelse(phoVers == 0,
                   'latest',
                   paste0('v_',phoVers))
  phoVers = paste0('PhotoProcs_',phoVers)

  if (exists('spectralResolution')){
    speReso = spectralResolution
    if(is.null(speReso) | !speReso  %in% c(0.1,1)) {
      speReso = REAC_DATA_default$spectralResolution
      ll = reacData()
      ll$spectralResolution = speReso
      reacData(ll)
    }
  } else {
    speReso = REAC_DATA_default$spectralResolution
    ll = reacData()
    ll$spectralResolution = speReso
    reacData(ll)
  }
  speReso = paste0(speReso,'nm')

  neuVers = neutralsVersion
  if(is.null(neuVers))
    neuVers = DB_DATA_default$neutralsVersion
  neuVers = ifelse(neuVers == 0,
                   'latest',
                   paste0('v_',neuVers))
  neuVers = paste0('Neutrals_',neuVers)

  ionVers = ionsVersion
  if(is.null(ionVers))
    ionVers = DB_DATA_default$ionsVersion
  ionVers = ifelse(ionVers == 0,
                   'latest',
                   paste0('v_',ionVers))
  ionVers = paste0('Ions_',ionVers)

  ui <- list(
    fixedRow(
      column(
        width = 3,
        selectInput(
          'phoVers',
          label    = 'PhotoProcs',
          choices  = allVersions[grepl('PhotoProcs',
                                       allVersions)],
          selected = phoVers,
          width    = '200px'
        ),
        selectInput(
          'neuVers',
          label    = 'Neutrals',
          choices  = allVersions[grepl('Neutrals',
                                       allVersions)],
          selected = neuVers,
          width    = '200px'
        ),
        selectInput(
          'ionVers',
          label    = 'Ions',
          choices  = allVersions[grepl('Ions',
                                       allVersions)],
          selected = ionVers,
          width    = '200px'
        )
      ),
      column(
        width = 3,
        selectInput(
          'speReso',
          label    = 'Spectral Resolution',
          choices  = allReso,
          selected = speReso,
          width    = '200px'
        )
      )
    )
  )
  ui
})
outputOptions(output, "chemDBVersions",suspendWhenHidden = FALSE)
output$checkChanges <- renderPrint({

  req(input$speReso)

  if(!is.null(reacData())) {

    numReso = as.numeric(sub('nm','',input$speReso))

    if(reacData()$spectralResolution != numReso) {
      ll = reacData()
      ll$spectralResolution = numReso
      reacData(ll)
    }

  }

})

# Generate Reacs ####
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
  nsp = length(rsp)

  width = 6

  # UI with 6 slots should be enough...
  nslots = 6
  if(nsp > nslots)
    showModal(modalDialog(
      title = ">>>> Too many initail species <<<< ",
      paste0('Your control data nb os species (',nsp,
             ') exceeds the capacity of the interface(',nslots,').\n',
             'If you need more, please raise an issue.'),
      easyClose = TRUE,
      footer = modalButton("OK"),
      size = 's'
    ))
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
  for (i in 1:nsp) {
    ui[[2+i]] = fixedRow(
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
  if (nsp < nslots) {
    # Additional empty lines
    for (i in (nsp+1):nslots) {
      ui[[2+i]] = fixedRow(
        column(
          width,
          textInput(
            paste0("sp_", i),
            label = NULL,
            value = NA,
            placeholder = "New sp."
          )
        ),
        column(
          width,
          numericInput(
            paste0("c0_", i),
            label = NULL,
            value = 0,
            min   = 0,
            max   = 1,
            step  = 0.1
          )
        )
      )
    }
  }
  ui
})
observeEvent(
  input$generateNetwork,
  {

    # Initial mixture

    ## Number of inputs
    nslots = 6

    ## Collect inputs
    spInit = cInit = c()
    for(sp in 1:nslots) {
      c0 = input[[paste0('c0_',sp)]]
      if(!is.finite(c0) | c0 == 0) next
      s0 = input[[paste0('sp_',sp)]]
      if(is.na(s0) | s0 == "") next
      cInit  = c(cInit,  c0)
      spInit = c(spInit, s0)
    }
    cInit = cInit / sum(cInit)

    ## Actuate parms list
    ll = reacData()
    ll$reactantsSpecies = paste0(spInit,collapse=',')
    ll$reactantsComposition = cInit
    reacData(ll)

    # Databases
    photoSourceDir    = file.path(chemDBDir(),input$phoVers,input$speReso)
    neutralsSourceDir = file.path(chemDBDir(),input$neuVers)
    ionsSourceDir     = file.path(chemDBDir(),input$ionVers)

   id = shiny::showNotification(
      h4('Generating chemistry, be patient...'),
      closeButton = FALSE,
      duration = 5,
      type = 'message'
    )
   # on.exit(removeNotification(id), add = TRUE)

    future({
      generateNetwork(
        spInit = spInit,
        photoSourceDir    = photoSourceDir,
        neutralsSourceDir = neutralsSourceDir,
        ionsSourceDir     = ionsSourceDir
      )
    }) %...>% reacScheme()
  })
observe({
  # Load reacScheme() from project if it exists
  req(projectDir())

  file = file.path(projectDir(),'Run','reacScheme.Rda')
  if(file.exists(file)) {
    load(file)     # Load RD list
    reacScheme(RS) # Populate reacScheme()
  }

})
observe({
  # Generate data files for reactor code when reacScheme() is ready
  req(projectDir())
  req(reacScheme())

  writeDL(
    outputDir = projectDir(),
    RS = reacScheme()
  )
})
observe({
  # Generate graph link matrices when reacScheme() is ready
  req(projectDir())
  req(reacScheme())

  # Save data reactiveVal
  RS = reacScheme()
  save(
    RS,
    file = file.path(projectDir(),'Run','reacScheme.Rda')
  )

  # Populate stoechList()
  stoechList(
    list(
      L = reacScheme()$L,
      R = reacScheme()$R
    )
  )

  # Populate graphList()
  graphsList(
    generateLinks(
      LR = stoechList()
    )
  )
})

output$summaryScheme <- renderPrint({
  if (is.null(reacScheme())) {
    cat('Please Generate Reactions...')
    return()
  }

  for (n in names(reacScheme()))
    assign(n, rlist::list.extract(reacScheme(), n))

  cat('Summary\n')
  cat('-------\n\n')
  cat('Nb species         = ', nbSpecies,'\n')
  cat('Nb reactions       = ', nbReac   ,'\n\n')

  # Nature of species
  types = c("neutrals","ions")
  chems = c("hydrocarbons","N-bearing","O-bearing")
  heavy = c("C0","C1","C2","C3","C4","C5","C6","Cmore")
  resu = matrix(0,
                nrow = 1 + length(chems),
                ncol = 2 + length(types))
  rownames(resu) = c(chems,'total')
  colnames(resu) = c(types,'radicals','total')
  for (type in types)
    resu['total',type] = sum(
      selectSpecies(species,c(type,chems,heavy)))
  resu['total','radicals'] = sum(
    selectSpecies(species,c(types,'radicals',chems,heavy)))

  for (chem in chems) {
    resu[chem,'total'] = sum(
      selectSpecies(species,c(types,chem,heavy)))
    resu[chem,'radicals'] = sum(
      selectSpecies(species,c(types,'radicals',chem,heavy)))
  }
  for (type in types)
    for (chem in chems)
      resu[chem,type] = sum(
        selectSpecies(species,c(type,chem,heavy)))

  cat('Compositions (excl. dummies)\n')
  cat('----------------------------\n\n')
  print(resu)
  cat('\n\n')

  # Volpert index analysis
  cat('Volpert analysis\n')
  cat('----------------\n\n')

  cat('Species distribution of Volpert index\n')
  print(table(vlpInd))
  cat('\n\n')

  maxVlpInd = max(vlpInd)
  # cat('Max. Volpert Index = ', maxVlpInd,'\n\n')
  for (i in 0:maxVlpInd)
    cat('VlpI = ',i,' / Species : ',names(vlpInd[vlpInd == i]),'\n\n\n')


  # Graph connectivity ####

  for (n in names(graphsList()))
    assign(n, rlist::list.extract(graphsList(), n))

  cat(' =============================\n',
      'Network connectivity analysis\n',
      '-----------------------------\n\n')

  g = simplify(
    graph_from_adjacency_matrix(
      linksR2,
      mode = "directed",
      weighted = TRUE
    )
  )
  nodeNames = rownames(linksR2)
  pmax = length(nodeNames)
  # cat(' Connectivity : ',igraph::vertex_connectivity(g),'\n',
  #     'Radius       : ',igraph::radius(g),'\n\n')

  deg    = degree(g)
  degOut = degree(g, mode = 'out') # Losses
  degIn  = degree(g, mode = 'in')  # Prods

  io = order(deg, decreasing = TRUE)
  cat(' Degree analysis\n',
      '---------------\n\n',
    'Node\t\t Out\t In\t Total\n')
  for (i in 1:min(pmax,20))
    cat(nodeNames[io][i], '\t\t',
        degOut[io][i],'\t',
        degIn[io][i],'\t',
        deg[io][i], '\n')

  cat('\n Species connectivity distribution\n',
         '---------------------------------\n')
  lSpecies = (nbReac + 1):(nbReac + nbSpecies)
  print(table(deg[lSpecies],
              dnn = 'Degree'
              )[1:min(length(unique(deg[lSpecies])),15)])

  cat('\n Reactions connectivity distribution\n',
         '-----------------------------------\n')
  print(table(deg[1:nbReac],
              dnn = 'Degree'
              )[1:min(length(unique(deg[1:nbReac])),15)])


})
# Sinks ####
output$quality   <- renderPrint({
  if (is.null(reacScheme())) {
    cat('Please Generate Reactions...')
    return()
  }

  for (n in names(reacScheme()))
    assign(n, rlist::list.extract(reacScheme(), n))

  cat(' Mass sorted sinks / Species with no loss:\n',
      '-----------------------------')
  cat('\n\n')
  sel = colSums(L) == 0
  if(sum(sel) != 0) {
    sp = species[sel]
    ms = mass[sel]
    io = order(ms)
    for (i in 1:length(sp)) {
      ind = io[i]
      s = sp[ind]
      prods = which(R[,s] != 0)
      cat(s,' (mass = ',round(ms[ind],digits=2),')\n  Productions :\n')
      for(j in 1:length(prods))
        cat('    ',reacTagFull[[prods[j]]],'\n')
      cat('\n')
    }
  } else {
    cat('None...')
  }
  cat('\n')
})
output$tabScheme <- renderDataTable({
  if (is.null(reacScheme())) {
    cat('Please Generate Reactions...')
    return()
  }
  id      = reacScheme()$reacTags
  # reacs   = reacScheme()$reacs
  species = reacScheme()$species
  params  = reacScheme()$params
  type    = reacScheme()$type
  L       = reacScheme()$L
  R       = reacScheme()$R

  formatReac <- function(i) {
    # Arrange reactants list
    sel = L[i, ] != 0
    react = colnames(L)[sel]
    stoec = L[i, ][sel]
    for (j in 1:length(react))
      if (stoec[j] != 1)
        react[j] = paste(rep(react[j], stoec[j]), collapse = ' + ')
    if (reacScheme()$type[[i]] == 'photo')
      react = c(react, 'hv')
    react = paste(react, collapse = ' + ')
    # Arrange products list
    sel = R[i, ] != 0
    prods = colnames(R)[sel]
    stoec = R[i, ][sel]
    for (j in 1:length(prods))
      prods[j] = paste(rep(prods[j], stoec[j]), collapse = ' + ')
    prods = paste(prods, collapse = ' + ')
    # Print reaction
    typ  = unlist(type[[i]])
    pars = unlist(params[[i]])
    if (typ != 'photo')
      pars = pars[-length(pars)]

    return(data.frame(
      Id        = id[i],
      Reactants = react,
      Products  = prods,
      Params    = paste0(pars, collapse = ', '),
      Type      = typ
    ))
  }

  dat = data.frame(
    Id = NA,
    Reactants = NA,
    Products = NA,
    Params = NA,
    Type = NA
  )
  if (is.na(input$targetSpecies) |
      input$targetSpecies == "") {
    for (i in 1:length(id))
      dat = rbind(dat, formatReac(i))

  } else {
    # Species-specific reaction list
    if (!input$targetSpecies %in% species)
      return(cat('Species not in list : ', input$targetSpecies))
    indx = which(species == input$targetSpecies)
    # Losses
    sel = which(L[, indx] != 0)
    for (i in sel)
      dat = rbind(dat, formatReac(i))
    # Productions
    sel = which(R[, indx] != 0)
    for (i in sel)
      dat = rbind(dat, formatReac(i))
  }

  return(dat[-1, ])
},
rownames = FALSE,
extensions = c('Buttons', 'Scroller'),
options = list(
  dom = 'Btip',
  buttons =
    list(
      'copy',
      list(
        extend = 'collection',
        buttons = c('csv', 'excel', 'pdf'),
        text = 'Download'
      )
    ),
  deferRender = TRUE,
  scrollY = 550,
  scroller = TRUE
))
# Network ####
output$plotScheme <- renderForceNetwork({
  if (is.null(reacScheme())) {
    cat('Please Generate Reactions...')
    return(NULL)
  }
  req(graphsList())

  for (n in names(reacScheme()))
    assign(n, rlist::list.extract(reacScheme(), n))
  for (n in names(graphsList()))
    assign(n, rlist::list.extract(graphsList(), n))
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  # Select max Volpert index to display
  selVp = vlpInd <= input$vlpMax

  if('digraph' %in% input$netCtrl) {

    linksR = linksR2
    helpNames = c(unlist(reacTagFull),species)
    nbReacs = nbReac

    if(sum(!selVp) > 0) {

      lR = matrix(R[,!selVp],ncol = sum(!selVp))
      lL = matrix(L[,!selVp],ncol = sum(!selVp))
      rsel = rowSums(lR) != 0 | rowSums(lL) != 0

      csel    = c(!rsel, selVp)
      linksR  = linksR[csel,csel]
      helpNames = helpNames[csel]

      nbReacs = sum(!rsel)
      species = species[selVp]
      nbSpecies = sum(selVp)
    }

    arrows = TRUE

  } else {

    species = species[selVp]
    linksR = linksR[selVp,selVp]

    arrows = FALSE

  }

  if(!is.matrix(linksR))
    return(NULL)

  g = simplify(
    graph_from_adjacency_matrix(
      linksR,
      mode = "directed",
      weighted = TRUE
    )
  )

  # Coloring
  compo = t(apply(as.matrix(species, ncol = 1), 1, get.atoms))
  colnames(compo) = elements
  charge = rep(0,length(species))
  ions   = grepl("\\+$",species)
  charge[ions] = 1

  charge[which(species == 'E')] = -1

  if (input$netColoring == 'volpert') {
    grp = vlpInd[selVp]
  } else if (input$netColoring == 'charge') {
    grp = charge
  } else if (input$netColoring == 'radicals') {
    radic = (apply(compo,1,numElec)-charge) %% 2
    radic[which(species == 'E')] = 1
    grp = radic
  } else if (input$netColoring == 'C/N/O') {
    azot = grepl("N",species)
    oxy  = grepl("O",species)
    grp  = rep('!N!O',length(species))
    sel  = azot & !oxy
    grp[sel] = 'N!O'
    grp[oxy] = 'O'
  }

  if('digraph' %in% input$netCtrl)
    grp = c(rep('reac',nbReacs),grp)

  graph_d3 <- igraph_to_networkD3(g, group = grp)
  graph_d3$nodes[['size']] = igraph::degree(g) * 5
  if('digraph'   %in% input$netCtrl)
    graph_d3$nodes[['size']][1:nbReacs] = 1
  if('digraph'   %in% input$netCtrl &
     'showNames' %in% input$netCtrl)
    graph_d3$nodes[['name']] = helpNames

  cs = JS("d3.scaleOrdinal(d3.schemeCategory10);")
  alert = 'alert(d.name);'

  lc = sprintf("%x",15-input$linkDensNet)
  lc = paste0('#',lc,lc,lc)

  networkD3::forceNetwork(
    Links   = graph_d3$links,
    Nodes   = graph_d3$nodes,
    Source  = 'source',
    Target  = 'target',
    Value   = 'value',
    NodeID  = 'name',
    Group   = 'group',
    Nodesize= 'size',
    bounded = FALSE,
    zoom    = TRUE,
    arrows  = arrows,
    opacity = 1.0,
    opacityNoHover = 0.4,
    colourScale = cs,
    clickAction = alert,
    linkColour = lc,
    legend   = 'legend' %in% input$netCtrl,
    charge   = input$forceNetCharge,
    fontSize = input$fontSizeNet
  )

})

# Sample ####
output$nMCButton <- renderUI({
  if (is.null(chemDBData())) {
    return(NULL)
  }

  # Get max number of samples
  photoSourceDir    = file.path(chemDBDir(),input$phoVers,input$speReso)
  neutralsSourceDir = file.path(chemDBDir(),input$neuVers)
  ionsSourceDir     = file.path(chemDBDir(),input$ionVers)

  maxNeutrals = length(
    list.files(path = neutralsSourceDir, pattern = '.csv'))
  maxIons = length(
    list.files(path = ionsSourceDir, pattern = '.csv'))
  files = list.files(path = photoSourceDir, pattern = '.dat')
  nFiles = length(files)
  maxPhoto = as.numeric(substr(files[nFiles],1,4))
  maxMC = min(maxNeutrals, maxIons, maxPhoto)

  nMC = chemDBData()$nMC
  if(is.null(nMC))
    nMC = DB_DATA_default$nMC

  ui <- list(
    selectInput(
      'nMC',
      label    = '# MC samples',
      choices  = c(10,seq(100,maxMC,by=100)),
      selected = nMC,
      width    = '200px'
    )
  )
  ui
})
observeEvent(input$sampleChem, {

    # Nb MC samples
    nMC = as.numeric(input$nMC)
    # Update chemDBData
    ll = chemDBData()
    ll$nMC = nMC
    chemDBData(ll)

    # Where to save files
    outputDir   = projectDir()
    targetMCDir = file.path(outputDir,'MC_Input')

    # Network params
    for (n in names(reacScheme()))
      assign(n, rlist::list.extract(reacScheme(), n))

    # Databases
    photoSourceDir    = file.path(chemDBDir(),input$phoVers,input$speReso)
    neutralsSourceDir = file.path(chemDBDir(),input$neuVers)
    ionsSourceDir     = file.path(chemDBDir(),input$ionVers)

    # # Generate data files for reactor code ###
    # sp_aux =paste(
    #   nbSpecies,'\n',
    #   paste(species,collapse = ' '), '\n',
    #   paste(mass,   collapse = ' '), '\n',
    #   paste(rep(0,nbSpecies), collapse = ' ')
    # )
    # writeLines(
    #   sp_aux,
    #   file.path(outputDir,'Run','species_aux.dat')
    # )

    # Split photodissociations and reactions
    photo = type == 'photo'
    Dphoto = matrix(D[ photo,],nrow=sum(photo), ncol=nbSpecies)
    D      = matrix(D[!photo,],nrow=sum(!photo),ncol=nbSpecies)
    Lphoto = matrix(L[ photo,],nrow=sum(photo), ncol=nbSpecies)
    L      = matrix(L[!photo,],nrow=sum(!photo),ncol=nbSpecies)

    # Treat Reactions ###

    if(nMC > 1) {
      # Clean target dir
      files = list.files(
        path = file.path(targetMCDir,'Reactions'),
        pattern = 'run_',
        full.names = TRUE
      )
      if(!is.null(files))
        file.remove(files)

      withProgress(message = 'Reactions', {
        sel = !photo
        origs = unique(unlist(orig[sel]))
        for (iMC in 0:nMC) {
          sampleFile = paste0('run_', sprintf('%04i', iMC), '.csv')
          data = ''
          for (io in origs) {
            sourceDir = io
            filename = file.path(sourceDir, sampleFile)
            # scheme   = read.csv(file = filename,
            #                    header = FALSE,
            #                    sep = ';')
            scheme   = as.data.frame(
              data.table::fread(file = filename,
                                header = FALSE,
                                sep = ';')
            )

            scheme  = t(apply(scheme, 1, function(x)
              gsub(" ", "", x)))
            paramsLoc = list()
            ii = 0
            for (i in 1:nrow(scheme)) {
              ii = ii + 1
              terms = scheme[i, 8:ncol(scheme)]
              paramsLoc[[ii]] = terms[!is.na(terms) & terms != '']
            }
            selOrig = which (orig == io)
            paramsLoc = paramsLoc[unlist(locnum[selOrig])]
            for (j in 1:length(paramsLoc))
              data = paste0(
                data,
                paste(paramsLoc[[j]], collapse = ' '),
                '\n')
          }
          writeLines(
            data,
            file.path(targetMCDir, 'Reactions', sampleFile)
          )

          incProgress(1/nMC, detail = paste('Sample', iMC))
        }
      })
    }

    # Treat Photo-processes ###

    if (dim(Dphoto)[1]!=0) {

      if(nMC > 1) {
        # Clean target dir
        files = list.files(
          path = file.path(targetMCDir,'Photoprocs'),
          pattern = '.dat',
          full.names = TRUE
        )
        if(!is.null(files))
          file.remove(files)

        withProgress(message = 'PhotoProcs', {
          for (iMC in 0:nMC) {
            prefix=paste0(sprintf('%04i',iMC),'_')
            for (i in (1:nbReac)[photo]) {
              sp=species[which(Lphoto[i,]!=0)]
              file = paste0(prefix,'se',sp,'.dat')
              fromFile = file.path(photoSourceDir, file)
              toFile   = file.path(targetMCDir, 'Photoprocs', file)
              file.copy(from = fromFile, to = toFile)
              if( params[[i]][1] != 0 ) {
                file = paste0(prefix,'qy',sp,'_',params[[i]][1],'.dat')
                fromFile = file.path(photoSourceDir, file)
                toFile   = file.path(targetMCDir, 'Photoprocs', file)
                file.copy(from = fromFile, to = toFile)
              }
            }
            incProgress(1/nMC, detail = paste('Sample', iMC))
          }
        })

        # Copy of nominal values in ./Photo for standalone test run
        # + Remove run prefix
        iMC = 0
        prefix=paste0(sprintf('%04i',iMC),'_')
        for (i in (1:nbReac)[photo]) {
          sp=species[which(Lphoto[i,]!=0)]
          file = paste0('se',sp,'.dat')
          fromFile = file.path(photoSourceDir, prefix, file)
          toFile   = file.path(outputDir, 'Run','Photo', file)
          file.copy(from = fromFile, to = toFile)
          if( params[[i]][1] != 0 ) {
            file = paste0('qy',sp,'_',params[[i]][1],'.dat')
            fromFile = file.path(photoSourceDir, prefix, file)
            toFile   = file.path(outputDir, 'Run','Photo', file)
            file.copy(from = fromFile, to = toFile)
          }
        }
      }
    }
  })

# Irradiation ####
getSpectrumReso = function(file) {
  # Get data and check resolution
  sp <- read.csv(file = file, header = FALSE, sep = "")
  reso = abs(sp[2,1]-sp[1,1])
  return(reso)
}
loadSpectrumFile = function(source, path, checkReso = TRUE) {

  if(checkReso) {
    # Check resolution
    speReso = reacData()$spectralResolution
    reso  = getSpectrumReso(path)
    check = abs((reso-speReso)/speReso) < 0.01
    if ( !check ) {
      showModal(modalDialog(
        title = ">>>> Spectrum problem <<<< ",
        paste0('File: ',source,
               ' is inconsistent with declared
                   spectral resolution.',
               '\nPlease choose another file...'),
        easyClose = TRUE,
        footer = modalButton("OK"),
        size = 's'
      ))
    }
    req(check)
  }

  fName  = basename(source)
  target = file.path(projectDir(),'Run',fName)
  if(path != target)
    file.copy(from = path, to = target, overwrite = TRUE)

  # Get spectrum data and store in in memory
  sp = read.csv(file = target, header = FALSE, sep = "")
  spectrumData(
    list(
      beamSpectrumFileName = fName,
      wavelength = sp[, 1],
      photonFlux = sp[, 2]
    )
  )

  # Update data
  ll = reacData()
  ll$beamSpectrumFile = fName
  reacData(ll)
}
output$irradUI <- renderUI({
  req(reacData())

  speReso = reacData()$spectralResolution
  if (is.null(speReso) | !speReso %in% c(0.1, 1))
    speReso = REAC_DATA_default$spectralResolution

  # Process spectrum file
  # & Check compatibility with resolution
  beamSpectrumFile = reacData()$beamSpectrumFile
  if(!is.null(beamSpectrumFile)) {
    bsf = file.path(projectDir(),'Run',beamSpectrumFile)
    if(file.exists(bsf)) {
      reso = getSpectrumReso(bsf)
      if ( abs((reso-speReso)/speReso) > 0.01) {
        showModal(modalDialog(
          title = ">>>> Spectrum problem <<<< ",
          paste0('Local file: ',beamSpectrumFile,
                 ' inconsistent with declared
                   spectral resolution.',
                 '\nPlease choose another file...'),
          easyClose = TRUE,
          footer = modalButton("OK"),
          size = 's'
        ))
        bsf = NULL
      } else {
        loadSpectrumFile(
          source    = bsf,
          path      = bsf,
          checkReso = FALSE)
      }
    } else {
      showModal(modalDialog(
        title = ">>>> Spectrum problem <<<< ",
        paste0('File: ',beamSpectrumFile,
               ' does not exist in local repertory.',
               '\nPlease add it or choose a file...'),
        easyClose = TRUE,
        footer = modalButton("OK"),
        size = 's'
      ))
      bsf = NULL
    }
  }

  # List of existing spectrum files with correct resolution
  source =  file.path( projectDir(),
    '..','..','ChemDBPublic','BeamSpectrumFiles',
    paste0(speReso,'nm'))
  files =  list.files(
    path = source,
    pattern = '.txt',
    full.names = TRUE
  )
  names(files)= file.path(paste0(speReso,'nm'),basename(files))

  tagList(
    fileInput(
      "beamSpectrumFileName",
      label = "Beam spectrum file",
      placeholder = ifelse(
        !is.null(bsf),
        basename(bsf),
        'Choose file...')
    ),
    selectizeInput(
      "bsfSelect",
      label = "Predefined files",
      choices =c(
        "Choose one:" = "",
        files)
    )
  )

})
outputOptions(output, "irradUI",suspendWhenHidden = FALSE)
observeEvent(input$beamSpectrumFileName, {
  source = input$beamSpectrumFileName
  loadSpectrumFile(
    source    = source$name,
    path      = source$datapath,
    checkReso = TRUE)
})
observeEvent(input$bsfSelect, {
  req(input$bsfSelect)
  source = input$bsfSelect
  loadSpectrumFile(
    source    = source,
    path      = source,
    checkReso = FALSE)
})
output$irradUI2 <- renderUI({
  req(reacData())
  req(spectrumData())

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))

  for (n in names(spectrumData()))
    assign(n, rlist::list.extract(spectrumData(), n))

  tagList(
    sliderInput(
      "spectrumRange",
      "Spectrum Range (nm)",
      min   = min(wavelength),
      max   = max(wavelength),
      value = range(wavelength),
      step  = 10
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
  req(reacData())
  req(spectrumData())
  req(input$beamIntensity)

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
  cat("Photons ***********************\n")
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
  cat("*******************************")

})
output$irradSpectrum <- renderPlot({
  req(reacData())
  req(spectrumData())
  req(input$beamIntensity) # Ensure that full UI is defined

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))
  for (n in names(spectrumData()))
    assign(n, rlist::list.extract(spectrumData(), n))
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  # reso = abs(wavelength[2]-wavelength[1])
  # if ( abs((reso-spectralResolution)/spectralResolution) > 0.01)
  #   return(NULL)
  # WIP ####

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
    main = paste0(
      "Photon flux :",
      basename(beamSpectrumFileName),
      " / ", spectralResolution, ' nm')
  )
  grid()
  box()
})

# Reactor ####
output$reactorUI <- renderUI({
  req( reacData() )

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))

  ui = list()
  for (n in listParsReac) {
    ui[[n]] = numericInput(
      n,
      label = paste0(n,' (',listParsReacUnits[n],')'),
      value = unlist(mget(n, ifnotfound = NA)),
      width = '200px'
    )
  }

  tagList(
    fixedRow(
      column(
        width = 4,
        ui[[1]],ui[[2]],ui[[3]],ui[[4]],
      ),
      column(
        width = 4,
        ui[[5]],ui[[6]],ui[[7]],ui[[8]],
      )
    )
  )
})
output$reactorParams <- renderPrint({
  req( input$reactorLength )

  for (n in listParsReac)
    cat(n, " = ", input[[n]], "\n")

  # Update data
  ll = reacData()
  for (n in listParsReac)
    ll[[n]] = input[[n]]
  reacData(ll)

})

# Bottom ####
