traceBack = function(sp1,
                     L,
                     R,
                     species,
                     spInit,
                     reacTags,
                     flMean,
                     vlpInd = NULL,
                     loopAvoid = TRUE,
                     maxVolpert = TRUE) {

  # Build (maximal flux / shortest) production pathway
  # from one of spInit to sp1

  nbSpecies = length(species)
  nbReacs = length(reacTags)

  Fl = matrix(flMean,
              nrow = nbReacs,
              ncol = nbSpecies,
              byrow = FALSE)
  KL = L * Fl
  KR = R * Fl

  if(maxVolpert)
    cat(paste('Shortest production path for ', sp1, '\n'))
  else
    cat(paste('Main production path for ', sp1, '\n'))

  start = sp1
  path = start
  stoppers =spInit # Primordial mixture species
  while (!(start %in% stoppers)) {

    # Define target species
    target = which(species == start)
    if (length(target) == 0)
      break
    vlpIndTarg = vlpInd[target]

    # Identify main productions
    fl = KR[, target]
    fl[!is.finite(fl)] = 0
    productions = order(fl, decreasing = TRUE) # Reactions producing target
    productions = productions[1:sum(fl !=0 )]

    ip = 1
    if(length(productions) > 1 & (loopAvoid | maxVolpert)) {
      # Search for reac implying at least one species
      # with Volpert's indices smaller than target's one
      vlpSum = c()
      for(ip in 1: length(productions)) {
        mainReac = productions[ip]
        sp = which(KL[mainReac,] != 0)
        vlpSum[ip] = sum(vlpInd[sp])
      }
      ip = which.min(vlpSum) # Largest flux with smallest vlpSum

      # for(ip in 1: length(productions)) {
      #   mainReac = productions[ip]
      #   sp = which(KL[mainReac,] != 0)
      #
      #   # Avoid loop
      #   if(any(names(sp) %in% path) & loopAvoid)
      #     next
      #
      #   # Select only reactants with smaller Volpert's index
      #   # (shortest path)
      #   if(all(vlpInd[sp] < vlpIndTarg) & maxVolpert)
      #     break
      #
      # }

      if(ip > length(productions)) # Failed. Might as well take main flux
        ip = 1

      # # Check if chosen reac produces a loop
      # if(length(KL[productions[ip], path]) != 0)
      #   while(any(KL[productions[ip], path] != 0) &
      #         ip < length(productions)) # Reac loops back
      #     ip = ip+1
    }
    mainReac = productions[ip]
    wgt = flMean[mainReac] / sum(flMean[productions])

    cat('  ', reacTags[mainReac], ' : ', signif(wgt, 2), '\n')

    # Prepare iteration

    if(!is.null(vlpInd) & maxVolpert) {
      # Select species with largest Volpert's index
      sp = which(KL[mainReac,] != 0)
      sel = which.max(vlpInd[sp])
      start = species[sp[sel]]
      # cat('Reactants:',sp,'; Volpert:',vlpInd[sp],'; Select:',start)

    } else {
      ## Select next target as max. flux reactant
      io = order(KL[mainReac,], decreasing = TRUE)
      start = species[io[1]]
      if (start %in% stoppers & KL[mainReac, io[2]] != 0)
        start = species[io[2]]

    }

    ## Stopping criteria
    start = start[!(start %in% stoppers)]
    if (length(start) == 0) # Reached primordial photolysis
      break

    if (start %in% path)    # Avoid loops
      break

    # Memorize path
    path = c(path, start)
  }
}
# Functions ####
getLR = function () {

  if(DEBUG)
    print('Start getLR')

  SP = readLines(
    con = file.path(ctrlPars$projectDir,'Run','species_aux.dat'),
    n = 4
  )
  nbSpecies = as.numeric(SP[[1]])

  ## get D and L matrices
  DL  = readLines(
    con = file.path(ctrlPars$projectDir,'Run','photo_DL.dat'),
    n   = 4
  )
  nnz = as.numeric(DL[1])
  Dsp = as.numeric(unlist(strsplit(trimws(DL[2]), split = ' ')))
  Dsp = matrix(Dsp, nrow = nnz, ncol = 3)
  nbPhoto = max(Dsp[,1])
  D = L = matrix(0,ncol=nbSpecies,nrow=nbPhoto)
  for (i in 1:nrow(Dsp))
    D[Dsp[i, 1], Dsp[i, 2]] = Dsp[i, 3]
  nnz = as.numeric(DL[3])
  Lsp = as.numeric(unlist(strsplit(trimws(DL[4]), split = ' ')))
  Lsp = matrix(Lsp, nrow = nnz, ncol = 3)
  for (i in 1:nrow(Lsp))
    L[Lsp[i, 1], Lsp[i, 2]] = Lsp[i, 3]

  DL  = readLines(
    con = file.path(ctrlPars$projectDir,'Run','reac_DL.dat'),
    n   = 4
  )
  nnz = as.numeric(DL[1])
  Dsp = as.numeric(unlist(strsplit(trimws(DL[2]), split = ' ')))
  Dsp = matrix(Dsp, nrow = nnz, ncol = 3)
  nbReac = max(Dsp[,1])
  Dr = Lr = matrix(0,ncol=nbSpecies,nrow=nbReac)
  for (i in 1:nrow(Dsp))
    Dr[Dsp[i, 1], Dsp[i, 2]] = Dsp[i, 3]
  nnz = as.numeric(DL[3])
  Lsp = as.numeric(unlist(strsplit(trimws(DL[4]), split = ' ')))
  Lsp = matrix(Lsp, nrow = nnz, ncol = 3)
  for (i in 1:nrow(Lsp))
    Lr[Lsp[i, 1], Lsp[i, 2]] = Lsp[i, 3]

  D = rbind(D,Dr)
  L = rbind(L,Lr)
  R = D + L

  if(DEBUG)
    print('Exit getLR')

  return(
    list(
      L = L,
      R = R
    )
  )

}
calcFluxes = function(coList,reList,stList) {

  if(DEBUG)
    print('Start calcFluxes')

  # Extract conc et al.
  for (n in names(coList))
    assign(n,rlist::list.extract(coList,n))
  # Extract rates et al.
  for (n in names(reList))
    assign(n,rlist::list.extract(reList,n))
  # Extract L & R mats
  for (n in names(stList))
    assign(n,rlist::list.extract(stList,n))

  # MC fluxes
  test = nfRates > 0
  if(!test)
    showNotification(
      h4('Not enough MC samples !'),
      duration = NULL,
      type = 'warning'
    )
  req(test)

  allRates   = cbind(photoRates,rates)
  reacNames  = colnames(allRates)
  nbReacs    = ncol(allRates)
  nMC        = nrow(allRates)
  nbSpecies  = nsp

  # Concentrations at stationary state (last snapshot)
  conc = matrix(conc[, nt, ], nrow = nMC, ncol = nbSpecies)

  ## Compute flux
  flux = matrix(NA, nrow = nMC, ncol = nbReacs)
  for (i in 1:nMC)
    flux[i, ] = allRates[i, ] *
    apply(L, 1, function(x) prod(conc[i, ] ^ x))

  logFlux = ifelse(flux == 0, NA, log(flux))

  # Summaries
  flMean  = apply(logFlux, 2, function(x)
    exp(mean(x, na.rm = TRUE)))
  flMed   = apply(logFlux, 2, function(x)
    exp(median(x, na.rm = TRUE)))
  flF     = apply(logFlux, 2, function(x)
    exp(sd(x, na.rm = TRUE)))
  flLow   = apply(logFlux, 2, function(x)
    exp(quantile(x, probs = 0.025, na.rm = TRUE)))
  flSup   = apply(logFlux, 2, function(x)
    exp(quantile(x, probs = 0.975, na.rm = TRUE)))

  names(flMean) = reacNames
  names(flMed)  = reacNames
  names(flF)    = reacNames
  names(flLow)  = reacNames
  names(flSup)  = reacNames

  if(DEBUG)
    print('Exit calcFluxes')

  return(
    list(
      flMean = flMean,
      flMed  = flMed,
      flF    = flF,
      flLow  = flLow,
      flSup  = flSup
    )
  )
}

# Interactive ####

observeEvent(
  ignoreInit = TRUE,
  input$calcFlux,
  {
    if(is.null(concList())) {
      id = shiny::showNotification(
        h4('Loading concentrations, be patient...'),
        closeButton = FALSE,
        duration = 5,
        type = 'message'
      )
      on.exit(removeNotification(id), add = TRUE)
      concList( getConc() )
      req(concList())
      if( !is.null(concList()$alert2) )
        id = shiny::showNotification(
          h4(concList()$alert2),
          closeButton = TRUE,
          duration = NULL,
          type = 'error'
        )
      req(is.null(concList()$alert2))
    }

    if(is.null(ratesList())){
      id = shiny::showNotification(
        h4('Loading rates, be patient...'),
        closeButton = FALSE,
        duration = 5,
        type = 'message'
      )
      on.exit(removeNotification(id), add = TRUE)
      ratesList( getRates() )
      req(ratesList())
      if( !is.null(ratesList()$alert) )
        id = shiny::showNotification(
          h4(ratesList()$alert),
          closeButton = TRUE,
          duration = NULL,
          type = 'error'
        )
      req(is.null(ratesList()$alert))
    }

    if(is.null(stoechList()))
      stoechList( getLR() )

    C = concList()
    R = ratesList()
    S = stoechList()

    future({
      calcFluxes(C,R,S)
    }) %...>% fluxesList()

  }
)

output$viewFlow <- renderPlot({
  req(concList(),
     ratesList(),
     fluxesList(),
     stoechList(),
     input$flSpec)

  # Extract data from lists
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))
  for (n in names(ratesList()))
    assign(n,rlist::list.extract(ratesList(),n))
  for (n in names(fluxesList()))
    assign(n,rlist::list.extract(fluxesList(),n))
  for (n in names(stoechList()))
    assign(n,rlist::list.extract(stoechList(),n))

  test = input$flSpec %in% species
  if(!test) # Rq: validate() would not display message !?!?!?
    showNotification(
      h4('Invalid species name !'),
      type = 'warning'
    )
  req(test)

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
    reacs    = c(
      paste0('Ph', 1:ncol(photoRates)),
      paste0('R',  1:ncol(rates))
    ),
    # reacs    = names(flMean),
    reacType = c(rep(1, ncol(photoRates)),
                 rep(2, ncol(rates))),
    reacTypeNames = reacTypeNames,
    flMean = flMean,
    spInit = spInit,
    topShow = 10^input$topShow,
    level = ifelse(input$level, 2, 1),
    showLegend = TRUE,
    curved = input$curvedArrow
  )

  layout = switch(
    input$fluxGraphAlgo,
    FR  = layout_with_fr(g),
    LGL = layout_with_lgl(g,
      root = input$flSpec,
      repulserad = vcount(g)^3 * 10^input$topShow),
    Bipartite = layout_as_bipartite(g,
      vgap = 100 * 10^input$topShow),
    layout_with_gem(g) # Default
  )
  plot(g, layout=layout)

})

# NetworkD3 ####
my_igraph_to_networkd3 = function (g, group, what = "both") {
  if (!("igraph" %in% class(g)))
    stop("g must be an igraph class object.", call. = FALSE)
  if (!(what %in% c("both", "links", "nodes")))
    stop("what must be either \"nodes\", \"links\", or \"both\".",
         call. = FALSE)

  temp_nodes <- V(g) %>% as.matrix %>% data.frame
  # temp_nodes$name <- row.names(temp_nodes)
  temp_nodes$name <- V(g)$label
  names(temp_nodes) <- c("id", "name")
  temp_nodes$id <- temp_nodes$id - 1
  nodes <- temp_nodes$name %>% data.frame %>% setNames("name")

  if (!missing(group)) {
    group <- as.matrix(group)
    if (nrow(nodes) != nrow(group))
      stop("group must have the same number of rows as the number of nodes in g.",
           call. = FALSE)
    nodes <- cbind(nodes, group)
  }
  row.names(nodes) <- NULL
  links <- as_data_frame(g, what = "edges")
  links <- merge(links, temp_nodes, by.x = "from", by.y = "name")
  links <- merge(links, temp_nodes, by.x = "to", by.y = "name")

  # print(str(links))

  # Depends on the edges structure generated in viewFLow()
  links <- links[, c(9:10, 4)] %>%
    setNames(c("source","target","value"))

  # if (ncol(links) == 5) {
  #   links <- links[, c(4:5, 3)] %>% setNames(c("source",
  #                                              "target", "value"))
  # }
  # else {
  #   links <- links[, c("id.x", "id.y")] %>% setNames(c("source",
  #                                                      "target"))
  # }

  if (what == "both") {
    return(list(links = links, nodes = nodes))
  }
  else if (what == "links") {
    return(links)
  }
  else if (what == "nodes") {
    return(nodes)
  }
}
viewFlow = function(sp1,
                    L,
                    R,
                    species,
                    reacs,
                    reacType,
                    reacTypeNames,
                    flMean,
                    spInit,
                    topShow = 0.1,
                    level = 1,
                    showLegend = TRUE,
                    PDF = FALSE,
                    curved = FALSE) {
  # Builds digraph of fluxes to and from sp1

  nbSpecies = length(species)
  nbReacs   = length(reacs)

  nedges = nbReacs + nbSpecies
  lSpecies = (nbReacs+1):nedges
  links = matrix(0, ncol = nedges, nrow = nedges)
  colnames(links) = c(reacs, species)
  rownames(links) = c(reacs, species)

  KL = L * matrix(flMean,
                  nrow = nbReacs,
                  ncol = nbSpecies,
                  byrow = FALSE)
  KR = R * matrix(flMean,
                  nrow = nbReacs,
                  ncol = nbSpecies,
                  byrow = FALSE)

  # Build bi-partite/digraph connectivity matrix
  for (i in 1:nbReacs) {
    links[lSpecies, i] = links[lSpecies, i] + KL[i, 1:nbSpecies]
    links[i, lSpecies] = links[i, lSpecies] + KR[i, 1:nbSpecies]
  }

  if(sp1 != '') {

  }
  # if(!is.null(sp1)) {
    # Build link matrix with first neighbors (reactants and products)
    # links = addLinks(sp1, links, species, KL, KR, nbReacs, nbSpecies)

    # if (level==2) {
    #   # Add 2nd neighnors to link matrix ### DOES NOT WORK AS IS !!!
    #   select = colSums(links[, lSpecies]) != 0 |
    #            rowSums(links[lSpecies, ]) != 0
    #   listSp = species[select]
    #   for (sp2 in listSp)
    #     if (!(sp2 %in% spInit))
    #       links = addLinks(sp2, links, species, KL, KR, nbReacs, nbSpecies,
    #                        wght = 1e-6)
    # }
  # }

  # Select most important links if there are too many
  linksThresh = links
  # if (sum(links != 0) > 50) {
  #   lnkThresh = quantile(abs(links)[abs(links) > 0],
  #                        probs = 1 - topShow,
  #                        na.rm = TRUE)
  #   linksThresh[abs(links) < lnkThresh] = 0
  # }
  # print(linksThresh)

  g = simplify(
    graph_from_adjacency_matrix(
      linksThresh,
      mode = "directed",
      weighted = TRUE)
  )

  cols=brewer.pal(9,"Set3")
  reacColor=c(cols[6:9])

  V(g)$label = c(reacs,species)
  V(g)$shape=c(rep("rectangle",nbReacs),rep("circle",nbSpecies))
  V(g)$size = c(rep(20,nbReacs),rep(20,nbSpecies))
  V(g)$size2 = 6
  V(g)$color = c(reacColor[reacType],rep("gold",nbSpecies))
  V(g)$label.cex = c(rep(1,nbReacs),rep(1,nbSpecies))
  V(g)$label.color = "black"
  V(g)$label.font = 2
  V(g)$type = c(rep(0,nbReacs),rep(1,nbSpecies))

  wid = abs(E(g)$weight)
  wid = wid^0.1 # Empirical transfo for better scaling...
  wid = 0.12 + (wid-min(wid))/(max(wid)-min(wid))
  E(g)$width = wid*5
  E(g)$color = ifelse(
    E(g)$weight > 0,
    col2tr('red',160),
    col2tr('blue',160))
  E(g)$arrow.size  = 0.25*max(E(g)$width)
  E(g)$arrow.width = 0.25*max(E(g)$width)
  E(g)$curved = curved

  loners = which(degree(g) == 0)
  g = delete.vertices(g, loners)

  return(g)
}
output$viewFlowD3 <- renderForceNetwork({
  if(is.null(concList()) |
     is.null(ratesList()) |
     is.null(fluxesList()) |
     is.null(stoechList())
  )
    return(NULL)

  # Extract data from lists
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))
  for (n in names(ratesList()))
    assign(n,rlist::list.extract(ratesList(),n))
  for (n in names(fluxesList()))
    assign(n,rlist::list.extract(fluxesList(),n))
  for (n in names(stoechList()))
    assign(n,rlist::list.extract(stoechList(),n))

  if(!is.null(input$flSpec) & input$flSpec != '' ) {
    test = input$flSpec %in% species
    if(!test) # Rq: validate() would not display message !?!?!?
      showNotification(
        h4('Invalid species name !'),
        type = 'warning'
      )
    req(test)
  }

  # Remove dummy species
  # sel  = ! species %in% spDummy
  sel = 1:length(species)

  reacTypeNames = c("Ph", "R")
  g = viewFlow(
    input$flSpec,
    L[, sel],
    R[, sel],
    species[sel],
    reacs    = c(
      paste0('Ph', 1:nPhotoRates),
      paste0('R',  1:nRates)
    ),
    # reacs = make.unique(trimws(names(flMean))),
    reacType = c(
      rep(1, ncol(photoRates)),
      rep(2, ncol(rates))
    ),
    reacTypeNames = reacTypeNames,
    flMean = flMean,
    spInit = spInit,
    topShow = 10^input$topShow,
    level = ifelse(input$level, 2, 1),
    showLegend = TRUE,
    curved = FALSE
  )

  grp = as.numeric(factor(V(g)$shape))-1

  graph_d3 = my_igraph_to_networkd3(g, group = grp)
  graph_d3$nodes[['size']]  = log10(igraph::strength(g))^0.7
  graph_d3$nodes[['size']][1:(nPhotoRates+nRates)] = 1
  # graph_d3$links[['value']] =
  #   pmax(1,10*abs(E(g)$weight)/max(abs(E(g)$weight)))

  # print(E(g)$weight)
  # print(graph_d3$nodes)
  # print(graph_d3$links, max = 1000)

  # Loss/prod links colors
  target = which(graph_d3$nodes[['name']] == input$flSpec)
  linkOut = graph_d3$links[['source']] == target-1
  for (i in which(linkOut) )
     linkOut = linkOut |
    graph_d3$links[['source']] == graph_d3$links[['target']][i]
  linkIn  = graph_d3$links[['target']] == target-1
  for (i in which(linkIn) )
    linkIn = linkIn |
    graph_d3$links[['target']] == graph_d3$links[['source']][i]
  linkColour = rep('#BBB',length(linkIn))
  linkColour[linkIn]  = gPars$cols[2]
  linkColour[linkOut] = gPars$cols[5]

  networkD3::forceNetwork(
    Links    = graph_d3$links,
    Nodes    = graph_d3$nodes,
    Source   = 'source',
    Target   = 'target',
    NodeID   = 'name',
    Group    = 'group',
    Value    = 'value',
    Nodesize = 'size',
    radiusCalculation = JS("d.nodesize + 3"),
    bounded  = FALSE,
    zoom     = TRUE,
    arrows   = TRUE,
    linkDistance = if(input$scaleDistD3)
      JS("function(d){return 50 / Math.pow(d.value,0.5) }")
    else
      60,
    opacity  = 0.9,
    colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"),
    clickAction = 'alert(d.name );',
    linkColour = linkColour,
    legend = FALSE,
    opacityNoHover = 0.9,
    charge = input$forceNetChargeFlux
  )

})
# Budget ####
output$viewBudget <- renderPrint({
  req(concList(),
     ratesList(),
     fluxesList(),
     stoechList())

  # Extract data from lists
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))
  for (n in names(ratesList()))
    assign(n,rlist::list.extract(ratesList(),n))
  for (n in names(fluxesList()))
    assign(n,rlist::list.extract(fluxesList(),n))
  for (n in names(stoechList()))
    assign(n,rlist::list.extract(stoechList(),n))

  if(!is.null(input$flSpec) & input$flSpec != '' ) {
    test = input$flSpec %in% species
    if(!test) # Rq: validate() would not display message !?!?!?
      showNotification(
        h4('Invalid species name !'),
        type = 'warning'
      )
    req(test)
  }

  # Remove dummy species
  sel  = ! species %in% spDummy
  spec = species[sel]
  LR   = L[,sel]
  RR   = R[,sel]

  budget(
    input$flSpec,
    LR,RR,spec,
    names(flMean),
    flMean,
    weightThresh=10^input$topShow
  )
})

output$viewTarget <- renderPrint({
  req(concList(),
      ratesList(),
      fluxesList(),
      stoechList())

  # Extract data from lists
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))
  for (n in names(ratesList()))
    assign(n,rlist::list.extract(ratesList(),n))
  for (n in names(fluxesList()))
    assign(n,rlist::list.extract(fluxesList(),n))
  for (n in names(stoechList()))
    assign(n,rlist::list.extract(stoechList(),n))

  if(!is.null(input$flSpec) & input$flSpec != '' ) {
    test = input$flSpec %in% species
    if(!test) # Rq: validate() would not display message !?!?!?
      showNotification(
        h4('Invalid species name !'),
        type = 'warning'
      )
    req(test)
  }

  # Remove dummy species
  sel  = ! species %in% spDummy
  spec = species[sel]
  LR   = L[,sel]
  RR   = R[,sel]

  traceBack(
    input$flSpec,
    LR, RR, spec, spInit,
    names(flMean),
    flMean,
    reacScheme()$vlpInd,
    loopAvoid = input$loopAvoid,
    maxVolpert = input$maxVolpert
  )
})
