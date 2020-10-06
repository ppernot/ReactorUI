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
  Dsp = as.numeric(unlist(strsplit(trimws(DL[2]), split = ' ')))
  Dsp = matrix(Dsp, nrow = nnz, ncol = 3)
  for (i in 1:nrow(Dsp))
    D[Dsp[i, 1], Dsp[i, 2]] = Dsp[i, 3]
  nnz = as.numeric(DL[3])
  Lsp = as.numeric(unlist(strsplit(trimws(DL[4]), split = ' ')))
  Lsp = matrix(Lsp, nrow = nnz, ncol = 3)
  for (i in 1:nrow(Lsp))
    L[Lsp[i, 1], Lsp[i, 2]] = Lsp[i, 3]
  nPh = ncol(photoRates)

  DL  = readLines(
    con = paste0(ctrlPars$projectDir,'/Run/reac_DL.dat'),
    n   = 4
  )
  nnz = as.numeric(DL[1])
  Dsp = as.numeric(unlist(strsplit(trimws(DL[2]), split = ' ')))
  Dsp = matrix(Dsp, nrow = nnz, ncol = 3)
  for (i in 1:nrow(Dsp))
    D[nPh + Dsp[i, 1], Dsp[i, 2]] = Dsp[i, 3]
  nnz = as.numeric(DL[3])
  Lsp = as.numeric(unlist(strsplit(trimws(DL[4]), split = ' ')))
  Lsp = matrix(Lsp, nrow = nnz, ncol = 3)
  for (i in 1:nrow(Lsp))
    L[nPh + Lsp[i, 1], Lsp[i, 2]] = Lsp[i, 3]

  R = D + L

  for (i in 1:nf)
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

      C = getConc()
      concList(C)
    }
    C = concList()

    test = dim(C$conc)[1] > 1
    if(!test) # Rq: validate() would not display message !?!?!?
      showNotification(
        h4('Flux analysis impossible for a single run !'),
        type = 'warning'
      )
    req(test)

    if(is.null(ratesList())){
      id = shiny::showNotification(
        h4('Loading rates, be patient...'),
        closeButton = FALSE,
        duration = 5,
        type = 'message'
      )
      on.exit(removeNotification(id), add = TRUE)

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
    Bipartite = layout_as_bipartite(g,
      vgap = 100 * input$topShow),
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
# NetworkD3 ####
my_igraph_to_networkd3 = function (g, group, what = "both") {
  if (!("igraph" %in% class(g)))
    stop("g must be an igraph class object.", call. = FALSE)
  if (!(what %in% c("both", "links", "nodes")))
    stop("what must be either \"nodes\", \"links\", or \"both\".",
         call. = FALSE)
  temp_nodes <- V(g) %>% as.matrix %>% data.frame
  temp_nodes$name <- row.names(temp_nodes)
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
output$viewFlowD3 <- renderForceNetwork({
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

  test = input$flSpec %in% species
  if(!test) # Rq: validate() would not display message !?!?!?
    showNotification(
      h4('Invalid species name !'),
      type = 'warning'
    )
  req(test)

  # Remove dummy species
  sel  = ! species %in% spDummy

  reacTypeNames = c("Ph", "R")
  g = viewFlow(
    input$flSpec,
    L[, sel],
    R[, sel],
    species[sel],
    # reacs    = c(
    #   paste0('Ph', 1:ncol(photoRates)),
    #   paste0('R',  1:ncol(rates))
    # ),
    reacs = names(flMean),
    reacType = c(
      rep(1, ncol(photoRates)),
      rep(2, ncol(rates))
    ),
    reacTypeNames = reacTypeNames,
    flMean = flMean,
    spInit = spInit,
    topShow = input$topShow,
    level = ifelse(input$level, 2, 1),
    showLegend = TRUE,
    curved = input$curvedArrow
  )

  grp = as.numeric(factor(V(g)$shape))-1

  graph_d3 = my_igraph_to_networkd3(g, group = grp)
  graph_d3$nodes[['size']]  = 2 * igraph::degree(g)
  # graph_d3$links[['value']] =
  #   pmax(1,10*abs(E(g)$weight)/max(abs(E(g)$weight)))

  print(E(g)$weight)
  print(graph_d3$nodes)
  print(graph_d3$links, max = 1000)

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
    # linkColour = "#BBB",
    linkColour = linkColour,
    legend = FALSE,
    opacityNoHover = 0.9,
    charge = input$forceNetChargeFlux
  )

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
