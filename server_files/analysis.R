# Functions ####
getConc  = function(concThresh = -50) {

  # Read Ctrl file to get spInit
  ctrlList = readCtrl(ctrlPars$projectDir)
  csp = ctrlList$REAC_DATA$reactantsComposition
  sel = which(is.finite(csp))
  rsp = ctrlList$REAC_DATA$reactantsSpecies
  rsp = strsplit(rsp,',')[[1]][sel]
  spInit = rsp

  # Load 1 sample files
  files = list.files(
    path = file.path(ctrlPars$projectDir,'MC_Output'),
    pattern = "fracmol_",
    full.names = TRUE
  )

  # x = read.table(files[1], header=TRUE, skip=0)
  x = as.data.frame(
    data.table::fread(files[1], header=TRUE, skip=0)
  )
  time=x[-1,1]

  # Species list
  line  = readLines(con=files[1], n=1)
  species = scan(text=line, what=character(),
                 strip.white=TRUE, quiet=TRUE)[-1]

  # Get all data
  nf   = length(files)
  nsp  = length(species)
  nt   = length(time)
  conc = moleFrac = array(0,dim=c(nf,nt,nsp))
  for(i in seq_along(files[1:nf]) ) {
    # tab  = read.table(files[i], header=TRUE, skip=0)
    tab  = as.data.frame(
      data.table::fread(files[i], header=TRUE, skip=0)
    )
    time = tab[-1,1]
    y    = as.matrix(tab[-1,-1])
    conc[i,1:nt,1:nsp]     = y[1:nt,1:nsp]
    moleFrac[i,1:nt,1:nsp] = y[1:nt,1:nsp] / rowSums(y)
  }

  # Pretreat moleFrac for plots
  moleFrac = ifelse(moleFrac == 0, NA, log10(moleFrac))
  # moleFrac = ifelse(moleFrac<=concThresh, NA, moleFrac)
  # --> Tresholding now done by plotting function

  yMean = apply(moleFrac,c(2,3),
                function(x) mean(x,na.rm=TRUE))
  ySd   = apply(moleFrac,c(2,3),
                function(x) sd(x,na.rm=TRUE))
  yLow95= apply(moleFrac,c(2,3),
                function(x) quantile(x,probs = 0.025,na.rm=TRUE))
  ySup95= apply(moleFrac,c(2,3),
                function(x) quantile(x,probs = 0.975,na.rm=TRUE))

  # Separate neutrals and ions
  charge = rep(0,length(species))
  ions   = grepl("\\+$",species)
  neus = !ions

  # Calculate mass
  compo = t(apply(as.matrix(species,ncol=1),1,get.atoms))
  colnames(compo)=elements
  mass  = apply(compo,1,massFormula)

  # Save final mole fractions to disk
  dirOut = file.path(ctrlPars$projectDir, 'Outputs')
  if(!dir.exists(dirOut))
    dir.create(dirOut, showWarnings = FALSE)
  cnv = log10(exp(1))
  mmol = 10^yMean[nt,]

  # Save Neutrals
  df = data.frame(
    mass  = mass[neus],
    mean  = mmol[neus],
    sd    = cnv * mmol[neus] * ySd[nt,neus],
    low95 = 10^yLow95[nt,neus],
    sup95 = 10^ySup95[nt,neus],
    row.names = species[neus]
  )
  io = order(mass[neus],na.last = TRUE)
  write.csv(
    df[io,],
    file = file.path(dirOut,'MoleFracNeutrals_tFinal.csv')
  )

  # Save Ions
  df = data.frame(
    mass  = mass[ions],
    mean  = mmol[ions],
    sd    = cnv * mmol[ions] * ySd[nt,ions],
    low95 = 10^yLow95[nt,ions],
    sup95 = 10^ySup95[nt,ions],
    row.names = species[ions]
  )
  io = order(mass[ions],na.last = TRUE)
  write.csv(
    df[io,],
    file = file.path(dirOut,'MoleFracIons_tFinal.csv')
  )

  return(
    list(
      nf       = nf,
      nt       = nt,
      nsp      = nsp,
      time     = time,
      conc     = conc,
      species  = species,
      spInit   = spInit,
      mfMean   = yMean,
      mfSd     = ySd,
      mfLow    = yLow95,
      mfSup    = ySup95
    )
  )
}
getRates = function() {
  # Load react. rates sample files
  files = list.files(
    path    = paste0(ctrlPars$projectDir,'/MC_Output'),
    pattern = 'reacs_rates_',
    full.names=TRUE)
  nf = length(files)

  ## Get MC rates
  x = as.data.frame(
    data.table::fread(files[1], header=FALSE, skip=0)
  )
  nRates=length(x)
  rates=matrix(NA,nrow=nf,ncol=nRates)
  for(i in seq_along(files) ) {
    tab = as.data.frame(
      data.table::fread(files[i], header=FALSE,
                     skip=0, colClasses='numeric')
    )
    rates[i,] = t(tab)
  }
  ## Get reac names
  reacs = readLines(
    file.path(ctrlPars$projectDir,'Run','reac_list.dat'))
  colnames(rates) = reacs

  # Load photorates sample files
  files = list.files(
    path=paste0(ctrlPars$projectDir,'/MC_Output'),
    pattern='photo_rates_',
    full.names=TRUE)
  nf = length(files)

  ## Get MC photo rates
  x = as.data.frame(
    data.table::fread(files[1], header=FALSE, skip=0)
  )
  nPhotoRates=length(x)
  photoRates=matrix(NA,nrow=nf,ncol=nPhotoRates)
  for(i in seq_along(files) ) {
    tab = as.data.frame(
      data.table::fread(files[i], header=FALSE,
                     skip=0,
                     colClasses=c('numeric')
      )
    )
    photoRates[i,] = t(tab)
  }

  ## Get reac names
  reacs = readLines(
    file.path(ctrlPars$projectDir,'Run','photo_list.dat'))
  colnames(photoRates)= reacs

  return(
    list(
      nfRates     = nf,
      nRates      = nRates,
      rates       = rates,
      nPhotoRates = nPhotoRates,
      photoRates  = photoRates
    )
  )
}
getIntegStats = function() {
  # Load integrator stats
  files = list.files(
    path    = paste0(ctrlPars$projectDir,'/MC_Output'),
    pattern = 'integ_stats_',
    full.names=TRUE)
  nf = length(files)

  ## Get MC values
  x = as.data.frame(
    data.table::fread(files[1], header=FALSE, skip=0)
  )
  time = x[,1]
  nstat = ncol(x)-1
  ntime = nrow(x)
  stats = array(0,dim=c(nf,ntime,nstat))
  for(i in seq_along(files) ) {
    tab = as.data.frame(
      data.table::fread(files[i], header=FALSE,
                        skip=0, colClasses='numeric')
    )
    stats[i,1:ntime,1:nstat] = as.matrix(tab[,-1])
  }

  yMean = apply(stats,c(2,3),
                function(x) mean(x,na.rm=TRUE))
  ySd   = apply(stats,c(2,3),
                function(x) sd(x,na.rm=TRUE))
  yLow95= apply(stats,c(2,3),
                function(x) quantile(x,probs = 0.025,na.rm=TRUE))
  ySup95= apply(stats,c(2,3),
                function(x) quantile(x,probs = 0.975,na.rm=TRUE))

  statNames = c('NFE', 'NFI', 'NSTEPS', 'NACCPT',
                'NREJCT', 'NFESIG', 'MAXM', 'SPRAD')[1:nstat]

  colnames(yMean) =
    colnames(ySd) =
    colnames(yLow95) =
    colnames(ySup95) =
    statNames

  return(
    list(
      nfStat   = nf,
      nt       = ntime,
      nStat    = nstat,
      time     = time,
      stats    = stats,
      statNames= statNames,
      statMean = yMean,
      statSd   = ySd,
      statLow  = yLow95,
      statSup  = ySup95
    )
  )
}
selectSpecies <- function(species, categs) {
  # Select species to plot from categsPlot

  compo   = t(apply(as.matrix(species,ncol=1),1,get.atoms))
  colnames(compo)=elements

  # Remove dummy species
  sel0 = ! species %in% spDummy

  if ('all' %in% categs) {
    return ( sel0)

  } else {

    # Charge
    charge = rep(0,length(species))
    ions   = grepl("\\+$",species)
    charge[ions] = 1
    neus = !ions
    sel  = rep(FALSE,length(species))
    for (cl in categs) {
      if (cl == "neutrals") {
        sel = sel | neus
      } else if (cl == "ions") {
        sel = sel | ions
      }
    }
    selCharge = sel

    # Radicals
    radic = (apply(compo,1,numElec)-charge)%%2
    for (cl in categs) {
      if (cl == "radicals") {
        sel = radic
      }
    }
    selRadic = sel

    # Composition Elementale
    azot = grepl("N",species)
    oxy  = grepl("O",species)
    sel  = rep(FALSE,length(species))
    for (cl in categs) {
      if (cl == "hydrocarbons") {
        sel = sel | (!azot & !oxy)
      } else if (cl == "N-bearing") {
        sel = sel | azot
      } else if (cl == "O-bearing") {
        sel = sel | oxy
      }
    }
    selElem = sel

    # Heavy elements
    nHeavy = nbHeavyAtoms(species)
    sel    = rep(FALSE,length(species))
    for (cl in categs) {
      if (cl == "C0") {
        sel = sel | nHeavy == 0
      } else if (cl == "C1") {
        sel = sel | nHeavy == 1
      } else if (cl == "C2") {
        sel = sel | nHeavy == 2
      } else if (cl == "C3") {
        sel = sel | nHeavy == 3
      } else if (cl == "C4") {
        sel = sel | nHeavy == 4
      } else if (cl == "C5") {
        sel = sel | nHeavy == 5
      } else if (cl == "C6") {
        sel = sel | nHeavy == 6
      } else if (cl == "Cmore") {
        sel = sel | nHeavy > 6
      }
    }
    selHeavy = sel

    return(sel0 & selCharge & selRadic & selElem & selHeavy)

  }

}
assignColorsToSpecies <- function(colSel, species, sel, nf,
                                  cols, col_tr, col_tr2,
                                  threshTransp = 50) {
  # Manage transparency wrt nb MC runs
  if(nf <= threshTransp)
    col_tr = col_tr2 # Darker

  if(colSel) {
    # Assign colors to species to avoid changes
    nrep     = ceiling(length(species)/length(cols))
    colsSp   = rep(cols  ,times = nrep)
    col_trSp = rep(col_tr,times = nrep)
  } else {
    # Use sequential colors for selection
    colsSp = col_trSp = rep(NA,length(species))
    nrep          = ceiling(length(species[sel])/length(cols))
    colsSp[sel]   = rep(cols  ,times = nrep)
    col_trSp[sel] = rep(col_tr,times = nrep)
  }

  return (
    list(
     colsSp   = colsSp,
     col_trSp = col_trSp
    )
  )
}
hsicMat <- function(C,S) {
  if(is.vector(C))
    C = matrix(C,ncol=1)
  if(is.vector(S))
    S = matrix(S,ncol=1)

  nC = ncol(C); nS = ncol(S)
  hMat = matrix(NA,nrow=nS,ncol=nC)
  sxx = future.apply::future_apply(
    X = C, MARGIN = 2,
    FUN = function(x) dHSIC::dhsic(x,x)$dHSIC
  )
  syy = future.apply::future_apply(
    X = S, MARGIN = 2,
    FUN = function(x) dHSIC::dhsic(x,x)$dHSIC
  )
  for(j in 1:nC) {
    x = C[,j]
    V = future.apply::future_apply(
      X = S, MARGIN = 2,
      FUN = function(y,x) dHSIC::dhsic(x,y)$dHSIC,
      x = x
    )
    # Filter problematic values
    V[!is.finite(V)] = 0
    # Reduced Indices
    V = sqrt(V / sqrt(sxx[j]*syy))
    # Normalize
    V = V / sum(V)
    hMat[,j] = V
  }
  # print('hsicMat OK')
  return(hMat)
}
dcorMat <- function(C,S) {
  if(is.vector(C))
    C = matrix(C,ncol=1)
  if(is.vector(S))
    S = matrix(S,ncol=1)

  nC = ncol(C); nS = ncol(S)
  hMat = matrix(NA,nrow=nS,ncol=nC)
  sxx = future.apply::future_apply(
    X = C, MARGIN = 2,
    FUN = function(x) fda.usc::dcor.xy(x,x,test=FALSE)
  )
  syy = future.apply::future_apply(
    X = S, MARGIN = 2,
    FUN = function(x) fda.usc::dcor.xy(x,x,test=FALSE)
  )
  for(j in 1:nC) {
    x = C[,j]
    V = future.apply::future_apply(
      X = S, MARGIN = 2,
      FUN = function(y,x) fda.usc::dcor.xy(x,y,test=FALSE),
      x = x
    )
    # Filter problematic values
    V[!is.finite(V)] = 0
    # Reduce
    V = sqrt(V / sqrt(sxx[j]*syy))
    # Normalize
    V = V / sum(V)
    hMat[,j] = V
  }
  # print('dcorMat OK')
  return(hMat)
}
# Interactive ####
observeEvent(
  ignoreInit = TRUE,
  input$loadMC, {

    id = shiny::showNotification(
      h4('Loading concentrations, be patient...'),
      closeButton = FALSE,
      duration = 5,
      type = 'message'
    )
    on.exit(removeNotification(id), add = TRUE)

    future({getConc()}) %...>% concList()
  }
)
observeEvent(
  ignoreInit = TRUE,
  input$loadMC,
  future({getRates()}) %...>% ratesList()
)
observeEvent(
  ignoreInit = TRUE,
  input$loadMC,
  future({getIntegStats()}) %...>% statsList()
)

output$loadMsg <- renderPrint({
  if(is.null(concList()) | is.null(ratesList()))
    return(NULL)

  # Extract conc et al.
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))

  cat('Read concentrations:\n',
      nf,' MC runs\n',
      nsp,' species\n',
      nt,' snapshots\n\n')

  # Extract rates
  for (n in names(ratesList()))
    assign(n,rlist::list.extract(ratesList(),n))

  cat('Read rates:\n',
      nfRates,' MC runs\n',
      nRates,' reactions\n',
      nPhotoRates,' photo-processes\n')

})

# Kinetics ####
rangesKinetics <- reactiveValues(x = NULL, y = NULL)
output$kinetics <- renderPlot({
  if(is.null(concList()))
    return(NULL)

  # Extract conc et al.
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))

  # Species list to plot
  sel = selectSpecies(species, input$categsPlot)
  if(sum(sel) == 0)
    showNotification(
      h4('Your selection is empty !'),
      type = 'warning'
    )
  req(sum(sel) > 0)

  # Define zoom range
  if (is.null(rangesKinetics$x)) {
    xlim <- range(time)# * 100 # expand for labels
  } else {
    xlim <- rangesKinetics$x
  }
  if (is.null(rangesKinetics$y)) {
    # ylim = range(c(mfLow[, sel], mfSup[, sel]), na.rm = TRUE)
    ylim = range(c(mfLow[, sel], 1), na.rm = TRUE)
    # ylim <- c(
    #   max(-30, min(mfLow[,sel])),
    #   min(0  , max(mfSup[,sel]))
    # )
  } else {
    ylim <- rangesKinetics$y
  }

  # Extract graphical params
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  # Generate colors per species
  colors = assignColorsToSpecies(
    input$colSel, species, sel, nf=1,
    cols, col_tr, col_tr2)

  # Show uncertainty bands ?
  showBands = (nf > 1) & input$mcPlot

  # Plot
  par(mfrow = c(1, 1),
      cex = cex, cex.main = cex, mar = c(mar[1:3],3.5),
      mgp = mgp, tcl = tcl, pty = pty, lwd = lwd)

  plot(time, time,
       type = 'n',
       log  = 'x',
       main = ifelse(!showBands,
                     'Mean values',
                     'Mean and 95 % Proba. Intervals'),
       xlab = 'Time / s',
       xlim = xlim,
       xaxs = 'i',
       ylab = 'Log10(mole fraction)',
       ylim = ylim,
       yaxs = 'i')
  grid()

  if(input$ppscale) {
    axis(side = 2,
         at = -3*(2:5),
         labels = FALSE,
         las = 1,
         adj = 1,
         line = 0,
         tck = 1,
         cex.lab = 0.8)
    mtext(at   = -3*(2:5),
          text = c("ppm", "ppb","ppt","ppq"),,
          side = 2,
          col  = 'gray30',
          las  = 1,
          adj  = 1,
          line = 0.2)
  }

  # CI
  if(showBands) {
    for(isp in (1:nsp)[sel])
      polygon(c(time     , rev(time      )),
              c(mfLow[,isp],rev(mfSup[,isp])),
              col = colors$col_trSp[isp],
              border = NA)
  }

  # Mean concentrations
  matplot(time, mfMean[,sel], type='l', lty = 1,
          col = colors$colsSp[sel],
          lwd = 1.2*lwd, add = TRUE)
  # text(x = time[nt], y= mfMean[nt,sel],
  #      labels = species[sel],
  #      col=colors$colsSp[sel],
  #      pos = 4, offset=0.2, cex=cex.leg)

  mtext(at  = mfMean[nt,sel],
       text = species[sel],
       side = 4,
       col  = colors$colsSp[sel],
       las  = 1,
       adj  = 0,
       line = 0.2)
  box()

})
outputOptions(output, "kinetics",
              suspendWhenHidden = FALSE)
observeEvent(
  input$kinetics_dblclick,
  {
    brush <- input$kinetics_brush
    if (!is.null(brush)) {
      rangesKinetics$x <- c(brush$xmin, brush$xmax)
      rangesKinetics$y <- c(brush$ymin, brush$ymax)
    } else {
      rangesKinetics$x <- NULL
      rangesKinetics$y <- NULL
    }
  }
)

# Pseudo MS ####
output$pseudoMS <- renderPlot({
  if(is.null(concList()))
    return(NULL)

  # Extract conc et al.
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))

  # Species list to plot
  selNeu = selectSpecies(species,
                         c("neutrals", input$categsPlotMS))
  selIon = selectSpecies(species,
                         c("ions",    input$categsPlotMS))

  # Stoechiometry & Mass
  compo = t(apply(as.matrix(species,ncol=1),1,get.atoms))
  colnames(compo)=elements
  mass  = apply(compo,1,massFormula)
  names(mass) = species

  # Extract graphical params
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  # Generate colors per species
  colorsNeu = assignColorsToSpecies(
    TRUE, species, selNeu,
    nf=1, cols, col_tr, col_tr2)
  colorsIon = assignColorsToSpecies(
    TRUE, species, selIon,
    nf=1, cols, col_tr, col_tr2)

  # Function to define bars width
  # (width decreases as mole fraction increases)
  lineWidth = function(y,threshMS) {
    ymin = threshMS[1]
    w    = ifelse(y < ymin, y/y, 2000/abs(y+35)^2)
    return( pmin(w,12) )
  }

  # MS sampling time
  t0 = max(min(time),min(max(time),10^input$timeMS))
  nt = which(time >= t0)[1]

  # Plot
  xlim = range(c(mass[selNeu],mass[selIon]), na.rm=TRUE)
  ylim = c(input$threshMS[1],input$threshMS[2]+2)

  par(mfrow = c(2, 1),
      cex = cex, cex.main = cex, mar = c(3.2, 3, 3, 3),
      mgp = mgp, tcl = tcl, pty = 'm', lwd = lwd, lend = 2)

  # Neutrals
  x    = mass[selNeu]
  y    = mfMean[nt, selNeu]
  yLow = mfLow[nt, selNeu]
  ySup = mfSup[nt, selNeu]
  w    = lineWidth(y,input$threshMS) * input$widthMS
  col  = colorsNeu$colsSp[selNeu]
  colt = colorsNeu$col_trSp[selNeu]

  plotMS( x, y, yLow, ySup, xlim, ylim, w, col, colt,
          species[selNeu],xlab ='mass',
          ppScale = input$ppScaleMS, errBar = input$mcPlotMS,
          main = paste0('Neutrals / Time = 10^',
                        round(log10(t0)),' s' ),
          text.cex = input$MStext.cex/5 )

  #Ions
  x    = mass[selIon]
  y    = mfMean[nt, selIon]
  yLow = mfLow[nt, selIon]
  ySup = mfSup[nt, selIon]
  w    = lineWidth(y,input$threshMS) * input$widthMS
  col  = colorsIon$colsSp[selIon]
  colt = colorsIon$col_trSp[selIon]


  plotMS( x, y, yLow, ySup, xlim, ylim,
          w, col, colt, species[selIon],
          ppScale = input$ppScaleMS, errBar = input$mcPlotMS,
          main = paste0('Ions / Time = 10^',
                        round(log10(t0)),' s' ),
          text.cex = input$MStext.cex/5
          )

})

# Sensitivity ####
MRList <- reactiveVal()
observeEvent(
  ignoreInit = TRUE,
  input$doSA,
  {
    if(is.null(concList()) |
       is.null(ratesList()) )
      return(NULL)

    MRList(NULL)

    # Extract conc et al.
    for (n in names(concList()))
      assign(n,rlist::list.extract(concList(),n))

    # Exit if single run
    test = dim(conc)[1] > 1
    if(!test) # Rq: validate() would not display message !?!?!?
      showNotification(
        h4('SA impossible for a single MC run !'),
        type = 'warning'
      )
    req(test)

    # Concentrations at stationary state (last snapshot)
    conc = conc[,nt,]

    # Filter out undesirable values
    conc = ifelse(conc==0, NA, log10(conc))
    sdc  = apply(conc,2,function(x) sd(x))
    selC = sdc!=0 & is.finite(sdc)
    C = conc[,selC]
    colnames(C) = species[selC]

    # Reaction rates
    for (n in names(ratesList()))
      assign(n,rlist::list.extract(ratesList(),n))

    rates = ifelse(rates==0, NA, log10(rates))
    sdc = apply(rates,2,function(x) sd(x))
    selR = sdc!=0 & is.finite(sdc)
    photoRates = ifelse(photoRates==0, NA, log10(photoRates))
    sdc = apply(photoRates,2,function(x) sd(x))
    selPR = sdc!=0 & is.finite(sdc)
    S = cbind(photoRates[,selPR],rates[,selR])

    SASpecies = input$SASpecies
    if(SASpecies != "") { # Check validity
      test = SASpecies %in% colnames(C)
      if(!test) # Rq: validate() would not display message !?!?!?
        showNotification(
          h4('Invalid species name, or insufficient data for SA !'),
          type = 'warning'
        )
      req(test)
    }

    ntop = min(20,ncol(S))

    if (input$anaType == 'spearman') {
      # CORR
      if (SASpecies == "") {
        main = 'TOP20 most influent reactions by Sum(RCC^2)'
        indRates = colSums(cor(C, S, method = "spearman") ^ 2)
        names(indRates) = colnames(S)
        MR = sort(indRates, decreasing = TRUE)[1:ntop]
      } else {
        main = paste0(SASpecies,
                      ' / TOP20 most correlated reactions by RCC')
        isp = which(colnames(C) == SASpecies)
        indRates = as.vector(cor(C[, isp], S, method = "spearman"))
        names(indRates) = colnames(S)
        indx = order(abs(indRates), decreasing = TRUE)[1:ntop]
        MR   = indRates[indx]
      }

    } else if (input$anaType == 'dcorr') {
      # fda.usc::dcor.xy
      id = shiny::showNotification(
        h4('Estimating dcor, be patient...'),
        closeButton = FALSE,
        duration = 5,
        type = 'message'
      )
      on.exit(removeNotification(id), add = TRUE)
      if(SASpecies == "") {
        main = 'TOP20 most influent reactions by Sum(dCor)'
        indRates = rowSums(dcorMat(C,S),na.rm = TRUE)
      } else {
        main = paste0(
          SASpecies,
          ' / TOP20 most influent reactions by dCor')
        isp = which(colnames(C) == SASpecies)
        indRates = as.vector(dcorMat(C[,isp],S))
      }
      names(indRates) = colnames(S)
      MR = sort(indRates,decreasing = TRUE)[1:ntop]

    } else if (input$anaType == 'hsic') {
      # HSIC
      id = shiny::showNotification(
        h4('Estimating HSIC, be patient...'),
        closeButton = FALSE,
        duration = 5,
        type = 'message'
      )
      on.exit(removeNotification(id), add = TRUE)
      if (SASpecies == "") {
        main = 'TOP20 most influent reactions by Sum(HSIC)'
        indRates = rowSums(hsicMat(C, S), na.rm = TRUE)
      } else {
        main = paste0(SASpecies,
                      ' / TOP20 most influent reactions by HSIC')
        isp = which(colnames(C) == SASpecies)
        indRates = as.vector(hsicMat(C[, isp], S))
      }
      names(indRates) = colnames(S)
      MR = sort(indRates, decreasing = TRUE)[1:ntop]

    } else {
      NULL
    }

    MRList(
      list(
        MR=MR,
        main=main,
        C = C,
        S = S
      )
    )

  })
output$sensitivity <- renderPlot({
  req(MRList())

  # Extract global lists
  for (n in names(MRList()))
    assign(n,rlist::list.extract(MRList(),n))
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  SASpecies  = isolate(input$SASpecies)

  if(input$SAPlotType == "scatterplot" &
     SASpecies  != ""           ) {
    # Scatterplots
    par(mfrow = c(4,5),
        cex = 0.7*cex, cex.main = 0.4*cex, mar = mar,
        mgp = mgp, tcl = tcl, pty = pty, lwd = lwd)
    isp = which(colnames(C) == SASpecies)[1]
    y = C[,isp]
    namesMR = names(MR)
    for (i in 1:length(MR)) {
      irp = which(colnames(S) == namesMR[i])[1]
      x = S[,irp]
      plot(
        x, y,
        xlab= 'log10(k)', xlim = range(x),
        ylab= paste0('log10[',SASpecies,']'), ylim = range(y),
        type ='p', pch=16,
        col = ifelse(cor(x,y) > 0, col_tr2[5], col_tr2[3]),
        main = paste0(namesMR[i])
      )
      a = signif(MR[i],2)
      legend(
        ifelse(cor(x,y) > 0, 'topleft', 'topright'),
        legend = '', bty='n',
        title = a,
        title.col = cols[2]
      )
    }

  } else {
    # Barplot
    par(mfrow = c(1,1),
        cex = 0.75*cex, cex.main = 0.8*cex, mar = c(3,20,2,1),
        mgp = mgp, tcl = tcl, pty = pty, lwd = lwd ,lend=2)

    MR = rev(MR) # To get largest at top of graph

    colbr = rep(col_tr2[5],length(MR))
    colbr[MR<0] = col_tr2[3]

    xlim = c(0,1.2*max(MR))
    if(sum(MR<0) != 0 )
      xlim = c(-1,1)

    barplot(MR,
            horiz     = TRUE, las=1,
            xlim      = xlim,
            beside    = FALSE,
            col       = colbr,
            border    = NA,
            main      = ''
    )
    grid(ny=0); box()
    mtext(main,side = 3, line = 0)
  }

})
# Sanity ####
rangesSanityInteg <- reactiveValues(x = NULL, y = NULL)
output$sanityInteg <- renderPlot({
  if(is.null(concList()))
    return(NULL)

  if(is.null(statsList()))
    statsList(getIntegStats())

  # Extract stats
  for (n in names(statsList()))
    assign(n,rlist::list.extract(statsList(),n))

  # Extract graphical params
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  sel = rep(TRUE, nStat)
  sel[3] = FALSE # do not show NSTEPS (=NACCPT+NREJCT)

  # Define zoom range
  if (is.null(rangesSanityInteg$x)) {
    xlim <- range(time)
  } else {
    xlim <- rangesSanityInteg$x
  }
  if (is.null(rangesSanityInteg$y)) {
    ylim = range(c(0, 1.2*statSup[, sel]), na.rm = TRUE)
  } else {
    ylim <- rangesSanityInteg$y
  }

  # Generate colors per stat
  colors = assignColorsToSpecies(
    input$colSel, statNames, sel, nf=1,
    cols, col_tr, col_tr2)

  # Show uncertainty bands ?
  showBands = (nfStat > 1)

  # Plot
  par(mfrow = c(1, 1),
      cex = cex, cex.main = cex, mar = mar,
      mgp = mgp, tcl = tcl, pty = pty, lwd = lwd)

  plot(time, time,
       type = 'n',
       log  = 'x',
       main = ifelse(!showBands,
                     'Mean values',
                     'Mean and 95 % Proba. Intervals'),
       xlab = 'Time / s',
       xlim = xlim,
       xaxs = 'i',
       ylab = 'Integrator statistics',
       ylim = ylim,
       yaxs = 'i')
  grid()


  # CI
  if(showBands) {
    for(isp in (1:nStat)[sel])
      polygon(c(time     , rev(time      )),
              c(statLow[,isp],rev(statSup[,isp])),
              col = colors$col_trSp[isp],
              border = NA)
  }

  # Mean concentrations
  matplot(time, statMean[,sel], type='l', lty = 1,
          col = colors$colsSp[sel],
          lwd = 1.2*lwd, add = TRUE)

  legend(
    'topleft', bty = 'n',
    legend = statNames[sel],
    col    = if(showBands) colors$col_trSp[sel] else colors$colsSp[sel],
    lty    = 1,
    lwd    = ifelse(showBands,20,4)
  )
  box()

})
outputOptions(output, "sanityInteg",
              suspendWhenHidden = FALSE)
observeEvent(
  input$sanity_dblclick,
  {
    brush <- input$sanity_brush
    if (!is.null(brush)) {
      rangesSanityInteg$x <- c(brush$xmin, brush$xmax)
      rangesSanityInteg$y <- c(brush$ymin, brush$ymax)
    } else {
      rangesSanityInteg$x <- NULL
      rangesSanityInteg$y <- NULL
    }
  }
)
# Spectral radius
rangesSanitySR <- reactiveValues(x = NULL, y = NULL)
output$sanitySR <- renderPlot({
  if(is.null(concList()))
    return(invisible())

  if(is.null(statsList()))
    statsList(getIntegStats())

  # Extract stats
  for (n in names(statsList()))
    assign(n,rlist::list.extract(statsList(),n))

  # Extract graphical params
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  sel = which( statNames == 'SPRAD')
  if(length(sel) == 0) {
    id = shiny::showNotification(
      h4('SPRAD unavailable'),
      closeButton = TRUE,
      # duration = 5,
      type = 'warning'
    )
    # on.exit(removeNotification(id), add = TRUE)
    return(invisible())
  }

  # Define zoom range
  if (is.null(rangesSanitySR$x)) {
    xlim <- range(time)
  } else {
    xlim <- rangesSanitySR$x
  }
  if (is.null(rangesSanitySR$y)) {
    ylim = range(c(0, 1.2*statSup[, sel]), na.rm = TRUE)
  } else {
    ylim <- rangesSanitySR$y
  }

  # Generate colors per stat
  colors = assignColorsToSpecies(
    input$colSel, statNames, sel, nf=1,
    cols, col_tr, col_tr2)

  # Show uncertainty bands ?
  showBands = (nfStat > 1)

  # Plot
  par(mfrow = c(1, 1),
      cex = cex, cex.main = cex, mar = mar,
      mgp = mgp, tcl = tcl, pty = pty, lwd = lwd)

  plot(time, time,
       type = 'n',
       log  = 'x',
       main = ifelse(!showBands,
                     'Mean value(s)',
                     'Mean and 95 % Proba. Interval(s)'),
       xlab = 'Time / s',
       xlim = xlim,
       xaxs = 'i',
       ylab = 'Spectral Radius',
       ylim = ylim,
       yaxs = 'i')
  grid()


  # CI
  if(showBands) {
    for(isp in (1:nStat)[sel])
      polygon(c(time     , rev(time      )),
              c(statLow[,isp],rev(statSup[,isp])),
              col = colors$col_trSp[isp],
              border = NA)
  }

  # Mean concentrations
  matplot(time, statMean[,sel], type='l', lty = 1,
          col = colors$colsSp[sel],
          lwd = 1.2*lwd, add = TRUE)

  legend(
    'topleft', bty = 'n',
    legend = statNames[sel],
    col    = if(showBands) colors$col_trSp[sel] else colors$colsSp[sel],
    lty    = 1,
    lwd    = ifelse(showBands,20,4)
  )
  box()

})
outputOptions(output, "sanitySR",
              suspendWhenHidden = FALSE)
observeEvent(
  input$sanitySR_dblclick,
  {
    brush <- input$sanitySR_brush
    if (!is.null(brush)) {
      rangesSanitySR$x <- c(brush$xmin, brush$xmax)
      rangesSanitySR$y <- c(brush$ymin, brush$ymax)
    } else {
      rangesSanitySR$x <- NULL
      rangesSanitySR$y <- NULL
    }
  }
)

output$sanityOutputs <- renderPrint({
  if(is.null(concList())   |
     is.null(ratesList())
  )
    return(NULL)

  # Extract data from lists
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))
  for (n in names(ratesList()))
    assign(n,rlist::list.extract(ratesList(),n))

  nMC = nf

  # Analyze final concentrations
  conc = conc[, nt, ]
  if(nMC == 1) {
    conc = matrix(conc,nrow=1)
    lconc = ifelse(conc <= 0, NA, log10(conc))
    sdc  = rep(1,ncol(conc))
  } else {
    lconc = ifelse(conc <= 0, NA, log10(conc))
    sdc  = apply(lconc, 2, function(x) sd(x,na.rm=TRUE))
  }
  sel  = which(sdc == 0 | !is.finite(sdc))


  if(length(sel) != 0) {
    cat(' Species with suspicious final concentrations\n',
        '---------------------------------------------\n\n')

    sd0 = nzero = ninf = c()
    for (ii in 1:length(sel)) {
      i = sel[ii]
      sd0[ii]   = ifelse(is.finite(sdc[i]),sdc[i] == 0,FALSE)
      nzero[ii] = paste0(sum(conc[,i] <= 0),'/',nMC)
      ninf[ii]  = paste0(sum(!is.finite(lconc[,i])),'/',nMC)
    }
    df = data.frame(
      'Name'  = species[sel],
      'Var=0' = sd0,
      'Nzero' = nzero,
      'Ninf'  = ninf
    )
    print(df)
  }

  sdc = apply(photoRates,2,function(x) sd(x))
  sel = which(sdc==0 & !is.finite(sdc))
  if(length(sel) != 0) {

    cat('\n\n')
    cat(' Photorates with suspicious samples\n',
        '-----------------------------------\n\n')
    sd0 = nzero = ninf = c()
    for (ii in 1:length(sel)) {
      i = sel[ii]
      sd0[ii]   = sdc[i] == 0
      nzero[ii] = paste0(sum(photoRates[,i] == 0),'/',nMC)
      ninf[ii]  = paste0(sum(!is.finite(photoRates[,i])),'/',nMC)
    }
    df = data.frame(
      'Name'  = colnames(photoRates)[sel],
      'Var=0' = sd0,
      'Nzero' = nzero,
      'Ninf'  = ninf
    )
    print(df)
  }


  sdc = apply(rates,2,function(x) sd(x))
  sel = which(sdc==0 & !is.finite(sdc))

  if(length(sel) != 0) {
    cat('\n\n')
    cat(' Reaction rates with suspicious samples\n',
        '---------------------------------------\n\n')

    sd0 = nzero = ninf = c()
    for (ii in 1:length(sel)) {
      i = sel[ii]
      sd0[ii]   = sdc[i] == 0
      nzero[ii] = paste0(sum(rates[,i] == 0),'/',nMC)
      ninf[ii]  = paste0(sum(!is.finite(rates[,i])),'/',nMC)
    }
    df = data.frame(
      'Name'  = colnames(rates)[sel],
      'Var=0' = sd0,
      'Nzero' = nzero,
      'Ninf'  = ninf
    )
    print(df)
  }

})
