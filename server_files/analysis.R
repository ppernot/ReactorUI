# Functions ####
getConc  = function(concThresh = -50) {
  # Load 1 sample files
  files = list.files(
    path=paste0(ctrlPars$projectDir,'/MC_Output'),
    pattern="fracmol_",
    full.names=TRUE
  )

  x = read.table(files[1], header=TRUE, skip=0)
  time=x[-1,1]

  # Species list
  line  = readLines(con=files[1], n=1)
  species = scan(text=line, what=character(),
                 strip.white=TRUE, quiet=TRUE)[-1]

  iNorm= which(species=="CH4")

  # Get all data
  nf   = length(files)
  nsp  = length(species)
  nt   = length(time)
  conc = moleFrac = array(0,dim=c(nf,nt,nsp))
  for(i in seq_along(files[1:nf]) ) {
    tab  = read.table(files[i], header=TRUE, skip=0)
    time = tab[-1,1]
    y    = as.matrix(tab[-1,-1])
    conc[i,1:nt,1:nsp]     = y[1:nt,1:nsp]
    moleFrac[i,1:nt,1:nsp] = y[1:nt,1:nsp] /y[1:nt,iNorm]
  }

  # Pretreat moleFrac for plots
  moleFrac = ifelse(moleFrac==0         , NA, log10(moleFrac))
  moleFrac = ifelse(moleFrac<=concThresh, NA, moleFrac)

  yMean = apply(moleFrac,c(2,3),
                function(x) mean(x,na.rm=TRUE))
  ySd   = apply(moleFrac,c(2,3),
                function(x) sd(x,na.rm=TRUE))
  yLow95= apply(moleFrac,c(2,3),
                function(x) quantile(x,probs = 0.025,na.rm=TRUE))
  ySup95= apply(moleFrac,c(2,3),
                function(x) quantile(x,probs = 0.975,na.rm=TRUE))

  return(
    list(
      nf       = nf,
      nt       = nt,
      nsp      = nsp,
      time     = time,
      conc     = conc,
      species  = species,
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
  x = read.table(files[1], header=FALSE, skip=0)
  nRates=length(x)
  rates=matrix(NA,nrow=nf,ncol=nRates)
  for(i in seq_along(files) ) {
    tab = read.table(files[i], header=FALSE,
                     skip=0, colClasses='numeric')
    rates[i,] = t(tab)
  }
  ## Get reac names
  reacs = readLines(
    paste0(ctrlPars$projectDir,'/Run/reac_list.dat'))
  colnames(rates) = reacs

  # Load photorates sample files
  files = list.files(
    path=paste0(ctrlPars$projectDir,'/MC_Output'),
    pattern='photo_rates_',
    full.names=TRUE)
  nf = length(files)

  ## Get MC photo rates
  x = read.table(files[1], header=FALSE, skip=0)
  nPhotoRates=length(x)
  photoRates=matrix(NA,nrow=nf,ncol=nPhotoRates)
  for(i in seq_along(files) ) {
    tab = read.table(files[i], header=FALSE,
                     skip=0,
                     colClasses=c('numeric')
    )
    photoRates[i,] = t(tab)
  }
  ## Get reac names
  reacs = readLines(
    paste0(ctrlPars$projectDir,'/Run/photo_list.dat'))
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
selectSpecies <- function(species, categs) {
  # Select species to plot from categsPlot

  compo   = t(apply(as.matrix(species,ncol=1),1,get.atoms))
  colnames(compo)=elements

  # Remove dummy species
  sel0 = ! species %in% spDummy

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
assignColorsToSpecies <- function(colSel, species, sel, nf,
                                  cols, col_tr, col_tr2,
                                  threshTransp = 50) {
  # Manage transparency wrt nb MC runs
  if(nf <= threshTransp)
    col_tr = col_tr2 # Darker

  if(colSel) {
    # Assign colors to species to avoid changes
    nrep = ceiling(length(species)/length(cols))
    colsSp   = rep(cols  ,times = nrep)
    col_trSp = rep(col_tr,times = nrep)
  } else {
    # Use sequential colors for selection
    colsSp = col_trSp = rep(NA,length(species))
    nrep = ceiling(length(species[sel])/length(cols))
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
    # Reduced Indices
    V = V / sqrt(sxx[j]*syy)
    # Normalize
    V[!is.finite(V)] = 0
    V = V / sum(V)
    hMat[,j] = V
  }

  return(hMat)
}
dcorMat <- function(C,S) {
  if(is.vector(C))
    C = matrix(C,ncol=1)
  if(is.vector(S))
    S = matrix(S,ncol=1)

  nC = ncol(C); nS = ncol(S)
  hMat = matrix(NA,nrow=nS,ncol=nC)
  for(j in 1:nC) {
    x = C[,j]
    V = future.apply::future_apply(
      X = S, MARGIN = 2,
      FUN = function(y,x) fda.usc::dcor.xy(x,y,test=FALSE),
      x = x
    )
    # Normalize
    V[!is.finite(V)] = 0
    V = V / sum(V)
    hMat[,j] = V
  }

  return(hMat)
}

testAna = function() {
  concList = getConc()
  for (n in names(concList))
    assign(n,rlist::list.extract(concList,n))

  ratesList = getRates()
  for (n in names(ratesList))
    assign(n,rlist::list.extract(ratesList,n))

  # Tests
  conc = conc[,nt,]

  # Filter out undesirable values
  conc = ifelse(conc==0, NA, log10(conc))
  sdc  = apply(conc,2,function(x) sd(x))
  selC = sdc!=0 & is.finite(sdc)

  rates = ifelse(rates==0, NA, log10(rates))
  sdc = apply(rates,2,function(x) sd(x))
  selR = sdc!=0 & is.finite(sdc)

  photoRates = ifelse(photoRates==0, NA, log10(photoRates))
  sdc = apply(photoRates,2,function(x) sd(x))
  selPR = sdc!=0 & is.finite(sdc)

  C = conc[,selC]
  S = cbind(photoRates[,selPR],rates[,selR])

  SASpecies = 'CH5+'
  isp = which(species[selC] == SASpecies)
  # indRates = as.vector(cor(C[,isp], S, method = "spearman")^2)
  # indRates = as.vector(hsicMat(C[,isp],S))
  indRates = as.vector(dcorMat(C[,isp],S))
  names(indRates) = colnames(S)
  indx = order(indRates,decreasing=TRUE)[1:9]
  print(indRates[indx])
  par(mfrow=c(3,3))
  for(i in indx) {
    plot(
      C[,isp],S[,i],
      main = names(indRates)[i])
    legend(
      'topleft', bty='n',
      legend = '',
      title = signif(indRates[i],2),
      title.col = 4
    )
  }

}
# Interactive ####

concList <- reactiveVal()
observeEvent(
  input$loadMC,
  future({getConc()}) %...>% concList()
)

ratesList <- reactiveVal()
observeEvent(
  input$loadMC,
  future({getRates()}) %...>% ratesList()
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

rangesKinetics <- reactiveValues(x = NULL, y = NULL)
output$kinetics <- renderPlot({
  if(is.null(concList()))
    return(NULL)

  # Extract conc et al.
  for (n in names(concList()))
    assign(n,rlist::list.extract(concList(),n))

  # Species list to plot
  sel = selectSpecies(species, input$categsPlot)

  # Define zoom range
  if (is.null(rangesKinetics$x)) {
    xlim <- range(time)*10 # expand for labels
  } else {
    xlim <- rangesKinetics$x
  }
  if (is.null(rangesKinetics$y)) {
    ylim <- c(
      max(-20, min(mfLow[,sel])),
      min(0  , max(mfSup[,sel]))
    )
  } else {
    ylim <- rangesKinetics$y
  }

  # Extract graphical params
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  # Generate colors per species
  colors = assignColorsToSpecies(
    input$colSel, species, sel, nf,
    cols, col_tr, col_tr2)

  # Plot
  par(mfrow = c(1, 1),
      cex = cex, cex.main = cex, mar = mar,
      mgp = mgp, tcl = tcl, pty = pty, lwd = lwd)
  plot(time, time, type='n', log='x',
       main = ifelse(!input$mcPlot, 'Mean values',
                     'Mean and 95 % Proba. Intervals'),
       xlab = 'Time / s'     , xlim = xlim,
       ylab = 'Log10(mole fraction)', ylim = ylim )
  # CI
  if(input$mcPlot) {
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
  text(x = time[nt], y= mfMean[nt,sel],
       labels = species[sel],
       col=colors$colsSp[sel],
       pos = 4, offset=0.2, cex=cex.leg)
  grid()
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

# Sensitivity ####
MRList <- reactiveVal()
observeEvent(
  input$doSA,
  {
    if(is.null(concList()) |
       is.null(ratesList()) )
      return(NULL)

    MRList(NULL)

    # Extract conc et al.
    for (n in names(concList()))
      assign(n,rlist::list.extract(concList(),n))
    for (n in names(ratesList()))
      assign(n,rlist::list.extract(ratesList(),n))

    # Concentrations at stationnary state (last snapshot)
    conc = conc[,nt,]
    colnames(conc) = species

    # Filter out undesirable values
    conc = ifelse(conc==0, NA, log10(conc))
    sdc  = apply(conc,2,function(x) sd(x))
    selC = sdc!=0 & is.finite(sdc)

    rates = ifelse(rates==0, NA, log10(rates))
    sdc = apply(rates,2,function(x) sd(x))
    selR = sdc!=0 & is.finite(sdc)

    photoRates = ifelse(photoRates==0, NA, log10(photoRates))
    sdc = apply(photoRates,2,function(x) sd(x))
    selPR = sdc!=0 & is.finite(sdc)

    C = conc[,selC]
    colnames(C) = species[selC]
    S = cbind(photoRates[,selPR],rates[,selR])

    SASpecies = input$SASpecies
    if(SASpecies != "") # Check validity
      validate(
        need(SASpecies %in% species[selC], 'Invalid species !')
      )

    if(input$anaType == 'spearman') {
      # CORR
      if(SASpecies == "") {
        main = 'TOP20 most influent reactions : Sum(RCC^2)'
        indRates = colSums(cor(C, S, method = "spearman")^2)
        names(indRates) = colnames(S)
        MR = sort(indRates,decreasing = TRUE)[1:20]
      } else {
        main = paste0(
          SASpecies,
          ' / TOP20 most correlated reactions : RCC')
        isp = which(species[selC]==SASpecies)
        indRates = as.vector(cor(C[,isp],S, method = "spearman"))
        names(indRates) = colnames(S)
        indx = order(abs(indRates),decreasing = TRUE)[1:20]
        MR   = indRates[indx]
      }

    } else if (input$anaType == 'dcorr') {
      # fda.usc::dcor.xy
      if(SASpecies == "") {
        main = 'TOP20 most influent reactions : Sum(dCor)'
        indRates = rowSums(dcorMat(C,S),na.rm = TRUE)
      } else {
        main = paste0(
          SASpecies,
          ' / TOP20 most influent reactions : dCor')
        isp = which(species == SASpecies)
        indRates = as.vector(dcorMat(C[,isp],S))
      }
      names(indRates) = colnames(S)
      MR = sort(indRates,decreasing = TRUE)[1:20]

    } else if (input$anaType == 'hsic') {
      # HSIC
      if(SASpecies == "") {
        main = 'TOP20 most influent reactions : Sum(HSIC)'
        indRates = rowSums(hsicMat(C,S),na.rm = TRUE)
      } else {
        main = paste0(
          SASpecies,
          ' / TOP20 most influent reactions : HSIC')
        isp = which(species == SASpecies)
        indRates = as.vector(hsicMat(C[,isp],S))
      }
      names(indRates) = colnames(S)
      MR = sort(indRates,decreasing = TRUE)[1:20]
    } else {
      NULL
    }

    MRList(list(MR=MR, main=main, C = C, S = S))
    # print('MRList updated !')
  })
output$sensitivity <- renderPlot({
  if(is.null(MRList()))
    return(NULL)

  # Extract data
  for (n in names(MRList()))
    assign(n,rlist::list.extract(MRList(),n))
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  SASpecies  = isolate(input$SASpecies)

  if(input$SAPlotType == "scatterplot" &
     SASpecies  != ""           ) {
    # Scatterplots
    par(mfrow = c(4,5),
        cex = cex, cex.main = 0.6*cex, mar = mar,
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
        type ='p', pch=20, col=cols[5],
        main = paste0(namesMR[i])
      )
      legend(
        'topleft', legend = '', bty='n',
        title = signif(MR[i],2),
        title.col = cols[2]
      )
    }

  } else {
    # Barplot
    par(mfrow = c(1,1),
        cex = cex, cex.main = cex, mar = c(3,20,2,1),
        mgp = mgp, tcl = tcl, pty = pty, lwd = lwd ,lend=2)
    MR = rev(MR)
    colbr = rep(cols[5],length(MR))
    colbr[MR<0] = cols[3]
    barplot(MR,
            horiz     = TRUE, las=1,
            beside    = FALSE,
            col       = colbr,
            names.arg = names(MR),
            main      = main
    )

  }

})

