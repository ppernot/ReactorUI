# enableBookmarking("server")
source("R/packages.R")
source("R/functions.R")
source("R/stoeFuns.R")
source("R/graphFuncs.R")

# Proportions of Side/Main Panels ####
sideWidth  <- 3
mainWidth  <- 12 - sideWidth
plotHeight <- 700

# Graphical parameters ####
gPars = list(
  cols    = rev(inlmisc::GetColors(8))[1:7],
  col_tr  = rev(inlmisc::GetColors(8,alpha=0.1))[1:7],
  col_tr2 = rev(inlmisc::GetColors(8,alpha=0.4))[1:7],
  pty     = 's',
  mar     = c(3.2,3,1.6,.2),
  mgp     = c(2,.75,0),
  tcl     = -0.5,
  lwd     = 2,
  cex     = 1,
  cex.leg = 0.8
)

# Control Parameters ####
ctrlPars = list(
  projectDir  = NULL
)

# Restore state
ctrlFile = 'ctrlParams.yaml'
if (file.exists(ctrlFile)) {
  ## Get from file
  lPars = rlist::list.load(ctrlFile)
  ## Expose parameters
  for (n in names(lPars))
    ctrlPars[[n]]= lPars[[n]]
}

# Define base directory
roots = c(wd='~')
if(!is.null(ctrlPars$projectDir)&
   is.character(ctrlPars$projectDir))
  roots = c( wd = dirname(ctrlPars$projectDir) )

