# Options ####
Sys.setlocale(
  category = "LC_NUMERIC",
  locale = "C")
options(
  shiny.maxRequestSize = 20 * 1024^2,
  width = 80,
  warn = 0,
  stringsAsFactors = FALSE,
  max.print = 10000,
  shiny.autoload.r = NULL)

# options(shiny.json.digits=32)
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
  col_tr  = rev(inlmisc::GetColors(8, alpha = 0.1))[1:7],
  col_tr2 = rev(inlmisc::GetColors(8, alpha = 0.4))[1:7],
  pty     = 's',
  mar     = c(3.2, 3, 1.6, .2),
  mgp     = c(2, .75, 0),
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
  ## Store in local list
  for (n in names(lPars))
    ctrlPars[[n]] = lPars[[n]]
}

# Define base directory
roots = c(wd = '~')
if (!is.null(ctrlPars$projectDir) &
    is.character(ctrlPars$projectDir))
  roots = c(wd = dirname(ctrlPars$projectDir))

# Defaults for MC-ChemDB
DB_DATA_default = list(
  nMC = 100,
  photoVersion = 0,
  neutralsVersion = 0,
  ionsVersion = 0
)

REAC_DATA_default = list(
  runId                = 'apsis',
  ifRestart            = FALSE,
  debug                = TRUE,
  reactorLength        = 50,
  reactorSection       = 105,
  beamSpectrumFile     = 'surf73.txt',
  spectralResolution   = 1,
  spectrumRange        = c(50,200),
  beamIntensity        = 5e+16,
  beamSection          = 0.78,
  gasTemperature       = 300,
  electronsTemperature = 300,
  totalPressure        = 700,
  reactantsPressure    = 700,
  reactantsFlux        = 7,
  reactantsSpecies     = c('N2,CH4'),
  reactantsComposition = c(0.9,0.1),
  reactionTime         = 3600,
  nbSnapshots          = 100
)

REAC_DATA_Titan = list(
  runId                = 'Titan',
  ifRestart            = FALSE,
  debug                = TRUE,
  reactorLength        = 5e7,
  reactorSection       = 96,
  beamSpectrumFile     = 'stellarflux.txt',
  spectralResolution   = 1,
  spectrumRange        = c(50,200),
  beamIntensity        = 1.1e11,
  beamSection          = 96,
  gasTemperature       = 150,
  electronsTemperature = 150,
  totalPressure        = 4e-4,
  reactantsPressure    = 4e-4,
  reactantsFlux        = 4e-4,
  reactantsSpecies     = c('N2,CH4'),
  reactantsComposition = c(0.90, 0.10),
  reactionTime         = 2e8,
  nbSnapshots          = 110
)

listParsReac <- c(
  "reactorLength",
  "reactorSection",
  "gasTemperature",
  "electronsTemperature",
  "totalPressure",
  "reactantsPressure",
  "reactantsFlux",
  "reactionTime"
)

listParsReacUnits <- c(
  reactorLength = 'cm',
  reactorSection = 'cm^2',
  gasTemperature = 'K',
  electronsTemperature = 'K',
  totalPressure = 'Pa',
  reactantsPressure = 'Pa',
  reactantsFlux = 'sccm',
  reactionTime = 's'
)

spectrumData_default = list(
  beamSpectrumFileName =
    '../../ChemDBPublic/BeamSpectrumFiles/1nm/stellarflux.txt'
)
