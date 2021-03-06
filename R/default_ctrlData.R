ctrlList = list(
  runId                = "apsis",
  ifRestart            = FALSE,
  debug                = FALSE,
  reactorLength        = 50,
  reactorSection       = 105,
  beamSpectrumFile     = "surf73.txt",
  spectrumRange        = c(50,200),
  beamIntensity        = 5e+16,
  beamSection          = 0.78,
  gasTemperature       = 300,
  electronsTemperature = 300,
  totalPressure        = 700,
  reactantsPressure    = 700,
  reactantsFlux        = 7,
  reactantsSpecies     = c("N2","CH4"),
  reactantsComposition = c(0.9,0.1),
  reactionTime         = 3600,
  nbSnapshots          = 100
)
