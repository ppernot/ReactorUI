# Function  ####
ascii_only <- function(file) {
  response <- what_ascii(file)

  if (length(response) > 0) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

what_ascii <- function(file) {
  response <- capture.output(tools::showNonASCIIfile(file))
  return(response)
}
buildVal <- function(textLine, lineNum, blckName) {
  #-----function appends nml list with new values-----
  # remove all text after comment string
  textLine <- strsplit(textLine, "!")[[1]][1]
  # Trim white spaces
  textLine <- trimws(textLine)
  # Remove trailing comma
  if (substr(textLine, nchar(textLine), nchar(textLine)) == ",") {
    textLine <- substr(textLine, 1, nchar(textLine) - 1)
  }

  if (!any(grep("=", textLine))) {
    stop(c(
      "no hanging lines allowed in .nml, used ", textLine,
      ".\nSee line number:", lineNum, ' in "&', blckName, '" section.'
    ))
  }
  params <- strsplit(textLine, "=") # break text at "="
  parNm <- params[[1]][1]
  parVl <- params[[1]][2]
  # figure out what parval is...if string, remove quotes and keep as string
  # ***for boolean text, use "indentical" so that 0!= FALSE
  # can be: string, number, comma-sep-numbers, or boolean

  # special case for date:
  if (is.na(parVl)) {
    stop('Empty values after "', textLine, '" on line ', lineNum,
      ". \nPerhaps the values are on the next line?",
      call. = FALSE
    )
  }
  if (nchar(parVl > 17) &
    substr(parVl, 14, 14) == ":" &
    substr(parVl, 17, 17) == ":") {
    parVl <-
      paste(c(
        substr(parVl, 1, 11), " ",
        substr(parVl, 12, nchar(parVl))
      ),
      collapse = ""
      )
  }
  if (any(grep("'", parVl))) {
    parVl <- gsub("'", "", parVl)
  } else if (any(grep("\"", parVl))) {
    parVl <- gsub("\"", "", parVl)
  } else if (isTRUE(grepl(".true.", parVl) || grepl(".false.", parVl))) {
    logicals <- unlist(strsplit(parVl, ","))
    parVl <- from.glm_boolean(logicals)
  } else if (any(grep(",", parVl))) { # comma-sep-nums
    parVl <- c(as.numeric(unlist(strsplit(parVl, ","))))
  } else { # test for number
    parVl <- as.numeric(parVl)
  }
  lineVal <- list(parVl)
  names(lineVal) <- parNm
  return(lineVal)
}
#' go from glm2.nml logical vectors to R logicals
#'
#' @param values a vector of strings containing either .false. or .true.
#' @return a logical vector
#' @keywords internal
from.glm_boolean <- function(values) {
  logicals <- sapply(values, FUN = function(x) {
    if (!isTRUE(grepl(".true.", x) || grepl(".false.", x))) {
      stop(x, " is not a .true. or .false.; conversion to TRUE or FALSE failed.",
        call. = FALSE
      )
    }
    return(ifelse(isTRUE(grepl(".true.", x)), TRUE, FALSE))
  })
  return(as.logical(logicals))
}
to.glm_boolean <- function(values) {
  val.logical <- values
  values[val.logical] <- ".true."
  values[!val.logical] <- ".false."
  return(values)
}
.nml <- function(list_obj) {
  nml <- list_obj
  class(nml) <- "nml"
  invisible(nml)
}
read_nml <- function(nml_file) {
  if (!ascii_only(nml_file)) {
    stop("non-ASCII characters found in nml file on line ", what_ascii(nml_file))
  }
  # skip all commented lines, return all variables and associated values
  # requires NO return line variables
  # (all variables must be completely defined on a single line)

  c <- file(nml_file, "r")
  fileLines <- readLines(c)
  close(c)

  lineStart <- substr(fileLines, 1, 1)
  # ignore comment lines or empty lines
  ignoreLn <- lineStart == "!" | fileLines == ""
  lineStart <- lineStart[!ignoreLn]
  fileLines <- fileLines[!ignoreLn]
  # find all lines which start with "&" * requires FIRST char to be value


  lineIdx <- seq(1, length(lineStart))
  blckOpen <- lineIdx[lineStart == "&"]
  blckClse <- lineIdx[lineStart == "/"]

  nml <- list()
  for (i in 1:length(blckOpen)) {
    blckName <- substr(fileLines[blckOpen[i]], 2, nchar(fileLines[blckOpen[i]]))
    blckName <- gsub("\\s", "", blckName)
    oldNms <- names(nml)
    nml[[i]] <- list()
    names(nml) <- c(oldNms, blckName)

    carryover <- ""
    for (j in (blckOpen[i] + 1):(blckClse[i] - 1)) {
      textLine <- paste0(carryover, gsub("\t", "", gsub(" ", "", fileLines[j])))
      if (substr(textLine, 1, 1) != "!") {

        # Add a check here, sometimes, if there is a hanging comma,
        # and only sometimes that means add next row
        if (substr(textLine, nchar(textLine), nchar(textLine)) == "," &&
          j + 1 <= length(fileLines) &&
          !any(grep("=", fileLines[j + 1])) &&
          !any(grep("/", fileLines[j + 1]))) {
          carryover <- textLine
          next
        } else {
          carryover <- ""
        }

        # else, line is commented out
        lineVal <- buildVal(textLine, lineNum = j, blckName)
        nml[[i]] <- c(nml[[i]], lineVal)
      }
    }
  }
  nml <- .nml(nml)
  return(nml)
}


# Interactive ####

reacData <- reactiveVal()

reacProc <- reactiveValues()
output$contentsNmlMsg <- renderPrint({
  if (is.null(input$projectDir)) {
    return(NULL)
  }

  # Count MC outputs
  mcf <- list.files(
    path = paste0(ctrlPars$projectDir, "/Run"),
    pattern = "control",
    full.names = TRUE
  )
  nmc <- length(mcf)

  if (nmc == 0) {
    cat("The project has not control file, using defaults !")
    source("R/default_ctrlData.R")
  } else {
    ctrlList <- read_nml(mcf[1])
  }

  # Update reactive value
  reacData(ctrlList$REAC_DATA)
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

output$chemistryParams <- renderPrint({
  if (is.null(reacData())) {
    return(NULL)
  }

  listPars <- c(
    "reactantsSpecies",
    "reactantsComposition"
  )

  for (n in names(reacData()))
    assign(n, rlist::list.extract(reacData(), n))

  for (n in listPars)
    cat(n, " = ", unlist(mget(n, ifnotfound = NA)), "\n")
})
