# Set global variables

elements    = c("H", "C", "N", "O", "S", "Ar")
numElecElem = c( 1 ,  6 ,  7 ,  8 , 16 ,  18 )
massElem    = mass(elements)
spDummy     = c('E', 'He',
                "SOOT",'SOOTC','SOOTN','SOOTC+','SOOTS','SOOTO',
                'CxHy+','CxHy','CxHyNz+','CxHyNz',
                'C4H2X','HC3NX','HC5NX',
                'Products')

# Functions
filterFormula <- function(sp) {
  # Filter out non-stoechimetric notations
  # to enable composition calculation in get.atoms()
  sp1 = sub("^l-", "", sub("^c-", "", sub("^t-", "", sp)))
  sp1 = sub("^l", "", sub("^c", "", sub("^t", "", sub("^i", "", sp1))))
  sp1 = sub("\\(.*\\)", "", sp1)
  sp1 = sub("^X", "", sp1)
  sp1 = sub("^1", "", sp1)
  sp1 = sub("^3", "", sp1)
  sp1 = sub("C3P", "C", sp1)
  sp1 = sub("N4S", "N", sp1)
  sp1 = sub("N2D", "N", sp1)
  sp1 = sub("N1D", "N", sp1)
  sp1 = sub("N2P", "N", sp1)
  sp1 = sub("N3P", "N", sp1)
  sp1 = sub("O1S", "O", sp1)
  sp1 = sub("O1D", "O", sp1)
  sp1 = sub("O3P", "O", sp1)
  sp1 = sub("S1D", "S", sp1)
  sp1 = sub("S3P", "S", sp1)
  sp1 = sub("\\+", "", sp1)
  return(sp1)
}
nbHeavyAtoms <- function(spList) {
  # Number of non-hydrogen atoms in formula
  compo = t(sapply(spList, get.atoms))
  nHeavy = rowSums(compo) - compo[, 1]
  return(nHeavy)
}
calcAtoms <- function(formula) {
  # Transform filtered formula to composition vector
  atoms = i2A(formula)
  compo = matrix(0, nrow = 1, ncol = length(elements))
  colnames(compo) = elements
  compo[1, colnames(atoms)] = atoms
  return(compo[1, ])
}
get.atoms <- function(sp) {
  # Transform formula to composition vector
  sp1 = filterFormula(sp)
  tryCatch(calcAtoms(sp1), error = function(x) rep(NA, length(elements)))
}
massFormula <- function(sto) {
  # Calculate mass of composition
  sum(sto * massElem)
}
numElec <- function(sto) {
  # Calculate electron number of composition
  sum(sto * numElecElem)
}
spCharge = function(species) {
  charge = rep(0,length(species))
  ions   = grepl("\\+$",species)
  charge[ions] = 1
  if(any(species == 'E'))
    charge[which(species == 'E')] = -1
  return(charge)
}
