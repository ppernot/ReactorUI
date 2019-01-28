#=================================================#
# Generate Monte Carlo chemistry data for reactor #
#=================================================#

# Load Stoechiomerty/Mass functions
source('./stoeFuns.R')

# Go to the main directory of project (one rank above 'Scripts')
setwd("..")

# Monte Carlo parameters (only for reactions at the moment)
nMC = 100

# Target directory for chemistry samples
targetMCDir = './MC_Input/'

# Source directories for reference chemistry data
reacSourceDirs = c('../ChemDBPublic/Neutrals/',
                   '../ChemDBPublic/Ions/')

# Scheme selection by initial composition of mixture
spInit = c("N2","CH4","CO2")   # Initial species

# Elements in initial species
elInit = names(CHNOSZ::count.elements(paste0(spInit,collapse='')))

# Misc functions ###########################################
lsp = function(sp) {
  lr = which(lapply(reactants,function(x) sp %in% x)==1)
  if(length(lr)!=0) {
    cat(sp,'\n')
    cat('As reactant:',lr,'\n')
    for(i in lr)
      print(paste0(paste(reactants[[i]],collapse=' + '),
                   ' >>> File: ',orig[[i]]))

    lp = which(lapply(products,function(x) sp %in% x)==1)
    if(length(lp)!=0) {
      cat('As product:',lp,'\n')
      for(i in lp)
        print(paste0(paste(reactants[[i]],collapse=' + '),
                     ' >>> File: ',orig[[i]]))
    }
    cat('\n')
  }
}

# Kinetic Parser ############################################################
nbReac=0
reactants = products = params = type = orig = locnum = list()

filename='../ChemDBPublic/Photoprocs/PhotoScheme.dat'
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
}

filename='../ChemDBPublic/Photoprocs/PhotoIonScheme.dat'
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
}

for (sourceDir in reacSourceDirs) {
  filename=paste0(sourceDir,'run_0000.csv')
  scheme  = read.csv(file=filename,header=FALSE,sep=';')
  scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
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
  }
}

# Build species list from reactants and products
species=levels(as.factor(unlist(c(reactants,products))))
nbSpecies=length(species)
print("Species List:")
print(species)

# Stoechiometry ######
compo = t(apply(as.matrix(species,ncol=1),1,get.atoms))
colnames(compo)=elements
mass  = apply(compo,1,massFormula)
names(mass) = species

## Attribute mass to species with undefined masse
dummySpecies = c('CxHy','CxHy+','CxHyNz','CxHyNz+',
                 'Products','E','SOOTC','SOOTC+',
                 'C4H2X','HC3NX','HC5NX','SOOTS','SOOTN')
dummyMass   = round(max(mass,na.rm=TRUE)+2)
mass[dummySpecies] = dummyMass

orderByMass=order(mass)
species1=species[orderByMass]
mass1=mass[orderByMass]

sink(file='speciesByMass.txt')
for (m in 1:round(max(mass1))) {
  sel = round(mass1) == m
  if(sum(sel)!=0) cat(m,species1[sel],'\n')
}
sink()

# Full R, L, D matrices #####
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
colnames(D)=species

# Volpert analysis #####
# Build species list and reactions list
# from initial species list by
# iterations over the reaction network
reacs=1:nbReac
spInit = spInit[order(mass[spInit])]
lReacs   = rowSums(L[,spInit] == 1) & type=='photo'
reacList = reacs[lReacs]
spProds  = species[colSums(R[lReacs,])!=0]
spProds  = spProds[! (spProds %in% spInit)]
spProds  = spProds[order(mass[spProds])]
while ( length(spProds)!=0 ) {
  spInit   = c(spInit,spProds)
  lReacsU  = rowSums(L[,spInit] == 1) & type=='photo'
  lReacsB  = rowSums(L[,! colnames(L) %in% spInit])==0 & type!='photo'
  lReacs   = lReacsU | lReacsB
  spProds  = species[colSums(R[lReacs,])!=0]
  spProds  = spProds[! (spProds %in% spInit)]
  spProds  = spProds[order(mass[spProds])]
}

# Orphan species: species without parents in scheme
## For info only. At the moment, they are not rejected from scheme
reject = apply(compo, 1, function(x)
  sum(x[elInit],na.rm=TRUE) != sum(x,na.rm=TRUE) )
orphans = species[!species %in% spInit &
                  !reject &
                  !species %in% dummySpecies ]
cat('Orphan Species:',orphans,'\n')
sink(file='orphanSpecies.txt')
for(sp in orphans) lsp(sp)
sink()

# Reduce reactions and species list to accessible species #####
species   = spInit
nbSpecies = length(species)
nbReac    = sum(lReacs)
L = L[lReacs,species]
D = D[lReacs,species]
params=params[lReacs]
type = type[lReacs]
locnum = locnum[lReacs]
orig = orig[lReacs]
mass = mass[species]

# Generate data files for reactor code #####
# Split photodissociations and reactions
photo = type == 'photo'
Dphoto = matrix(D[ photo,],nrow=sum(photo),ncol=nbSpecies)
D      = matrix(D[!photo,],nrow=sum(!photo),ncol=nbSpecies)
Lphoto = matrix(L[ photo,],nrow=sum(photo),ncol=nbSpecies)
L      = matrix(L[!photo,],nrow=sum(!photo),ncol=nbSpecies)

# Treat Reactions #####
# 1/ Convert D & L to triplets (sparse matrices)
Dsp = cbind(which(D!=0,arr.ind=TRUE),D[D!=0])
Lsp = cbind(which(L!=0,arr.ind=TRUE),L[L!=0])

sink(file ='./Run/reac_DL.dat')
cat(nrow(Dsp))
cat('\n')
cat(Dsp)
cat('\n')
cat(nrow(Lsp))
cat('\n')
cat(Lsp)
sink()

# 2/ generate parameters files

# Single output file from nominal data
sink(file ='./Run/reac_params.dat')
for (i in (1:nbReac)[!photo])
    cat(params[[i]],'\n')
sink()

if(nMC > 1) {
  # Clean target dir
  command = paste0('rm -rf ',targetMCDir,'Reactions/*')
  dummy   = system(command)

  # MC samples for reduced reactions list
  sel = !photo
  origs = unique(unlist(orig[sel]))
  for (iMC in 0:nMC) {
    sampleFile=paste0('run_',sprintf('%04i',iMC),'.csv')
    sink(file =paste0(targetMCDir,'Reactions/',sampleFile))
    for (io in origs) {
      sourceDir = io
      filename=paste0(sourceDir,sampleFile)
      scheme  = read.csv(file=filename,header=FALSE,sep=';')
      scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
      paramsLoc = list()
      ii = 0
      for (i in 1:nrow(scheme)) {
        ii = ii + 1
        terms=scheme[i,8:ncol(scheme)]
        paramsLoc[[ii]]    = terms[!is.na(terms) & terms!=""]
      }
      selOrig = which (orig == io)
      paramsLoc = paramsLoc[unlist(locnum[selOrig])]
      for (j in 1:length(paramsLoc)) cat(paramsLoc[[j]],'\n')
    }
    sink()
  }
}


# Treat Nominal Photo-processes #####

if (dim(Dphoto)[1]!=0) {
  # Convert D & L to triplets (sparse matrices)
  Dsp = cbind(which(Dphoto!=0,arr.ind=TRUE),Dphoto[Dphoto!=0])
  Lsp = cbind(which(Lphoto!=0,arr.ind=TRUE),Lphoto[Lphoto!=0])

  sink(file ='./Run/photo_DL.dat')
  cat(nrow(Dsp))
  cat('\n')
  cat(Dsp)
  cat('\n')
  cat(nrow(Lsp))
  cat('\n')
  cat(Lsp)
  sink()

  sink(file ='./Run/photo_params.dat')
  for (i in (1:nbReac)[photo]) {
    sp=species[which(Lphoto[i,]!=0)]
    if( params[[i]][1] == 0 )
      cat(paste0('se',sp,'.dat'),
          '\n')
    else
      cat(paste0('se',sp,'.dat'),
          paste0('qy',sp,'_',params[[i]][1],'.dat'),
          '\n')
  }
  sink()
}

sink(file='./Run/species_aux.dat')
cat(nbSpecies)
cat('\n')
cat(species)
cat('\n')
cat(mass)
cat('\n')
cat(rep(0,nbSpecies))
cat('\n')
sink()
