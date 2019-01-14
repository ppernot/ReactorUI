# require('compiler',quietly=TRUE); enableJIT(3)
require(RColorBrewer,quietly=TRUE)
library('CHNOSZ')
data('thermo')
source('Scripts/graphFuncs.R')

# Start #####
PRINT=TRUE
MASS=FALSE

concThresh = -50
maxMass = 120

# Load conc sample files ####
files = list.files(pattern='fracmol_',path='MC_Output',full.names=TRUE)
nf = length(files)

# Species list in output files
line  = readLines(con=files[1], n=1)
spOut = scan(text=line, what=character(), 
             strip.white=TRUE, quiet=TRUE)[-1]
nSpOut= length(spOut)

# Get mole fractions at stationary state (final time)
x = read.table(files[1], header=TRUE, skip=0)
time=x[-1,1]
nt = length(time)
conc=matrix(NA,nrow=nf,ncol=nSpOut)
colnames(conc)=spOut
for(i in seq_along(files) ) {
  tab = as.vector(read.table(files[i], header=TRUE, 
                             skip=nt, colClasses='numeric'))
  tab = unlist(tab[-1])
  conc[i,] = tab   
}


# Load react. rates sample files ####
files = list.files(pattern='reacs_rates_',path='MC_Output',full.names=TRUE)
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
reacs = readLines('./Run/reac_list.dat')
colnames(rates) = reacs

# Load photorates sample files ####
files = list.files(pattern='photo_rates_',path='MC_Output',full.names=TRUE)
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
reacs = readLines('./Run/photo_list.dat')
colnames(photoRates)= reacs


# MC fluxes ####
allRates   = cbind(photoRates,rates)
reacNames  = colnames(allRates)
nbReacs    = ncol(allRates)
nf         = nrow(allRates)
flux       = matrix(NA,nrow=nf,ncol=nbReacs)
species    = spOut
nbSpecies  = nSpOut

## get D and L matrices
D = L = matrix(0,ncol=nbSpecies,nrow=nbReacs)
DL  = readLines(con='Run/photo_DL.dat',n=4)
nnz = as.numeric(DL[1])
Dsp   = as.numeric(unlist(strsplit(DL[2],split=' ')))
Dsp   = matrix(Dsp,nrow=nnz, ncol=3)
for (i in 1:nrow(Dsp))
  D[Dsp[i,1],Dsp[i,2]]=Dsp[i,3]
nnz = as.numeric(DL[3])
Lsp   = as.numeric(unlist(strsplit(DL[4],split=' ')))
Lsp   = matrix(Lsp,nrow=nnz, ncol=3)
for (i in 1:nrow(Lsp))
  L[Lsp[i,1],Lsp[i,2]]=Lsp[i,3]
nPh = ncol(photoRates)

DL  = readLines(con='Run/reac_DL.dat',n=4)
nnz = as.numeric(DL[1])
Dsp   = as.numeric(unlist(strsplit(DL[2],split=' ')))
Dsp   = matrix(Dsp,nrow=nnz, ncol=3)
for (i in 1:nrow(Dsp))
  D[nPh+Dsp[i,1],Dsp[i,2]]=Dsp[i,3]
nnz = as.numeric(DL[3])
Lsp   = as.numeric(unlist(strsplit(DL[4],split=' ')))
Lsp   = matrix(Lsp,nrow=nnz, ncol=3)
for (i in 1:nrow(Lsp))
  L[nPh+Lsp[i,1],Lsp[i,2]]=Lsp[i,3]

R = D + L

for (i in 1:nf)
  flux[i,] = allRates[i,] * apply(L,1,function(x) prod(conc[i,]^x))

logFlux = ifelse(flux==0, NA, log10(flux))
flMean = apply(logFlux,2,function(x) 10^mean(x,na.rm=TRUE))
flMed  = apply(logFlux,2,function(x) 10^median(x,na.rm=TRUE))
flF    = apply(logFlux,2,function(x) 10^sd(x,na.rm=TRUE))
names(flMean) = reacNames
names(flMed)  = reacNames
names(flF)    = reacNames

# Results #####
# Remove useless species
spDel= c('HV','dummy','E','CxHyNz+','Z+','SOOT','C4H2X')
sel  = ! species %in% spDel
spec = species[sel]
LR   = L[,sel]
RR   = R[,sel]


# for ( sp1 in spec )
  # sp1='N2'
  viewFlow('CH4',LR,RR,spec,
           reacNames,
           reacType=c(rep(1,ncol(photoRates)),rep(2,ncol(rates))),
           reacTypeNames=c("Ph","R"),
           flMean,
           topShow=0.05,
           level=1,
           showLegend=TRUE,PDF=TRUE)
  viewFlow('HCN',LR,RR,spec,
           reacNames,
           reacType=c(rep(1,ncol(photoRates)),rep(2,ncol(rates))),
           reacTypeNames=c("Ph","R"),
           flMean,
           topShow=0.1,
           level=1,
           showLegend=TRUE,PDF=TRUE)
  
  
  
  # budget('CH4' ,LR,RR,spec,reacNames,flMean,weightThresh=0.001)

budget('CH4',LR,RR,spec,reacNames,flMean,weightThresh=0.01)

traceBack('C2H2',LR,RR,spec,reacNames,flMean)
 
traceBack('HCN',LR,RR,spec,reacNames,flMean) 
traceBack('H2CN',LR,RR,spec,reacNames,flMean) 
traceBack('C4H8',LR,RR,spec,reacNames,flMean)
