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


# Filter out undesirable values
conc = ifelse(conc==0, NA, log10(conc))
conc = ifelse(conc<=concThresh,NA,conc)
sdc = apply(conc,2,function(x) sd(x))
selC = sdc!=0 & !is.na(sdc)

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

rates = ifelse(rates==0, NA, log10(rates))

sdc = apply(rates,2,function(x) sd(x))
selR = sdc!=0 & !is.na(sdc)

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

photoRates = ifelse(photoRates==0, NA, log10(photoRates))

sdPR = apply(photoRates,2,function(x) sd(x))
selPR = sdPR!=0 & !is.na(sdPR)

# Correlation Analysis ####
S = cbind(photoRates[,selPR],rates[,selR])
cMat = cor(conc[,selC],S)

ncolors=256
library('fields')
mypalette=two.colors(ncolors,start='blue',middle='white',end='red')
hM=heatmap(cMat,symm=FALSE,col=mypalette,
           breaks=seq(-1,1,by=2/ncolors),
           cexRow=0.5,cexCol=0.5,margins=c(4,4),
           Rowv=NULL,Colv=NULL,scale="none",
           # labRow=species,labCol=reacTags,
           keep.dendro=FALSE)

# # Get major processes ####
indRates=colSums(cMat^2)

## All processes
MR = sort(indRates,decreasing = TRUE)[1:20]
nn = names(MR)
print(data.frame(nn,MR),row.names=FALSE)

## HSIC SA ### PARALLELISER !!!!!!!!!!!!!
C = conc[,selC]
hMat = matrix(NA,nrow=ncol(S),ncol=ncol(C))
syy = sxx = c()
for (i in 1:ncol(S)) {
  y = S[,i]
  syy[i] = dHSIC::dhsic(y,y)$dHSIC
}
for(j in 1:ncol(C)) {
  x = C[,j]
  sxx[j] = dHSIC::dhsic(x,x)$dHSIC
}
for(j in 1:ncol(C)) {
  x = C[,j]
  for (i in 1:ncol(S)) {
    y = S[,i]
    hMat[i,j] = dHSIC::dhsic(x,y)$dHSIC /sqrt( sxx[j] * syy[i] )
  }
  sel = !is.finite(hMat[,j])
  hMat[sel,j] = 0
  hMat[,j] = hMat[,j] / sum(hMat[,j])
}
colnames(hMat) = colnames(C)
rownames(hMat) = colnames(S)

hM=heatmap(hMat,symm=FALSE,
           # col=mypalette,
           # breaks=seq(-1,1,by=2/ncolors),
           cexRow=0.5,cexCol=0.5,margins=c(4,4),
           Rowv=NA,Colv=NA,scale="none",
           # labRow=species,labCol=reacTags,
           keep.dendro=FALSE)

indRates = rowSums(hMat,na.rm = TRUE)
names(indRates) = colnames(S)
sel = is.finite(indRates)
MR = sort(indRates[sel],decreasing = TRUE)[1:20]
nn = names(MR)
print(data.frame(nn,MR),row.names=FALSE)

# CH4
head(sort(hMat[,1],decreasing = TRUE))









# ## Photo processes
MR = sort(indRates[1:sum(selPR)],decreasing=TRUE)
nn = names(MR)
print(data.frame(nn,MR),row.names=FALSE)

## Contributions to CH4
# Compute Fluxes #####

# MC fluxes
reacNames = c(colnames(photoRates),colnames(rates))
nbReac     = ncol(rates)+ncol(photoRates)
nf         = nrow(rates)
flux       = matrix(NA,nrow=nbReac,ncol=nf)
species    = spOut
nbSpecies  = nSpOut
D = L = matrix(0,ncol=nbSpecies,nrow=nbReac)

# get D and L matrices
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


rates=cbind(photoRates,rates)
for (i in 1:nf)
  flux[,i] = rates[,i] * apply(L,1,function(x) prod(conc[,i]^x))

logFlux = ifelse(flux==0, NA, log10(flux))
flMean = apply(logFlux,1,function(x) 10^mean(x,na.rm=TRUE))
flMed  = apply(logFlux,1,function(x) 10^median(x,na.rm=TRUE))
flF    = apply(logFlux,1,function(x) 10^sd(x,na.rm=TRUE))
names(flMean)=reacNames
names(flMed)=reacNames
names(flF)  =reacNames
sort(flMean,decreasing=TRUE)[1:10]
sort(flMed,decreasing=TRUE)[1:10]

# save(rMean,rF,reacType,reacTags,L,R,flMean,flF,
#      file=paste('reacStats_',case,'.Rda',sep=''))

# Build graph ######
# Remove useless species
spDel= c('HV','dummy','E','CxHyNz+','Z+','SOOT','C4H2X')
sel  = ! species %in% spDel
spec = species[sel]
LR   = L[,sel]
RR   = D[,sel] + L[,sel]
# yMean=yMean[spec]
# logConc=logConc[spec,]
# save.image(file='networks.Rda')


# for ( sp1 in spec )
  # sp1='N2'
  # viewFlow(sp1,LR,RR,spec,
  #          reacNames,
  #          reacType=c(rep(1,ncol(photoRates)),rep(2,ncol(rates))),
  #          reacTypeNames=c("Ph","R"),
  #          flMean,
  #          topShow=0.1,level=1,
  #          showLegend=TRUE,PDF=TRUE)
  #

  budget('CH4' ,LR,RR,spec,reacNames,flMean,weightThresh=0.001)
  budget('C2H2',LR,RR,spec,reacNames,flMean,weightThresh=0.01)

  traceBack('HCN',LR,RR,spec,reacNames,flMean)


