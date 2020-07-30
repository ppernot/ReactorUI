library(RColorBrewer) # Temporary

col2tr =function(x,alpha=80)
  rgb(unlist(t(col2rgb(x))),alpha=alpha,maxColorValue = 255)

# Misc network plotting funcs

hiveFlowGraph <- function (L,R,species,nbReac,reacType,
                           alpha=50,PDF=FALSE) {
  require('HiveR')
  require('grid')

  vplayout <- function(x,y,bkgnd="black"){
    viewport(layout.pos.row=x,
             layout.pos.col=y,
             gp = gpar(fill=bkgnd))
  }
  compo = t(sapply(species,get.atoms))
  mass  = apply(compo,1,massCxHyNz)

  bkgnd = "black"
  if(PDF) bkgnd = "white"
  transp=col2tr(bkgnd,0)
  cols=brewer.pal(9,"Set1")
  edgeColor= c(cols[3], cols[4] )
  spColor=c(cols[1],cols[2],cols[5])
  reacColor=c(cols[6:9])

  # Losses
  h1 = adj2HPD(L,desc="Hive Network",type="2D",axis.cols = bkgnd)

  # Define axes
  h1$nodes$size=0.5

  # Reacs
  sel = grepl("^R",h1$nodes$lab)
  h1$nodes$axis  [sel] = as.integer(1)
  h1$nodes$radius[sel] = 1:nbReac
  h1$nodes$color [sel] = reacColor[reacType]
  h1$nodes$size  [sel] = 0.25

  # Species
  sel = !sel
  h1$nodes$axis  [sel] = as.integer(2)
  h1$nodes$radius[sel] = as.numeric(mass[species])
  h1$nodes$color [sel] = spColor[1]

  # Ions: discriminate pos and neg
  sel = grepl("\\+$",h1$nodes$lab)
  h1$nodes$axis  [sel] = as.integer(3)
  h1$nodes$color [sel] = spColor[2]

  h1$edges$color=col2tr(edgeColor[1],alpha)

  # Productions
  h2 = adj2HPD(R, desc = "Hive Network",type="2D", axis.cols = bkgnd)
  h2$edges$color=col2tr(edgeColor[2],alpha)

  #Merge prod and loss edges
  h1$edges=rbind(h1$edges,h2$edges)

  # Rescale axes
  method="scale"
  action=c(1/nbReac,1/max(unlist(mass)),1/max(unlist(mass)))
  h1 = manipAxis(h1, method=method, action=action)

  axLab.gpar = gpar(col=c(do.call(rgb,as.list(1-col2rgb(bkgnd)/255)),
                          spColor[1:2]),fontsize = 16)
  axLabs=c("Reactions","Neutrals / mass","Ions / mass")

  if(PDF) pdf(file="hivePlotLRTitan.pdf",width=12,height=12)
  plotHive(h1, ch=0.05, bkgnd = bkgnd,axLabs=axLabs, axLab.gpar=axLab.gpar)
  if(PDF) dev.off()
}
netGraph <- function (links, species, PDF=FALSE) {
  # Display graph
  g = graph.adjacency(linksThresh,
                      mode = "directed",
                      weighted = TRUE)
  g = simplify(g)
  V(g)$label = spec
  V(g)$shape="circle"
  str=graph.strength(g,mode="in")-graph.strength(g,mode="out")
  wid=abs(str)
  wid=wid^0.1
  wid=0.01+(wid-min(wid))/(max(wid)-min(wid))
  V(g)$size = wid*10
  V(g)$color = ifelse(str>0,'gold','green')
  V(g)$label.cex = 2
  V(g)$label.color = "red"
  V(g)$label.font = 1


  wid=E(g)$weight
  wid=wid^0.1
  wid=0.01+(wid-min(wid))/(max(wid)-min(wid))
  E(g)$width = wid*20
  E(g)$color = col2tr('blue',wid*80)
  E(g)$arrow.size = 1
  E(g)$curved = TRUE

  loners=which(degree(g)==0)
  g=delete.vertices(g,loners)

  set.seed(23451)
  if (PDF) pdf(file="graphFlow.pdf",width=30,height=30)
  plot(g, layout=layout.graphopt(g,niter=10000,charge=0.5,spring.length=1))
  if (PDF) dev.off()
}
# Target graph functions
spCoord= function(i,count,species,range=50) {
  if(length(species)==0)
    return(c(0,0))
  else
  {radius = 0.5 + count
   dangle = pi/144
   angle  = pi + mass[[species[i]]]*dangle*ifelse(grepl("\\+$",species[i]),-1,1)
   radius*c(cos(angle),sin(angle))
  }
}
targetPlot = function(g,layout,species,reacTypeNames) {

  labels=V(g)$label
  V(g)$label = ""
  V(g)$label.degree = -pi/8
  V(g)$label.dist = 0.4
  V(g)$label[which(labels == "N2")]="N2"
  V(g)$label[which(labels == "CH4")]="CH4"

  #   V(g)$color[which(grepl("N",labels))]="green"
  V(g)$size = 2
  #   V(g)$size[which(grepl("N",labels))]=1
  #   V(g)$shape[which(grepl("N",labels))]="square"

  scale=1/6
  plot(g, layout=scale*layout, rescale=FALSE)

  rmin=scale*1.5
  rmax=scale*6.5

  for (i in -8:8) lines(x=c(rmin*cos(pi-i*pi/8),rmax*cos(pi-i*pi/8)),
                        y=c(rmin*sin(pi-i*pi/8),rmax*sin(pi-i*pi/8)),
                        col="brown",lwd=4,lty=2)
  for (i in c(0,8)) lines(x=c(rmin*cos(pi-i*pi/8),rmax*cos(pi-i*pi/8)),
                          y=c(rmin*sin(pi-i*pi/8),rmax*sin(pi-i*pi/8)),
                          col="brown",lwd=4,lty=1)

  rtext=rmax+0.05
  for (i in 0:8) text(x=rtext*cos(pi-i*pi/8),
                      y=rtext*sin(pi-i*pi/8),
                      labels= i*144/8,
                      col="red",cex=3,
                      srt=90-i*180/8)
  for (i in 1) text(x=1.1*rtext*cos(pi-i*pi/8),
                    y=1.1*rtext*sin(pi-i*pi/8),
                    labels= "mass",
                    col="red",cex=4,
                    srt=90-i*180/8)
  draw.arc(0, 0, 1.1*rtext, deg1=164, deg2=170, col="red",lwd=3)
  draw.arc(0, 0, 1.1*rtext, deg1=140, deg2=151, col="red",lwd=3)
  arrows(x0=1.1*rtext*cos(140*pi/180),
         y0=1.1*rtext*sin(140*pi/180),
         x1=1.1*rtext*cos(138*pi/180),
         y1=1.1*rtext*sin(138*pi/180),
         col="red",lwd=3
  )

  for (i in 1:7) text(x=rtext*cos(pi-i*pi/8),
                      y=-rtext*sin(pi-i*pi/8),
                      labels= i*144/8,
                      col="blue",cex=3,
                      srt=-90+i*180/8)
  for (i in 1:1) text(x=1.1*rtext*cos(pi-i*pi/8),
                      y=-1.1*rtext*sin(pi-i*pi/8),
                      labels= "mass",
                      col="blue",cex=4,
                      srt=-90+i*180/8)
  draw.arc(0, 0, 1.1*rtext, deg1=-164, deg2=-170, col="blue",lwd=3)
  draw.arc(0, 0, 1.1*rtext, deg1=-140, deg2=-151, col="blue",lwd=3)
  arrows(x0=1.1*rtext*cos(-140*pi/180),
         y0=1.1*rtext*sin(-140*pi/180),
         x1=1.1*rtext*cos(-138*pi/180),
         y1=1.1*rtext*sin(-138*pi/180),
         col="blue",lwd=3
  )

  for (i in 1:6) draw.circle(0,0,radius = scale*(0.5+i),
                             border="gray50",lwd=4)
  for (i in 1:6) text(x=scale*(0.65+i)*cos(-pi/16),
                      y=scale*(0.65+i)*sin(-pi/16),
                      labels=i-1,
                      col="gray50",
                      cex=3)
  text(x=scale*7.0*cos(-pi/16),
       y=scale*6.7*sin(-pi/16),
       labels="Generations",
       col="gray50",
       cex=3,pos=4)

  text(x=scale*(-7.5),y=scale*5,label="Ions",col="red",cex=6)
  text(x=scale*(-7.5),y=scale*(-5),label="Neutrals",col="blue",cex=6)
  #   text(x=scale*(7),y=scale*5,label="(Azotés)",col="green",cex=6)

  plot(g, layout=scale*layout, rescale=FALSE, add=TRUE)
  V(g)$label=labels
  legend('topright',legend=reacTypeNames,cex=3,bty='n',
         lty=1,lwd=4,col=1:length(reacTypeNames))
}
targetGraph = function (links,L,R,species,reacTags,
                        reacTypeNames,PDF=FALSE) {
  g = graph.adjacency(linksThresh, mode = "directed", weighted=TRUE)
  g = simplify(g)

  nedges=length(species)
  compo = t(sapply(species,get.atoms))
  mass  = apply(compo,1,massCxHyNz)

  V(g)$label = species
  V(g)$color = rep("blue",nedges)
  V(g)$shape = rep("circle",nedges)
  ions = grepl("\\+$",species)
  V(g)$color[ions]="red"

  V(g)$size = 0
  V(g)$label.cex = 2.0
  V(g)$label.color = V(g)$color
  V(g)$label.font = 1
  V(g)$frame.color = NA

  wid=E(g)$weight
  wid=wid^0.1  #log10(wid)
  wid=0.01+(wid-min(wid))/(max(wid)-min(wid))
  #   hist(wid)

  E(g)$color = col2tr('brown',wid*80)
  E(g)$width = wid*20
  E(g)$arrow.size = 0
  E(g)$curved = FALSE

  # _ 2_0 - Create the world
  labels = V(g)$label
  layout = matrix(NA,ncol=2,nrow=length(labels))

  if(PDF) pdf(file="graphGrowth.pdf",width=30,height=20)

  count = 1
  spInit = c("N2","CH4")
  spInit = spInit[order(mass[spInit])]
  for (i in 1:length(spInit))
    layout[which(labels==spInit[i]),] = spCoord(i,count,species=spInit)

  # png(file=paste('graphGrowth-',count,'.png',sep=''),width=1800,height=1200)
  targetPlot(g,layout,species,reacTypeNames)
  # dev.off()

  lReacs   = rowSums(L[,spInit] == 1) & reacType==1
  reacList = reacTags[lReacs]
  spProds  = species[colSums(R[lReacs,])!=0]
  spProds  = spProds[! (spProds %in% spInit)]
  spProds  = spProds[order(mass[spProds])]

  while ( count <= 7 ) {
    count=count+1
    for (i in 1:length(spProds))
      layout[which(labels==spProds[i]),] = spCoord(i,count,spProds)

    #   png(file=paste('graphGrowth-',count,'.png',sep=''),width=1800,height=1200)
    targetPlot(g,layout,spec,reacTypeNames)
    #   dev.off()

    spInit   = c(spInit,spProds)
    lReacsU  = (reacType==1 | reacType==4) & rowSums(L[,spInit] != 0 )
    lReacsB  = (reacType==2 | reacType==3) &
      rowSums(L[,! colnames(L) %in% spInit])==0
    lReacs   = lReacsU | lReacsB
    spProds  = species[colSums(R[lReacs,])!=0]
    spProds  = spProds[! (spProds %in% spInit)]
    spProds  = spProds[order(mass[spProds])]
  }

  if(PDF) dev.off()
}
viewFlowOld = function(sp1,L,R,species,reacs,reacType,reacTypeNames,
                    flMean,topShow=0.5,level=1,showLegend=TRUE,PDF=FALSE) {

  nbSpecies=length(species)
  nbReacs=length(reacs)

  nedges=nbReacs+nbSpecies
  links=matrix(0,ncol=nedges,nrow=nedges)
  dimnames(links)[[1]]=c(reacs,species)
  dimnames(links)[[2]]=c(reacs,species)

  n1 = which (species == sp1)
  for (n2 in 1:nbSpecies) {
    sp2=species[n2]
    for (m in 1:nbReacs) {
      if(L[m,n1]!=0 & L[m,n2]!=0) {
        links[nbReacs+n1,m]=links[nbReacs+n1,m]-flMean[m]
        links[nbReacs+n2,m]=links[nbReacs+n2,m]-flMean[m]
      }
      if(L[m,n1]!=0 & R[m,n2]!=0) {
        links[nbReacs+n1,m]=links[nbReacs+n1,m]-flMean[m]
        links[m,nbReacs+n2]=links[m,nbReacs+n2]-flMean[m]*R[m,n2]
      }
      if(R[m,n1]!=0 & L[m,n2]!=0) {
        links[m,nbReacs+n1]=links[m,nbReacs+n1]+flMean[m]*R[m,n1]
        links[nbReacs+n2,m]=links[nbReacs+n2,m]+flMean[m]
      }
      if(R[m,n1]!=0 & R[m,n2]!=0) {
        links[m,nbReacs+n1]=links[m,nbReacs+n1]+flMean[m]*R[m,n1]
        links[m,nbReacs+n2]=links[m,nbReacs+n2]+flMean[m]*R[m,n2]
      }
    }
  }

  if ( level==2 ) {
    select=colSums(links[,(nbReacs+1):nedges]) !=0 |
      rowSums(links[(nbReacs+1):nedges,]) !=0
    listSp = species[select]
    for (sp1 in listSp) {
      n1 = which (species == sp1)
      for (n2 in 1:nbSpecies) {
        sp2=species[n2]
        for (m in 1:nbReacs) {
          if(L[m,n1]!=0 & L[m,n2]!=0) {
            links[nbReacs+n1,m]=links[nbReacs+n1,m]-flMean[m]
            links[nbReacs+n2,m]=links[nbReacs+n2,m]-flMean[m]
          }
          if(L[m,n1]!=0 & R[m,n2]!=0) {
            links[nbReacs+n1,m]=links[nbReacs+n1,m]-flMean[m]
            links[m,nbReacs+n2]=links[m,nbReacs+n2]-flMean[m]*R[m,n2]
          }
          if(R[m,n1]!=0 & L[m,n2]!=0) {
            links[m,nbReacs+n1]=links[m,nbReacs+n1]+flMean[m]*R[m,n1]
            links[nbReacs+n2,m]=links[nbReacs+n2,m]+flMean[m]
          }
          if(R[m,n1]!=0 & R[m,n2]!=0) {
            links[m,nbReacs+n1]=links[m,nbReacs+n1]+flMean[m]*R[m,n1]
            links[m,nbReacs+n2]=links[m,nbReacs+n2]+flMean[m]*R[m,n2]
          }
        }
      }
    }
  }

  linksThresh=links
  if(sum(links!=0) > 50) {
    lnkThresh=quantile(abs(links)[abs(links)>0],probs=1-topShow,na.rm=TRUE)
    linksThresh[abs(links)<lnkThresh]=0
  }

  g = graph.adjacency(linksThresh, mode = "directed",weighted=TRUE)
  g = simplify(g)

  cols=brewer.pal(9,"Set3")
  reacColor=c(cols[6:9])

  reacLab=paste0('R',1:nbReacs)
  V(g)$label = c(reacLab,species)
  V(g)$shape=c(rep("rectangle",nbReacs),rep("circle",nbSpecies))
  V(g)$size = c(rep(10,nbReacs),rep(10,nbSpecies))
  V(g)$size2 = 20
  V(g)$color = c(reacColor[reacType],rep("gold",nbSpecies))
  V(g)$label.cex = c(rep(1,nbReacs),rep(1.5,nbSpecies))
  V(g)$label.color = "black"
  V(g)$label.font = 1

  wid=abs(E(g)$weight)
  wid=wid^0.1
  wid=0.01+(wid-min(wid))/(max(wid)-min(wid))
  E(g)$width = wid*20
  E(g)$color = ifelse(E(g)$weight>0,col2tr('red',80),col2tr('blue',80))
  E(g)$arrow.size = 3
  E(g)$curved = TRUE

  loners=which(degree(g)==0)
  g=delete.vertices(g,loners)

  set.seed(23451)
  if (PDF) pdf(file=paste0('Graphs/',sp1,'_graphFlow.pdf'),width=30,height=30)
  plot(g, layout=layout.auto(g))
  if(showLegend) legend('topright',legend=reacTypeNames,cex=3,bty='n',
                        lty=1,lwd=10,col=reacColor[1:length(reacTypeNames)])
  if (PDF) dev.off()

}

addLinks = function (sp1,
                     links,
                     species,
                     KL,
                     KR,
                     nbReacs,
                     nbSpecies,
                     wght = 1) {
  # Update flow digraph
  lSpecies = (nbReacs + 1):(nbReacs + nbSpecies)

  # Losses as negative links
  n1 = which (species == sp1)
  reacSel = which(KL[, n1] != 0)
  if (length(reacSel) != 0) {
    # Loss Reactions
    links[lSpecies, reacSel] =
      links[lSpecies, reacSel] - wght * t(KL)[1:nbSpecies, reacSel]
    # Products of loss reactions
    links[reacSel, lSpecies] =
      links[reacSel, lSpecies] - wght * KR[reacSel, 1:nbSpecies]
  }

  # Productions as positive links
  reacSel = which(KR[, n1] != 0)
  if (length(reacSel) != 0) {
    # Prod Reactions
    links[lSpecies, reacSel] =
      links[lSpecies, reacSel] + wght * t(KL)[1:nbSpecies, reacSel]
    # Reactants of prod. reacs
    links[reacSel, lSpecies] =
      links[reacSel, lSpecies] + wght * KR[reacSel, 1:nbSpecies]
  }
  return(links)
}

viewFlow = function(sp1,
                    L,
                    R,
                    species,
                    reacs,
                    reacType,
                    reacTypeNames,
                    flMean,
                    spInit,
                    topShow = 0.5,
                    level = 1,
                    showLegend = TRUE,
                    PDF = FALSE,
                    curved = FALSE) {
  # Builds digraph of fluxes to and from sp1

  nbSpecies = length(species)
  nbReacs   = length(reacs)

  nedges = nbReacs + nbSpecies
  lspecies = (nbReacs+1):nedges
  links = matrix(0, ncol = nedges, nrow = nedges)
  dimnames(links)[[1]] = c(reacs, species)
  dimnames(links)[[2]] = c(reacs, species)

  KL = L * matrix(flMean,
                  nrow = nbReacs,
                  ncol = nbSpecies,
                  byrow = FALSE)
  KR = R * matrix(flMean,
                  nrow = nbReacs,
                  ncol = nbSpecies,
                  byrow = FALSE)

  # Build link matrix with first neighbors (reactants and products)
  links = addLinks(sp1, links, species, KL, KR, nbReacs, nbSpecies)

  if (level==2) {
    # Add 2nd neighnors to link matrix ### DOES NOT WORK AS IS !!!
    select = colSums(links[, lspecies]) != 0 |
             rowSums(links[lspecies, ]) != 0
    listSp = species[select]
    for (sp2 in listSp)
      if (!(sp2 %in% spInit))
        links = addLinks(sp2, links, species, KL, KR, nbReacs, nbSpecies,
                         wght = 1e-6)
  }

  # Select most important links if there are too many
  linksThresh = links
  if (sum(links != 0) > 50) {
    lnkThresh = quantile(abs(links)[abs(links) > 0],
                         probs = 1 - topShow,
                         na.rm = TRUE)
    linksThresh[abs(links) < lnkThresh] = 0
  }

  g = graph.adjacency(linksThresh, mode = "directed", weighted=TRUE)
  g = simplify(g)

  cols=brewer.pal(9,"Set3")
  reacColor=c(cols[6:9])

  reacLab=paste0('R',1:nbReacs)
  V(g)$label = c(reacLab,species)
  V(g)$shape=c(rep("rectangle",nbReacs),rep("circle",nbSpecies))
  V(g)$size = c(rep(20,nbReacs),rep(20,nbSpecies))
  V(g)$size2 = 6
  V(g)$color = c(reacColor[reacType],rep("gold",nbSpecies))
  V(g)$label.cex = c(rep(1,nbReacs),rep(1,nbSpecies))
  V(g)$label.color = "black"
  V(g)$label.font = 2
  V(g)$type = c(rep(0,nbReacs),rep(1,nbSpecies))

  wid = abs(E(g)$weight)
  wid = wid^0.1 # Empirical transfo for better scaling...
  wid = 0.12 + (wid-min(wid))/(max(wid)-min(wid))
  E(g)$width = wid*5
  E(g)$color = ifelse(E(g)$weight>0,col2tr('red',160),col2tr('blue',160))
  E(g)$arrow.size  = 0.25*max(E(g)$width)
  E(g)$arrow.width = 0.25*max(E(g)$width)
  E(g)$curved = curved

  loners=which(degree(g)==0)
  g=delete.vertices(g,loners)

  return(g)
}

volpertGraph = function (linkMat,specList,vlpInd,topShow=0.5,PDF=FALSE,
                         curves=NULL) {

  # Calculate net flux between species (upperTri-lowerTri)
  nedges=ncol(linkMat)
  for (n1 in 1:(nedges-1) ) {
    for (n2 in (n1+1):nedges) {
      del=linkMat[n1,n2]-linkMat[n2,n1]
      if( del >=0) {
        linkMat[n1,n2]=del
        linkMat[n2,n1]=0
      } else {
        linkMat[n1,n2]=0
        linkMat[n2,n1]=-del
      }
    }
  }

  linksThresh=linkMat
  if(sum(linkMat!=0) > 50) {
    lnkThresh=quantile(abs(linkMat)[abs(linkMat)>0],probs=1-topShow,na.rm=TRUE)
    linksThresh[abs(linkMat)<lnkThresh]=0
  }



  nedges=length(specList)
  compo = t(sapply(specList,get.atoms))
  mass  = apply(compo,1,massCxHyNz)

  V(g)$label = specList
  V(g)$shape ='circle'
  V(g)$size = 1.8
  V(g)$color = 'gold'
  V(g)$label.cex = 1.8
  V(g)$label.color = 'black'
  V(g)$label.font = 1

  wid=E(g)$weight
  wid=wid^0.1
  wid=0.01+(wid-min(wid))/(max(wid)-min(wid))
  E(g)$width = wid*10
  E(g)$color = col2tr('blue',wid*80)
  E(g)$arrow.size = 2
  E(g)$curved = TRUE

  labels = V(g)$label
  layout = matrix(NA,ncol=2,nrow=length(labels))
  rownames(layout)=labels
  for ( vI in min(vlpInd):max(vlpInd) ) {
    spList = names(vlpInd[which(vlpInd==vI)])
    spList = labels[labels %in% spList]
    ions=grepl('\\+$',spList)
    if(sum(ions) != 0) {
      listNeus=spList[!ions]
      listNeus = listNeus[order(mass[listNeus])]
      listIons=spList[ions]
      listIons = listIons[order(mass[listIons])]
      spList=c(listNeus,listIons)
    } else {
      spList = spList[order(mass[spList])]
    }
    del=10/(length(spList)+1)
    for (sp in spList) {
      layout[sp,1] = vI
      layout[sp,2] = which(spList==sp)*del
    }
  }

  if(is.null(curves)) curves = seq(-0.5,0.5,len=nedges) #autocurve.edges(g,start=0.5)
  if(PDF) pdf(file="matrixFlowGraph.pdf",width=30,height=30)
  plot.igraph(g,layout=0.1*layout,edge.curved=curves)
  if(PDF) dev.off()
  curves
}

volpertLayout <- function (labels, vlpInd, mass, yMass=FALSE) {
  layout = matrix(NA,ncol=2,nrow=length(labels))
  rownames(layout)=labels
  if (yMass)
    rankMass=rank(mass,ties.method='first')
  for ( vI in min(vlpInd):max(vlpInd) ) {
    spL = names(vlpInd[which(vlpInd==vI)])
    spL = labels[labels %in% spL]
    ions= grepl('\\+$',spL)
    if(sum(ions) != 0) {
      listNeus = spL[!ions]
      listNeus = listNeus[order(mass[listNeus])]
      listIons = spL[ions]
      listIons = listIons[order(mass[listIons])]
      spL=c(listNeus,listIons)
    } else {
      spL = spL[order(mass[spL])]
    }
    del=10/(length(spL)+1)
    for (sp in spL) {
      layout[sp,1] = vI
      layout[sp,2] = ifelse(yMass,
                            rankMass[sp],
                            which(spL==sp)*del)
    }
  }
  return(layout)
}
volpertCircleLayout <- function (labels, vlpInd, mass) {
  layout = matrix(NA,ncol=2,nrow=length(labels))
  rownames(layout)=labels
  scale = 5
  y0  = 1 + scale*max(vlpInd)
  for ( vI in min(vlpInd):max(vlpInd) ) {
    spL = names(vlpInd[which(vlpInd==vI)])
    spL = labels[labels %in% spL]
    spL = spL[order(mass[spL])]
    R   = 1 + scale*vI
    del = 0
    if(length(spL)!=1)
      del = 0.5*pi/(length(spL)-1)
    for (i in 1:length(spL)) {
      p = rnorm(2) * scale / 6
      layout[spL[i], 1] =      R * cos((i - 1) * del) + p[1]
      layout[spL[i], 2] = y0 - R * sin((i - 1) * del) + p[2]
    }
  }
  return(layout)
}
animVolpertGraph <- function (spList,L,R,species,
                              reacs,reacType,reacTypeNames,
                              yMean,logConc,flMean,flux,vlpInd,
                              topShow=1,
                              fileTag=0, tag='',
                              fullConnectPlot=TRUE,
                              meanFluxPlot=TRUE,
                              animMeanFluxPlot=TRUE,
                              mcFluxPlot=TRUE,
                              output='screen',yMass=FALSE
                              ) {
  # Build hard net
  nbSpecies=length(species)
  nbReacs=length(reacs)

  compo = t(sapply(species,get.atoms))
  mass  = apply(compo,1,massCxHyNz)

  linksR=matrix(0,ncol=nbSpecies,nrow=nbSpecies)
  dimnames(linksR)[[1]]=species
  dimnames(linksR)[[2]]=species

  KL=L*matrix(flMean,nrow=nbReacs,ncol=nbSpecies,byrow=FALSE)
  KR=R*matrix(flMean,nrow=nbReacs,ncol=nbSpecies,byrow=FALSE)
  KD=KR-KL

  if(is.null(spList)) {
    for (n1 in 1:(nbSpecies-1) ) {
      for (n2 in (n1+1):nbSpecies) {
        reacSelLoss=which(L[,n1]*R[,n2]!=0)
        reacSelProd=which(L[,n2]*R[,n1]!=0)
        netFlow=sum(KD[reacSelLoss,n1])+sum(KD[reacSelProd,n1])
        if(netFlow<0)
          linksR[n1,n2]= -netFlow
        else
          linksR[n2,n1]= netFlow
      }
    }

  } else {
    for (sp in spList ) {
      n1 = which(species==sp)
      for (n2 in (1:nbSpecies)[-n1]) {
        reacSelLoss=which(L[,n1]*R[,n2]!=0)
        reacSelProd=which(L[,n2]*R[,n1]!=0)
        netFlow=sum(KD[reacSelLoss,n1])+sum(KD[reacSelProd,n1])
        if(netFlow<0)
          linksR[n1,n2]= -netFlow
        else
          linksR[n2,n1]= netFlow
      }
    }

  }

  refMat = linksR!=0

  g = graph.adjacency(linksR, mode = "directed", weighted=TRUE)
  size=2

  V(g)$label = species
  V(g)$shape ='circle'
  V(g)$size = size*1.7
  V(g)$color = 'brown'
  ions = grepl("\\+$",species)
  V(g)$color[ions]='darkgreen'
  V(g)$label.cex = size*1
  V(g)$label.color = ifelse(output=='screen','white','red')
  V(g)$label.font = 2
  V(g)$label.dist = size*2.5
  V(g)$label.degree = 0

  E(g)$width = 2
  E(g)$color = ifelse(output=='screen','gold','blue')
  E(g)$arrow.size = size*0.5
  E(g)$curved = TRUE

  labels = V(g)$label
  layout = volpertLayout(labels,vlpInd,mass,yMass=yMass)

  tab=as.matrix(get.edgelist(g))
  withIons=apply(tab,1,function(x) sum(grepl("\\+",x)) )

  if(fullConnectPlot) {
    # All connections with equal weight
    E(g)$color = ifelse(withIons==0,
                        col2tr(ifelse(output=='screen','gold','blue'),80),
                        col2tr('purple',80))

    outFile = paste0("./connectFlux.png")
    size=   200*nbSpecies^0.5 # 2000
    png(outFile, width=size, height=size,
        bg=ifelse(output!='screen','white','black'))
    vlp=vlpInd-1
    layout[,1]=(layout[,1]-1)*nbSpecies/max(vlp)
    par(cex.lab=1.5)
    plot(g,layout=layout,rescale=FALSE,
         axes=FALSE,
         xlim=c(0,nbSpecies),xlab='Volpert Index',
         ylim=c(0,nbSpecies),ylab='Molecular Mass Rank'
    )
    axis(side=1,at=as.numeric(levels(factor(vlp)))*nbSpecies/max(vlp),
         labels=levels(factor(vlp)))
    axis(side=2,at=c(0,nbSpecies),labels=FALSE,tick=TRUE)
    dev.off()
  }

  if(meanFluxPlot) {
    # Edges Weigted by mean flux &
    # Labels weighted by mean conc
    wid=log10(yMean[species])
    wThresh=-15
    wid[wid<wThresh]=wThresh
    wid['N2']=wThresh
    wid['CH4']=wThresh
    wid=0.01+(wid-min(wid))/(max(wid)-min(wid))
    wid[!is.finite(wid)]=0.01
    wid['N2']=1.01
    wid['CH4']=1.01
    V(g)$color = ifelse(ions,
                        col2tr('green',wid*255/1.01),
                        col2tr('red',wid*255/1.01))
    wid= E(g)$weight
    wThresh=quantile(wid[wid>0],probs=1-topShow,na.rm=TRUE)
    wid[wid<wThresh]=0
    wid=wid^0.2
    wid=0.01+(wid-min(wid))/(max(wid)-min(wid))
    E(g)$color = ifelse(withIons==0,
                        col2tr('gold',wid*255/1.01),
                        col2tr('purple',wid*255/1.01))

    outFile = paste0("./meanFlux.png")
    png(outFile, width=2000, height=2000, bg='black')
    plot(g,layout=layout)
    dev.off()
  }

  if(animMeanFluxPlot) {
    # Hide labels
    V(g)$label = ''

    # Edges Weigted by mean flux &
    # Labels weighted by mean conc
    wid=log10(yMean[species])
    wThresh=-20
    wid[wid<wThresh]=wThresh
    wid['N2']=wThresh
    wid['CH4']=wThresh
    wid=0.01+(wid-min(wid))/(max(wid)-min(wid))
    wid[!is.finite(wid)]=0.01
    wid['N2']=1.01
    wid['CH4']=1.01
    V(g)$color = ifelse(ions,
                        col2tr('green',wid*255/1.01),
                        col2tr('red',wid*255/1.01))
    wid= E(g)$weight
    wThresh=quantile(wid[wid>0],probs=1-topShow,na.rm=TRUE)
    wid[wid<wThresh]=0
    wid=wid^0.2
    wid=0.01+(wid-min(wid))/(max(wid)-min(wid))
    E(g)$color = ifelse(withIons==0,
                        col2tr('gold',wid*255/1.01),
                        col2tr('purple',wid*255/1.01))

    outFile = paste0("./Anim/meanFlux_",sprintf('%04i',fileTag),".png")
    png(outFile, width=1000, height=1000, bg='black')
    plot(g,layout=layout)
    title(main=tag,col.main='white')
    dev.off()
  }

  if(mcFluxPlot) {
    for (i in 1:nf) { #1:nf
      fl1=flux[,i]
      KL=L*matrix(fl1,nrow=nbReacs,ncol=nbSpecies,byrow=FALSE)
      KR=R*matrix(fl1,nrow=nbReacs,ncol=nbSpecies,byrow=FALSE)
      KD=KR-KL
      links=linksR*0
      if(is.null(spList)) {
        for (n1 in 1:(nbSpecies-1) ) {
          for (n2 in (n1+1):nbSpecies) {
            reacSelLoss=which(L[,n1]*R[,n2]!=0)
            reacSelProd=which(L[,n2]*R[,n1]!=0)
            netFlow=sum(KD[reacSelLoss,n1])+sum(KD[reacSelProd,n1])
            if(netFlow<0)
              links[n1,n2]= -netFlow
            else
              links[n2,n1]= netFlow
          }
        }

      } else {
        for (sp in spList ) {
          n1 = which(species==sp)
          for (n2 in (1:nbSpecies)[-n1]) {
            reacSelLoss=which(L[,n1]*R[,n2]!=0)
            reacSelProd=which(L[,n2]*R[,n1]!=0)
            netFlow=sum(KD[reacSelLoss,n1])+sum(KD[reacSelProd,n1])
            if(netFlow<0)
              links[n1,n2]= -netFlow
            else
              links[n2,n1]= netFlow
          }
        }

      }

      # Hide labels
      V(g)$label = ''

      # Weight vertices
      wid=logConc[species,i]
      wThresh=-15
      noShow=is.na(wid)
      wid[wid<wThresh]=wThresh
      wid[noShow]=wThresh
      wid['N2']=wThresh
      wid['CH4']=wThresh
      wid=0.1+(wid-min(wid))/(max(wid)-min(wid))
      wid[noShow]=0
      wid['N2']=1.1
      wid['CH4']=1.1
      V(g)$color = ifelse(ions,
                          col2tr('green',wid*255/1.1),
                          col2tr('red',wid*255/1.1))

      # Weight edges
      wid= as.vector(t(links)[t(refMat)])
      wThresh=quantile(wid[wid>0],probs=1-topShow,na.rm=TRUE)
      wid[wid<wThresh]=0
      wid=wid^0.2
      wid=(wid-min(wid))/(max(wid)-min(wid))
      E(g)$color = ifelse(withIons==0,
                          col2tr('gold',wid*255),
                          col2tr('purple',wid*255))

      outFile = paste0("./Anim/img_",sprintf('%03i',i),".png")
      png(outFile, width=1000, height=1000, bg='black')
      plot(g,layout=layout)
      dev.off()

    }
  }
}

viewFlowClust = function(clust,L,R,species,colSp,
                         reacs,reacType,reacTypeNames,colRe,
                         flMean,topShow=0.5,level=1,
                         showLegend=TRUE,keepAll=FALSE) {

  # Builds digraph of fluxes within species cluster

  nbSpecies=length(species)
  nbReacs=length(reacs)

  nedges=nbReacs+nbSpecies
  links=matrix(0,ncol=nedges,nrow=nedges)
  dimnames(links)[[1]]=c(reacs,species)
  dimnames(links)[[2]]=c(reacs,species)

  KL=L*matrix(flMean,nrow=nbReacs,ncol=nbSpecies,byrow=FALSE)
  KR=R*matrix(flMean,nrow=nbReacs,ncol=nbSpecies,byrow=FALSE)

  lSpecies=(nbReacs+1):(nbReacs+nbSpecies)
  for (isp1 in 1:(length(clust)-1) ) {
    n1=which(species==clust[isp1])
    for (isp2 in (isp1+1):length(clust) ) {
      n2=which(species==clust[isp2])

      reacSel = which(KL[,n1]*KR[,n2]!=0)
      if(length(reacSel)!=0) {
        links[lSpecies,reacSel]=
          links[lSpecies,reacSel]+t(KL)[,reacSel]
        links[reacSel,lSpecies]=
          links[reacSel,lSpecies]+KR[reacSel,]
      }
#       print(reacSel)

      reacSel=which(KL[,n2]*KR[,n1]!=0)
      if(length(reacSel)!=0) {
        links[lSpecies,reacSel]=
          links[lSpecies,reacSel]+t(KL)[,reacSel]
        links[reacSel,lSpecies]=
          links[reacSel,lSpecies]+KR[reacSel,]
      }
#       print(reacSel)
    }
  }

  linksThresh=links
  if(sum(links!=0) > 50) {
    lnkThresh=quantile(abs(links)[abs(links)>0],probs=1-topShow,na.rm=TRUE)
    linksThresh[abs(links)<lnkThresh]=0
  }

  g = graph.adjacency(linksThresh, mode = "directed",weighted=TRUE)

  cols=brewer.pal(9,"Set3")
  reacColor=c(cols[6:9])

  reacLab=paste0('R',1:nbReacs)
  V(g)$label = c(reacLab,species)
  V(g)$shape=c(rep("rectangle",nbReacs),rep("circle",nbSpecies))
  V(g)$size = c(rep(10,nbReacs),rep(10,nbSpecies))
  V(g)$size2 = 5
  V(g)$color = c(reacColor[reacType],rep("gold",nbSpecies))
  V(g)$label.cex = c(rep(1.5,nbReacs),rep(1.5,nbSpecies))
  V(g)$label.color = "black"
  V(g)$label.font = 2
  #   V(g)$label.dist = 1
  #   V(g)$label.degree = 0

  wid=abs(E(g)$weight)
  wid=wid^0.1
  wid=0.01+(wid-min(wid))/(max(wid)-min(wid))
  E(g)$width = 2 #wid*10
  E(g)$color = ifelse(E(g)$weight>0,col2tr('red',180),col2tr('blue',180))
  E(g)$arrow.size = 2
  E(g)$curved = FALSE

  loners=which(degree(g)==0)
  if(keepAll)
    loners=which(degree(g)==0 & !(V(g)$label %in% clust) )
  g=delete.vertices(g,loners)

  set.seed(23451)
#   plot(g, layout=layout.graphopt(g,charge=0.1,spring.length=1))
#   return(g)
  outFile = paste0("./digraph.png")
  size=   200*(nbReacs+nbSpecies)^0.5 # 2000
  png(outFile, width=size, height=size)
#   layout=layout.graphopt(g,charge=0.3,spring.length=1)
  compo = t(sapply(species,get.atoms))
  mass  = apply(compo,1,massCxHyNz)
  rankMass=rank(mass,ties.method='first')
  xmax=max(nbSpecies,sum(reacType==1),sum(reacType==2))
  dx1=xmax/sum(reacType==1)
  dx2=xmax/sum(reacType==2)
  layout=matrix(0,nrow=nedges,ncol=2)
  for(i in 1:nbReacs) {
    layout[i,2]=reacType[i]
    layout[i,1]=ifelse(reacType[i]==1,
                       i*dx1,
                       (i-sum(reacType==1))*dx2)
  }
  for(i in (nbReacs+1):nedges) {
    layout[i,2]= 1.5
    layout[i,1]= rankMass[i-nbReacs]
  }
  plot(g,layout=layout,rescale=TRUE)
  dev.off()
}

traceBack = function(sp1,L,R,species,spInit,reacTags,flMean) {

  # Builds digraph of fluxes to and from sp1

  nbSpecies=length(species)
  nbReacs=length(reacTags)

  Fl = matrix(flMean,nrow=nbReacs,ncol=nbSpecies,byrow=FALSE)
  KL = L * Fl
  KR = R * Fl

  cat(paste('Main prod path for ', sp1, '\n'))
  start = sp1
  path = start
  stoppers = spInit
  while (!(start %in% stoppers)) {
    center = which(species == start)
    if (center == 0)
      break
    prods  = which(KR[, center] != 0)
    mainReac = prods[which.max(flMean[prods])]
    wgt = flMean[mainReac] / sum(flMean[prods])
    cat(paste0('  ', reacTags[mainReac], ' : ', signif(wgt, 2), '\n'))
    io = order(KL[mainReac, ], decreasing = TRUE)
    start = species[io[1]]
    if(start %in% stoppers)
      start = species[io[2]]
    start = start[!(start %in% stoppers)]

    if (length(start) == 0)
      break
    if (start %in% path)
      break
    path = c(path, start)
  }
}
budget = function(sp1,L,R,species,reacTags,flMean,weightThresh=0.01) {

  # Builds digraph of fluxes to and from sp1

  nbSpecies=length(species)
  nbReacs=length(reacTags)

  nedges=nbReacs+nbSpecies
  links=matrix(0,ncol=nedges,nrow=nedges)
  dimnames(links)[[1]]=c(reacTags,species)
  dimnames(links)[[2]]=c(reacTags,species)

  KL=L*matrix(flMean,nrow=nbReacs,ncol=nbSpecies,byrow=FALSE)
  KR=R*matrix(flMean,nrow=nbReacs,ncol=nbSpecies,byrow=FALSE)

  links = addLinks(sp1,links,species,KL,KR,nbReacs,nbSpecies)

  g = graph.adjacency(links, mode = "directed",weighted=TRUE)

  reacLab=paste0('R',1:nbReacs)
  V(g)$label = c(reacLab,species)

  center=which(V(g)$label==sp1)
  for (mode in c('in','out')) {
    reac=incident(g,v=center,mode=mode)
    tags=sapply(reac, function(x)
      get.edge(g,x)[ifelse(mode=='in',1,2)])
    if(length(tags)==0) next
    name=reacTags[tags]
    wgt=abs(E(g)[reac]$weight)
    wgt=wgt/sum(wgt)
    iord= order(wgt,decreasing=TRUE)
    wgt=wgt[iord]
    name=name[iord]
    len=sum(wgt>weightThresh)
    tab = cbind(name[1:len],signif(wgt[1:len],digits=2))
    colnames(tab)=c('Reaction','Weight')
    rownames(tab)=rep('',len)
    cat('\n')
    if (mode == 'in')
      cat(paste('Main productions for ',sp1,'\n'))
    else
      cat(paste('Main losses for ',sp1,'\n'))
    print(as.table(tab))
  }

}
volpertIndex <- function (rootSpecies, species, L, R) {
  vlpInd  = rep(NA,length(species))
  names(vlpInd)= species
  spInit  = rootSpecies
  vlpInd[spInit] = 1

  count   = 2
  lReacs  = rowSums(L[,spInit]) >= 2
  spProds = species[colSums(R[lReacs,])!=0]
  spProds = spProds[!(spProds %in% spInit)]
  vlpInd[spProds] = count

  while ( length(spProds) !=0) {
    count=count+1
    spInit   = c(spInit,spProds)
    lReacs  = rowSums(L[,spInit]) >= 2
    spProds = species[colSums(R[lReacs,])!=0]
    spProds = spProds[!(spProds %in% spInit)]
    vlpInd[spProds] = count
  }
  vlpInd
}

plotMS <- function(x, y, yLow, ySup, xlim, ylim,
                   w, col, colt=col, species, xlab= 'm/z',
                   main='', ppScale=FALSE, errBar = FALSE) {
  plot(x, y, type = 'n',
       xlim = xlim, xlab = xlab, cex.axis=0.75,
       ylim = ylim, ylab = 'Mole fraction', yaxs ='i',
       main = main)
  ti = seq(floor(xlim[1]), ceiling(xlim[2]), by=1)
  axis(side = 1, at = ti, labels = ti, cex.axis=0.75)
  grid()
  if(ppScale)
    axis(side = 4,
         labels=c("ppm", "ppb","ppt","ppq"),
         at = -3*(2:5), las = 1, tck = 1)
  for(i in 1:2) # Increase color density
    segments(x, y, x, rep(-30, length(y)),
             lwd = w,
             col = colt)
  # CI
  if(errBar)
    Hmisc::errbar(x, y, ySup, yLow, cap = 0.02,
                  add = TRUE, lty = 1, type = 'p',
                  pch = 16, lwd = 1.2,
                  col = col, errbar.col = col)
  text( x, y, labels = species, col = col,
    pos = 4, offset = 0.25, font = 2)
  box()
}
