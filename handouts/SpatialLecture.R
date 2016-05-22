library(knitr)
opts_chunk$set(fig.path='figure/beamer-',fig.align='center',fig.show='hold',size='footnotesize')
options(width=60)  # make the printing fit on the page
# some setup
library( fields)
load("~/Home/Teaching/RCTutorial/COMAMExample.rda")
ls()
 data(COmonthlyMet)
 
  x<- CO.loc
  y<- CO.tmin.MAM.climate
  elev<- CO.elev
  good<- !is.na( y)
  x<- x[good,]
  y<- y[good]
  elev<- elev[good]
 # interpolate elevations to a useful a grid (will use these later)
  NGRID <- 50 
  COGrid<- fields.x.to.grid( x, nx=NGRID, ny=NGRID)
  names(COGrid)<- c("x","y")
  data( RMelevation) # elevations for the ROcky mountain area.
  elevGrid<- interp.surface.grid( RMelevation, COGrid )
quilt.plot( x, y) 
US( add=TRUE)
title("CO MAM average minimum temperatures")
fields.style()
plot( elev, y, xlab="elevation", ylab="MAM minimum (C) ")
 X<- cbind( x, elev) 
 lmObj<- lm( y ~ lon + lat + elev, data=X )  
 summary( lmObj)
 quilt.plot( x, lmObj$residuals) 
 US( add=TRUE, col="grey", lwd=2)
 fit1<- spatialProcess( x,y)
 print( fit1)
set.panel( 2,2)
plot( fit1)
surface( fit1)
US(add= TRUE)
 fit1E<- spatialProcess( x,y, Z= elev, ngrid=20)
 surFull<- predictSurface( fit1E, grid.list= COGrid,
                          ZGrid= elevGrid)
 surSmooth<-predictSurface( fit1E, grid.list= COGrid, drop.Z=TRUE)
 set.panel( 1,2)
par( mar=c(1,1,1,4))
image.plot( surFull, col=terrain.colors(256),  
            axes=FALSE,xlab="", ylab="")
US(add=TRUE, col="grey")
image.plot( surSmooth, col=two.colors(256),
            axes=FALSE,xlab="", ylab="")
US(add=TRUE, col="grey")
# unroll the grids in locations and vectors
 COGridPoints<- make.surface.grid( COGrid)
 COGridElevations<- c(elevGrid$z)
 dim( COGRidPoints)
 SEout<- sim.Krig( fit1E,
             xp=COGridPoints,  Z= COGridElevations, M= 10, drop.Z=TRUE)
 dim( SEout)
set.panel( 3,3)
  ctab<- rainbow( 256)
  image.plot( surSmooth, zlim=c(3.5,14.5), col=ctab)
  US(add=TRUE)
  par( mar=c(3,3,1,1))
title("Best estimate", adj=0)
  for( k in 1:8){
      image( as.surface( COGridPoints, SEout[k,]),
               zlim =c(3.5,14.5), axes=FALSE, col=ctab)
       US(add=TRUE)
      title( paste(k), adj=0)
  }
