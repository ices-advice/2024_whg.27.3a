## Functions
#' Calculate index from given data
#'
#' @param dat Input data as DATRASraw object
#' @param output_folder Where the output should be saved
#'
#' @return Invisible NULL
#'
#' @details Saves a file with the model fit (m0.Rdata) and one with the model summary (modelSummary.txt).
#' @export
#'
#' @examples
calc_Index <- function(dat, output_folder) {
  if (output_folder != ".") system(paste("mkdir", output_folder))
  owd <- setwd(output_folder)
  on.exit({
    setwd(owd)
    print(paste("Changed directory to", getwd()))
  })
  print(paste("Changed to directory", getwd()))
  load(dat)
  print(paste("Loaded data file", dat))
  grid <- getGrid(d, nLon=40)
  # plot(grid)

  tt = median(d$timeOfYear[d$Quarter=="1"])

  fm = "Gear + te(lon,lat,bs=c('ds','ds'),k=c(10,10),m=c(1,0))+te(timeOfYear,lon,lat,bs=c('cc','ds','ds'),k=c(6,5,5),m=c(1,0)) + te(ctime,lon,lat,bs=c('ds','ds','ds'),k=c(24,4,4),m=c(1,0))+s(Depth,bs='ds',k=5,m=c(1,0))+s(ShipG,bs='re',by=dum)+offset(log(HaulDur))"

  knots=list(timeOfYear=seq(0,1,length=6)) ##,TimeShotHour=seq(0,24,length=6))
  ages=1
  
  system.time( m0 <- getSurveyIdx(d,ages,myids=grid[[3]],cutOff=500,fam="Tweedie",
                                  mc.cores=1,modelP=fm,knotsP=knots,nBoot=0,
                                  predfix=list(timeOfYear=tt,ShipG="77AR:GOV"), # to do
                                  control=list(trace=TRUE,maxit=15) ) )

  save(m0,file="m0.RData")

  par(mfrow=c(2,1))
  resi = residuals(m0)
  resi[which(!is.finite(resi))] <- max(resi[is.finite(resi)], na.rm = TRUE)
  qqnorm(resi)
  abline(0,1)

  sink("modelSummary.txt")
  summary(m0$pModels[[1]])
  sink()
  print("Done.")
  invisible(NULL)

}

makeSurveyIndexPlot <- function(mqs, add = FALSE, col = 1) {
  x <- mqs[[1]]
  a <- 1
  ys=sort(unique(x$yearNum))
  idx = x$idx
  lo = x$lo
  up = x$up
  idx[x$idx <= 0] = NA
  lo[x$idx <= 0] = NA
  up[x$idx <= 0] = NA

  par(mar=c(3.5,3.5,1,1))
  if (! add) {
    plot(ys, ys, ylim=c(0,2.5), type="n", xlab="", ylab = "", yaxs="i")
    title(xlab = "Year", ylab="Standardised biomass index", line = 2)
  }
  ysp <- ys
  ysp[1] <- ysp[1] - 0.3
  ysp[length(ysp)] <- ysp[length(ysp)] + 0.3
  polygon(c(ysp, rev(ysp)),
          c(lo[, a]/mean(idx[, a], na.rm = TRUE), rev(up[, a]/mean(idx[, a], na.rm = TRUE))), col = adjustcolor(col, 0.3), border = FALSE)
  lines(ys, idx[, a]/mean(idx[, a], na.rm = TRUE),lwd = 2, pch = 20, type ='b', col = col)
  box()

}

make_all_plots <- function(dat, wd = ".", quarters = 1:4, years = 1983:2019) {
  setwd(wd)
  stopifnot(file.exists(dat))
  stopifnot(file.exists("m0.RData"))
  load("m0.RData")
  load(dat)

  png("entireAreaBubbles.png",width=1200,height=800)
  bubblePlot(d)
  dev.off()

  grid <- getGrid(d, nLon=40)
  # plot(grid)

  factorplot<-function(x, name, invlink=function(xx)xx, ... ){
    sel<- grep( name, names(coef(x)))
    est <- coef(x)[ sel ]
    mynames <- names(est)
    if( substr( mynames[1],0,2)=="s(") mynames <- levels(x$var.summary[[name]])
    sds <- sqrt(diag( vcov(x) )[ sel ])
    lo <- invlink(est - 2*sds)
    hi <- invlink(est + 2*sds)
    ylims <- range(c(lo,hi))
    xs <- 1:length(est)

    plot( xs, invlink(est), ylim=ylims,...,xaxt="n",xlab="", ylab="Estimate")
    axis(1,labels=mynames,at=xs)
    arrows( xs,lo,y1=hi,angle=90,code=3,length=0.1)
    abline(h=invlink(0))
  }

  pdf("gearEffects.pdf")
  factorplot(m0$pModels[[1]],"Gear",invlink=exp)
  dev.off()

  pdf("depthEffect.pdf")
  surveyIdxPlots(m0,d,select="4",par=list(mfrow=c(1,1)),xlim=c(0,150),ylim=c(-1.3,1),main="")
  dev.off()

  qts = 0:3/4+1/8
  mqs = list()
  for(i in quarters) mqs[[i]] = redoSurveyIndex(d,m0,myids=grid[[3]],predfix=list(timeOfYear=qts[i],ShipG="77AR:GOV"),nBoot=400,mc.cores=1)

  ##surveyIdxPlots(m0,d,myids=grid[[3]],select="absolutemap",year=years,colors=rev(heat.colors(10)),par=list(mfrow=c(4,7),mar=c(0,0,2,0)),legend=TRUE,legend.signif=2,map.cex=1)
  surveyIdxPlots(mqs[[1]],d,myids=grid[[3]],select="index",year=years,colors=rev(heat.colors(10)),par=list(mfrow=c(1,1),mar=c(4,3,3,3)),legend=FALSE,legend.signif=2,map.cex=1,ylim=c(0,3),main="")

  for(qq in quarters){
    png(paste0("absolutemapsQ",qq,".png"),width=1200,height=800)
    surveyIdxPlots(mqs[[qq]],d,myids=grid[[3]],select="absolutemap",year=years,colors=rev(heat.colors(10)),par=list(mfrow=c(6,7),mar=c(0,0,2,0)),legend=TRUE,legend.signif=2,map.cex=1)
    dev.off()
  }

  for(qq in quarters){
    png(paste0("concentrationmapsQ",qq,".png"),width=1200,height=800)
    par(mfrow=c(6,7),mar=c(0,0,2,0))
    for(yy in years){
      surveyIdxPlots(mqs[[qq]],d,myids=grid[[3]],select="map",year=yy,colors=rev(heat.colors(10)),par=list(),legend=FALSE,map.cex=1,main="",axes=FALSE)
      box()
      title(yy, line = 1)
    }
    dev.off()
  }


  ###########
  # surveyIdxPlots(m0,d,select="residuals",par=list(mfrow=c(1,1)))
  # surveyIdxPlots(m0,d,select="resVsYear",par=list(mfrow=c(1,1)))
  # surveyIdxPlots(m0,d,select="fitVsRes",par=list(mfrow=c(1,1),mar=c(4,3,3,1)))


  resi = residuals(m0)
  if(all(is.finite(resi))) {
    pdf("qqplot.pdf")
    qqnorm(resi)
    abline(0,1)
    dev.off()
  } else {
    warning("Non finite residuals.")
    pdf("qqplot.pdf")
    qqnorm(resi[which(is.finite(resi))])
    abline(0,1)
    dev.off()
  }


  lastyr <- max(as.numeric(as.character(d$Year)))
  starty=seq(1983,lastyr,by=4)


  for(sy in starty){
    png(paste0("spatialResid",sy,".png"),width=1000,height=1000)
    par(mfrow=c(4,4),mar=c(0,0,2,0))
    for(yy in sy:(sy+3)){
      for(qq in 1:4){
        scale <- 3
        sel <- which(d[[2]]$Year == as.character(yy) & d[[2]]$Quarter==qq)
        plot(d$lon, d$lat, type = "n",
             xlab = "Longitude", ylab = "Latitude", main = "",axes=FALSE)
        maps::map("world", fill = TRUE, plot = TRUE,
                  add = TRUE, col = grey(0.5))
        positive = resi[sel] > 0
        points(d$lon[sel][positive], d$lat[sel][positive],
               pch = 1, cex = scale * sqrt(resi[sel][positive]),
               col = "blue")
        points(d$lon[sel][!positive], d$lat[sel][!positive],
               pch = 1, cex = scale * sqrt(-resi[sel][!positive]),
               col = "red")
        box()
        title(paste0(yy," Q",qq), line = 1)
      }
    }
    dev.off()
    DATRASdataReport(x = d,figure.type="png")
  }


  ############## Export ###############
  mm = mean(mqs[[1]]$idx)
  if(1 %in% quarters)
    resQ1=data.frame(Index=mqs[[1]]$idx[,1]/mm,Year=rownames(mqs[[1]]$idx),Quarter=1, sdlogI=(log(mqs[[1]]$up)-log(mqs[[1]]$lo))/4 )

  if(2 %in% quarters)
    resQ2=data.frame(Index=mqs[[2]]$idx[,1]/mm,Year=rownames(mqs[[2]]$idx),Quarter=1, sdlogI=(log(mqs[[2]]$up)-log(mqs[[2]]$lo))/4 )

  if(3 %in% quarters)
    resQ3=data.frame(Index=mqs[[3]]$idx[,1]/mm,Year=rownames(mqs[[3]]$idx),Quarter=1, sdlogI=(log(mqs[[3]]$up)-log(mqs[[3]]$lo))/4 )

  if(4 %in% quarters)
    resQ4=data.frame(Index=mqs[[4]]$idx[,1]/mm,Year=rownames(mqs[[4]]$idx),Quarter=1, sdlogI=(log(mqs[[4]]$up)-log(mqs[[4]]$lo))/4 )

  ys=years
  plot(ys,resQ1$Index,ylim=c(0,2), type = "b")
  if(2 %in% quarters)
    points(ys+0.25,resQ2$Index,col=2, type= "b")
  if(3 %in% quarters)
    points(ys+0.5,resQ3$Index,col=3, type = "b")
  if(4 %in% quarters)
    points(ys+0.75,resQ4$Index,col=4, type = "b")

  rownames(resQ1)<-NULL
  write.csv(resQ1,file="indexQ1.csv")

  ##png("indexQ1.png", pointsize = 15, width = 640)
  pdf("indexQ1.pdf")
  makeSurveyIndexPlot(mqs = mqs)
  dev.off()

  #######################
  pdf("residSGQY%03d.pdf",onefile=FALSE)
  d$SGQ = factor(paste(d$Ship,d$Gear,d$Quarter,sep=":"))
  nlevels(d$SGQ)
  sum(table(d$SGQ)>=30)
  resi = residuals(m0)
  par(mfrow=c(5,3),mar=c(2,1,3,1))
  for(ss in levels(d$SGQ)){
    sel = which(d$SGQ==ss)
    if(length(sel)>=30){
      plot(d$Year[sel],resi[sel],ylim=c(-3,3),main=ss)
      abline(h=0,col=2)
    }
  }
  dev.off()
  save(dat, mqs, grid, file = "dat_mqs_grid.RData")
}

run_retro_lo <- function(dat, model, wd) {
  owd <- setwd(wd)
  stopifnot(file.exists(dat))
  load(dat)
  load(model)

  retro3a = retro.surveyIdx(mqs[[1]],d,grid,npeels=3,control=list(trace=TRUE,maxit=10))

  save(retro3a, mqs, file = "retro_results.Rdata")

  pdf("retro.pdf",width=10,height=7)
  ## png("retro.png", width=640, pointsize=15)
  plot(retro3a,mqs[[1]],rescale=TRUE,main="Retrospective analysis")
  dev.off()

  rm(retro3a)

  d$Survey2 = as.character(d$Survey)
  d$Survey2[ d$Survey2 %in% c("TN","TOR") ] <- "TN/TOR"
  d$Survey2=factor(d$Survey2)

  lo3a = leaveout.surveyIdx(mqs[[1]],d,grid,fac=d$Survey2,control=list(trace=TRUE,maxit=10))
  save(lo3a, mqs, file = "leave_one_out_results.Rdata")

  pdf("leaveout.pdf",width=10,height=7)
  #png("leaveout.png",width=640, pointsize=15)
  plot(lo3a,mqs[[1]],rescale=TRUE,main="Leave-one-out")
  dev.off()
  setwd(owd)
}


## Data overview for DATRAS data
## Tables:
## Number of hauls by -- survey, year x quarter, year x ship, ship x gear
## Gear plot space
## Bubble plots
## Time of year ( by ship )
## Length distribution (by gear)
DATRASdataReport<-function(x, type="latex",file=paste0("tables.",ifelse(type=="latex","tex","html")),maxBubble=4, figure.type="pdf"){
  requireNamespace("xtable")
  nsu = nlevels(x$Survey)
  nq = nlevels(x$Quarter)
  ns = nlevels(x$Ship)
  ng = nlevels(x$Gear)

  sink(file)
  if(nsu>1){
    print(xtable::xtable(xtabs(~Survey,data=x[[2]]),digits=0,caption="Number of hauls survey"),type=type, caption.placement = "top")

    print(xtable::xtable(xtabs(~Year+Survey,data=x[[2]]),digits=0,caption="Number of hauls by year and survey"),type=type, caption.placement = "top")
  }

  if(nq>1){
    print(xtable::xtable(xtabs(~Year+Quarter,data=x[[2]]),digits=0,caption="Number of hauls by year and quarter"),type=type, caption.placement = "top")
  }

  if(ng>1){
    print(xtable::xtable(xtabs(~Gear,data=x[[2]]),digits=0,caption="Number of hauls by gear"),type=type, caption.placement = "top")
  }

  if(ns>1){
    print(xtable::xtable(xtabs(~Ship+Survey,data=x[[2]]),digits=0,caption="Number of hauls by ship and survey"),type=type, caption.placement = "top")

    if(ng>1){
      print(xtable::xtable(xtabs(~Ship+Gear,data=x[[2]]),digits=0,caption="Number of hauls by ship and gear"),type=type, caption.placement = "top")
    }
  }
  sink()

  if(figure.type=="pdf") pdf(filename="dataplot%03d.pdf",onefile=FALSE,width=10,height=10) else png(filename="dataplot%03d.png",width=1024,height=1024,pointsize=10/480*1024)
  par(mfrow=c(1,1))

  ## Bubble plot all
  x = addSpectrum(x,by=1)
  x = addWeightByHaul(x,FALSE)
  myscale=maxBubble/max(sqrt(x$HaulWgt))
  bubblePlot(x,scale=myscale,pch.zero=".",rim=TRUE)
  legkg = c( 0, round((max(sqrt(x$HaulWgt))/maxBubble)^2/1000,1),
             round((max(sqrt(x$HaulWgt)))^2/1000,1))
  legend("topright",pch=c(46,16,16),col=c(2,1,1),pt.cex=c(1,1,maxBubble),legend=paste(legkg,"kg"),bg="white")
  ## Bubbles by survey
  if(nsu>1){
    par(mfrow=n2mfrow(nlevels(x$Survey)),mar=c(4,3,1,1))
    for(ss in levels(x$Survey)){
      tmp = subset(x,Survey==ss)
      myscale=maxBubble/max(sqrt(tmp$HaulWgt),na.rm=TRUE)
      bubblePlot(tmp,scale=myscale,pch.zero=".",rim=TRUE)
      title(ss)
      legkg = c( 0,
                 round((max(sqrt(tmp$HaulWgt))/maxBubble)^2/1000,1),
                 round((max(sqrt(tmp$HaulWgt)))^2/1000,1))
      legend("topright",pch=c(46,16,16),col=c(2,1,1),pt.cex=c(1,1,maxBubble),legend=paste(legkg,"kg"),bg="white")
    }
  }


  ## Gear plot
  cols = rainbow(nlevels(x$Gear))
  set.seed(123456789); cols = sample(cols,replace=FALSE)
  par(mfrow=c(1,1))
  plot(x,col=cols[as.numeric(x$Gear)],pch=as.numeric(x$Gear),cex=0.4,plot.response=FALSE)
  legend("topright",legend=levels(x$Gear),col=cols,pch=1:nlevels(x$Gear))


  ## length distributions by gear type
  par(mfrow=n2mfrow(nlevels(x$Gear)),mar=c(4,3,1,1))
  for(gg in levels(x$Gear)){
    tmp = subset(x, Gear==gg)
    spec = colSums(tmp$N)
    spec = spec/sum(spec)
    plot(spec,type="h",main=gg,xlim=c(0,70),xlab="Length",ylab="Frequency")
  }

  ## survey time of year
  par(mfrow=c(nlevels(x$Survey),1),mar=c(4,3,1,1))
  for(ss in levels(x$Survey)){
    tmp = subset(x,Survey==ss)
    hist(tmp$timeOfYear,xlim=c(0,1),breaks=seq(0,1,length=52),main=ss,xlab="Time of year")
  }

  dev.off()

}

change2commaInZip <- function(fn, outdir = dirname(fn)) {
  fn <- normalizePath(fn)
  csvfn <- unzip(fn, list = TRUE)$Name
  stopifnot(length(csvfn) == 1)
  con <- unz(fn, csvfn)
  on.exit(close(con))
  dat <- readLines(con)
  dat2 <- gsub(";", ",", dat)
  tempcsv <- file.path(tempdir(), csvfn)
  writeLines(dat2, tempcsv, sep = "\n")
  zip(file.path(outdir, basename(fn)), tempcsv, flags = "-j9")
}

fname <- function(d) {
  surname <- levels(d[[2]]$Survey)
  yrs <- paste(range(levels(d[[2]]$Year)), collapse = "-")
  paste0(surname,"_", yrs, ".Rds")
}
