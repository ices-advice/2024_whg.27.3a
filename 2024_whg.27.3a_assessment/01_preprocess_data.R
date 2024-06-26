## pre-process data ahead of estimation of survey indices
#library(DATRAS)

lastyear <- max(years)
load(file = paste0("../data/WhitingData_1983-", lastyear,".RData"))

geartab = sort(table(d$Gear),TRUE)
d = subset(d, Gear %in% names(geartab[geartab>100]))

## Remove long and short hauls 
d = subset(d, HaulDur >= 15) #  short hauls
d = subset(d, HaulDur <= 100) # long hauls

d <- addSpectrum(d,cm.breaks=0:100)
d <- addWeightByHaul(d,to1min=FALSE)
d$ctime = as.numeric(as.character(d$Year)) + round(d$timeOfYear,1)

d$Nage = matrix(d$HaulWgt, nrow=nrow(d$N),ncol=1)
colnames(d$Nage) <- 1

d$ShipG = factor(paste(d$Ship,d$Gear,sep=":"))

save(d, file = "WhitingData_subset_noCorrection.RData")