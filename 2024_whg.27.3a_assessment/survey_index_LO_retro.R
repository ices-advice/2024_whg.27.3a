#######################################################
## Retro + LO
#######################################################
retro3a = retro.surveyIdx(mqs[[1]],d,grid,npeels=3,control=list(trace=TRUE,maxit=10))

save(retro3a, mqs, file = "retro_results.Rdata")

  pdf("retro.pdf",width=10,height=7)
plot(retro3a,mqs[[1]],rescale=TRUE,main="Retrospective analysis")
dev.off()



d$Survey2 = as.character(d$Survey)
d$Survey2[ d$Survey2 %in% c("TN","TOR") ] <- "TN/TOR"
d$Survey2=factor(d$Survey2)


lo3a = leaveout.surveyIdx(mqs[[1]],d,grid,fac=d$Survey2,control=list(trace=TRUE,maxit=10))
save(lo3a, mqs, file = "leave_one_out_results.Rdata")

pdf("leaveout.pdf",width=10,height=7)
plot(lo3a,mqs[[1]],rescale=TRUE,main="Leave-one-out",basename="All data")
dev.off()
