filepath <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/ComponentTemporalStability/CompAutoCorrFinal_20Sept2013.csv"

data <- read.csv(paste(filepath, sep = ""), header = T, sep = "\t")
data$overlap <- data$intersection.size / data$limited.union.size

#-- plot out correlations by lag --#
bb <- subset(data, Pop == "BLACKBUTTE")
rb <- subset(data, Pop == "REDBIRD")
we <- subset(data, Pop == "WENAHA")
im <- subset(data, Pop == "IMNAHA")

#-- plot v1 --#
plot(data$overlap ~ data$Lag, col = data$Pop, xlab = "Lag (years)", ylab =
		 "Within-component Autocorrelation")
lines(lowess(bb$overlap ~ bb$Lag), col = 1)
lines(lowess(rb$overlap ~ rb$Lag), col = 3)
lines(lowess(we$overlap ~ we$Lag), col = 4)
lines(lowess(im$overlap ~ im$Lag), col = 2)


#-- plot v2 --#
new.cols <- c(rgb( red = 0, green = 0, blue = 0, alpha = .5),
							rgb(red = 1, green = 0, blue = 0, alpha = .5),
							rgb(red = 0, green = 1, blue = 0, alpha = .5),
							rgb(red = 0, green = 0, blue = 1, alpha = .5))

#-- ID components with only one stationarity points --@
data$duration <- rep(NA, dim(data)[1])
for(i in 1:dim(data)[1]){
	k <- subset(data, as.character(StableCompoLabel) ==
							as.character(data$StableCompoLabel)[i])
	data$duration[i] <- max(k$Lag) 
}
multiyr.dat <- subset(data, duration >= 2)
bb.mu <- subset(multiyr.dat, Pop == "BLACKBUTTE")
rb.mu <- subset(multiyr.dat, Pop == "REDBIRD")
we.mu <- subset(multiyr.dat, Pop == "WENAHA")
im.mu <- subset(multiyr.dat, Pop == "IMNAHA")

write.path <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Plots/RevisedPlots_19Sept2013/CompoTempAutoCorr_20Sept2013.svg"
svg(paste(write.path), height = 5, width = 7)
boxplot(data$overlap ~ data$Pop + data$Lag, col = new.cols[rep(1:4, 10)], xlab
				= "Lag (years)", ylab = "Within-component Autocorrelation", xaxt = "n",
				xlim = c(0, 30))
lines(lowess(bb$overlap ~ (bb$Lag - 1) * 4  +1) , col = 1, lwd = 2)
lines(lowess(rb$overlap ~ (rb$Lag - 1) * 4 + 3) , col = 3, lwd = 2)
lines(lowess(we$overlap ~ (we$Lag - 1) * 4 + 4), col = 4, lwd = 2)
lines(lowess(im$overlap ~ (im$Lag - 1) * 4 + 2), col = 2, lwd = 2)
legend(x = 24, y = 1, levels(data$Pop), bty = "n",  fill = new.cols[1:4])
dev.off()

boxplot(multiyr.dat$overlap ~ multiyr.dat$Pop + multiyr.dat$Lag, col = new.cols[rep(1:4, 10)], xlab
				= "Lag (years)", ylab = "Within-component Autocorrelation", xaxt = "n",
				xlim = c(0, 30))
lines(lowess(bb.mu$overlap ~ (bb.mu$Lag - 1) * 4  +1) , col = 1, lwd = 2)
lines(lowess(rb.mu$overlap ~ (rb.mu$Lag - 1) * 4 + 3) , col = 3, lwd = 2)
lines(lowess(we.mu$overlap ~ (we.mu$Lag - 1) * 4 + 4), col = 4, lwd = 2)
lines(lowess(im.mu$overlap ~ (im.mu$Lag - 1) * 4 + 2), col = 2, lwd = 2)
legend(x = 24, y = 1, levels(multiyr.dat$Pop), bty = "n",  fill = new.cols[1:4])

#-- Stationarity --#
pop.stationarity <- tapply(data$overlap, data$Pop, mean)
compo.stationarity <- tapply(data$overlap, data$StableCompoLabel, mean)

#-- what I really need is summer lamb survival within each component --#
lambdat.path <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/LambSurvDatWithLastYearCompCovs/LambDataWithLastYearCompCovs_19Sept2013.csv"
lambdat <- read.csv(paste(lambdat.path, sep = ""), header = T)

#-- map lambdat$ThisYearLambMortOrNoLamb to component.names to --#
#-- compo.stationarity --#
data$Y1CompoLambMortOrNoLamb <- data$Y1PNYear <- rep(NA, dim(data)[1])
data$Y2CompoLambMortOrNoLamb <- data$Y2PnYear <- rep(NA, dim(data)[1])

for(i in 1:dim(data)[1]){
	k <- subset(lambdat, as.character(component.ind) ==
										 as.character(data$Y1component.ind[i]))
	data$Y1CompoLambMortOrNoLamb[i] <- k$ThisYearCompLambDiedOrNoLamb[1]
	data$Y1PNYear[i] <- k$PNYear[1]
	l <- subset(lambdat, as.character(component.ind) ==
							as.character(data$Y2component.ind)[i])
	data$Y2CompoLambMortOrNoLamb[i] <- l$ThisYearCompLambDiedOrNoLamb[1]
	data$Y2PNYear[i] <- k$PNYear[1]
}

data$DiffInLambMort <- data$Y1CompoLambMortOrNoLamb -
data$Y2CompoLambMortOrNoLamb

compo. <- subset(data,  is.na(Y1CompoLambMortOrNoLamb) ==
												F)
compo.summary <- data.frame(rep(NA, length(levels(compo.$StableCompoLabel))))
compo.summary$duration <- compo.summary$PNyrtot <- rep(NA, length(levels(compo.$StableCompoLabel)))
	for(i in 1:dim(compo.summary)[1]){
		k <- subset(data, as.character(StableCompoLabel) ==
								as.character(levels(data$StableCompoLabel))[i])
		compo.summary$duration[i] <- max(k$Lag)
		compo.summary$PNyrtot[i] <- table(k$Y1PNYear == 1)["TRUE"]
	}
compo.summary$PropPNyrs <- compo.summary$duration / compo.summary$PNyrtot

compo.summary$stationarity <- tapply(compo.$overlap,
																		compo.$StableCompoLabel,
																		mean)
compo.summary$lambmort <- tapply(compo.$Y1CompoLambMortOrNoLamb,
																compo.$StableCompoLabel, mean)
compo.summary$lambmort.sd <- tapply(compo.$Y1CompoLambMortOrNoLamb,
																compo.$StableCompoLabel, sd)
compo.summary$diffinlambmort <- tapply(compo.$DiffInLambMort,
																			compo.$StableCompoLabel, mean)
compo.summary$size <- table(compo.$StableCompoLabel)
compo.summary$pop <- factor(c(rep("BlackButte", 8), rep("Imnaha", 22-8), rep("Redbird", 45 -
																															22),
							 rep("Wenaha", 7))) 
compo.summary$StableCompoLabel <-
	as.character(levels(compo.$StableCompoLabel))
compo.summary <- compo.summary[, -1]

#-- boxplots of component duration by population --#
install.packages("vcd")
require(vcd)

no.comps.tab <- table(compo.summary$pop)
mosaic()
mosaic(art, gp = shading_max, split_vertical = TRUE, main = "Arthritis:
			 [Treatment]")

compo.geq2yr <- subset(compo.summary, duration >= 2)

new.col <- c(rgb(red = 0, green = 0, blue = 0, alpha = .5),
						 rgb(red = 1, green = 0, blue = 0, alpha = .5),
						 rgb(red = 0, green = 1, blue = 0, alpha = .5),
						 rgb(red = 0, green = 0, blue = 1, alpha = .5))
plot(compo.geq2yr$lambmort ~ compo.geq2yr$stationarity, pch = 16, ylim = c(0,
																																					 1),
		 xlim = c(0, 1), cex = compo.geq2yr$size * .5, col =
		 new.col[compo.geq2yr$pop])

plot(compo.geq2yr$lambmort ~ compo.geq2yr$stationarity, pch = 16, ylim = c(0,
																																					 1),
		 xlim = c(0, 1), cex = 3, col = new.col[compo.geq2yr$pop])

compo.highPNprop <- subset(compo.summary, PropPNyrs == 1)
plot(1 - compo.highPNprop$lambmort ~
		 compo.highPNprop$stationarity, pch = 16, cex =
		 compo.highPNprop$size, xlab = "Component stationarity", ylab = "Mean
		 Component reproductive success rate")

compo.geq2yr.highPNprop <- subset(compo.geq2yr, PropPNyrs == 1)

plot(compo.geq2yr.highPNprop$lambmort ~
		 compo.geq2yr.highPNprop$stationarity, pch = 16, cex =
		 compo.geq2yr.highPNprop$size * .5)
lines(lowess(compo.geq2yr.highPNprop$lambmort ~
						 compo.geq2yr.highPNprop$stationarity, f = .9))

fit <- lm(compo.geq2yr.highPNprop$lambmort ~
					compo.geq2yr.highPNprop$stationarity +
					I(compo.geq2yr.highPNprop$stationarity ^ 2))

#-- read in compiled data to compare stationarity with SLS --#
compd.path <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/compiled_data_summary_0124_efc_withlambfirstbirth.csv"
compd.data <- read.csv(paste(compd.path, sep = ""), header = T) 
popstouse <- c("BlackButte", "Redbird", "Imnaha", "Wenaha")
pnyears <- subset(compd.data, (CLASS == "ALL_AGE"| CLASS == "ALL_AGE_INITIAL" |
									CLASS == "LAMBS" ) & Pop %in% popstouse) 
pnyears <- subset(pnyears, is.na(SumLambSurv) == F)
pnyears$Pop <- factor(pnyears$Pop)

pop.SLS <- tapply(pnyears$SumLambSurv, pnyears$Pop, mean)

plot(pop.SLS ~ pop.stationarity, pch = 16, xlim = c(0, 1), ylim = c(0, 1))


