#-----------------------------------------------------------#
#-- code for fitting cox models with decomposed variances --#
#-----------------------------------------------------------#

filepath =
"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/CleanLambSurvData/"

#-- path data with all edges --#
#allrelocs <- read.csv(paste(filepath, "FullEweRelocsLambSurvDat_allewerelocs_17Sept2013.csv", sep = ""), header=T)

#-- path to data with edgeweights >= 0.1 --#
allrelocs <- read.csv(paste(filepath, "FullEweRelocsLambSurvDat_MinEdge.1_allewerelocs_18Sept2013.csv", sep = ""), header=T)

#install.packages("coxme")
#install.packages("plotrix")
require(coxme)
require(plotrix)

#-- rename dataset for analysis "lamb.data" --#
lamb.data <- allrelocs

#-- build popyear.ind covariate by pasting each lamb's --#
#-- Pop and Year covariates --#
lamb.data$popyear.ind <- rep(NA, dim(lamb.data)[1])
for(i in 1:dim(lamb.data)[1]){
	lamb.data$popyear.ind[i] <- paste(lamb.data$Pop[i], "_", lamb.data$Year[i], sep = "")
}

#--subset datasets to include only pn years, with multiple components --#
pnyears <- subset(lamb.data, PNYear == 1)
mult.comps.pn <- subset(pnyears, no.components >= 2)
mult.comps.pn$component.ind2 <- rep(NA, dim(mult.comps.pn)[1])

#-- note: event should be coded 0 = alive, 1 = dead
#-- fit nested shared fraitly models using coxme --#
pnyear.pop.coxfit <- coxme(Surv(SURV_DAYS, CENSOR2) ~ (1 | factor(Pop)), data = mult.comps.pn)
pnyear.pop.py.coxfit <- coxme(Surv(SURV_DAYS, CENSOR2) ~ (1 | factor(Pop) / factor(popyear.ind)), data = mult.comps.pn)
pnyear.pop.py.compo.coxfit <- coxme(Surv(SURV_DAYS, CENSOR2) ~ (1 | factor(Pop) / factor(popyear.ind) / factor(component.ind)), data = mult.comps.pn)
#
#pnyear.pop.coxfit <- coxme(Surv(SURV_DAYS, CENSOR2) ~ (1 | factor(Pop)), data =
#													 pnyears)
#pnyear.pop.py.coxfit <- coxme(Surv(SURV_DAYS, CENSOR2) ~ (1 | factor(Pop) /
#																													factor(popyear.ind)),
#															data = pnyears)
#pnyear.pop.py.compo.coxfit <- coxme(Surv(SURV_DAYS, CENSOR2) ~ (1 | factor(Pop)
#																																/
#																																factor(popyear.ind)
#																																/
#																																factor(component.ind)),
#																		data = pnyears) 
#-- ANOVA comparison of model fits --#
finestanova <- anova(pnyear.pop.py.compo.coxfit, pnyear.pop.py.coxfit)
midanova    <- anova(pnyear.pop.py.coxfit, pnyear.pop.coxfit)
coarseanova <- anova(pnyear.pop.coxfit)

#-- code for triangular grid of RE hierarchical histograms --#
k <- k1 <- c(-0.9, hist(pnyear.pop.py.compo.coxfit$frail$`factor.Pop./factor.popyear.ind./factor.component.ind.`, breaks = 15)$breaks)
main1    = "Population REs"
main2    = "Year in population REs"
main3    = "Component in year in population REs"
x.lab    = ""
y.lab    = "Pop/Popyear/Component"
x.lim    = c(-1.5, 1.5)
filename = "REHestsNestedModel_19Mar2013.svg"
col      = "grey90"

#svg(filename = filename, width = 7, height = 3)
par(mfrow=c(1, 3), oma=c(5, 2, 0, 0), mex=.7, mar=c(3, 5, 3, 3))
hist(pnyear.pop.py.compo.coxfit$frail$factor.Pop., main = main1, xlab = x.lab,
		 ylab = y.lab, ylim = c(0, 10), xlim = c(-1.5, 1.5), breaks = k, col = col, cex.lab = 1.4)
hist(pnyear.pop.py.compo.coxfit$frail$'factor.Pop./factor.popyear.ind.', main =
		 main2, xlab = "", ylim = c(0, 3), xlim = c(-1.5, 1.5), breaks = k, col = col, ylab = y.lab)
hist(pnyear.pop.py.compo.coxfit$frail$`factor.Pop./factor.popyear.ind./factor.component.ind.`,
		 main = main3, xlab = "", ylim = c(0, 2), xlim = c(-1.5, 1.5), breaks = k, col = "grey90", ylab = "")
mtext(side = 1, outer = T, line = .5, "Model is PN years only, Lamb Surv ~1|Pop/Popyear/Component", cex = .9)
mtext(side = 1, outer = T, line = 2.5, "Fit using coxme (2012); networks based only on all ewe observations from May 1 to Oct 1.", cex = .9)
#dev.off()

#-- code to write out a table of variance components for the saturated model --#
component.var <- VarCorr(pnyear.pop.py.compo.coxfit)$`factor.Pop./factor.popyear.ind./factor.component.ind.`
popyear.var   <- VarCorr(pnyear.pop.py.compo.coxfit)$`factor.Pop./factor.popyear.ind.`
pop.var       <- VarCorr(pnyear.pop.py.compo.coxfit)$factor.Pop.
var.vec       <- c(component.var, popyear.var, pop.var)

component.sd <- sqrt(component.var)
popyear.sd   <- sqrt(popyear.var)
pop.sd       <- sqrt(pop.var)
sd.vec       <- c(component.sd, popyear.sd, pop.sd)

component.sig <- anova(pnyear.pop.py.compo.coxfit, pnyear.pop.py.coxfit)$P[2]
popyear.sig   <- anova(pnyear.pop.py.coxfit, pnyear.pop.coxfit)$P[2]
anova.p.vec   <- c(component.sig, popyear.sig, NA)

model.level <- c("Component", "Population-Year", "Population")
vardecomp.out <- data.frame(cbind(model.level, var.vec, sd.vec, anova.p.vec))
names(vardecomp.out) <- c("ModelLevel", "RE Variance (Nested)", "RE Standard deviation (Nested)", "Anova P-value (Sequentially Nested)")

#-- same variance decomposition for healthy years --#
healthy <- subset(lamb.data, PNYear == 0)
mult.comps.he <- subset(healthy, no.components >= 2)
mult.comps.he$component.ind2 <- rep(NA, dim(mult.comps.he)[1])

#-- fit mixed effects models for health years --#
heyear.pop.coxfit <- coxme(Surv(SURV_DAYS, CENSOR2) ~ (1|factor(Pop)), data=mult.comps.he)
heyear.pop.py.coxfit <- coxme(Surv(SURV_DAYS, CENSOR2)~(1|factor(Pop)/factor(popyear.ind)),data=mult.comps.he)
heyear.pop.py.compo.coxfit <- coxme(Surv(SURV_DAYS,CENSOR2) ~ (1|factor(Pop)/factor(popyear.ind)/factor(component.ind)), data=mult.comps.he)

#-- ANOVA for healthy years --#
finestanova.he <- anova(heyear.pop.py.compo.coxfit,
												heyear.pop.py.coxfit)
midanova.he <- anova(heyear.pop.py.coxfit, heyear.pop.coxfit)
coarseanova.he <- anova(heyear.pop.coxfit)

#-- variance in REs --#
component.var.he <- VarCorr(heyear.pop.py.compo.coxfit)$`factor.Pop./factor.popyear.ind./factor.component.ind.`
popyear.var.he   <- VarCorr(heyear.pop.py.compo.coxfit)$`factor.Pop./factor.popyear.ind.`
pop.var.he       <- VarCorr(heyear.pop.py.compo.coxfit)$factor.Pop.
var.vec.he       <- c(component.var.he, popyear.var.he, pop.var.he)
sig.vec.he <- sqrt(var.vec.he)

#-- Variance decomposition density plots for both PN and HE --#
#-- first extract all densities --#
pop.dens <- density(pnyear.pop.py.compo.coxfit$frail$factor.Pop.)
py.dens <-
	density(pnyear.pop.py.compo.coxfit$frail$'factor.Pop./factor.popyear.ind.')
compo.dens <-
	density(pnyear.pop.py.compo.coxfit$frail$'factor.Pop./factor.popyear.ind./factor.component.ind.')

pop.dens.he <- density(heyear.pop.py.compo.coxfit$frail$factor.Pop., from =
											 -1.5, to = 1.5)
py.dens.he <-
	density(heyear.pop.py.compo.coxfit$frail$'factor.Pop./factor.popyear.ind.',
					from = -1.5, to = 1.5)
compo.dens.he <-
	density(heyear.pop.py.compo.coxfit$frail$'factor.Pop./factor.popyear.ind./factor.component.ind.')

leg.text <- c("Population", "Population-year", "Component")

write.path <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Plots/RevisedPlots_19Sept2013/VarianceDecompDensities_19Sept2013.svg"

svg(write.path, width = 12, height = 7)
par(mfrow = c(1, 2), oma = c(3, 0, 0, 0),  bty = "n")
#axis.break()
#axis.break(axis = 2, breakpos = c(2.25, 7.75), style = "zigzag", breakcol = "black", brw = 0.05)
gap.plot(0, 0, gap = c(2.25, 7.75), gap.axis = "y",  pch = "", xlim = c(-1.5, 1.5),
				 ylim = c(0, 9), breakcol = "red", brw = .05, ytics = c(0, 1, 2,
																																 
																																 8, 9), xtics =
				 seq(-1.5, 1.5, by = .5), ylab = "Density", main = "PN Year Random Effects")
gap.plot(pop.dens$x[300:356], pop.dens$y[300:356], gap = c(2.25, 7.75), add =
				 T, type = "l", col = "orange", lwd = 2)
lines(pop.dens$x[1:98], pop.dens$y[1:98], col = "orange", lwd = 2)
lines(pop.dens$x[427:511], pop.dens$y[427:511], col = "orange", lwd = 2)
lines(py.dens, col = "blue", lwd = 2)
lines(compo.dens, col = "black", lwd = 2)
mtext(side = 1, line = 3, outer = F, "Random effect values")
legend(x = 0.25, y = 2, bty = "n", leg.text, col = c("orange", "blue", "black"), lwd
			 = rep(2, 3), lty = rep(1, 3))

gap.plot(0, 0, gap = c(2.25, 7.75), gap.axis = "y", pch = "", xlim = c(-1.5,
				 1.5), ylim = c(0, 9), breakcol = "red", brw = 0.05, ytics = c(0, 1, 2, 8,
				9), xtics = seq(-1.5, 1.5, by = .5), ylab = "Density", main = "Healthy Year Random Effects")
lines(pop.dens.he, col = "orange", lwd = 2)
lines(py.dens.he, col = "blue", lwd = 2)
lines(compo.dens.he, col = "black", lwd = 2)
mtext(side = 1, line = 3, outer = F, "Random effect values")
legend(x = 0.25, y = 2, bty = "n", leg.text, col = c("orange", "blue", "black"), lwd
			 = rep(2, 3), lty = rep(1, 3))
dev.off()

#-- inset of component mort rates in healthy vs. pn years --#
#-- bring in new dataset: Sept 19 version of netowrk covariate lamb data --#
#-- this dataset has info on component-level mort rates for each component
read.path <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/LambSurvDatWithLastYearCompCovs/"
comp.covs.data <- read.csv(paste(read.path,
																 "LambDataWithLastYearCompCovs_19Sept2013.csv",
																 sep = ""), header = T)

compo.data <- vector("list", length = length(levels(lamb.data$component.ind)))
for(i in 1:length(compo.data)){
	k <- subset(comp.covs.data, as.character(component.ind) ==
							levels(lamb.data$component.ind)[i])
	compo.data[[i]] <- k[1, ]
}

compo.dat <- do.call(rbind, compo.data)
no.comps <- table(compo.dat$PNYear)
compo.pn <- subset(compo.dat, PNYear == 1)
compo.he <- subset(compo.dat, PNYear == 0)
compo.sdest <- c(sd(na.omit(compo.pn$ThisYearCompLambDiedOrNoLamb)),
									sd(na.omit(compo.he$ThisYearCompLambDiedOrNoLamb)))

#-- for Pop-year --#
popyr.data <- vector("list", length(levels(lamb.data$popyear.ind))) 
for(i in 1:length(popyr.dat)){
	k <- subset(comp.covs.data, as.character(popyear.ind) ==
							levels(factor(lamb.data$popyear.ind))[i])
	popyr.data[[i]] <- k[1, ]
}
popyr.dat <- do.call(rbind, popyr.data)
no.popyrs <- table(popyr.dat$PNYear)
popyr.pn <- subset(popyr.dat, PNYear == 1)
popyr.he <- subset(popyr.dat, PNYear == 0)
popyr.sdest <- c(sd(na.omit(popyr.pn$ThisPopyrLambDiedOrNoLamb)),
									sd(na.omit(popyr.he$ThisPopyrLambDiedOrNoLamb)))

write.path <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Plots/RevisedPlots_19Sept2013/CompoLambSurvBoxplots_19Sept2013.svg"

#-- 
svg(paste(write.path))
par(mfrow = c(1, 2), mex = .8)
boxplot(1 - compo.dat$ThisYearCompLambDiedOrNoLamb ~ compo.dat$PNYear, xaxt = "n",
				ylab = "Component reproductive success through weaning", col =
				"blue", ylim = c(0, 1))
text(1, 0.8, paste("N = ", no.comps[1], sep = ""), col = "lightblue", cex = .8)
text(2, 0.2, paste("N = ", no.comps[2], sep = ""), col = "lightblue", cex = .8)
text(1, 0.95, paste("sd = ", round(compo.sdest[1], 2), sep = ""), col =
		 "lightblue", cex = .8)
text(2, 0.6, paste("sd = ", round(compo.sdest[2], 2), sep = ""), col =
		 "lightblue", cex = .8)
axis(side = 1, at = c(1, 2), labels = c("Healthy", "Pneumonia"))
mtext(side = 1, outer = T, line = 1, "Population-Year's lamb health classification")

boxplot(1 - popyr.dat$ThisPopyrLambDiedOrNoLamb ~ popyr.dat$PNYear, xaxt = "n",
				ylab = "Population Year reproductive success through weaning", col =
				"orange", ylim = c(0, 1))
text(1, 0.76, paste("N = ", no.popyrs[1], sep = ""), cex = .8)
text(2, 0.5, paste("N = ", no.popyrs[2], sep = ""), cex = .8)
text(1, 0.83, paste("sd = ", round(popyr.sdest[1], 2), sep = ""), cex = .8)
text(2, 0.63, paste("sd = ", round(popyr.sdest[2], 2), sep = ""), cex = .8)
axis(side = 1, at = c(1, 2), labels = c("Healthy", "Pneumonia"))
mtext(side = 1, outer = T, line = 1, "Population-Year's lamb health classification")

dev.off()
