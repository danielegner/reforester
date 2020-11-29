# FOURTH YEAR PROJECT
# Last Modified: 03/04/2020 15:32
# main.R

# This script runs repeated simulations based on user selected parameters.
# It calls "functions.R". For my final results, I ran approximately four
# simulations per day to reach a total of 50, all using the following options:
# - "full (541 / 358.9 Myrs)"
# - "diversity-dependent"

# ----------------------------------------------------------------------
# Copyright © 2020 Daniel Egner

# This file is part of Reforester.

# Reforester is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Reforester is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Reforester. If not, see <https://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------



# Prepare R workspace ------------------------------------

graphics.off()
rm(list = ls(all = TRUE))

library(caper)   # For growTree()
library(phytools)   # For bind.tip() and drop.tip()
library(paleotree)   # For timePaleoPhy()

source("functions.R")   # For all custom functions



# Set simulation parameters ------------------------------------

createplots <<- TRUE   # Set global variable to trigger figure generation
samplingmodels <- c("t_occ", "t_fm", "t_loc", "m_occ", "m_fm", "m_loc")

# Ask user which factor to scale-down simulation time by
scalechoice <- NA
while (is.na(scalechoice)) {
	scalechoice <- menu(c("full (541 / 358.9 Myrs)", "one tenth (54.1 / 35.89 Myrs)", "one hundredth (5.41 / 3.589 Myrs)"), title="Which simulation time-scale?")
}

if (scalechoice == 1) {   # Full length simulations for gathering final results
	marstart <<- 541
	terstart <<- 358.9
} else if (scalechoice == 2) {   # Medium length simulations for debugging
	marstart <<- 54.1
	terstart <<- 35.89
} else if (scalechoice == 3) {   # Short length simulations for debugging
	marstart <<- 5.41
	terstart <<- 3.589
}

# Ask user which diversificaiton model to use and how many simulation runs to perform
# Note: exponential diversification models are not feasible when running full or medium length simulations
divmodelchoice <- NA
while (is.na(divmodelchoice)) {
	divmodelchoice <- menu(c("diversity-dependent", "exponential"), title="Which diversification model?")
}

numberofsimruns <- NA
while (is.na(numberofsimruns)) {
	numberofsimruns <- as.numeric(readline("How many simulation runs? "))   # Careful, this user input is unsanitised
}



# Run appropriate simulation loops ------------------------------------

if (divmodelchoice == 1) {
	# Run N simulations for each sampling model using diversity-dependent diversification model
	allsims.slopediffs.divdep.abs <- matrix(data = NA, ncol = length(samplingmodels), nrow = numberofsimruns)
	allsims.slopediffs.divdep.rel <- matrix(data = NA, ncol = length(samplingmodels), nrow = numberofsimruns)

	successfulruns.divdep <- 0
	# Retry loop if error occurs
	while (successfulruns.divdep < numberofsimruns) {
		if (is.null(attr(try(outputs <- RunSimulations(divmodel = "diversitydependent", allsamplemodels = samplingmodels)), "class"))) {

			allsims.slopediffs.divdep.abs[successfulruns.divdep + 1, ] <- outputs[[1]]
			allsims.slopediffs.divdep.rel[successfulruns.divdep + 1, ] <- outputs[[2]]

			successfulruns.divdep <- successfulruns.divdep + 1
			cat("Simulation run", successfulruns.divdep, "complete\n", sep = " ")   # To keep track of progress
		}
	}

} else if (divmodelchoice == 2) {
	# Run N simulations for each sampling model using exponential diversification model
	allsims.slopediffs.exp.abs <- matrix(data = NA, ncol = length(samplingmodels), nrow = numberofsimruns)
	allsims.slopediffs.exp.rel <- matrix(data = NA, ncol = length(samplingmodels), nrow = numberofsimruns)

	successfulruns.exp <- 0
	# Retry loop if error occurs
	while (successfulruns.exp < numberofsimruns) {
		if (is.null(attr(try(outputs <- RunSimulations(divmodel = "exponential", allsamplemodels = samplingmodels)), "class"))) {

			allsims.slopediffs.exp.abs[successfulruns.exp + 1, ] <- outputs[[1]]
			allsims.slopediffs.exp.rel[successfulruns.exp + 1, ] <- outputs[[2]]

			successfulruns.exp <- successfulruns.exp + 1
			cat("Simulation run", successfulruns.exp, "complete\n", sep = " ")   # To keep track of progress
		}
	}
}



# Plot all simulations' slopedifferences ------------------------------------
# Note: only supports plotting diversity-dependent simulations, not exponential ones

dev.new(height = 6, width = 8)
close.screen(all.screens = TRUE)

par(mar = c(4, 4, 1, 1))
plot(x = matrix(data = 1, ncol = 1, nrow = numberofsimruns),
     y = allsims.slopediffs.divdep.rel[, 1],
     bty = "l",
     cex.lab = 1,
     cex.axis = 1,
     pch = 4,
     col = "blue",
     xaxt = "n",
     xlab = "samplingmodels",
     ylab = "re-scaled reconstructed slope - true slope",
     xaxs = "i",
     yaxs = "i",
     xlim = c(0, 7),
     ylim = c(-0.01 + min(allsims.slopediffs.divdep.rel), 0.01 + max(allsims.slopediffs.divdep.rel)))

for (i in 2:length(samplingmodels)) {
	points(x = matrix(data = i, ncol = 1, nrow = numberofsimruns),
		 y = allsims.slopediffs.divdep.rel[, i],
		 pch = 4,
		 col = "blue")
}

axis(side = 1, at = 1:6, labels = samplingmodels)
title(main = paste("REL: All Sampling Models,", numberofsimruns, "Simulation Runs, Each", marstart, "Myrs Long", sep = " "), cex.main = 1)


# Alternative plot
dev.new(height = 6, width = 8)
close.screen(all.screens = TRUE)
split.screen(c(2, 3))

for (i in 1:length(samplingmodels)) {
	screen(i)
	par(mar = c(4, 4, 1, 1))
	plot(x = matrix(data = 1, ncol = 1, nrow = numberofsimruns),
		 y = allsims.slopediffs.divdep.rel[, i],
		 bty = "l",
		 cex.lab = 0.5,
		 cex.axis = 0.7,
		 col = "blue",
		 xlab = "",
		 ylab = "re-scaled reconstructed slope - true slope",
		 xaxs = "i",
		 yaxs = "i",
		 xlim = c(0, 2),
		 ylim = c(-0.01 + min(allsims.slopediffs.divdep.rel[, i]), 0.01 + max(allsims.slopediffs.divdep.rel[, i])))

	title(main = paste("Sim runs: ", numberofsimruns, "; Sampling model: ", samplingmodels[i], sep = ""), cex.main = 0.75)
}
