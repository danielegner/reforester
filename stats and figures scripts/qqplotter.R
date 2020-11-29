# FOURTH YEAR PROJECT
# Last Modified: 03/05/2020 18:05
# qqplotter.R

# Script for plotting qq-plots of all six proxy-derived sampling models'
# "relative slope difference" results, as well as for the two uniform
# sampling models' results.

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



graphics.off()


# Set number of displayed significant figures
options(scipen = 999)

ORIGsamplingmodels <- c("t_occ", "t_fm", "t_loc", "m_occ", "m_fm", "m_loc")
USsamplingmodels <- c("t_uniform", "m_uniform")


# Initialise figure
dev.new(height = 6, width = 10)
close.screen(all.screens = TRUE)
split.screen(c(2, 5))

ORIGplacer <- c(1, 2, 3, 6, 7, 8)
USplacer <- c(5, 10)


# Plot qq-plots of the "relative slope difference" results for all six sampling models
for (i in 1:length(ORIGonefiftyREL[1, ])) {

	screen(ORIGplacer[i])
	par(mar = c(4, 4, 1, 1))

	qqnorm(ORIGonefiftyREL[, i], main = ORIGsamplingmodels[i], xlab = "", ylab = "", col = "#66CDAA", bty = "l", cex.axis = 0.9)
	qqline(ORIGonefiftyREL[, i])


	if (i == 1) {
		title(ylab = "Sample quantiles                                                                  ", line = 2.5, cex.lab = 0.8)
	}

	if (i == 6) {
		title(xlab = "Theoretical quantiles", cex.lab = 0.8)
	}

}


# Plot qq-plots of the "relative slope difference" results for both uniform sampling models
for (i in 1:length(USonefiftyREL[1, ])) {

	screen(USplacer[i])
	par(mar = c(4, 4, 1, 1))

	qqnorm(USonefiftyREL[, i], main = USsamplingmodels[i], xlab = "", ylab = "", col = "#F88379", bty = "l", cex.lab = 0.8, cex.axis = 0.9)
	qqline(USonefiftyREL[, i])

}


# Export to pdf
dev.print(pdf, "qqplots.pdf", height = 6, width = 10)
