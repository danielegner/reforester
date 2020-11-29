# FOURTH YEAR PROJECT
# Last Modified: 06/05/2020 12:48
# histplotter.R

# Script for plotting histograms for all eight relative slopedifference
# results sets (sixe proxy-derived sampling model sets and two uniform
# sampling model sets).

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



# Initialise figure
ORIGsamplingmodels <- c("t_occ", "t_fm", "t_loc", "m_occ", "m_fm", "m_loc")
USsamplingmodels <- c("t_uniform", "m_uniform")

dev.new(height = 6, width = 10)
close.screen(all.screens = TRUE)
split.screen(c(2, 5))

# Set which sub-plot positions to plot each histogram in
ORIGplacer <- c(1, 2, 3, 6, 7, 8)
USplacer <- c(5, 10)

# Set title
title(ylab = "Frequency                                                                 ", cex.lab = 0.8)


# Plot all six proxy-derived sampling model histograms
for (i in 1:length(ORIGonefiftyREL[1, ])) {

	screen(ORIGplacer[i])
	par(mar = c(4, 4, 1, 1))

	# Plot the histogram
	h <- hist(ORIGonefiftyREL[, i], density = 25, col = "#66CDAA", xlab = "",  ylab = "", main = ORIGsamplingmodels[i], xlim = c(min(0, ORIGonefiftyREL[, i]) - 0.001, max(ORIGonefiftyREL[, i]) + 0.001), ylim = c(0, 20), bty = "l", cex.axis = 0.9)
	
	# Superimpose a normal distribution curve centred on the histogram
	xfit <- seq(min(0, ORIGonefiftyREL[, i]) - 0.001, max(ORIGonefiftyREL[, i]) + 0.001, length = 100) 
	yfit <- dnorm(xfit, mean = mean(ORIGonefiftyREL[, i]), sd = sd(ORIGonefiftyREL[, i])) 
	yfit <- yfit * diff(h$mids[1:2]) * length(ORIGonefiftyREL[, i]) 
	lines(xfit, yfit, col = "black", lwd = 1)

	# Plot axis label under the last proxy-derived histogram
	if (i == 6) {
		title(xlab = "Slope difference", line = 2.5, cex.lab = 0.8)
	}

	# Plot blue vertical dashed line at slopedifference of zero
	lines(x = c(0, 0), y = c(0, 20), lty = 5, col = "darkblue", lwd = 1)

}

# Plot both uniform sampling model histograms
for (i in 1:length(USonefiftyREL[1, ])) {

	screen(USplacer[i])
	par(mar = c(4, 4, 1, 1))

	# Plot the histogram
	h <- hist(USonefiftyREL[, i], density = 25, col = "#F88379", xlab = "",  ylab = "", main = USsamplingmodels[i], xlim = c(min(0, USonefiftyREL[, i]) - 0.001, max(USonefiftyREL[, i]) + 0.001), ylim = c(0, 20), bty = "l", cex.axis = 0.9)
	
	# Superimpose a normal distribution curve centred on the histogram
	xfit <- seq(min(0, USonefiftyREL[, i]) - 0.001, max(USonefiftyREL[, i]) + 0.001, length = 100) 
	yfit <- dnorm(xfit, mean = mean(USonefiftyREL[, i]), sd = sd(USonefiftyREL[, i])) 
	yfit <- yfit * diff(h$mids[1:2]) * length(USonefiftyREL[, i]) 
	lines(xfit, yfit, col = "black", lwd = 1)

	# Plot blue vertical dashed line at slopedifference of zero
	lines(x = c(0, 0), y = c(0, 20), lty = 5, col = "darkblue", lwd = 1)

}


# Export to pdf
dev.print(pdf, "histplots.pdf", height = 6, width = 10)
