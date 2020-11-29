# FOURTH YEAR PROJECT
# Last Modified: 06/05/2020 12:43
# boxplotter.R

# Script for plotting boxplots of all calculated relative slopedifferences.
# Plots all eight boxplots (three from the terrestrial tetrapod proxies,
# three from the marine eumetazoa proxies, followed by the two uniform
# sampling models).

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

library(scales)



# Initialise plot with correct axis limits, margins, axis labels, etc.
ylimit <- c(min(ORIGonefiftyREL, USonefiftyREL) - 0.0001, max(ORIGonefiftyREL, USonefiftyREL) + 0.0001)

dev.new(height = 6, width = 8)
close.screen(all.screens = TRUE)
par(mar = c(6, 4, 1, 1))

plot(x = 0, y = 0, xlim = c(0.5, 9.5), ylim = ylimit,  type = "n", xaxt = "n", ylab = "Relative slope difference", xlab = "", bty = "l", cex.lab = 0.8 , cex.axis = 0.9 )
title(xlab = "Sampling model", line = 4, cex.lab = 0.8)
axis(side = 1, at = seq(from = 1, to = 9, by = 1), labels = c("t_occ", "t_fm", "t_loc", "m_occ", "m_fm", "m_loc", "", "t_uniform", "m_uniform"), tick = FALSE, las = 2, line = -0.5)

# Plot horizontal lines across whole graph to aid readability of the vertical axis
lines(x = c(0, 9.5), y = c(0, 0), lty = 5, col = "darkblue")   # Draw a dashed blue hoizontal line at relative slope difference of zero
lines(x = c(0, 9.5), y = c(-0.001, -0.001), lty = 1, col = alpha("darkgrey", 0.2))
lines(x = c(0, 9.5), y = c(0.001, 0.001), lty = 1, col = alpha("darkgrey", 0.2))
lines(x = c(0, 9.5), y = c(0.002, 0.002), lty = 1, col = alpha("darkgrey", 0.2))
lines(x = c(0, 9.5), y = c(0.003, 0.003), lty = 1, col = alpha("darkgrey", 0.2))
lines(x = c(0, 9.5), y = c(0.004, 0.004), lty = 1, col = alpha("darkgrey", 0.2))


# Plot all six proxy-derived sampling model boxplots
for (i in 1:length(ORIGonefiftyREL[1, ])) {

	boxplot(ORIGonefiftyREL[, i], add = TRUE, at = i, axes = FALSE, col = "#66CDAA", outline = FALSE, whisklty = 1)
	points(x = matrix(data = i + 0.3, nrow = length(ORIGonefiftyREL[, i]), ncol = 1), y = ORIGonefiftyREL[, i], pch = ".")

}

# Plot both uniform sampling model boxplots
for (i in 1:length(USonefiftyREL[1, ])) {

	boxplot(USonefiftyREL[, i], add = TRUE, at = i + 7, axes = FALSE, col = "#F88379", outline = FALSE, whisklty = 1)
	points(x = matrix(data = i + 7 + 0.3, nrow = length(USonefiftyREL[, i]), ncol = 1), y = USonefiftyREL[, i], pch = ".")

}

# Export to pdf
dev.print(pdf, "boxplots.pdf", height = 6, width = 8)
