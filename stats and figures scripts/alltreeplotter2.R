# FOURTH YEAR PROJECT
# Last Modified: 06/05/2020 16:48
# alltreeplotter2.R

# Script for plotting a single figure, containing six columns and four
# rows of sub-figures. Row one contains all six reconstructed trees,
# row two contains all six corresponding rescaled-reconstructed trees,
# row three contains all six corresponding lineage-through-time curves,
# and row four contains all six corresponding "slopedifferences".

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

library(caper)
library(phytools)
library(paleotree)



# Initialise figure
ORIGsamplingmodels <- c("t_occ", "t_fm", "t_loc", "m_occ", "m_fm", "m_loc")
marstart <- 541
terstart <- 358.9

dev.new(height = 18, width = 25)
split.screen(c(4, 6))


for (i in 1:6) {

	# Plot all six reconstructed (but not yet rescaled) trees from a single simulation run in a row ---------------------------------------
	screen(i)
	par(mar = c(4, 4, 3, 1))

	thisrectree <- tdl_allsix[[i]]
	plot(thisrectree, show.tip.label = FALSE, edge.width = 0.001, edge.color = "#FFD300")
	title(ORIGsamplingmodels[i])

	if (i %in% 1:3) {
		axis(side = 1, at = terstart - abs(c(0, 50, 100, 150, 200, 250, 300, 350, terstart)), labels = c(0, 50, 100, 150, 200, 250, 300, 350, terstart))

	} else if (i %in% 4:6) {
		axis(side = 1, at = marstart - abs(c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, marstart)), labels = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, marstart))

	}

	# Plot axis labels on the first tree of the row
	if (i == 1) {
		title(xlab = "Time before end of simulation / Myrs", cex.lab = 0.8, line = 2)
		title(ylab = "Intermediate Reconstructed Tree", cex.lab = 1, font.lab = 2)
	}


	# Plot the corresponding rescaled-reconstructed trees on the second row ---------------------------------------
	screen(i + 6)
	par(mar = c(4, 4, 1, 1))
	thisresrectree <- trr_allsix[[i]]
	plot(thisresrectree, show.tip.label = FALSE, edge.width = 0.001, edge.color = "#C71585")

	if (i %in% 1:3) {
		axis(side = 1, at = terstart - abs(c(0, 50, 100, 150, 200, 250, 300, 350, terstart)), labels = c(0, 50, 100, 150, 200, 250, 300, 350, terstart))

	} else if (i %in% 4:6) {
		axis(side = 1, at = marstart - abs(c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, marstart)), labels = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, marstart))

	}

	# Plot axis labels on the first tree of the row
	if (i == 1) {
		title(xlab = "Time before end of simulation / Myrs", cex.lab = 0.8, line = 2)
		title(ylab = "Final Reconstructed Tree", cex.lab = 1, font.lab = 2)
	}


	# Plot the corresponding lineages-through-time curves on the third row ---------------------------------------
	screen(i + 12)
	par(mar = c(4, 4, 1, 1))
	sampmodel <- ORIGsamplingmodels[i]
	trueltt <- ltt_allsix[[i]][[1]]
	reconltt <- ltt_allsix[[i]][[2]]
	resrecltt <- ltt_allsix[[i]][[3]]

	# Plot the three terrestrial tetrapod-derived trees' lineage-through-time curves
	if (sampmodel %in% c("t_occ", "t_fm", "t_loc")) {
		plot(x = trueltt[, 2] - (marstart - terstart),
			     y = trueltt[, 1],
			     type = "l",
			     col = "black",
			     xlab = "",
			     ylab = "",
			     xaxs = "i",
			     yaxs = "i",
			     xlim = c(terstart, 0),
			     ylim = c(0, 2 + max(c(trueltt[, 1], reconltt[, 1], resrecltt[, 1]))),
			     bty = "l",
			     lwd = 0.0001,
			     cex.axis = 0.9)

		lines(x = reconltt[, 2] - (marstart - terstart), y = reconltt[, 1], col = "#FFD300", lwd = 0.0001)
		lines(x = resrecltt[, 2] - (marstart - terstart), y = resrecltt[, 1], col = "	#C71585", lwd = 0.0001)

	# Plot the three marine eumetazoa-derived trees' lineage-through-time curves
	} else {
		plot(x = trueltt[, 2],
			     y = trueltt[, 1],
			     type = "l",
			     col = "black",
			     xlab = "",
			     ylab = "",
			     xaxs = "i",
			     yaxs = "i",
			     xlim = c(marstart, 0),
			     ylim = c(0, 2 + max(c(trueltt[, 1], reconltt[, 1], resrecltt[, 1]))),
			     bty = "l",
			     lwd = 0.0001,
			     cex.axis = 0.9)

		lines(x = reconltt[, 2], y = reconltt[, 1], col = "#FFD300", lwd = 0.0001)
		lines(x = resrecltt[, 2], y = resrecltt[, 1], col = "	#C71585", lwd = 0.0001)

	}
	
	# Plot axis labels on the first graph of the row
	if (i == 1) {
		title(xlab = "Time before end of simulation / Myrs", cex.lab = 0.8, line = 2)
		mtext("Number of lineages present", side = 2, line = 2, cex = 0.8)
		mtext("Diversity Patterns", side = 2, line = 3, font = 2)
	}


	# Display the corresponding "slopedifferences" on the fourth row ---------------------------------------
	screen(i + 18)
	slopediffs <- c(0.001597446, 0.00258289, 0.001858116, 0.0001697627, 0.0003667578, 0.0002488586)
	
	par(mar = c(4, 4, 1, 1))
	text(x = 0.5, y = 0.5, labels = slopediffs[i])

	# Draw axis labels on the first sub-plot of the row
	if (i == 1) {
		title(ylab = "Relative Slope Difference", cex.lab = 1, font.lab = 2)
	}

}


# Export figure to pdf
dev.print(pdf, "onepagetrees.pdf", height = 18, width = 25)
