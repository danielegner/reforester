# FOURTH YEAR PROJECT
# Last Modified: 04/05/2020 19:13
# lttplotter.R

# Script for plotting an example set of lineage-through-time curves,
# demonstrating my method of trimming the start and end sections off the
# curves, then rescaling the remaining middle sections so that they all end on
# some arbitrary number (here, three). The gradients of the linear regressions
# of the trimmed, rescaled true and final-reconstructed ltt curves are
# subtracted from each other to give the "relative slope difference" for that
# particular sampling model in this particular simulation run.

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



library(scales)

graphics.off()



# Plot all six lineage-through-time curves ----------------------------------

# Initialise figure
ORIGsamplingmodels <- c("t_occ", "t_fm", "t_loc", "m_occ", "m_fm", "m_loc")
marstart <- 541
terstart <- 358.9

dev.new(height = 6, width = 28)
close.screen(all.screens = TRUE)
split.screen(c(2, 3))


for (i in 1:6) {

	screen(i)
	par(mar = c(4, 4, 1, 1))

	sampmodel <- ORIGsamplingmodels[i]

	trueltt <- ltt_allsix[[i]][[1]]
	reconltt <- ltt_allsix[[i]][[2]]
	resrecltt <- ltt_allsix[[i]][[3]]


	# Plot all three terrestrial tetrapod-derived ltt curves
	if (sampmodel %in% c("t_occ", "t_fm", "t_loc")) {
		plot(x = trueltt[, 2] - (marstart - terstart),
			     y = trueltt[, 1],
			     type = "l",
			     col = "black",
			     xlab = "time before end of simulation",
			     ylab = "number of lineages present",
			     xaxs = "i",
			     yaxs = "i",
			     xlim = c(terstart, 0),
			     ylim = c(0, 2 + max(c(trueltt[, 1], reconltt[, 1], resrecltt[, 1]))),
			     bty = "l",
			     lwd = 0.0001,
			     cex.axis = 0.9)

		lines(x = reconltt[, 2] - (marstart - terstart), y = reconltt[, 1], col = "#FFD300", lwd = 0.0001)
		lines(x = resrecltt[, 2] - (marstart - terstart), y = resrecltt[, 1], col = "	#C71585", lwd = 0.0001)

	# Plot all three marine eumetazoa-derived ltt curves
	} else {
		plot(x = trueltt[, 2],
			     y = trueltt[, 1],
			     type = "l",
			     col = "black",
			     xlab = "time before end of simulation",
			     ylab = "number of lineages present",
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

	title(main = sampmodel)

}

# Export to pdf
dev.print(pdf, "lttplots.pdf", height = 6, width = 18)



# Separately plot the step-by-step lineage-through-time explanation figure ----------------------------------

# Initialise figure
dev.new(height = 10, width = 28)
close.screen(all.screens = TRUE)
split.screen(c(3, 2))

screen(1)

# Plot two columns, one with t_occ and one with m_occ ltts
for (i in c(1, 4)) {

	# Plot true, reconstructed, and rescaled-reconstructd ltts on same graph --------

	# Set placer variable to determine which sub-plot to plot in
	if (i == 1) {
		placer <- 1
	} else if (i == 4) {
		placer <- 2
	}

	screen(placer)
	par(mar = c(4, 4, 1, 1))

	sampmodel <- ORIGsamplingmodels[i]

	trueltt <- ltt_allsix[[i]][[1]]
	reconltt <- ltt_allsix[[i]][[2]]
	resrecltt <- ltt_allsix[[i]][[3]]


	# Plot all three terrestrial tetrapod-derived ltt curves
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

	# Plot all three marine eumetazoa-derived ltt curves
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

	# Plot vertical lines at the start and end cut-off points
	if (i == 1) {
		lines(x = c(20, 20), y = c(0, 40), col = "darkgrey", lwd = 3)
		lines(x = c(terstart - 20, terstart - 20), y = c(0, 40), col = "darkgrey", lwd = 3)

	} else if (i == 4) {
		lines(x = c(20, 20), y = c(0, 40), col = "darkgrey", lwd = 3)
		lines(x = c(marstart - 20, marstart - 20), y = c(0, 40), col = "darkgrey", lwd = 3)
	}

	title(main = sampmodel)


	# On next row, plot only the trimmed true and rescaled-reconstructed ltts --------

	screen(placer + 2)
	par(mar = c(4, 4, 1, 1))

	sampmodel <- ORIGsamplingmodels[i]

	if (i == 1) {
		trim_trueltt <- trueltt[-1 * which(trueltt[, 2] > (marstart - 20)), ]
		trim_trueltt <- trim_trueltt[-1 * which(trim_trueltt[, 2] < ((marstart - terstart) + 20)), ]

		trim_reconltt <- reconltt[-1 * which(reconltt[, 2] > (marstart - 20)), ]
		trim_reconltt <- trim_reconltt[-1 * which(trim_reconltt[, 2] < ((marstart - terstart) + 20)), ]

		trim_resrecltt <- resrecltt[-1 * which(resrecltt[, 2] > (marstart - 20)), ]
		trim_resrecltt <- trim_resrecltt[-1 * which(trim_resrecltt[, 2] < ((marstart - terstart) + 20)), ]

	} else if (i == 4) {
		trim_trueltt <- trueltt[-1 * which(trueltt[, 2] > (marstart - 20)), ]
		trim_trueltt <- trim_trueltt[-1 * which(trim_trueltt[, 2] < 20), ]

		trim_reconltt <- reconltt[-1 * which(reconltt[, 2] > (marstart - 20)), ]
		trim_reconltt <- trim_reconltt[-1 * which(trim_reconltt[, 2] < 20), ]

		trim_resrecltt <- resrecltt[-1 * which(resrecltt[, 2] > (marstart - 20)), ]
		trim_resrecltt <- trim_resrecltt[-1 * which(trim_resrecltt[, 2] < 20), ]

	}

	if (sampmodel %in% c("t_occ", "t_fm", "t_loc")) {
		plot(x = trim_trueltt[, 2] - (marstart - terstart),
			     y = trim_trueltt[, 1],
			     type = "l",
			     col = alpha("black", 0.5),
			     xlab = "",
			     ylab = "Number of lineages present",
			     xaxs = "i",
			     yaxs = "i",
			     xlim = c(terstart, 0),
			     ylim = c(0, 2 + max(c(trim_trueltt[, 1], trim_reconltt[, 1], trim_resrecltt[, 1]))),
			     bty = "l",
			     lwd = 0.0001,
			     cex.axis = 0.9,
			     cex.lab = 0.8)

		lines(x = trim_resrecltt[, 2] - (marstart - terstart), y = trim_resrecltt[, 1], col = alpha("#C71585", 0.5), lwd = 0.0001)

	} else {
		plot(x = trim_trueltt[, 2],
			     y = trim_trueltt[, 1],
			     type = "l",
			     col = alpha("black", 0.5),
			     xlab = "",
			     ylab = "",
			     xaxs = "i",
			     yaxs = "i",
			     xlim = c(marstart, 0),
			     ylim = c(0, 2 + max(c(trim_trueltt[, 1], trim_reconltt[, 1], trim_resrecltt[, 1]))),
			     bty = "l",
			     lwd = 0.0001,
			     cex.axis = 0.9)

		lines(x = trim_resrecltt[, 2], y = trim_resrecltt[, 1], col = alpha("#C71585", 0.5), lwd = 0.0001)

	}


	# Perform linear regression on the trimmed true ltts --------

	xdata <- trim_trueltt[, 2]
	ydata <- trim_trueltt[, 1]

	# Check xdata and ydata are same length
	if (length(xdata) != length(ydata)) {
		cat("Error: Lengths of regression xdata and ydata are not equal.")
		quit(save = "ask")
	}

	# Calculate mean of x data and mean of y data
	xmean <- sum(xdata) / length(xdata)
	ymean <- sum(ydata) / length(ydata)

	# Calculate mean of xy
	eachxy <- matrix(data = NA, ncol = 1, nrow = length(xdata))
	for (j in 1:length(xdata)) {
		eachxy[j] <- xdata[j] * ydata[j]
	}
	xymean <- sum(eachxy) / length(eachxy)

	# Calculate mean of x^2
	eachx2 <- matrix(data = NA, ncol = 1, nrow = length(xdata))
	for (j in 1:length(xdata)) {
		eachx2[j] <- xdata[j]^2
	}
	x2mean <- sum(eachx2) / length(xdata)

	# Calculate slope of linear regression line
	m <- (xymean - xmean * ymean) / (x2mean - xmean^2)

	# Find intercept and therefore calculate y values at the start and end of the trimmed sections, using "y = m*x + c"
	intercept <- ymean  - m * xmean

	linearmodel_true_1 <- m * 20 + intercept
	linearmodel_true_2 <- m * (terstart - 20) + intercept
	

	# Perform linear regression on the trimmed rescaled-reconstructed ltts --------

	xdata <- trim_resrecltt[, 2]
	ydata <- trim_resrecltt[, 1]

	# Check xdata and ydata are same length
	if (length(xdata) != length(ydata)) {
		cat("Error: Lengths of regression xdata and ydata are not equal.")
		quit(save = "ask")
	}

	# Calculate mean of x data and mean of y data
	xmean <- sum(xdata) / length(xdata)
	ymean <- sum(ydata) / length(ydata)

	# Calculate mean of xy
	eachxy <- matrix(data = NA, ncol = 1, nrow = length(xdata))
	for (j in 1:length(xdata)) {
		eachxy[j] <- xdata[j] * ydata[j]
	}
	xymean <- sum(eachxy) / length(eachxy)

	# Calculate mean of x^2
	eachx2 <- matrix(data = NA, ncol = 1, nrow = length(xdata))
	for (j in 1:length(xdata)) {
		eachx2[j] <- xdata[j]^2
	}
	x2mean <- sum(eachx2) / length(xdata)

	# Calculate slope of linear regression line
	m <- (xymean - xmean * ymean) / (x2mean - xmean^2)

	# Find intercept and therefore calculate y values at the start and end of the trimmed sections, using "y = m*x + c"
	intercept <- ymean  - m * xmean

	linearmodel_resrec_1 <- m * 20 + intercept
	linearmodel_resrec_2 <- m * (marstart - 20) + intercept


	# Superimpose the vertical trimming lines and the linear regression lines --------

	if (i == 1) {
		lines(x = c(20, 20), y = c(0, 40), col = "darkgrey", lwd = 3)
		lines(x = c(terstart - 20, terstart - 20), y = c(0, 40), col = "darkgrey", lwd = 3)

		lines(x = c(20, terstart - 20), y = c(linearmodel_true_1, linearmodel_true_2), col = "black")
		lines(x = c(20, terstart - 20), y = c(linearmodel_resrec_1, linearmodel_resrec_2), col = "#C71585")

	} else if (i == 4) {
		lines(x = c(20, 20), y = c(0, 40), col = "darkgrey", lwd = 3)
		lines(x = c(marstart - 20, marstart - 20), y = c(0, 40), col = "darkgrey", lwd = 3)

		lines(x = c(20, marstart - 20), y = c(linearmodel_true_1, linearmodel_true_2), col = "black")
		lines(x = c(20, marstart - 20), y = c(linearmodel_resrec_1, linearmodel_resrec_2), col = "#C71585")

	}


	# Re-scale the trimmed ltts, making them all end on an arbitrary number, in this case three  --------

	x_true <- trim_trueltt[length(trim_trueltt[, 1]), 1] / 3
	trim_trueltt[, 1] <- trim_trueltt[, 1] / x_true

	x_resrec <- trim_resrecltt[length(trim_resrecltt[, 1]), 1] / 3
	trim_resrecltt[, 1] <- trim_resrecltt[, 1] / x_resrec


	# On third row, plot these re-scaled and trimmed ltts --------

	screen(placer + 4)
	par(mar = c(4, 4, 1, 1))

	if (sampmodel %in% c("t_occ", "t_fm", "t_loc")) {
		plot(x = trim_trueltt[, 2] - (marstart - terstart),
			     y = trim_trueltt[, 1],
			     type = "l",
			     col = alpha("black", 0.5),
			     xlab = "Time before end of simulation / Myrs",
			     ylab = "",
			     xaxs = "i",
			     yaxs = "i",
			     xlim = c(terstart, 0),
			     ylim = c(0, 6),
			     bty = "l",
			     lwd = 0.0001,
			     cex.axis = 0.9,
			     cex.lab = 0.8)

		lines(x = trim_resrecltt[, 2] - (marstart - terstart), y = trim_resrecltt[, 1], col = alpha("#C71585", 0.5), lwd = 0.0001)

	} else {
		plot(x = trim_trueltt[, 2],
			     y = trim_trueltt[, 1],
			     type = "l",
			     col = alpha("black", 0.5),
			     xlab = "Time before end of simulation / Myrs",
			     ylab = "",
			     xaxs = "i",
			     yaxs = "i",
			     xlim = c(marstart, 0),
			     ylim = c(0, 6),
			     bty = "l",
			     lwd = 0.0001,
			     cex.axis = 0.9,
			     cex.lab = 0.8)

		lines(x = trim_resrecltt[, 2], y = trim_resrecltt[, 1], col = alpha("#C71585", 0.5), lwd = 0.0001)

	}


	# Perform linear regression on the re-scaled trimmed true ltts --------

	xdata <- trim_trueltt[, 2]
	ydata <- trim_trueltt[, 1]

	# Check xdata and ydata are same length
	if (length(xdata) != length(ydata)) {
		cat("Error: Lengths of regression xdata and ydata are not equal.")
		quit(save = "ask")
	}

	# Calculate mean of x data and mean of y data
	xmean <- sum(xdata) / length(xdata)
	ymean <- sum(ydata) / length(ydata)

	# Calculate mean of xy
	eachxy <- matrix(data = NA, ncol = 1, nrow = length(xdata))
	for (j in 1:length(xdata)) {
		eachxy[j] <- xdata[j] * ydata[j]
	}
	xymean <- sum(eachxy) / length(eachxy)

	# Calculate mean of x^2
	eachx2 <- matrix(data = NA, ncol = 1, nrow = length(xdata))
	for (j in 1:length(xdata)) {
		eachx2[j] <- xdata[j]^2
	}
	x2mean <- sum(eachx2) / length(xdata)

	# Calculate slope of linear regression line
	m <- (xymean - xmean * ymean) / (x2mean - xmean^2)

	# Find intercept and therefore calculate y values at the start and end of the trimmed sections, using "y = m*x + c"
	intercept <- ymean  - m * xmean

	linearmodel_true_1 <- m * 20 + intercept
	linearmodel_true_2 <- m * (terstart - 20) + intercept
	

	# Perform linear regression on the re-scaled trimmed rescaled-reconstructed ltts --------

	xdata <- trim_resrecltt[, 2]
	ydata <- trim_resrecltt[, 1]

	# Check xdata and ydata are same length
	if (length(xdata) != length(ydata)) {
		cat("Error: Lengths of regression xdata and ydata are not equal.")
		quit(save = "ask")
	}

	# Calculate mean of x data and mean of y data
	xmean <- sum(xdata) / length(xdata)
	ymean <- sum(ydata) / length(ydata)

	# Calculate mean of xy
	eachxy <- matrix(data = NA, ncol = 1, nrow = length(xdata))
	for (j in 1:length(xdata)) {
		eachxy[j] <- xdata[j] * ydata[j]
	}
	xymean <- sum(eachxy) / length(eachxy)

	# Calculate mean of x^2
	eachx2 <- matrix(data = NA, ncol = 1, nrow = length(xdata))
	for (j in 1:length(xdata)) {
		eachx2[j] <- xdata[j]^2
	}
	x2mean <- sum(eachx2) / length(xdata)

	# Calculate slope of linear regression line
	m <- (xymean - xmean * ymean) / (x2mean - xmean^2)

	# Find intercept and therefore calculate y values at the start and end of the trimmed sections, using "y = m*x + c"
	intercept <- ymean  - m * xmean

	linearmodel_resrec_1 <- m * 20 + intercept
	linearmodel_resrec_2 <- m * (marstart - 20) + intercept


	# Superimpose the vertical trimming lines and the linear regression lines --------

	if (i == 1) {
		lines(x = c(20, 20), y = c(0, 6), col = "darkgrey", lwd = 3)
		lines(x = c(terstart - 20, terstart - 20), y = c(0, 6), col = "darkgrey", lwd = 3)

		lines(x = c(20, terstart - 20), y = c(linearmodel_true_1, linearmodel_true_2), col = "black")
		lines(x = c(20, terstart - 20), y = c(linearmodel_resrec_1, linearmodel_resrec_2), col = "#C71585")

	} else if (i == 4) {
		lines(x = c(20, 20), y = c(0, 6), col = "darkgrey", lwd = 3)
		lines(x = c(marstart - 20, marstart - 20), y = c(0, 6), col = "darkgrey", lwd = 3)

		lines(x = c(20, marstart - 20), y = c(linearmodel_true_1, linearmodel_true_2), col = "black")
		lines(x = c(20, marstart - 20), y = c(linearmodel_resrec_1, linearmodel_resrec_2), col = "#C71585")

	}

}

# Export to pdf
dev.print(pdf, "lttplots_stepbystep.pdf", height = 10, width = 18)
