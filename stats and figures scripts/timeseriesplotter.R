# FOURTH YEAR PROJECT
# Last Modified: 07/05/2020 12:15
# timeseriesplotter.R

# Script for plotting the raw terrestrial tetrapod and marine eumetazoa
# time series data, as well as their log-transformed versions and log-
# transformed versions that were re-scaled to make their means equal to 0.5.
# This demonstrates the process I used to take the raw fossil sampling rate
# proxy data and turn it into a form useable as sampling rate models in
# my simulations.

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
dev.new(height = 6, width = 10)
close.screen(all.screens = TRUE)
split.screen(c(3, 2))

importfilenames <- c("OUTPUT_terrestrial_tetrapods.csv", "OUTPUT_marine_eumetazoa.csv", "OUTPUT_terrestrial_tetrapods.csv", "OUTPUT_marine_eumetazoa.csv", "OUTPUT_terrestrial_tetrapods.csv", "OUTPUT_marine_eumetazoa.csv")
marstart <- 541
terstart <- 358.9


# Two columns (one for terrestrial tetrapod time series and one for marine eumetazoa time series), three rows
for (outerloop in 1:6) {

	screen(outerloop)
	par(mar = c(4, 4, 1, 1))

	importfilename <- importfilenames[outerloop]
	colours_ter <- c("#7B3F00", "#228B22", "#FFD300")
	colours_mar <- c("#0087BD", "#15F4EE", "#FF1493")


	# For first row, plot raw fossil data time series ------
	if (outerloop %in% c(1, 2)) {

		for (j in 1:3) {

			if (j == 1) {
				plot(x = 0, y = 0, type = "n", xlim = c(0, 541), ylim = c(0, 40000), bty = "l", xlab = "", ylab = "", cex.axis = 0.9, xaxt = "n")

				if (outerloop == 1) {
					title(main = "Terrestrial Tetrapod Time Series")
				} else if (outerloop == 2) {
					title(main = "Marine Eumetazoa Time Series")
				}

			}

			if (outerloop == 1) {
				mtext(text = "Raw Fossil Data", side = 2, line = 3, font = 2)

				if (j == 1) {
					title(ylab = "Sampling proxy", cex.lab = 0.8, line = 2)
					legend(x = "topleft", legend = c("t_occ", "t_fm", "t_loc", "m_occ", "m_fm", "m_loc"), lty = 1, col = c("#7B3F00", "#228B22", "#FFD300", "#0087BD", "#15F4EE", "#FF1493"), cex = 0.7, lwd = 2, bty = "n")
				}
			}

			whichcol <- j + 3

			importeddata <- as.data.frame(read.csv(file = paste(importfilename), header = TRUE, sep = ","))
			samplingproxy <- importeddata[[whichcol]]

			# Fill first two columns of a 3-column matrix with interval start and end times respectively
			intervalsinfo <- matrix(data = NA, nrow = length(importeddata[, 1]), ncol = 3)
			intervalsinfo[, 1] <- importeddata[, 2]
			intervalsinfo[, 2] <- importeddata[, 3]

			# Set each time bin's end time as equal to the start time of the following bin.
			# This closes the temporal gaps between fossil data time bins (these gaps are errors and shouldn't exist). 
			for (i in 1:(length(intervalsinfo[, 1]) - 1)) {
				intervalsinfo[i, 2] <- intervalsinfo[i + 1, 1]
			}

			# Fill third column with real-data-derived sampling rates
			intervalsinfo[, 3] <- samplingproxy

			plotinfo <- matrix(data = NA, ncol = 2, nrow = 2 * nrow(intervalsinfo))
			plotinfo[seq(from = 1, to = nrow(plotinfo), by = 2), ] <- intervalsinfo[, c(1, 3)]
			plotinfo[seq(from = 2, to = nrow(plotinfo), by = 2), ] <- intervalsinfo[, c(2, 3)]

			# Plot axis label and legend on first column (terrestrial tetrapod) graph, plot axes on both columns
			if (outerloop == 1) {
				lines(x = abs(plotinfo[, 1] - 541), y = plotinfo[, 2], col = colours_ter[j])
				if (j == 1) {
					axis(side = 1, at = 182.1 + c(-182.1, -141.1, -91.1, -41.1, 0, 8.9, 58.9, 108.9, 158.9, 208.9, 258.9, 308.9, 358.9), labels = c(541, 500, 450, 400, 358.9, 350, 300, 250, 200, 150, 100, 50, 0), cex.axis = 0.9)
				}
			} else if (outerloop == 2) {
				lines(x = abs(plotinfo[, 1] - 541), y = plotinfo[, 2], col = colours_mar[j])
				if (j == 1) {
					axis(side = 1, at = c(0, 41, 91, 141, 191, 241, 291, 341, 391, 441, 491, 541), labels = c(541, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0), cex.axis = 0.9)
				}
			}

		}

	# For second row, plot log-transformed time series ------
	} else if (outerloop %in% c(3, 4)) {

		# Plot axes and some labels
		for (j in 1:3) {

			if (j == 1) {
				plot(x = 0, y = 0, type = "n", xlim = c(0, 541), ylim = c(0, 5), bty = "l", xlab = "", ylab = "", cex.axis = 0.9, xaxt = "n")
			}

			if (outerloop == 3) {
				mtext(text = "Log-Transformed", side = 2, line = 3, font = 2)

				if (j == 1) {
					title(ylab = "Sampling proxy", cex.lab = 0.8, line = 2)
				}
			}

			whichcol <- j + 3

			importeddata <- as.data.frame(read.csv(file = paste(importfilename), header = TRUE, sep = ","))
			samplingproxy <- importeddata[[whichcol]]

			# Log-transform the time series
			samplingrates <- log10(samplingproxy)

			# Fill first two columns of a 3-column matrix with interval start and end times respectively
			intervalsinfo <- matrix(data = NA, nrow = length(importeddata[, 1]), ncol = 3)
			intervalsinfo[, 1] <- importeddata[, 2]
			intervalsinfo[, 2] <- importeddata[, 3]

			# Set each time bin's end time as equal to the start time of the following bin.
			# This closes the temporal gaps between fossil data time bins (these gaps are errors and shouldn't exist). 
			for (i in 1:(length(intervalsinfo[, 1]) - 1)) {
				intervalsinfo[i, 2] <- intervalsinfo[i + 1, 1]
			}

			# Fill third column with real-data-derived sampling rates
			intervalsinfo[, 3] <- samplingrates

			plotinfo <- matrix(data = NA, ncol = 2, nrow = 2 * nrow(intervalsinfo))
			plotinfo[seq(from = 1, to = nrow(plotinfo), by = 2), ] <- intervalsinfo[, c(1, 3)]
			plotinfo[seq(from = 2, to = nrow(plotinfo), by = 2), ] <- intervalsinfo[, c(2, 3)]

			# Plot axis label and legend on first column (terrestrial tetrapod) graph, plot axes on both columns
			if (outerloop == 3) {
				lines(x = abs(plotinfo[, 1] - 541), y = plotinfo[, 2], col = colours_ter[j])
				if (j == 1) {
					axis(side = 1, at = 182.1 + c(-182.1, -141.1, -91.1, -41.1, 0, 8.9, 58.9, 108.9, 158.9, 208.9, 258.9, 308.9, 358.9), labels = c(541, 500, 450, 400, 358.9, 350, 300, 250, 200, 150, 100, 50, 0), cex.axis = 0.9)
				}
			} else if (outerloop == 4) {
				lines(x = abs(plotinfo[, 1] - 541), y = plotinfo[, 2], col = colours_mar[j])
				if (j == 1) {
					axis(side = 1, at = c(0, 41, 91, 141, 191, 241, 291, 341, 391, 441, 491, 541), labels = c(541, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0), cex.axis = 0.9)
				}
			}

		}

	# For third row, plot re-scaled log-transformed time series ------
	} else if (outerloop %in% c(5, 6)) {

		for (j in 1:3) {

			# Plot axes and some labels
			if (j == 1) {
				plot(x = 0, y = 0, type = "n", xlim = c(0, 541), ylim = c(0, 1), bty = "l", xlab = "Time / Ma", ylab = "", cex.axis = 0.9, cex.lab = 0.8, xaxt = "n")

			}

			if (outerloop == 5) {
				mtext(text = "Final (Mean Set To 0.5)", side = 2, line = 3, font = 2)

				if (j == 1) {
					title(ylab = "Sampling rate", cex.lab = 0.8, line = 2)
				}
			}

			whichcol <- j + 3

			importeddata <- as.data.frame(read.csv(file = paste(importfilename), header = TRUE, sep = ","))
			samplingproxy <- importeddata[[whichcol]]

			# Log-transform the time series
			samplingrates <- log10(samplingproxy)

			# Re-scale the log-transformed time series such that its mean == 0.5 and its max value <= 1
			x <- mean(samplingrates) / 0.5
			if (x >= max(samplingrates)) {
				samplingrates <- samplingrates / x
			} else {
				cat("Error: Sampling rate in one or more time bins is greater than 1.")
				quit(save = "ask")
			}

			# Fill first two columns of a 3-column matrix with interval start and end times respectively
			intervalsinfo <- matrix(data = NA, nrow = length(importeddata[, 1]), ncol = 3)
			intervalsinfo[, 1] <- importeddata[, 2]
			intervalsinfo[, 2] <- importeddata[, 3]

			# Set each time bin's end time as equal to the start time of the following bin.
			# This closes the temporal gaps between fossil data time bins (these gaps are errors and shouldn't exist). 
			for (i in 1:(length(intervalsinfo[, 1]) - 1)) {
				intervalsinfo[i, 2] <- intervalsinfo[i + 1, 1]
			}

			# Fill third column with real-data-derived sampling rates
			intervalsinfo[, 3] <- samplingrates

			plotinfo <- matrix(data = NA, ncol = 2, nrow = 2 * nrow(intervalsinfo))
			plotinfo[seq(from = 1, to = nrow(plotinfo), by = 2), ] <- intervalsinfo[, c(1, 3)]
			plotinfo[seq(from = 2, to = nrow(plotinfo), by = 2), ] <- intervalsinfo[, c(2, 3)]

			# Plot axis label and legend on first column (terrestrial tetrapod) graph, plot axes on both columns
			if (outerloop == 5) {
				lines(x = abs(plotinfo[, 1] - 541), y = plotinfo[, 2], col = colours_ter[j])
				if (j == 1) {
					axis(side = 1, at = 182.1 + c(-182.1, -141.1, -91.1, -41.1, 0, 8.9, 58.9, 108.9, 158.9, 208.9, 258.9, 308.9, 358.9), labels = c(541, 500, 450, 400, 358.9, 350, 300, 250, 200, 150, 100, 50, 0), cex.axis = 0.9)
				}
			} else if (outerloop == 6) {
				lines(x = abs(plotinfo[, 1] - 541), y = plotinfo[, 2], col = colours_mar[j])
				if (j == 1) {
					axis(side = 1, at = c(0, 41, 91, 141, 191, 241, 291, 341, 391, 441, 491, 541), labels = c(541, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0), cex.axis = 0.9)
				}
			}

		}
	}

}


# Export to pdf
dev.print(pdf, "timeseriesplots.pdf", height = 6, width = 10)
