# FOURTH YEAR PROJECT
# Last Modified: 20/04/2020 15:03
# statstester_uniformsampling.R

# Script for calculating various statistical parameters of the two uniform
# sampling model-derived "relative slope difference" results sets. These
# include checking for normality by plotting histograms, plotting Q-Q graphs,
# and using the Shapiro-Wilk test. Also runs one-sample t-tests and calculates
# means, medians, standard deviaitons, SEMs and interquartile ranges for each set.

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



# Close all graphics windows
graphics.off()


# Assign raw results data
absresults <- onefiftyABS   # Not currently used, only looking at the relative slope differences, not the absolute differences
relresults <- onefiftyREL

samplingmodels <- c("t_occ", "m_occ")


# For each sampling model's results, check for normality by plotting a histogram with a normal distribution superimposed
dev.new(height = 6, width = 8)
close.screen(all.screens = TRUE)
split.screen(c(2, 3))
for (i in 1:length(relresults[1, ])) {

	screen(i)
	par(mar = c(4, 4, 1, 1))

	h <- hist(relresults[, i], breaks = 10, density = 10, col = "lightgray", xlab = "",  ylab = "", main = paste("samp model:", samplingmodels[i], sep = " "))
	xfit <- seq(min(relresults[, i]), max(relresults[, i]), length = 100) 
	yfit <- dnorm(xfit, mean = mean(relresults[, i]), sd = sd(relresults[, i])) 
	yfit <- yfit * diff(h$mids[1:2]) * length(relresults[, i]) 
	lines(xfit, yfit, col = "black", lwd = 2)
	title(xlab = "Slope difference", ylab = "Frequency", line = 2.5, cex.lab = 1)

}


# For each sampling model's results, check for normality by plotting Q-Q graphs
dev.new(height = 6, width = 8)
close.screen(all.screens = TRUE)
split.screen(c(2, 3))
for (i in 1:length(relresults[1, ])) {

	screen(i)
	par(mar = c(4, 4, 1, 1))

	qqnorm(relresults[, i], main = paste("samp model:", samplingmodels[i], sep = " "), xlab = "", ylab = "")
	qqline(relresults[, i])
	title(xlab = "Theoretical quantiles", ylab = "Sample quantiles", line = 2.5, cex.lab = 1)

}


# For each sampling model's results, check for normality with S-W test and perform a one-sample t-test
shapwilks <- list(NA, NA, NA, NA, NA, NA)
swpvalue <- matrix(data = NA, ncol = 1, nrow = length(relresults[1, ]))
tvalues <- matrix(data = NA, ncol = 1, nrow = length(relresults[1, ]))
pvalues <- matrix(data = NA, ncol = 1, nrow = length(relresults[1, ]))
for (i in 1:length(relresults[1, ])) {

	shapwilks[[i]] <- shapiro.test(relresults[, i])
	swpvalue[i] <- shapwilks[[i]]$p.value

	tvalues[i] <- (mean(relresults[, i]) - 0) / (sd(relresults[, i]) / sqrt(length(relresults[, i])))
	pvalues[i] <- 2 * pt(q = -abs(tvalues[i]), df = length(relresults[, i]) - 1)

	# The line below gives the same t- and p-values as the manual calculations above
	# current_ttest <- t.test(relresults[, i], alternative = "two.sided", mu = 0, conf.level = 0.99)

}


# Calculate other descriptive statistics for each set of results
means <- matrix(data = NA, ncol = 1, nrow = length(relresults[1, ]))
medians <- matrix(data = NA, ncol = 1, nrow = length(relresults[1, ]))
stddevs <- matrix(data = NA, ncol = 1, nrow = length(relresults[1, ]))
sems <- matrix(data = NA, ncol = 1, nrow = length(relresults[1, ]))
iqrs <- matrix(data = NA, ncol = 1, nrow = length(relresults[1, ]))
for (i in 1:length(relresults[1, ])) {

	means[i] <- sum(relresults[, i]) / length(relresults[, i])
	medians[i] <- median(relresults[, i])
	stddevs[i] <- sd(relresults[, i])
	sems[i] <- stddevs[i] / sqrt(length(relresults[, i]))
	iqrs[i] <- IQR(relresults[, i])

}


# Save raw results to output file
outputfilename <- "usrelresults.csv"
outputfile <- file(outputfilename, "w")
write.table(relresults, outputfile, quote = FALSE, row.names = FALSE, sep = " & ")
close(outputfile)


# Collate and save results to output file
collateddata <- data.frame(mean = means, median = medians, sd = stddevs, SEM = sems, IQR = iqrs, swtestpvalue = swpvalue, tvalue = tvalues, ttestpvalue = pvalues)
collateddata <- t(collateddata)

outputfilename <- "usanalysis.csv"
outputfile <- file(outputfilename, "w")
write.table(formatC(collateddata, format = "g", digits = 5), outputfile, quote = FALSE, row.names = FALSE, sep = " & ")
close(outputfile)
