# FOURTH YEAR PROJECT
# Last Modified: 20/03/2020 13:29
# dataimporter.R

# This script takes the fossil database files "terrestrial_tetrapods.csv"
# and "marine_eumetazoa.csv", extracts the required fossil occurrence data
# from them (divided by time bin), and saves the resulting filtered data
# as "OUTPUT_terrestrial_tetrapods.csv" and "OUTPUT_marine_eumetazoa.csv".

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



# Functions -----------------------------------------------------

# Get vector of unique time bin labels
GetUniqueTimeBins <- function(dat, bincol, maxagecol, minagecol) {

	uniquebins <- levels(dat[, bincol])

	# Get mode (most common) bin start and end times (in order to ignore typos in the database)
	binstarts <- matrix(data = NA, ncol = 1, nrow = length(uniquebins))
	binends <- matrix(data = NA, ncol = 1, nrow = length(uniquebins))
	for (i in 1:length(uniquebins)) {
		thisbin.starttimes <- matrix(data = NA, ncol = 1, nrow = length(dat[, maxagecol]))
		thisbin.endtimes <- matrix(data = NA, ncol = 1, nrow = length(dat[, minagecol]))
		for (j in 1:length(dat[, bincol])) {
			if (uniquebins[i] == paste(dat[j, bincol])) {
				thisbin.starttimes[j] <- dat[j, maxagecol]
				thisbin.endtimes[j] <- dat[j, minagecol]
			}
		}

		# Get mode start and end time for the bin
		thisbin.starttimes <- thisbin.starttimes[!is.na(thisbin.starttimes)]
		u <- unique(thisbin.starttimes)
		binstarts[i] <- u[which.max(tabulate(match(thisbin.starttimes, u)))]

		thisbin.endtimes <- thisbin.endtimes[!is.na(thisbin.endtimes)]
		u <- unique(thisbin.endtimes)
		binends[i] <- u[which.max(tabulate(match(thisbin.endtimes, u)))]
	}

	# Collate info and re-order by descending start time (i.e. chronological time order)
	timebins <- data.frame(binlabel = uniquebins, starttime = binstarts, endtime = binends)
	timebins <- timebins[order(timebins[["starttime"]], decreasing = TRUE), ]

	return(timebins)

}


# Get number of occurences per unique time bin
CountBinOccs <- function(bins, dat, bincol) {

	binocccounts <- matrix(data = 0, ncol = 1, nrow = length(bins))

	for (i in 1:length(dat[, bincol])) {
		binocccounts[which(bins == paste(dat[i, bincol]))] <- binocccounts[which(bins == paste(dat[i, bincol]))] + 1
	}

	return(binocccounts)

}


# Get number of formations per unique time bin
CountBinFms <- function(bins, dat, bincol, fmcol, datafilename) {

	# Get list of lists of each unique formation name in each time bin
	fmsbybin <- vector("list", length(bins))
	for (i in 1:length(bins)) {
		for (j in 1:length(dat[, fmcol])) {
			if (bins[i] == paste(dat[j, bincol])) {
				if (paste(dat[j, fmcol]) != "NA") {
					if (paste(dat[j, fmcol]) %in% fmsbybin[[i]] == FALSE) {
						fmsbybin[[i]] <- c(fmsbybin[[i]], paste(dat[j, fmcol]))
					}
				}
			}
		}
	}

	# Order each bin's fomrations alphabetically for easier manual typo counting
	for (i in 1:length(fmsbybin)) {
		fmsbybin[[i]] <- fmsbybin[[i]][order(fmsbybin[[i]])]
	}

	# Count number of unique formations in each bin
	binfmcounts <- matrix(data = NA, ncol = 1, nrow = length(fmsbybin))
	for (i in 1:length(fmsbybin)) {
		binfmcounts[i] <- length(fmsbybin[[i]])
	}

	# Remove typo-caused duplicates from each bin's formation count using manually-counted values
	if (datafilename == "terrestrial_tetrapods.csv") {
		fmtypos <- c(0, 0, 0, 0, 0, 6, 0, 3, 0, 0,
				 5, 0, 1, 1, 0, 0, 0, 0, 0, 0,
				 0, 0, 0, 2, 3, 1, 0, 5, 4, 2,
				 9, 3, 9, 2, 3, 2, 14, 15, 11)
		binfmcounts <- binfmcounts - fmtypos
	} else if (datafilename == "marine_eumetazoa.csv") {
		fmtypos <- c(1, 19, 0, 1, 17, 1, 32, 6, 41, 28,
				 85, 29, 9, 35, 11, 10, 7, 13, 2, 7,
				 21, 20, 18, 24, 14, 13, 21, 9, 12, 7,
				 12, 10, 9, 27, 13, 4, 5, 2, 15, 2,
				 22, 7, 23, 9, 10, 37, 29, 25, 35)
		binfmcounts <- binfmcounts - fmtypos
	}

	return(binfmcounts)

}


# Get number of localities (i.e. unique collection no.s) per unique time bin
CountBinLocs <- function(bins, dat, bincol, loccol) {

	locsbybin <- vector("list", length(bins))

	for (i in 1:length(bins)) {
		for (j in 1:length(dat[, loccol])) {
			if (bins[i] == paste(dat[j, bincol])) {
				if (paste(dat[j, loccol]) != "NA") {
					if (paste(dat[j, loccol]) %in% locsbybin[[i]] == FALSE) {
						locsbybin[[i]] <- c(locsbybin[[i]], paste(dat[j, loccol]))
					}
				}
			}
		}
	}
	
	binloccounts <- matrix(data = NA, ncol = 1, nrow = length(locsbybin))
	for (i in 1:length(locsbybin)) {
		binloccounts[i] <- length(locsbybin[[i]])
	}

	return(binloccounts)

}



# Main -----------------------------------------------------

datafilenames <- c("terrestrial_tetrapods.csv", "marine_eumetazoa.csv")

for (i in 1:length(datafilenames)) {

	# Check which file data is coming from to get data from correct columns
	if (datafilenames[i] == "terrestrial_tetrapods.csv") {
		bincolumn <- 2
		maxagecolumn <- 18
		minagecolumn <- 19
		fmcolumn <- 85
		loccolumn <- 1
	} else if (datafilenames[i] == "marine_eumetazoa.csv") {
		bincolumn <- 3
		maxagecolumn <- 15
		minagecolumn <- 16
		fmcolumn <- 81
		loccolumn <- 1
	}

	# Import data csv file
	data <- read.csv(datafilenames[i], header = TRUE)

	# Run data collection functions
	timebins <- GetUniqueTimeBins(dat = data, bincol = bincolumn, maxagecol = maxagecolumn, minagecol = minagecolumn)

	occnums <- CountBinOccs(bins = timebins[["binlabel"]], dat = data, bincol = bincolumn)

	fmnums <- CountBinFms(bins = timebins[["binlabel"]], dat = data, bincol = bincolumn, fmcol = fmcolumn, datafilename = datafilenames[i])

	locnums <- CountBinLocs(bins = timebins[["binlabel"]], dat = data, bincol = bincolumn, loccol = loccolumn)

	# Gather all data into one dataframe
	collateddata <- data.frame(timebins, numofoccurrences = occnums, numofformations = fmnums, numoflocalities = locnums)

	# Initialise, write, and close output file
	outputfilename <- paste("OUTPUT_", datafilenames[i], sep = "")
	outputfile <- file(outputfilename, "w")
	write.csv(collateddata, outputfile, quote = FALSE, row.names = FALSE)
	close(outputfile)

}
