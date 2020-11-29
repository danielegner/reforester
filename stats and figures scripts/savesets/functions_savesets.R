# FOURTH YEAR PROJECT
# Last Modified: 05/04/2020 12:16
# functions_savesets.R

# This script is almost the same as "functions.R", with some modifications
# to make it display the true ltt curve alongside all six reconstructed
# and all six rescaled-reconstructed ltt curves, for a figure in the
# dissertation. It is called by "main_savesets.R".

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



GrowTree <- function(diversificationmodel) {

	# Grow tree based on given diversification model
	if (diversificationmodel == "diversitydependent") {
		# The halt arguement used excludes root edge from simulation growth time.
		# However, the root node is never exactly at 541 Ma due to simulation time-stepping.
		# Thus, must cut-off small amounts of time at the end of each simulation (see trimlength).
		
		# Note: growTree()'s 'halt' arguement cannot accept marstart for some reason
		if (marstart == 541) {
			tree_orig <- growTree(b = expression(1 - nExtantTip / 60),
						    d = expression(nExtantTip / 60),
						    grain = Inf,
						    halt = expression(clade.age - lin.age[1] >= 541))
		} else if (marstart == 54.1) {
			tree_orig <- growTree(b = expression(1 - nExtantTip / 60),
						    d = expression(nExtantTip / 60),
						    grain = Inf,
						    halt = expression(clade.age - lin.age[1] >= 54.1))
		} else if (marstart == 5.41) {
			tree_orig <- growTree(b = expression(1 - nExtantTip / 60),
						    d = expression(nExtantTip / 60),
						    grain = Inf,
						    halt = expression(clade.age - lin.age[1] >= 5.41))
		}

	} else if (diversificationmodel == "exponential") {
		# Currently unsupported, infeasible to run exponential simulation for 541 Myrs, so will cause errors
		tree_orig <- growTree(b = 0.3,
					    d = 0.05,
					    grain = Inf,
					    halt = expression(nExtantTip >= 20))

	} else {
		cat("Error: Neither 'exponential' nor 'diversitydependent' diversification model set.")
		quit(save = "ask")
	}

	# Ladderize and "kick" tree (i.e. add and remove a 0-length edge to force new node/edge/tip labeling convention)
	tree_lad <- ladderize(tree_orig[["phy"]], right = FALSE)

	tree_addedtip <- bind.tip(tree_lad,
					  tip.label = "sometip",
					  edge.length = 0,
					  where = min(tree_lad[["node.label"]]),
					  position = 0)

	tree_kicked <- drop.tip(phy = tree_addedtip, tip = which(tree_addedtip[["tip.label"]] == "sometip"))

	# Add root.edge component to kicked tree
	tree_kicked[["root.edge"]] <- tree_orig[["phy"]][["root.edge"]]

	# Plot tree_kicked
	if (createplots == TRUE) {
		dev.new()
		plot(tree_kicked, show.tip.label = FALSE)

		# Plot both sets of axes (i.e. for both the marine and terrestrial datasets)
		axis(side = 1, at = seq(0, marstart, length.out = 50), labels = round(seq(marstart, 0, length.out = 50), 2))
		axis(side = 1, at = seq(0, terstart, length.out = 50), labels = round(seq(terstart, 0, length.out = 50), 2), line = 3)

		title("tree_kicked")
	}

	outputs <- list(tree_orig, tree_kicked)
	return(outputs)

}


# This function relies on tree$node.data
# Ladderize(tree) requires removal of tree$node.data, so below function can only be used on the original unaltered tree_orig
VerifyTreeLengths <- function(tree) {

	# Copy tree$phy$edge
	edges <- tree[["phy"]]["edge"]

	# Copy root length
	rootlength <- tree[["phy"]]["root.edge"][[1]]

	# Copy node numbers and associated death.time values (i.e. node ages) from node.data
	nodedata <- tree[["node.data"]]   # Doesn't include root node
	nodeages <- data.frame(nodedata[1], nodedata[4], row.names = NULL)
	colnames(nodeages) <- c("node", "age")

	# Find earliest node from node.data, use its birth.time as the root node's age, and append root node to nodeages.
	# This is not necessary, but can be useful for debugging.
	secondnodenum <- min(nodeages[, 1], na.rm = TRUE)
	rootnodenum <- secondnodenum - 1
	append <- data.frame(rootnodenum, nodedata[which(nodedata["node"] == secondnodenum), "birth.time"])
	colnames(append) <- c("node", "age")
	nodeages <- rbind(nodeages, append)   # Now includes root node

	# Calculate total simulation length...
	finalnodenum <- max(nodeages[, 1], na.rm = TRUE)   # In unaltered trees, the greatest node number is always given to the clostest-to-extant node
	totallength <- nodeages["age"][which(nodeages["node"] == finalnodenum), 1]

	# ... Ensuring that chosen tip is extant
	finaldaughteredges <- which(edges[[1]][, 1] == finalnodenum)
	picktipnum <- edges[[1]][finaldaughteredges[1], 2]
	if (tree[["data"]]["extinct"][picktipnum, 1] == "FALSE") {
		addlength <- tree[["data"]]["lin.age"][picktipnum, 1]
		totallength <- totallength + addlength
	} else {
		picktipnum <- edges[[1]][finaldaughteredges[2], 2]
		if (tree[["data"]]["extinct"][picktipnum, 1] == "FALSE") {
			addlength <- tree[["data"]]["lin.age"][picktipnum, 1]
			totallength <- totallength + addlength
		} else {
			cat("Error: Neither of finalnode's daughter edges are extant, can't calculate total tree length.")
			quit(save = "ask")
		}
	}

	outputs <- list(totallength, rootlength)   # totallength includes the root length, but is NOT trimmed!
	return(outputs)

}


# Does not rely on tree$node.data so can be used on all modified trees to get edge times
GetEdgeTimes <- function(phy, totlength) {

	# Copy phy$edge and make empty copy
	edges <- phy[["edge"]]
	edgetimes <- matrix(data = NA, nrow = nrow(edges), ncol = ncol(edges))

	# Get edge lengths
	edgelengths <- phy[["edge.length"]]

	# Get root edge's length and root node's number
	rootlength <- phy[["root.edge"]]
	rootnodenum <- length(phy[["tip.label"]]) + 1

	# Set root node's corresponding start edgetime to length of root (the undrawn initial edge leading into the root node)
	edgetimes[edges[, 1] == rootnodenum, 1] <- rootlength

	# Loop through all nodes, filling-in the corresponding two daughter edges' start and end edgetimes
	allnodenums <- rootnodenum:((rootnodenum - 1) + phy[["Nnode"]])
	for (currentnodenum in allnodenums) {
		currentedges.nums <- which(edges[, 1] == currentnodenum)
		currentedges.starttimes <- edgetimes[currentedges.nums, 1]
		currentedges.lengths <- edgelengths[currentedges.nums]
		currentedges.endtimes <- currentedges.starttimes + currentedges.lengths
		edgetimes[currentedges.nums, 2] <- currentedges.endtimes

		# Check that each of the two endnodenums is actually a node, not a tip
		currentedges.endnodenums <- edges[currentedges.nums, 2]
		# For each of the two current edges' end nodes...
		for (i in 1:length(currentedges.endnodenums)) {
			# ... If it's a node (not a tip), set the start edgetimes of that node's two daughter edges to the end edgetime of...
			# ... the parent edge (i.e. the current edge).
			if (currentedges.endnodenums[i] > length(phy[["tip.label"]])) {
				edgetimes[edges[, 1] == currentedges.endnodenums[i], 1] <- currentedges.endtimes[i]
			}
		}
	}

	# "Reverse" the ages using total length calculated with a confirmed-extant tip.
	# Must use VerifyTreeLengths() on an unladderized verison of the tree to get total length via a confirmed extant tip...
	# ... because ladderize() removes tree$data, wherein lies extinct/extant data for each edge.
	if (isTRUE(all.equal(max(edgetimes[, 2]), totlength))) {
		edgetimes <- abs(edgetimes - max(edgetimes[, 2]))
		# Note: gives very small errors (e.g. 4e-16 instead of 0)
		# Almost certainly rounding errors from calculating edge times by adding-up preceding lengths; inconsequential
	} else {
		cat("Error: Total length of tree calculated using an extant tip != maximum death time of any edge, can't reverse edge times.")
		quit(save = "ask")
	}

	# Trim-off a final section of the tree by subtracting "trimlength" from all edge times.
	# This allows all trees, which cannot be grown to an exact length excluding the root edge,...
	# ... to have the same length (from the root node at 541 Ma up to 0 Ma).
	trimlength <- max(edgetimes) - marstart
	edgetimes.trim <- edgetimes - trimlength
	trimmededges.starts <- vector()
	trimmededges.ends <- vector()
	for (i in 1 :length(edgetimes.trim[, 1])) {
		if (edgetimes.trim[i, 1] < 0) {
			edgetimes.trim[i, 1] <- 0
			trimmededges.starts <- c(trimmededges.starts, i)   # Keep record of which edges' starts were trimmed
		}
	}
	for (i in 1 :length(edgetimes.trim[, 2])) {
		if (edgetimes.trim[i, 2] < 0) {
			edgetimes.trim[i, 2] <- 0
			trimmededges.ends <- c(trimmededges.ends, i)   # Keep record of which edges' ends were trimmed
		}
	}
	# The max value in edgetimes.trim (i.e birth time of root node) should now equal 541

	outputs <- list(edgetimes.trim, trimlength, trimmededges.starts, trimmededges.ends)
	return(outputs)

}


# As above, but specifically for use on "reconstructed" trees, re-scaled or otherwise
GetReconEdgeTimes <- function(phy, rootnodenum, rootnodeage) {

	edges <- phy[["edge"]]
	edgetimes <- matrix(data = NA, nrow = nrow(edges), ncol = ncol(edges))

	edgelengths <- phy[["edge.length"]]

	edgetimes[edges[, 1] == rootnodenum, 1] <- rootnodeage

	allnodenums <- rootnodenum:((rootnodenum - 1) + phy[["Nnode"]])
	for (currentnodenum in allnodenums) {
		currentedges.nums <- which(edges[, 1] == currentnodenum)
		currentedges.starttimes <- edgetimes[currentedges.nums, 1]
		currentedges.lengths <- edgelengths[currentedges.nums]
		currentedges.endtimes <- currentedges.starttimes - currentedges.lengths   # Minus because using already-reversed rootnodeage
		edgetimes[currentedges.nums, 2] <- currentedges.endtimes

		currentedges.endnodenums <- edges[currentedges.nums, 2]
		for (i in 1:length(currentedges.endnodenums)) {
			if (currentedges.endnodenums[i] > length(phy[["tip.label"]])) {
				edgetimes[edges[, 1] == currentedges.endnodenums[i], 1] <- currentedges.endtimes[i]
			}
		}
	}

	return(edgetimes)
	# Note: again, very small errors (e.g. 4e-16 instead of 0) due to inherited rounding errors; again, inconsequential

}


GetIntervalsSamplingRates <- function(importfilename, whichcol) {

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

	# In not running full-length simulations, set ratio of simulation time units:sampling proxy time units, dependent on user choice
	if (marstart == 54.1 & terstart == 35.89) {
		intervalsinfo[, 1] <- intervalsinfo[, 1] / 10
		intervalsinfo[, 2] <- intervalsinfo[, 2] / 10
	} else if (marstart == 5.41 & terstart == 3.589) {
		intervalsinfo[, 1] <- intervalsinfo[, 1] / 100
		intervalsinfo[, 2] <- intervalsinfo[, 2] / 100
	}

	return(intervalsinfo)

}


SampleEdges <- function(intsinfo, etimes) {

	# Create outer list
	multipleints.whensampled <- vector("list", nrow(intsinfo))

	# Iterating over each row (i.e. each sampling interval) in the intsinfo matrix
	for (i in 1:nrow(intsinfo)) {

		# Find all lineages present (at least partially) within the interval
		presence <- which((etimes[, 1] >= intsinfo[i, 2]) & (etimes[, 2] <= intsinfo[i, 1]))

		# Only sample if some lineage(s) is/are present
		if (length(presence) > 0) {
			# Record length of time each lineage was present within the interval for
			# Also record new matrix containing linintervals start and end times for each sampled (section of) lineage
			linintervals <- matrix(data = 0, nrow = nrow(etimes), ncol = 1)
			linint.etimes <- matrix(data = NA, nrow = nrow(etimes), ncol = 2)
			for (j in 1:length(presence)) {
				if (etimes[presence[j], 1] >= intsinfo[i, 1]) {
					linintervals[presence[j]] <- intsinfo[i, 1] - max(etimes[presence[j], 2], intsinfo[i, 2])
					linint.etimes[presence[j], 1] <- intsinfo[i, 1]
					linint.etimes[presence[j], 2] <- max(etimes[presence[j], 2], intsinfo[i, 2])
				} else if (etimes[presence[j], 1] < intsinfo[i, 1]) {
					linintervals[presence[j]] <- etimes[presence[j], 1] - max(etimes[presence[j], 2], intsinfo[i, 2])
					linint.etimes[presence[j], 1] <- etimes[presence[j], 1]
					linint.etimes[presence[j], 2] <- max(etimes[presence[j], 2], intsinfo[i, 2])
				}
			}

			# Use Poisson distribution to generate values for number of times "present-in-interval" lineages were sampled
			occslinsampled <- rpois(n = length(linintervals), lambda = linintervals * intsinfo[i, 3])

			# Use unifrom distribution to generate values for when each sample of each lineage occurred
			whensampled <- vector("list", length(occslinsampled))
			for (k in 1:length(occslinsampled)) {
				whensampled[[k]] <- runif(n = occslinsampled[k], max = linint.etimes[k, 1], min = linint.etimes[k, 2])
			}

			# Save this i's whensampled as new element in outer list
			multipleints.whensampled[[i]] <- whensampled
		}
	}

	# Reduce multipleints.whensampled into list of lineages containing only first and last samples
	maxminsampletimes <- vector("list", length(etimes[, 1]))
	for (j in 1:length(etimes[, 1])) {
		for (tbin in 1:length(intsinfo[, 1])) {
			maxminsampletimes[[j]] <- c(maxminsampletimes[[j]], multipleints.whensampled[[tbin]][[j]])
		}
	}
	for (j in 1:length(maxminsampletimes)) {
		if (length(maxminsampletimes[[j]]) == 0) {
			maxminsampletimes[[j]] <- NA
		}
	}

	firstlastsamptimes <- matrix(data = NA, nrow = length(maxminsampletimes), ncol = 2)
	for (j in 1:length(maxminsampletimes)) {
		if (!anyNA(maxminsampletimes[[j]])) {
			firstlastsamptimes[j, 1] <- max(maxminsampletimes[[j]], na.rm = TRUE)
			firstlastsamptimes[j, 2] <- min(maxminsampletimes[[j]], na.rm = TRUE)
		}
	}

	return(firstlastsamptimes)

}


BindToNodes <- function(kickedtree, flstimes, etimes, trmlength, trmedgeends) {

	tree_nodeloop <- ladderize(kickedtree, right = FALSE)
	tree_nodeloop <- makeNodeLabel(tree_nodeloop)
	original.numofnodes <- length(tree_nodeloop[["node.label"]])
	original.numoftips <- length(tree_nodeloop[["tip.label"]])
	original.edgeendnodenums <- tree_nodeloop[["edge"]][, 2]

	tree_nodeloop.original <- tree_nodeloop

	addedtipcounter <- 0

	# Rename each original tip, necessary for BindToTips
	for (i in 1:length(tree_nodeloop[["tip.label"]])) {
		tree_nodeloop[["tip.label"]][i] <- paste("Tip", i, sep = "")
	}

	# Bind new tips to each original node, excluding the root node
	for (currentnode in (original.numofnodes + original.numoftips):(original.numoftips + 2)) {
		currentnode.name <- paste("Node", currentnode - original.numoftips, sep = "")
		currentedge <- which(original.edgeendnodenums == currentnode)

		if (!is.na(flstimes[currentedge, 1] - flstimes[currentedge, 2]) &&
		    !is.na(flstimes[currentedge, 2] - etimes[currentedge, 2])) {

			# Check that the (untrimmed) positon argument > 0, if <= 0 then no new node will be created and the counter gets out of sync
			if (flstimes[currentedge, 1] - etimes[currentedge, 2] <= 0) {
				cat("Error: Position argument <= 0, counter out of sync.")
				quit(save = "ask")
			}

			# Add trimlength to position arguement if currentedge's tail portion was trimmed
			if (currentedge %in% trmedgeends) {
				# The below 'where' argument only works when binding to a node; it fails if it tries to bind to tips because they have no custom labels yet
				tree_nodeloop <- bind.tip(tree_nodeloop,
								  tip.label = paste("tipforedge", currentedge, sep = ""),
								  edge.length = flstimes[currentedge, 1] - flstimes[currentedge, 2],
								  where = addedtipcounter + original.numoftips + which(tree_nodeloop.original[["node.label"]] == currentnode.name),
								  position = trmlength + flstimes[currentedge, 1] - etimes[currentedge, 2])
				addedtipcounter <- addedtipcounter + 1
				# This is really counting the number of new nodes added (as long as position > 0, each added tip adds a new node too)

			} else {
				# The below 'where' argument only works when binding to a node; it fails if it tries to bind to tips because they have no custom labels yet
				tree_nodeloop <- bind.tip(tree_nodeloop,
								  tip.label = paste("tipforedge", currentedge, sep = ""),
								  edge.length = flstimes[currentedge, 1] - flstimes[currentedge, 2],
								  where = addedtipcounter + original.numoftips + which(tree_nodeloop.original[["node.label"]] == currentnode.name),
								  position = flstimes[currentedge, 1] - etimes[currentedge, 2])
				addedtipcounter <- addedtipcounter + 1
				# This is really counting the number of new nodes added (as long as position > 0, each added tip adds a new node too)

			}
		}
	}

	outputs <- list(tree_nodeloop, original.numoftips, original.edgeendnodenums)
	return(outputs)

}


BindToTips <- function(nlooptree, orignumoftips, origedgeendnodenums, flstimes, etimes, trmlength, trmedgeends) {

	tree_tiploop <- nlooptree

	# Bind new tips to each original tip
	for (currenttip in orignumoftips:1) {
		currenttipname <- paste("Tip", currenttip, sep = "")
		currentedge <- which(origedgeendnodenums == currenttip)

		if (!is.na(flstimes[currentedge, 1] - flstimes[currentedge, 2] ) == TRUE &&
			!is.na(flstimes[currentedge, 2] - etimes[currentedge, 2]) == TRUE) {
			if (!startsWith(tree_tiploop[["tip.label"]][which(tree_tiploop[["tip.label"]] == currenttipname)], "tipforedge")) {

				# Add trimlength to position arguement if currentedge's tail portion was trimmed
				if (currentedge %in% trmedgeends) {
					tree_tiploop <- bind.tip(tree_tiploop,
									 tip.label = paste("tipforedge", currentedge, sep = ""),
									 edge.length = flstimes[currentedge, 1] - flstimes[currentedge, 2],
									 where = which(tree_tiploop[["tip.label"]] == currenttipname),
									 position = trmlength + flstimes[currentedge, 1] - etimes[currentedge, 2])
				} else {
					tree_tiploop <- bind.tip(tree_tiploop,
									 tip.label = paste("tipforedge", currentedge, sep = ""),
									 edge.length = flstimes[currentedge, 1] - flstimes[currentedge, 2],
									 where = which(tree_tiploop[["tip.label"]] == currenttipname),
									 position = flstimes[currentedge, 1] - etimes[currentedge, 2])
				}
			}
		}
	}

	return(tree_tiploop)

}


DropOldTips <- function(tlooptree, orignumoftips) {

	tree_droploop <- tlooptree

	# Remove all old/original tips
	for (currentoldtip in 1:orignumoftips) {
		tree_droploop <- drop.tip(tree_droploop, tip = which(tree_droploop[["tip.label"]] == paste("Tip", currentoldtip, sep = "")))
	}
	tree_droploop <- ladderize(tree_droploop, right = FALSE)

	return(tree_droploop)

}


# Produces diversity through time curve for the true tree
GetTrueDiversityCurve <- function(trueetimes) {

	# Counts number of lineages just before and just after each entry in etimes
	linsthrutime <- matrix(data = NA, nrow = 2 * length(trueetimes), ncol = 2)
	for (i in 1:length(trueetimes)) {
		pre_time <- trueetimes[i] + 0.001
		post_time <- trueetimes[i] - 0.001

		linsthrutime[2 * i - 1, 1] <- sum(trueetimes[, 1] >= pre_time & trueetimes[, 2] <= pre_time, na.rm = TRUE)
		linsthrutime[2 * i - 1, 2] <- pre_time

		linsthrutime[2 * i, 1] <- sum(trueetimes[, 1] >= post_time & trueetimes[, 2] <= post_time, na.rm = TRUE)
		linsthrutime[2 * i, 2] <- post_time
	}

	# Re-order linsthrutime (by time/column two) into chronological order
	linsthrutime <- linsthrutime[order(linsthrutime[, 2], decreasing = TRUE), ]

	# Remove the two pre_time entries of the root node (sum will == 0 when it should == 1, because root edge is not in etimes)
	linsthrutime <- linsthrutime[c(-1, -2), ]

	# Remove any entries with negative time values (they are post_times of extant lineages, i.e. times after simulation stopped)
	if (sum(linsthrutime[, 2] < 0) > 0) {
		linsthrutime <- linsthrutime[-1 * which(linsthrutime[, 2] < 0), ]
	}
	
	# Append lineage count of 2 at the exact true root's time.
	# Then append lineage count equal to the final lineage count at exactly time 0.
	# This ensures the all the ltt curves start and end at the same times.
	truerootage <- max(trueetimes)
	if (max(linsthrutime[, 2]) < truerootage) {
		unsampledtime.initial <- c(2, truerootage)
		linsthrutime <- rbind(unsampledtime.initial, linsthrutime, deparse.level = 0)
	}
	if (min(linsthrutime[, 2]) > 0) {
		append <- c(linsthrutime[which(linsthrutime[, 2] == min(linsthrutime[, 2])), 1][1], 0)
		linsthrutime <- rbind(append, linsthrutime, deparse.level = 0)
	}

	# Re-order linsthrutime again
	linsthrutime <- linsthrutime[order(linsthrutime[, 2], decreasing = TRUE), ]

	return(linsthrutime)

}


# Produces diversity through time curve for a reconstructed tree, re-scaled or otherwise
GetReconstructionDiversityCurve <- function(truerootage, trueetimes, truetree, thisrecontree, origrecontree) {

	# Identify the oldest "NodeX" in thisrecontree, i.e. the oldest original node remaining
	remainingorig.nodenames <- thisrecontree[["node.label"]][startsWith(origrecontree[["node.label"]], "Node")]
	currentmaxnodeage <- 0
	currentnodename <- remainingorig.nodenames[1]
	truetree <- makeNodeLabel(truetree)
	for (i in 1:length(remainingorig.nodenames)) {
		currentnodeage <- trueetimes[min(which(truetree[["edge"]][, 1] == (length(truetree[["tip.label"]]) +
									 which(truetree[["node.label"]] == remainingorig.nodenames[i])))), 1]
		if (currentnodeage > currentmaxnodeage) {
			currentmaxnodeage <- currentnodeage
			currentnodename <- remainingorig.nodenames[i]
		}
	}

	# Get the node number of thisrecontree's oldest surviving "NodeX"
	recrootnodenum <- length(thisrecontree[["tip.label"]]) + which(thisrecontree[["node.label"]] == currentnodename)

	# Check whether any bind.tip()-created "NA" nodes are in fact older than the oldest surviving "NodeX"
	if (length(which(thisrecontree[["edge"]][, 2] == recrootnodenum)) > 0 & recrootnodenum != length(thisrecontree[["tip.label"]]) + 1) {
		recrootnodenum.new <- length(thisrecontree[["tip.label"]]) + 1
		sumaddlengths <- 0
		for (i in recrootnodenum:(recrootnodenum.new + 1)) {
			addedge <- which(thisrecontree[["edge"]][, 2] == i)
			addlength <- thisrecontree[["edge.length"]][addedge]
			sumaddlengths <- sumaddlengths + addlength
		}
		recrootnodeage.new <- currentmaxnodeage + sumaddlengths

		# Get start and end times for all edges, including the internal non-sampled branches
		reconetimes <- GetReconEdgeTimes(phy = thisrecontree, rootnodenum = recrootnodenum.new, rootnodeage = recrootnodeage.new)

	} else {
		# Get start and end times for all edges, including the internal non-sampled branches
		reconetimes <- GetReconEdgeTimes(phy = thisrecontree, rootnodenum = recrootnodenum, rootnodeage = currentmaxnodeage)

	}

	# Counts number of lineages just before and just after each entry in etimes
	linsthrutime <- matrix(data = NA, nrow = 2 * length(reconetimes), ncol = 2)
	for (i in 1:length(reconetimes)) {
		pre_time <- reconetimes[i] + 0.001
		post_time <- reconetimes[i] - 0.001

		linsthrutime[2 * i - 1, 1] <- sum(reconetimes[, 1] >= pre_time & reconetimes[, 2] <= pre_time, na.rm = TRUE)
		linsthrutime[2 * i - 1, 2] <- pre_time

		linsthrutime[2 * i, 1] <- sum(reconetimes[, 1] >= post_time & reconetimes[, 2] <= post_time, na.rm = TRUE)
		linsthrutime[2 * i, 2] <- post_time
	}

	# Re-order linsthrutime (by time/column two) into chronological order
	linsthrutime <- linsthrutime[order(linsthrutime[, 2], decreasing = TRUE), ]

	# Remove the any entries with counts of zero (including pre_time entries of the root node)
	linsthrutime <- linsthrutime[-1 * which(linsthrutime[, 1] == 0), ]

	# Remove any entries with negative time values (they are post_times of extant lineages, i.e. times after simulation stopped).
	# This is very unlikely for a reconstructed tree.
	if (sum(linsthrutime[, 2] < 0) > 0) {
		linsthrutime <- linsthrutime[-1 * which(linsthrutime[, 2] < 0), ]
	}

	# Add exact-time starting entry with a lineage count equal to the previous first count, then re-order
	unsampledtime.atreconroot <- c(linsthrutime[which(linsthrutime[, 2] == max(linsthrutime[, 2])), 1][1], linsthrutime[which(linsthrutime[, 2] == max(linsthrutime[, 2])), 2][1] + 0.001)
	linsthrutime <- rbind(unsampledtime.atreconroot, linsthrutime, deparse.level = 0)
	linsthrutime <- linsthrutime[order(linsthrutime[, 2], decreasing = TRUE), ]

	# If recontree's root node's age is < true root age, add count of 0 lineages at true root age.
	# Also add a count of zero at recontree's root node's pre_time.
	if (max(linsthrutime[, 2]) < truerootage) {
		unsampledtime.prereconroot <- c(0, linsthrutime[1, 2] + 0.001)
		linsthrutime <- rbind(unsampledtime.prereconroot, linsthrutime, deparse.level = 0)

		unsampledtime.attrueroot <- c(0, truerootage)
		linsthrutime <- rbind(unsampledtime.attrueroot, linsthrutime, deparse.level = 0)
		linsthrutime <- linsthrutime[order(linsthrutime[, 2], decreasing = TRUE), ]
	}

	# If recontree's final sample was at time > 0, add final entry at time 0 with a lineage count of 0.
	# This should be the case for almost every reconstructed tree.
	# Also add a count of zero at recontree's final edge's exact time.
	if (min(linsthrutime[, 2]) > 0) {
		unsampledtime.atfinaltip <- c(0, min(linsthrutime[, 2]) - 0.001)
		linsthrutime <- rbind(unsampledtime.atfinaltip, linsthrutime, deparse.level = 0)

		unsampledtime.attimezero <- c(0, 0)
		linsthrutime <- rbind(unsampledtime.attimezero, linsthrutime, deparse.level = 0)
		linsthrutime <- linsthrutime[order(linsthrutime[, 2], decreasing = TRUE), ]
	}

	# Re-order linsthrutime again
	linsthrutime <- linsthrutime[order(linsthrutime[, 2], decreasing = TRUE), ]

	outputs <- list(linsthrutime, reconetimes)
	return(outputs)

}


RescaleReconstructedTree <- function(origrecontree, flstimes) {

	# Generate timedata from the first and last sample times, setting each row's name to the corresponding "tipforedgeX"
	timedata <- NULL
	for (i in 1:length(which(!is.na(flstimes[, 1])))) {
		thisrow <- matrix(data = NA, nrow = 1, ncol = 2)
		thisrow[1, 1] <- flstimes[, 1][which(!is.na(flstimes[, 1]))][i]
		thisrow[1, 2] <- flstimes[, 2][which(!is.na(flstimes[, 2]))][i]
		rownames(thisrow) <- paste("tipforedge", which(flstimes[, 1] == thisrow[1, 1] & flstimes[, 2] == thisrow[1, 2]), sep = "")

		timedata <- rbind(timedata, thisrow)
	}

	# Produce a re-scaled reconstructed tree
	rescaledrecontree <- timePaleoPhy(tree = origrecontree, timeData = timedata, type = "basic", dateTreatment = "firstLast", add.term = TRUE)

	return(rescaledrecontree)

}


DoLinearRegression <- function(xdata, ydata) {

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
	for (i in 1:length(xdata)) {
		eachxy[i] <- xdata[i] * ydata[i]
	}
	xymean <- sum(eachxy) / length(eachxy)

	# Calculate mean of x^2
	eachx2 <- matrix(data = NA, ncol = 1, nrow = length(xdata))
	for (i in 1:length(xdata)) {
		eachx2[i] <- xdata[i]^2
	}
	x2mean <- sum(eachx2) / length(xdata)

	# Calculate slope of linear regression line
	m <- (xymean - xmean * ymean) / (x2mean - xmean^2)

	return(m)

}


CalculateSlopeDifference <- function(kickedtree, etimes, origrectree, flstimes, sampmodel) {

	# Get age of true tree's root node
	truertage <- max(etimes)   # Should be 541 due to "trimming"

	# Get ltt for the true tree and the non-re-scaled reconstructed tree
	trueltt <- GetTrueDiversityCurve(trueetimes = etimes)
	outputs <- GetReconstructionDiversityCurve(truerootage = truertage,
								 trueetimes = etimes,
								 truetree = kickedtree,
								 thisrecontree = origrectree,
								 origrecontree = origrectree)
	reconltt <- outputs[[1]]
	reconetimes <- outputs[[2]]

	# Plot reconstructed tree (i.e. tree_droploop), had to wait until here as reconetimes is needed for x-axis
	if (createplots == TRUE) {
		dev.new()
		plot(origrectree, show.tip.label = FALSE)
		reconrootage <- max(reconetimes)   # Can be no older than 541 Ma
		reconfinalage <- min(reconetimes)   # Should equal the youngest sampling event's time

		# If terrestrial sampling data was used, set correct start and end times for x-axis
		if (sampmodel %in% c("t_occ", "t_fm", "t_loc")) {
			reconrootage <- reconrootage - (marstart - terstart)
			reconfinalage <- reconfinalage - (marstart - terstart)
		}
		axis(side = 1, at = seq(reconfinalage, reconrootage, length.out = 50), labels = round(seq(reconrootage, reconfinalage, length.out = 50), 2))
		title("tree_droploop")
	}

	# Rescale reconstructed tree
	resrectree <- RescaleReconstructedTree(origrecontree = origrectree, flstimes = flstimes)

	# Get ltt for the re-scaled reconstructed tree
	outputs <- GetReconstructionDiversityCurve(truerootage = truertage,
								 trueetimes = etimes,
								 truetree = kickedtree,
								 thisrecontree = resrectree,
								 origrecontree = origrectree)
	resrecltt <- outputs[[1]]
	resrecetimes <- outputs[[2]]   # Unnecessary, but useful for debugging

	# Save all the lineage-through-time curves to allow the true curve to be plotted against all 6 reconstructed and all 6 rescaled-reconstructed curves
	thisltts <- list(trueltt, reconltt, resrecltt)
	
	# Plot all diversity curves on same graph
	if (createplots == TRUE) {
		dev.new()

		if (sampmodel %in% c("t_occ", "t_fm", "t_loc")) {
			plot(x = trueltt[, 2] - (marstart - terstart),
			     y = trueltt[, 1],
			     type = "l",
			     col = "black",
			     xlab = "time before end of simulation",
			     ylab = "number of lineages present",
			     xaxs = "i",
			     yaxs = "i",
			     xlim = c(truertage - (marstart - terstart), 0),
			     ylim = c(0, 2 + max(c(trueltt[, 1], reconltt[, 1], resrecltt[, 1]))))

			lines(x = reconltt[, 2] - (marstart - terstart), y = reconltt[, 1], col = "green")
			lines(x = resrecltt[, 2] - (marstart - terstart), y = resrecltt[, 1], col = "red")

		} else {
			plot(x = trueltt[, 2],
			     y = trueltt[, 1],
			     type = "l",
			     col = "black",
			     xlab = "time before end of simulation",
			     ylab = "number of lineages present",
			     xaxs = "i",
			     yaxs = "i",
			     xlim = c(truertage, 0),
			     ylim = c(0, 2 + max(c(trueltt[, 1], reconltt[, 1], resrecltt[, 1]))))

			lines(x = reconltt[, 2], y = reconltt[, 1], col = "green")
			lines(x = resrecltt[, 2], y = resrecltt[, 1], col = "red")

		}

		title(main = "true, reconstructed, and re-scaled reconstructed ltts")
		legend(x = "topleft", c("true tree", "reconstructed tree", "re-scaled reconstructed tree"),
			 lty = c(1, 1, 1),
			 lwd = c(1, 1, 1),
			 col = c("black", "green", "red"))
	}

	# Cut-off start and end portions of the ltts to remove any edge effect and capture only the in-equilibrium portion of the true tree.
	# Note: this will break if using the short or medium length simulation options, as those simulations (at least the terrestrial parts)...
	# ... are shorter than the lengths cut-off.
	if (sampmodel %in% c("t_occ", "t_fm", "t_loc")) {
		# For terrestrial sampling data iterations:
		# Cut-off 20 Myr start portion and 20 Myr end portion (i.e. removes all ltt counts before 338.9 Ma and after 20 Ma, in "terrestrial time").
		# Remember, terrestrial iterations have their start times scaled to 541 Ma, so must use "marine time" to determine when to cut-off.
		trueltt <- trueltt[-1 * which(trueltt[, 2] > (marstart - 20)), ]
		trueltt <- trueltt[-1 * which(trueltt[, 2] < ((marstart - terstart) + 20)), ]

		reconltt <- reconltt[-1 * which(reconltt[, 2] > (marstart - 20)), ]
		reconltt <- reconltt[-1 * which(reconltt[, 2] < ((marstart - terstart) + 20)), ]

		resrecltt <- resrecltt[-1 * which(resrecltt[, 2] > (marstart - 20)), ]
		resrecltt <- resrecltt[-1 * which(resrecltt[, 2] < ((marstart - terstart) + 20)), ]

	} else {
		# For marine sampling data iterations:
		# Cut-off 20 Myr start portion (before 521 Ma) and end portion (after 20 Ma) of all ltts
		trueltt <- trueltt[-1 * which(trueltt[, 2] > (marstart - 20)), ]
		trueltt <- trueltt[-1 * which(trueltt[, 2] < 20), ]

		reconltt <- reconltt[-1 * which(reconltt[, 2] > (marstart - 20)), ]
		reconltt <- reconltt[-1 * which(reconltt[, 2] < 20), ]

		resrecltt <- resrecltt[-1 * which(resrecltt[, 2] > (marstart - 20)), ]
		resrecltt <- resrecltt[-1 * which(resrecltt[, 2] < 20), ]

	}

	# Save absolute slope difference ----------------------------------------------------------------------------------------
	trueslope <- -1 * DoLinearRegression(xdata = trueltt[, 2], ydata = trueltt[, 1])
	reconslope <- -1 * DoLinearRegression(xdata = reconltt[, 2], ydata = reconltt[, 1])
	resrecslope <- -1 * DoLinearRegression(xdata = resrecltt[, 2], ydata = resrecltt[, 1])

	slopedifference.abs <- resrecslope - trueslope

	# Plot all linear regression slopes
	if (createplots == TRUE) {
		dev.new()
		plot(x = 1,
		     y = trueslope,
		     pch = 4,
		     col = "black",
		     xlab = "tree type (1 = true, 2 = recon and resrec)",
		     ylab = "linear regression gradient",
		     xaxs = "i",
		     yaxs = "i",
		     xlim = c(0, 3),
		     ylim = c(min(c(trueslope, reconslope, resrecslope)) - 0.01, max(c(trueslope, reconslope, resrecslope)) + 0.01))

		points(x = 2, y = reconslope, pch = 4, col = "green")
		points(x = 2, y = resrecslope, pch = 4, col = "red")

		title(main = "all trees' linear regression gradients (absolute)")
		legend(x = "topleft", c("true tree", "reconstructed tree", "re-scaled reconstructed tree"), pch = c(4, 4, 4), col = c("black", "green", "red"))
	}

	# Save relative slope difference ----------------------------------------------------------------------------------------
	# Re-scale all ltts to end on same value, e.g. 3, so as to calculate relative slope difference (not absolute).
	x <- trueltt[length(trueltt[, 1]), 1] / 3
	trueltt[, 1] <- trueltt[, 1] / x

	x <- reconltt[length(reconltt[, 1]), 1] / 3
	reconltt[, 1] <- reconltt[, 1] / x

	x <- resrecltt[length(resrecltt[, 1]), 1] / 3
	resrecltt[, 1] <- resrecltt[, 1] / x

	# Do linear regression for every ltt curve.
	# Note: geologcial time counts down towards zero, so seemingly positive slopes will actually produce negative linear regression gradients.
	# However, for ease of interpretaion, I multiply all linear regression gradients by -1.
	# i.e. if diversity is increasing as time moves towards zero, the slope will be reported as positive.
	trueslope <- -1 * DoLinearRegression(xdata = trueltt[, 2], ydata = trueltt[, 1])
	reconslope <- -1 * DoLinearRegression(xdata = reconltt[, 2], ydata = reconltt[, 1])
	resrecslope <- -1 * DoLinearRegression(xdata = resrecltt[, 2], ydata = resrecltt[, 1])

	# Plot all linear regression slopes
	if (createplots == TRUE) {
		dev.new()
		plot(x = 1,
		     y = trueslope,
		     pch = 4,
		     col = "black",
		     xlab = "tree type (1 = true, 2 = recon and resrec)",
		     ylab = "linear regression gradient",
		     xaxs = "i",
		     yaxs = "i",
		     xlim = c(0, 3),
		     ylim = c(min(c(trueslope, reconslope, resrecslope)) - 0.01, max(c(trueslope, reconslope, resrecslope)) + 0.01))

		points(x = 2, y = reconslope, pch = 4, col = "green")
		points(x = 2, y = resrecslope, pch = 4, col = "red")

		title(main = "all trees' linear regression gradients (relative)")
		legend(x = "topleft", c("true tree", "reconstructed tree", "re-scaled reconstructed tree"), pch = c(4, 4, 4), col = c("black", "green", "red"))
	}

	slopedifference.rel <- resrecslope - trueslope

	outputs <- list(slopedifference.abs, slopedifference.rel, thisltts, resrectree)
	return(outputs)

}


# Runs simulations for one tree, sampled using each sampling model
RunSimulations <- function(divmodel, allsamplemodels) {

	graphics.off()

	outputs <- GrowTree(diversificationmodel = divmodel)
	tree_orig <- outputs[[1]]
	tree_kicked <- outputs[[2]]

	ttsave <- tree_kicked

	outputs <- VerifyTreeLengths(tree = tree_orig)
	totallength <- outputs[[1]]
	rootlength <- outputs[[2]][[1]]

	outputs <- GetEdgeTimes(phy = tree_kicked, totlength = totallength)
	edgetimes <- outputs[[1]]   # Remember, these are "trimmed" edgetimes
	trimlength <- outputs[[2]]
	trimmededges.starts <- outputs[[3]]   # Unnecessary because if an edge's start has been trimmed, its end will have too and it will thus never be sampled
	trimmededges.ends <- outputs[[4]]   

	treelength.tr <- totallength - rootlength - trimlength   # Should == 541

	allsamplemodels.slopediffs.abs <- matrix(data = NA, ncol = 1, nrow = length(allsamplemodels))
	allsamplemodels.slopediffs.rel <- matrix(data = NA, ncol = 1, nrow = length(allsamplemodels))

	tdlsaves <- list(NA, NA, NA, NA, NA, NA)
	lttsaves <- list(NA, NA, NA, NA, NA, NA)
	resrectreesaves <- list(NA, NA, NA, NA, NA, NA)

	# Cycle through all 3 data columns in both file names
	progbar = txtProgressBar(min = 0, max = length(allsamplemodels), initial = 0)   # Initialise progress bar
	for (i in 1:length(allsamplemodels)) {
		current.samplemodel <- allsamplemodels[i]
		if (current.samplemodel == "t_occ"){ targetfile <- "OUTPUT_terrestrial_tetrapods.csv"; targetcol <- 4
		} else if (current.samplemodel == "t_fm"){ targetfile <- "OUTPUT_terrestrial_tetrapods.csv"; targetcol <- 5
		} else if (current.samplemodel == "t_loc") { targetfile <- "OUTPUT_terrestrial_tetrapods.csv"; targetcol <- 6
		} else if (current.samplemodel == "m_occ") { targetfile <- "OUTPUT_marine_eumetazoa.csv"; targetcol <- 4
		} else if (current.samplemodel == "m_fm") { targetfile <- "OUTPUT_marine_eumetazoa.csv"; targetcol <- 5
		} else if (current.samplemodel == "m_loc") { targetfile <- "OUTPUT_marine_eumetazoa.csv"; targetcol <- 6
		}

		intervalsinfo <- GetIntervalsSamplingRates(importfilename = targetfile, whichcol = targetcol)

		# The terrestrial tetrapod sampling data only reach back to 358.9 Ma.
		# The difference between this and the marine data's 541 Ma is 182.1 Ma.
		# So, to ensure terrestrial sampling begins at the true tree's root node's time (541 Ma)...
		# ... I must add 182.1 Myrs to all terrestrial time bins' start and end times.
		if (current.samplemodel %in% c("t_occ", "t_fm", "t_loc")) {
			intervalsinfo[, 1] <- intervalsinfo[, 1] + (marstart - terstart)
			intervalsinfo[, 2] <- intervalsinfo[, 2] + (marstart - terstart)
		}
		# Note: no sampling of the true tree will occur after 182.1117 Ma in the iterations using terrestrial sampling data...
		# ... because the last terrestrial data time bin ends at 0.0117 Ma.

		firstlastsampletimes <- SampleEdges(intsinfo = intervalsinfo, etimes = edgetimes)

		outputs <- BindToNodes(kickedtree = tree_kicked,
					     flstimes = firstlastsampletimes,
					     etimes = edgetimes,
					     trmedgeends = trimmededges.ends,
					     trmlength = trimlength)
		tree_nodeloop <- outputs[[1]]
		original.numoftips <- outputs[[2]]
		original.edgeendnodenums <- outputs[[3]]

		tree_tiploop <- BindToTips(nlooptree = tree_nodeloop,
					   orignumoftips = original.numoftips,
					   origedgeendnodenums = original.edgeendnodenums,
					   flstimes = firstlastsampletimes,
					   etimes = edgetimes,
					   trmedgeends = trimmededges.ends,
					   trmlength = trimlength)

		tree_droploop <- DropOldTips(tlooptree = tree_tiploop, orignumoftips = original.numoftips)
		tdlsaves[[i]] <- tree_droploop

		outputs <- CalculateSlopeDifference(kickedtree = tree_kicked,
								etimes = edgetimes,
								origrectree = tree_droploop,
								flstimes = firstlastsampletimes,
								sampmodel = current.samplemodel)
		allsamplemodels.slopediffs.abs[i] <- outputs[[1]]
		allsamplemodels.slopediffs.rel[i] <- outputs[[2]]
		lttsaves[[i]] <- outputs[[3]]
		resrectreesaves[[i]] <- outputs[[4]]

		setTxtProgressBar(progbar, i)   # Update progress bar
	}

	outputs <- list(allsamplemodels.slopediffs.abs, allsamplemodels.slopediffs.rel, ttsave, tdlsaves, lttsaves, resrectreesaves)
	return(outputs)

}
