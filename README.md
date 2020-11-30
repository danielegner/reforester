# Reforester

Reforester, an R program that:
1. Uses a logistic birth-death model to generate simulated phylogenetic trees, referred to as "true" trees.
2. Samples the true trees according to various empirically-derived fossil record sampling rate proxies.
3. Uses the phylogenetic diversity estimation (PDE) method to produce reconstructed trees based on the sampled branch sections.
4. Computes species diversity through time for both the true and reconstructed trees, also computing the difference in gradient between these diversification histories.



## Contents

  - [Research Summary](#research-summary)
  - [Getting Started](#getting-started)
  - [Authors](#authors)
  - [License](#license)
  - [Acknowledgements](#acknowledgements)



## Research Summary

### Introduction and aims
Did life diversify exponentially or logistically from the "Cambrian explosion" to the Recent? To answer this we must count the number of species present on Earth through geological time, as preserved in the fossil record. However, the incompleteness of the fossil record is heterogeneous and biased through geological time due to factors such as the inverse relationship between the age of a sedimentary deposit and its probability of having survived to the present (thereby having its fossil assemblages recorded). Phylogenetic diversity estimates (PDEs) attempt to partially correct for unpreserved sections of lineages' histories by combining inferred topological relationships with first and last appearance times from the fossil record. PDE studies are common and are often thought to be reliable estimates of past diversity.

The aim of this program is to demonstrate whether it is possible for inherently biased fossil sampling rate patterns to introduce bias into the phylogenies and diversity curves that are reconstructed using said fossil data, i.e. as in the PDE method. In other words, *is it theoretically possible that the increase in fossil sampling rate towards the Recent could make a logistic true diversity history appear to have been exponential when viewed through the lens of a phylogeny reconstructed using the fossil record/the PDE method?*

### Methods
I repeated the simulation 50 times, reconstructing each true tree using six time-heterogeneous sampling rate models (derived from global terrestrial tetrapod and marine eumetazoa fossil occurrence data). I compared these results with those produced by 50 simulation runs using a sampling rate that remains constant through time. In this way, I isolated the effects of biases in the sampling rate model.

### Results, conclusions, and implications
My dissertation, results, and statistical analyses are not published in this repository. In summary, they indicate that it *is* possible for such a bias to be introduced.

There are imperfections in my method, most-notably that my simulated true phylogenies are generated using birth-death parameters that give an equilibrium diversity of 30 species. The total number of extant terrestrial tetrapods is closer to 30,000, yet computational limits prevented me from running these simulations with such vast numbers of lineages over the full 541 Myr simulation time. Similarly, computational restrictions prevented me from generating exponentially-diversifying true trees, and so I was only able to demonstrate the effect on logistic true trees.

Nonetheless, if Phanerozoic life diversified logistically, my results show it is at least possible that such a diversification history could be misinterpreted as exponential due simply to biased fossil record sampling rate patterns through time. Due to the ubiquity of PDEs and the long-held assumption that they sufficiently account for fossil record incompleteness, my conclusions have significant implications for many past PDE studies. Further research is required to better understand this systematic bias in the PDE methodology, and perhaps subsequently to reinterpret post-Cambrian global diversity patterns.



## Getting Started

This project was *not* intended to be easily re-purposed or modified to use different input data, particularly because the simulation relies on durations and time bin start and end points that are hard-coded into many different functions. However, it is certainly possible to re-use or re-purpose the scripts, with some modification. I adhered to the [Google R style guide](https://google.github.io/styleguide/Rguide.html) throughout, and followed the [tidyverse R style guide](https://style.tidyverse.org) when further clarification was needed. This should make understanding the scripts easier. The meat and bones of the project, that is *main.R* and *functions.R*, will hopefully be the most useful parts of the project for re-use.

If simply testing the project and not modifying it, skip to [Run main.R](#run-main.R).

### Prerequisites
- R Windows x64 v3.5.2 (not tested in later releases), with working directory set as Reforester's base directory
- R packages:
    - caper v1.0.1
    - phytools v0.7-20
    - paleotree v3.3.25

### (Optional) Modify and run dataimporter.R
*dataimporter.R* takes as input two .csv files ("terrestrial_tetrapods.csv" and "marine_eumetazoa.csv") containing fossil occurrence data from The Paleobiology Database. These must be placed in the base directory along with *dataimporter.R*. The script extracts relevant data from the databases and collates this into the six sampling rate proxy time series.

Unfortunately I cannot include these .csv files in this repository. Editing *dataimporter.R* so that it is compatible with your own database files is likely more trouble than it's worth because my code is specific to the exact database format used in my input files. Instead of editing *dataimporter.R*, I recommend extracting your own time-binned sampling rate proxies separately and then feeding them into my simulation in place of *OUTPUT_terrestrial_tetrapods.csv* and *OUTPUT_marine_eumetazoa.csv*, the outputs of *dataimporter.R*. These output files are included and will give you an idea of the format that your time-binned sampling proxies will need to adhere to.

### Run main.R
If left unmodified, the simulations will use *OUTPUT_terrestrial_tetrapods.csv* and *OUTPUT_marine_eumetazoa.csv* as the sampling rate time proxies. These are what I used in my simulations.

To begin the simulation(s), **run main.R**.

In the console, answer the prompts as they appear. These will set the simulated duration and diversification model of each tree, as well as the total number of simulation runs. To generate my final results, I used the following options:

    - Choose "full (541 / 358.9 Myrs)" simulation length.
    - Choose "diversity-dependent" diversification model.
    - Set desired number of simulation runs (please note I could only complete a maximum of ~5 per working day).

### Run main_uniformsampling.R
To run simulations that use time-constant sampling rates (rather than the proxy-derived heterogeneous sampling models), **run main_uniformsampling.R** (located in the *uniformsampling* folder).

Answer the console prompts as above.

### Additional scripts
For completeness I have included all the additional scripts used to analyse my data and produce figures for the dissertation. Many sections of these additional scripts are copied across from the original core scripts with slight modifications.



## Authors

  - **Daniel Egner** - [GitHub](https://github.com/danielegner)



## License

Reforester is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Reforester is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Reforester (see the [COPYING.txt](COPYING.txt) file for details). If not, see <https://www.gnu.org/licenses/>.



## Acknowledgements

  - **R. Benson** - Project supervisor and originator of the project's aim
  - **R. Close** - Provided filtered fossil datasets, screening-out non-marine eumetazoa and marine tetrapods from the marine eumetazoa and terrestrial tetrapod databases respectively
  - **The Paleobiology Database** - The original source of all the fossil data used - [paleobiodb.org](https://paleobiodb.org)
  - **Creators and maintainers of R packages *caper*, *phytools*, and *paleotree*** - Integral to the operation of this program
  - **B. Thompson (PurpleBooth)** - Provided README template - [GitHub](https://github.com/PurpleBooth)
