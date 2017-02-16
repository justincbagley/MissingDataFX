# MissingData
Testing effects of missing data on phylogenetic inferences

## LICENSE

All code within the MissingData v0.1.0 repository is available "AS IS" under a generous GNU license. See the [LICENSE](LICENSE) file for more information.

## CITATION

If you use scripts from this repository as part of your published research, I require that you cite the repository as follows (also see DOI information below): 
  
- Bagley, J.C. 2017. MissingData v0.1.0. GitHub repository, Available at: http://github.com/justincbagley/MissingData.

Alternatively, please provide the following link to this software repository in your manuscript:

- https://github.com/justincbagley/MissingData

## DOI

The DOI for MissingData v0.1.0, via [Zenodo](https://zenodo.org), is coming soon... stay tuned.

## CONTENTS

- [Introduction](https://github.com/justincbagley/MissingData#introduction)
- [Getting Started](https://github.com/justincbagley/MissingData#getting-started)
- [Troubleshooting](https://github.com/justincbagley/MissingData#troubleshooting)
- [Acknowledgements](https://github.com/justincbagley/MissingData#acknowledgements)
- [References](https://github.com/justincbagley/MissingData#references)
- [TODO List](https://github.com/justincbagley/MissingData#todo-list)


## INTRODUCTION

> *"As fossil taxa will necessarily introduce a large proportion of missing data in any combined data set (Wiens 2009), the impact of these absent characters on branch length estimation and support may be substantial. In particular, the interaction between among partition rate variation and missing data may negatively affect phylogenetic inference (Lemmon et al. 2009). Additionally, heterogeneity in evolutionary rates between the morphological characters may also affect estimated branch lengths (Clarke and Middleton 2008)."* - Pyron (2011)

This repository automates exploratory analyses calculating the amount of missing data in phylogenetic datasets (NEXUS character partitions), as well as testing the form and potential significance of relationships between missing characters and phylogenetic tree parameters. Specifically, the current version of this software focuses on developing MissingData.sh, a shell script that 1) characterizes the contents of data blocks in a NEXUS input file, including proportions of data versus missing data for each partition, and 2) conducts a customized R analysis to a) extract tree parameters (terminal branch lengths, terminal/subtending branch heights, posterior support) from BEAST or MrBayes trees generated from the input NEXUS, and then b) use appropriate standard or nonparametric tests for correlations between missing data proportions and the tree parameters, and plot relationships among variables. More details are given below. This software was inspired by, and recreates, several correlational/exploratory analyses of the effects of missing data on phylogenies used previously by Wiens et al. (2005) and Pyron (2011).

As in the case of the author's software package for phylogenetic and phylogeographic inference, [PIrANHA](https://github.com/justincbagley/PIrANHA), the MissingData package is fully command line-based and is available as open-source software according to the license. 

## GETTING STARTED

### Dependencies

![alt tag](https://cloud.githubusercontent.com/assets/10584087/21243724/b12a94d2-c2df-11e6-9d20-5cf06877ad94.png)  

Code in the MissingData repository depends on R and a series of R packages, as follows:
- The [R](https://www.r-project.org) software environment
  * available from download mirrors from The Comprehensive R Archive Network (CRAN), such as https://cloud.r-project.org
- tools - a default R package with "base" priority
- [ape](https://cran.r-project.org/web/packages/ape/index.html)
- [ips](https://cran.r-project.org/web/packages/ips/index.html)
- [phytools](https://cran.r-project.org/web/packages/phytools/index.html)
- [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html)

### Installation
:computer: MissingData uses UNIX shell and R scripts and was developed on Mac OSX, thus runs on a variety of operating systems, but especially UNIX/LINUX-like systems. MissingData code should run "out-of-the-box" from most any folder on your machine. To 'install' MissingData, download the repository, move into the repository folder and enter ```$ chmod u+x ./*.sh``` into the command line interface (e.g. Terminal app on Mac). This changes file mode bits in the .sh files to allow the user permission to access and execute them, which is an important prerequisite for the MissingData code to run.

### Run directory structure
:warning: MissingData code assumes that you are running within a sub-folder specific to your analysis, located within the MissingData-master distro folder. Thus, we assume that the current working directory for any particular run contains **ONLY** 1) the MissingData.sh script (which you copy and paste into the dir), 2) the input NEXUS file, 3) the input BEAST or MrBayes tree file that you wish to analyze, and/or 4) an *optional* 'drop file' listing taxa to be pruned from the tree prior to analysis. **Importantly**, as a result of this directory structure, the "R" folder from the distro, which contains a modified function for ips, will be located at the relative path "../R/". See README and in-script commenting for further details on file types and instructions for running basic analyses. Ideally, users will run MissingData on multiple NEXUS-tree file combinations, for example 1) mtDNA only, 2) nuclear DNA only, 3) morphological characters only, and 4) all data combined--or a 'combined data' (mtDNA + nuclear) or 'total-evidence' (sequence + morphology) matrix and tree. Each of these analyses would be run within a separate sub-folder in the master distro folder.

### Input files
:warning: The main input files for MissingData are NEXUS files, tree files, and 'drop' files. For the NEXUS input file, we assume that there are no spaces or special characters in the filename (though underscores are OK), and that the file contains only DNA sequences or morphological data in simplified NEXUS format--i.e. header followed by a matrix block, and no subsequent data/info blocks (e.g. sets or MrBayes blocks; these must be removed). Prior to running MissingData, users should check their input NEXUS filenames and contents to ensure that they meet these assumptions.

Regarding tree files, let me say I'm very much a Bayesian, and my empirical research focuses frequently entails making model-based inferences in phylogenetics and phylogeography based on Bayesian or approximate Bayesian computation (e.g. ABC, hABC) methods (see my [Research](http://www.justinbagley.org/research), [Publications](http://www.justinbagley.org/publications), and [ResearchGate](https://www.researchgate.net/profile/Justin_Bagley2) pages). Thus, it's no surprise that the initial development offering of MissingData focuses **Only** on Bayesian phylogenetic trees. Allowed tree files include BEAST maximum clade credibility (MCC) trees annotated with posterior parameter distributions in TreeAnnotator (these files are output directly from TreeAnnotator and have '.out' extensions), and MrBayes consensus tree files (with '.con.tre' extensions, generated by running sumt).

'Drop' files are files the user creates in a text editor that contain a list of taxa (with names exactly matching taxon labels in the NEXUS and tree files), with one taxon name per line, and saved with '.drop' file extensions. Drop files specify taxa to be pruned from the phylogeny in R prior to further analyses. Drop files provide additional flexibility to the analysis, and they handy when the user would like to exclude certain taxa, for example a) extinct taxa or b) taxa with low or high amounts of missing data.

### What happens in R?
Following several steps summarizing the NEXUS input file using operations in the shell, MissingData.sh creates and runs a customized R script that loads the 'Interfaces to Phylogenetic Software in R' or 'ips' R package (plus my fixed version of one of its functions) and related phylogenetics packages, and then reads in the tree, plots the tree (and saves to file), extracts the node and branch length annotations into a matrix (also saved to file). Next, Shapiro-Wilk tests are conducted to evaluate whether the log-transformed data and tree parameters meet normality criteria for subsequent analyses. Then MissingData tests for correlations between the proportions of data or missing data and (y-axis/independent var.:) 1) posterior support for terminal taxa and 2) length of terminal branch (same as height estimate for terminal taxon tmrca/node, in case of BEAST chronograms). If the data are normal (Shapiro-Wilk tests were non-significant at the alpha = 0.05 level), regular Pearson correlations and generalized linear modeling analyses are conducted; however, if the data are non-normal (Shapiro-Wilk tests were significant), then correlations are conducted using nonparametric Spearman's rank correlation coefficient and linear relationships among log-transformed variables are plotted, but not given trendlines. Results are output to file in text or PDF (graphics) files. 

### Usage
It's so easy to use, the MissingData.sh script doesn't display any sophisticated Usage or help flag info yet, because it doesn't need to. Everything you need to know is given here in the README! Assuming that you have installed the dependencies and the repo and followed guidelines for input files above, you may run MissingData by simply entering the script name at the command line and hitting return!
```
$ ./MissingData.sh
```
#### Screen output example:
Here is an example of what output to screen would look like for a normal MissingData analysis including a NEXUS file, a drop file, and a MrBayes tree file:
```
$ ./MissingData.sh
##########################################################################################
#                           MissingData v0.1.0, February 2017                            #
##########################################################################################

INFO      | Thu Feb 16 12:25:37 CST 2017 | Starting MissingData analysis... 
INFO      | Thu Feb 16 12:25:37 CST 2017 | STEP #1: SETUP AND USER INPUT. 
INFO      | Thu Feb 16 12:25:38 CST 2017 |          Setting working directory to: /Users/justinbagley/Documents/GitHub/MissingData/Example_run 
INFO      | Thu Feb 16 12:25:38 CST 2017 |          Reading in input NEXUS file(s)... 
INFO      | Thu Feb 16 12:25:38 CST 2017 | STEP #2: PROCESSING INPUT NEXUS, SPLITTING TAXON LABELS AND DATA BLOCKS INTO SEPARATE FILES. 
INFO      | Thu Feb 16 12:25:38 CST 2017 |          4lociplusMorpho_n84_simple.NEX 
INFO      | Thu Feb 16 12:25:38 CST 2017 | STEP #3: USE TAXON LABEL AND CONCATENATED SEQUENCE FILES GENERATED DURING PREVIOUS STEP TO CREATE ONE 
INFO      | Thu Feb 16 12:25:38 CST 2017 |          FILE WITH MISSING DATA COUNTS AND PROPORTIONS FOR EACH INDIVIDUAL. 
INFO      | Thu Feb 16 12:25:38 CST 2017 | STEP #4: PRE-PROCESSING MRBAYES CONSENSUS TREE INPUT FILE, IF PRESENT: SPLIT .con.tre FILE, EXTRACT 
INFO      | Thu Feb 16 12:25:38 CST 2017 |          TERMINAL BRANCH LENGTHS & THEIR CONFIDENCE INTERVALS, AND CREATE BRANCH LENGTH SUMMARY TABLE. 
INFO      | Thu Feb 16 12:25:38 CST 2017 |          If MrBayes consensus tree file ('.con.tre' extension) present, split file, extract terminal 
INFO      | Thu Feb 16 12:25:38 CST 2017 |          branch lengths (term_brL), and tabulate branch length data including CIs. 
INFO      | Thu Feb 16 12:25:38 CST 2017 |          Found '.con.tre' file. Reading in MrBayes tree from current directory... 
INFO      | Thu Feb 16 12:25:38 CST 2017 |          Running mbTreeStatMiner script... 
INFO      | Thu Feb 16 12:25:38 CST 2017 | STEP #5: MAKE R SCRIPT THAT A) EXTRACTS PARAMETER ESTIMATES FROM BEAST TREES OR READS IN MRBAYES DATA
INFO      | Thu Feb 16 12:25:38 CST 2017 |          TABLE (STEP #3) IN WORKING DIR THEN B) TESTS FOR IMPACT OF MISSING DATA ON PHYLO SUPPORT AND BRANCH LENGTHS. 
INFO      | Thu Feb 16 12:25:38 CST 2017 | STEP #6: RUN THE R SCRIPT (WHICH ALSO SAVES RESULTS TO FILE). 
INFO      | Thu Feb 16 12:25:43 CST 2017 | STEP #7: CLEANUP: ORGANIZE RESULTS, REMOVE UNNECESSARY FILES. 
INFO      | Thu Feb 16 12:25:43 CST 2017 | Done analyzing the amount and potential effects of missing data on phylogenetic support and branch 
INFO      | Thu Feb 16 12:25:43 CST 2017 | lengths using MissingData. 
INFO      | Thu Feb 16 12:25:43 CST 2017 | Bye.

```

## TROUBLESHOOTING
How to troubleshoot some potentially common problems:

(**1**) One problem you may run into is the 'permission denied' error returned by the shell when attempting to execute the shell scripts; for example, bash might return "-bash: ./MissingData.sh: Permission denied". This indicates that permissions for .sh files were not correctly set during [Installation](https://github.com/justincbagley/MissingData#installation). Fix this by moving into the repository master folder and entering ```$ chmod u+x ./*.sh``` at the command line, then using the shell script(s) in the repo folder from that point forward.

## ACKNOWLEDGEMENTS
During the development of this software, J.C.B. received stipend support from a Ciência Sem Fronteiras (Science Without Borders) postdoctoral fellowship from the Brazilian Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq; Processo 314724/2014-1). Lab and computer space was also supplied by The University of Alabama, during an internship in the Lozier Lab in the UA Department of Biological Sciences.

## REFERENCES
- Bagley et al. (in revision) 
- Clarke JA, Middleton K (2008) Mosaicism, modules, and the evolution of birds: results from a Bayesian approach to the study of morphological evolution using discrete character data. Systematic Biology, 57, 185-201.
- Lemmon AR., Brown JM, Stanger-Hall K, Lemmon EM (2009) The effect of missing data on phylogenetic estimates obtained by maximum-likelihood and Bayesian inference. Systematic Biology, 58, 130-145.
- Pyron RA (2011) Divergence time estimation using fossils as terminal taxa and the origins of Lissamphibia. Systematic Biology, 60, 466-481.
- Wiens JJ, Fetzner JW, Parkinson CL, Reeder TW (2005) Hylid frog phylogeny and sampling strategies for speciose clades. Systematic Biology, 54, 719-748.

## TODO LIST
Current:

- Provide more graphical output options, including summaries of data/missing data proportions for different data blocks.
- Modify to accomodate and extract parameters from annotated maximum-likelihood trees from programs like RAxML.
- Add screenshots of R analyses, as well as example text and graphical output produced by the software.

Recently finished/fixed:

- Solve two input file problem. **DONE!** :white_check_mark:

February 16, 2017
Justin C. Bagley, Tuscaloosa, AL, USA
