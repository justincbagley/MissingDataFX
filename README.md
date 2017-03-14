# MissingDataFX
Testing effects of missing data on phylogenetic inferences

## LICENSE

All code within the MissingDataFX v0.1.0 repository is available "AS IS" under a generous GNU license. See the [LICENSE](LICENSE) file for more information.

## CITATION

If you use scripts from this repository as part of your published research, I require that you cite the repository as follows (also see DOI information below): 
  
- Bagley, J.C. 2017. MissingDataFX v0.1.0. GitHub repository, Available at: http://github.com/justincbagley/MissingDataFX.

Alternatively, please provide the following link to this software repository in your manuscript:

- https://github.com/justincbagley/MissingDataFX

## DOI

The DOI for MissingDataFX v0.1.0, via [Zenodo](https://zenodo.org), is coming soon... stay tuned.

## CONTENTS

- [Introduction](https://github.com/justincbagley/MissingDataFX#introduction)
- [Getting Started](https://github.com/justincbagley/MissingDataFX#getting-started)
- [Troubleshooting](https://github.com/justincbagley/MissingDataFX#troubleshooting)
- [Acknowledgements](https://github.com/justincbagley/MissingDataFX#acknowledgements)
- [References](https://github.com/justincbagley/MissingDataFX#references)
- [TODO List](https://github.com/justincbagley/MissingDataFX#todo-list)


## INTRODUCTION

> *"The long term goal of phylogenetics, both neontological and paleontological, is to reconstruct an accurate phylogeny for all species of living and fossil organisms. The problem of missing data has been considered to be the major obstacle to accurately reconstructing the phylogeny of fossil taxa and their relationships to living taxa. Recent simulation studies show that there is not a single missing data problem... there are potentially two problems... Adding incomplete taxa ...[and] adding incomplete characters ... Identifying the mechanisms that may cause incomplete taxa and characters to be problematic is an important step in devising effective solutions."* - Wiens (2003)

> *"As fossil taxa will necessarily introduce a large proportion of missing data in any combined data set (Wiens 2009), the impact of these absent characters on branch length estimation and support may be substantial. In particular, the interaction between among partition rate variation and missing data may negatively affect phylogenetic inference (Lemmon et al. 2009). Additionally, heterogeneity in evolutionary rates between the morphological characters may also affect estimated branch lengths (Clarke and Middleton 2008)."* - Pyron (2011)

Missing data is an important consideration in the theory and practice of phylogenetic systematics, as well as the design of phylogenetics studies (Wiens 2003; Wiens et al. 2005; Wiens 2006; Pyron 2011). This software automates exploratory analyses calculating the amount of missing data in phylogenetic datasets (NEXUS character partitions), as well as testing the form and potential significance of relationships between missing characters and phylogenetic tree parameters. Parts of MissingDataFX were inspired by, and recreate, correlational/exploratory analyses of the effects of missing data on phylogenies used previously by Wiens et al. (2005) and Pyron (2011). I wrote MissingDataFX code to help me automate these and related analyses for a recent project on molecular phylogenetics and Bayesian total-evidence dating, as well as the effects of missing data on phylogenetic results, in 'sucker' fishes of the family Catostomidae (Bagley et al., in revision).

The current version of this software focuses on developing **MissingDataFX.sh**, a shell script that looks at potential effects of missing data on phylogenetic analyses in two basic ways: 1) by characterizing the contents of data blocks in a NEXUS input file, including proportions of data versus missing data for each partition, and 2) by performing a customized R analysis. Three main operations are performed in the R environment: a) extracting tree parameters (terminal branch lengths, terminal/subtending branch heights, posterior support) from BEAST or MrBayes trees generated from the input NEXUS; b) using appropriate standard or nonparametric tests for correlations between missing data proportions and the tree parameters; and c) plotting relationships among variables. More details are given below. 

As in the case of the author's software package for phylogenetic and phylogeographic inference, [PIrANHA](https://github.com/justincbagley/PIrANHA), the MissingDataFX package is fully command line-based and is available as open-source software according to the license. 

## GETTING STARTED

### Dependencies

![alt tag](https://cloud.githubusercontent.com/assets/10584087/21243724/b12a94d2-c2df-11e6-9d20-5cf06877ad94.png)  

Code in the MissingDataFX repository depends on R and a series of R packages, as follows:
- The [R](https://www.r-project.org) software environment
  * available from download mirrors from The Comprehensive R Archive Network (CRAN), such as https://cloud.r-project.org
- tools - a default R package with "base" priority
- [ape](https://cran.r-project.org/web/packages/ape/index.html)
- [ips](https://cran.r-project.org/web/packages/ips/index.html)
- [phytools](https://cran.r-project.org/web/packages/phytools/index.html)
- [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html)

### Installation
:computer: MissingDataFX uses UNIX shell and R scripts, thus runs on a variety of operating systems, but was especially designed with UNIX/LINUX-like systems in mind. MissingDataFX code should run "out-of-the-box" from most any folder on your machine. To 'install' MissingDataFX, download the repository, move into the repository folder and enter ```$ chmod u+x ./*.sh``` into the command line interface (e.g. Terminal app on Mac). This changes file mode bits in the .sh files to allow the user permission to access and execute them, which is an important prerequisite for the MissingDataFX code to run.

### Run directory structure
:warning: MissingDataFX code assumes that you are running within a sub-folder specific to your analysis, located within the MissingDataFX-master distro folder. **Importantly**, as a result of this directory structure, the "R" folder from the distro, which contains a modified function for ips, will be located at the relative path "../R/". See README and in-script commenting for further details on file types and instructions for running basic analyses. Ideally, users will run MissingDataFX on multiple NEXUS-tree file combinations, for example 1) mtDNA only, 2) nuclear DNA only, 3) morphological characters only, and 4) all data combined--or a 'combined data' (mtDNA + nuclear) or 'total-evidence' (sequence + morphology) matrix and tree. Each of these analyses would be run within a separate sub-folder in the master distro folder.

### Input files and filenames
:warning: The main input files for MissingDataFX are NEXUS files, tree files, and 'drop' files. **This software assumes that the current working directory (i.e. sub-folder) for any particular run contains _ONLY_ the following files\:** 

1) **MissingDataFX.sh script (copy/paste into dir)**
2) **Input NEXUS file**
3) **Input BEAST or MrBayes tree file**
4) ***Optional* 'drop file' (one per tree file)**

Please follow the guidelines below when constructing input files for analysis.

#### NEXUS file
For the NEXUS input file, we assume that there are no spaces or special characters in the filename (though underscores are OK), and that the file contains only DNA sequences or morphological data in simplified NEXUS format--i.e. header followed by a matrix block, and no subsequent data/info blocks (e.g. sets or MrBayes blocks; these must be removed). The input NEXUS may combine morphological and DNA sequence characters and it doesn't matter how these are specified after 'FORMAT DATATYPE' (so long as DATATYPE input sits on a single line); however, if one or more morphological data blocks are included in the NEXUS, then users **must NOT include a 'CharStateLabels' or 'CHARACTERS' block** in the file, or there will be issues. Prior to running MissingDataFX, users should check their input NEXUS filenames and contents to ensure that they meet these assumptions.

#### Tree file
Regarding tree files, let me say I'm very much a Bayesian, and my empirical research focuses frequently entails making model-based inferences in phylogenetics and phylogeography based on Bayesian or approximate Bayesian computation (e.g. ABC, hABC) methods (see my [Research](http://www.justinbagley.org/research), [Publications](http://www.justinbagley.org/publications), and [ResearchGate](https://www.researchgate.net/profile/Justin_Bagley2) pages). Thus, it's no surprise that the initial development offering of MissingDataFX **Only** provides support for Bayesian phylogenetic trees. Accepted tree files include BEAST maximum clade credibility (MCC) trees annotated with posterior parameter distributions in TreeAnnotator, and MrBayes consensus tree files generated by running sumt command (only tested with version 3.2+). **Multiple tree file extensions are supported; however, ...**
- **You may use any MrBayes consensus tree file with the default '.con.tre' extension, or you may change the extension to '.tree'.** 
- **By contrast, users *must* rename any BEAST tree file with '.out' (TreeAnnotator) extension to have the extension '.tree' instead.**

#### Drop file
*Optional* 'drop files' are files the user creates in a text editor that contain a list of taxa (with names exactly matching taxon labels in the NEXUS and tree files), with one taxon name per line followed by an empty line, and saved with '.drop' file extensions. Drop files specify taxa to be pruned from the phylogeny in R prior to further analyses. Drop files provide additional flexibility to the analysis, and they are handy when the user would like to exclude certain taxa, for example a) extinct taxa or b) taxa with low or high amounts of missing data. Why extinct taxa? Some BEAST or MrBayes trees may result from analyses including extinct taxa as tips, e.g. tip-dating analyses (Pyron 2011) or Bayesian total-evidence dating using fossilized birth-death (FBD) models (Heath et al. 2014). In such cases, the extinct taxa usually have tip dates older than a few hundred to thousands of years ago (>300-1000 yr BP), and this produces a set of non-contemporaneous tips in resulting phylogenies--enough to majorly effect the distribution of branch lengths (especially with more than just a 1-5 extinct taxa). **Remember:**
- **Drop files are optional.** MissingDataFX will run just fine whether they are included or not.
- **At present, there can only be *one* '.drop' file and *one* tree file per input NEXUS analyzed in MissingDataFX.** So, if you want to analyze two NEXUS files with only one corresponding Bayesian tree file (e.g. you're running different alignments from a partitioned Bayesian analysis that yielded a single MrBayes or BEAST tree), then the tree file MUST be duplicated and given separate names matching each NEXUS, and sets of files for each analysis should be run in separate sub-folders, as per 'Run directory structure' section above.

#### Guidelines for input file NAMES
For simplicity, all filenames supplied to MissingDataFX for a given analysis should have the same basename. A suitable set of files for analysis might look like the following example, where MY_BASENAME is the basename applied to all of the files:  

| Input file type        | Example filename                                                    |
| :--------------------- |:--------------------------------------------------------------------|
| NEXUS file             | MY_BASENAME.nex                                                     |
| Tree file              | MY_BASENAME.tree (BEAST, MrBayes), *OR* MY_BASENAME.con.tre (MrBayes) |
| Drop file              | MY_BASENAME.drop                                                    |

We use this naming convention so that the tree filenames can be linked to the original NEXUS input file(s) without conflicting with the NEXUS filenames used in other procedures used in the shell script.

### What happens in R?
Following several steps summarizing the NEXUS input file using operations in the shell, **MissingDataFX.sh** creates and runs a customized R script that loads the 'Interfaces to Phylogenetic Software in R' or 'ips' R package (plus my fixed version of one of its functions) and related phylogenetics packages, and then reads in the tree, plots the tree (and saves to file), extracts the node and branch length annotations into a matrix (also saved to file). Next, [Shapiro-Wilk tests](https://en.wikipedia.org/wiki/Shapiro–Wilk_test) are conducted to evaluate whether the log-transformed data and tree parameters meet normality criteria for subsequent analyses. Then MissingDataFX tests for correlations between the proportions of data or missing data and (y-axis/independent var.:) 1) posterior support for terminal taxa and 2) length of terminal branch (same as height estimate for terminal taxon tmrca/node, in case of BEAST chronograms). If the data are normal (Shapiro-Wilk tests were non-significant at the alpha = 0.05 level), regular [Pearson correlations](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) and [generalized linear modeling](https://en.wikipedia.org/wiki/Generalized_linear_model) analyses are conducted; however, if the data are non-normal (Shapiro-Wilk tests were significant), then correlations are conducted using nonparametric [Spearman's rank correlation coefficient](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient) and linear relationships among log-transformed variables are plotted, but not given trendlines. Results are output to file in text or PDF (graphics) files. 

### Usage
It's so easy to use, the MissingDataFX.sh script doesn't display any sophisticated Usage or help flag info yet, because it doesn't need to. Everything you need to know is given here in the README! Assuming that you have installed the dependencies and the repo and followed guidelines for input files above, you may run MissingDataFX by simply entering the script name at the command line and hitting return!
```
$ ./MissingDataFX.sh
```
#### Screen output example:
Here is an example of what output to screen would look like for a normal MissingDataFX analysis including a NEXUS file and a MrBayes tree file with '.tree' extension, and no drop file:
```
$ ./MissingDataFX.sh

##########################################################################################
#                            MissingDataFX v0.1.0, March 2017                            #
##########################################################################################

INFO      | Tue Mar 14 15:35:18 CDT 2017 | Starting MissingDataFX analysis... 
INFO      | Tue Mar 14 15:35:18 CDT 2017 | STEP #1: SETUP AND USER INPUT. 
INFO      | Tue Mar 14 15:35:18 CDT 2017 |          Setting working directory to: /Users/justinbagley/Documents/GitHub/MissingDataFX-master-Feb17/4loci_new_50mil_long_run1-PARTIAL46mill_MrBayes_tree2 
INFO      | Tue Mar 14 15:35:18 CDT 2017 |          Reading in input NEXUS file(s)... 
INFO      | Tue Mar 14 15:35:18 CDT 2017 | STEP #2: PROCESSING INPUT NEXUS, SPLITTING TAXON LABELS AND DATA BLOCKS INTO SEPARATE FILES. 
INFO      | Tue Mar 14 15:35:18 CDT 2017 |          Catostomidae_4loci_new_forMrBayes.NEX 
INFO      | Tue Mar 14 15:35:18 CDT 2017 | STEP #3: USING TAXON LABEL AND CONCATENATED SEQUENCE FILES GENERATED DURING PREVIOUS STEP TO CREATE ONE 
INFO      | Tue Mar 14 15:35:18 CDT 2017 |          FILE WITH MISSING DATA COUNTS AND PROPORTIONS FOR EACH INDIVIDUAL. 
INFO      | Tue Mar 14 15:35:19 CDT 2017 | STEP #4: PRE-PROCESSING MRBAYES CONSENSUS TREE INPUT FILE, IF PRESENT: SPLIT .con.tre FILE, EXTRACT 
INFO      | Tue Mar 14 15:35:19 CDT 2017 |          TERMINAL BRANCH LENGTHS & THEIR CONFIDENCE INTERVALS, AND CREATE BRANCH LENGTH SUMMARY TABLE. 
INFO      | Tue Mar 14 15:35:19 CDT 2017 |          Encountered one or more '.tree' files. Testing them... 
INFO      | Tue Mar 14 15:35:19 CDT 2017 |          No '.con.tre' tree file in current working directory. Checking for '.tree' file... 
INFO      | Tue Mar 14 15:35:19 CDT 2017 |          Your tree file looks like it is from MrBayes. Assuming MrBayes tree(s) available hereafter... 
INFO      | Tue Mar 14 15:35:19 CDT 2017 |          Running mbTreeStatMiner script modified for '.tree' files... 
INFO      | Tue Mar 14 15:35:19 CDT 2017 | STEP #5: MAKE R SCRIPT THAT A) EXTRACTS PARAMETER ESTIMATES FROM BEAST TREES OR READS IN MRBAYES DATA
INFO      | Tue Mar 14 15:35:19 CDT 2017 |          TABLE (STEP #3) IN WORKING DIR THEN B) TESTS FOR IMPACT OF MISSING DATA ON PHYLO SUPPORT AND BRANCH LENGTHS. 
INFO      | Tue Mar 14 15:35:19 CDT 2017 | STEP #6: RUN THE R SCRIPT (WHICH ALSO SAVES RESULTS TO FILE). 
INFO      | Tue Mar 14 15:35:25 CDT 2017 | STEP #7: CLEANUP: ORGANIZE RESULTS, REMOVE UNNECESSARY FILES. 
INFO      | Tue Mar 14 15:35:25 CDT 2017 |          R workspace file confirmed in dir: /Users/justinbagley/Documents/GitHub/MissingDataFX-master-Feb17/4loci_new_50mil_long_run1-PARTIAL46mill_MrBayes_tree2 
INFO      | Tue Mar 14 15:35:25 CDT 2017 | Done analyzing the amount and potential effects of missing data on phylogenetic support and branch 
INFO      | Tue Mar 14 15:35:25 CDT 2017 | lengths using MissingDataFX. 
INFO      | Tue Mar 14 15:35:25 CDT 2017 | Bye.
```

## TROUBLESHOOTING
How to troubleshoot some potentially common problems:

(**1**) One problem you may run into is the 'permission denied' error returned by the shell when attempting to execute the shell scripts; for example, bash might return "-bash: ./MissingDataFX.sh: Permission denied". This indicates that permissions for .sh files were not correctly set during [Installation](https://github.com/justincbagley/MissingDataFX#installation). Fix this by moving into the repository master folder and entering ```$ chmod u+x ./*.sh``` at the command line, then using the shell script(s) in the repo folder from that point forward.

## ACKNOWLEDGEMENTS
During the development of this software, J.C.B. received stipend support from a Ciência Sem Fronteiras (Science Without Borders) postdoctoral fellowship from the Brazilian Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq; Processo 314724/2014-1). Lab and computer space was also supplied by The University of Alabama, during an internship in the Lozier Lab in the UA Department of Biological Sciences.

## REFERENCES
- Bagley JC, Mayden RL, Harris PM (in revision) Phylogeny, temporal diversification, and conflicting relationships of suckers (Cypriniformes: Catostomidae) inferred from Bayesian total-evidence dating. Molecular Phylogenetics and Evolution.
- Clarke JA, Middleton K (2008) Mosaicism, modules, and the evolution of birds: results from a Bayesian approach to the study of morphological evolution using discrete character data. Systematic Biology, 57, 185-201.
- Lemmon AR., Brown JM, Stanger-Hall K, Lemmon EM (2009) The effect of missing data on phylogenetic estimates obtained by maximum-likelihood and Bayesian inference. Systematic Biology, 58, 130-145.
- Pyron RA (2011) Divergence time estimation using fossils as terminal taxa and the origins of Lissamphibia. Systematic Biology, 60, 466-481.
- Wiens JJ (2003) Incomplete taxa, incomplete characters, and phylogenetic accuracy: is there a missing data problem? Journal of Vertebrate Paleontology, 23, 297-310.
- Wiens JJ (2006) Missing data and the design of phylogenetic analyses. Journal of Biomedical Informatics, 39, 34-42.
- Wiens JJ, Fetzner JW, Parkinson CL, Reeder TW (2005) Hylid frog phylogeny and sampling strategies for speciose clades. Systematic Biology, 54, 719-748.

## TODO LIST
**Current:**

- Provide more graphical output options:
 o Summaries of data/missing data proportions for different data blocks (though it's difficult to carry names of blocks/partitions into R).
 o Linear relationships between posterior support and data characteristics, at the level of blocks/partitions rather than individual nodes. 
- Modify to accomodate and extract parameters from annotated maximum-likelihood trees from programs like RAxML (check extent and format of branch/node labels on ML trees from different programs).
- Add screenshots of R analyses, as well as example text and graphical output produced by the software.
- Allow MrBayes tree files renamed with extension '.tree' (currently automatically converts '.con.tre' to '.tree').

**Recently finished/fixed:**

- Solve two input file problem--accomodate MrBayes trees, in addition to BEAST MCC trees. **DONE!** :white_check_mark:
- Specifically accomodate MrBayes tree files with '.con.tre' or '.tree' extensions. **DONE!** :white_check_mark:
- Correct issue in R code with 'ALL data' dataframe name for the non-drop file case. **DONE!** :white_check_mark:

March 14, 2017
Justin C. Bagley, Tuscaloosa, AL, USA
