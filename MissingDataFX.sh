#!/bin/sh

##########################################################################################
#                            MissingDataFX v0.1.1, March 2017                            #
#   SHELL SCRIPT FOR CONDUCTING TESTS OF THE EFFECTS OF MISSING DATA ON PHYLOGENETIC     #
#   ANALYSES (NODAL SUPPORT, BRANCH LENGTHS)                                             #
#   Copyright (c)2017 Justin C. Bagley, Universidade de Brasília, Brasília, DF, Brazil.  #
#   See the README and license files on GitHub (http://github.com/justincbagley) for     #
#   further information. Last update: March 15, 2017. For questions, please email        #
#   jcbagley@unb.br.                                                                     #
##########################################################################################

echo "
##########################################################################################
#                            MissingDataFX v0.1.1, March 2017                            #
##########################################################################################
"

######################################## START ###########################################
echo "INFO      | $(date) | Starting MissingDataFX analysis... "
echo "INFO      | $(date) | STEP #1: SETUP AND USER INPUT. "
###### Set paths and filetypes as different environmental variables:
	MY_PATH=`pwd -P`				## This script assumes it is being run in a sub-folder
									## of the MissingDataFX master directory specific to the 
									## current analysis. Hence, the R directory in the 
									## distro is relatively located at path "../R/".
echo "INFO      | $(date) |          Setting working directory to: $MY_PATH "
	CR=$(printf '\r');
	calc () {
	bc -l <<< "$@" 
}

echo "INFO      | $(date) |          Reading in input NEXUS file(s)... "
###### Read in the NEXUS file(s), by getting filenames matching .nex or .NEX patterns.

	MY_NEXUS_FILES="$(find . \( -name "*.nex" -o -name "*.NEX" -o -name "*.Nex" \) -type f | sed 's/\.\///g')"

echo "INFO      | $(date) | STEP #2: PROCESSING INPUT NEXUS, SPLITTING TAXON LABELS AND DATA BLOCKS INTO SEPARATE FILES. "
###### Several things to do here. First, 1) is the data interleaved? Extract info on whether 
##--file contains interleaved data or not from data block at start of NEXUS file. 2) Isolate 
##--the taxon labels (which are in the same exact order across all data blocks) and send
##--them to a separate file named "taxonLabels.txt". 3) Isolate sequence data block(s)
##--and send them each to a block named "seqDataBlock0n.txt", where n ranges from 1 to the
##--total number of blocks. Do this using seqDataSplitter script, as function:

################################## seqDataSplitter.sh ####################################

	seqDataSplitter () {

	(
		for i in $MY_NEXUS_FILES; do
			echo "INFO      | $(date) |          $i "

			##--This is the base name of the original nexus file, so you have it. This WILL work 
			##--regardless of whether the NEXUS filename extension is written in lowercase or in 
			##--all caps, ".NEX".
			MY_NEXUS_BASENAME="$(echo $i | sed 's/\.\///g; s/\.[A-Za-z]\{3\}$//g')"


			###### CHECK WHETHER NEXUS FILE IS IN A) INTERLEAVED FORMAT OR B) NON-INTERLEAVED
			###### FORMAT, AND MODIFY FILE ACCORDINGLY TO SPLIT TAXON LABELS AND DATA BLOCKS
			###### INTO SEPARATE FILES, THEN GET ONE FINAL SEQUENCE FILE.
			## MY_NEXUS_INTERL_STATUS="$(grep -o 'interleave\=[A-Za-z]\{2,\}' $i | sed 's/.*\=//g; s/\ //g')"
			MY_NEXUS_INTERL_STATUS="$(grep -o 'interleave\=[A-Za-z]\{2\}\|INTERLEAVE\=[A-Za-z]\{2\}\|Interleave\=[A-Za-z]\{2\}' $i | sed 's/.*\=//g; s/\ //g')"
			
			if [[ "$MY_NEXUS_INTERL_STATUS" = "ye" ]] || [[ "$MY_NEXUS_INTERL_STATUS" = "Ye" ]] || [[ "$MY_NEXUS_INTERL_STATUS" = "YE" ]]; then

			## SECTION A. INTERLEAVED FILES:
				##--Gather useful information for splitting file. Count number of lines in header (from ^#NEXUS to ^MATRIX) and whole file:
				NUM_HEADER_LINES="$(grep -n "matrix\|MATRIX" $i | sed 's/:.*$//g')"
				MY_EOF_LINE="$(wc -l $i | sed 's/\.\/.*$//g')"
			 			
				##--Remove header from NEXUS file:
				tail -n+"$(calc $NUM_HEADER_LINES+1)" $i > ${MY_NEXUS_BASENAME}_headless.nex

				##--Count total number of sequence data block lines in headless file, and 
				##--save sequence data to separate file:
				NUM_SEQ_DATABLOCK_LINES="$(calc $(grep -n '\;' ${MY_NEXUS_BASENAME}_headless.nex | \
				head -n1 | sed 's/:.*$//g')-1)"
				head -n$NUM_SEQ_DATABLOCK_LINES ${MY_NEXUS_BASENAME}_headless.nex > ${MY_NEXUS_BASENAME}_seqData.tmp
				
				##--Count total number of blank line occurrences, which mark the end of
				##--each data block
				HEADLESS_BLANK_LINES="$(grep -n '^$' ${MY_NEXUS_BASENAME}_headless.nex)"
				NUM_HEADLESS_BLANK_LINES="$(calc $(grep -n '^$' ${MY_NEXUS_BASENAME}_headless.nex | \
				wc -l | sed 's/\ //g')-1)"
				
					##--Make START line, and then add it to beginning of temporary seqData files
					##--to make a new seqData file. Then edit the blank lines so that they
					##--say START instead and save to new, final seqData file.
					echo "START" > START.tmp
					sed 's/^$/START/g' ${MY_NEXUS_BASENAME}_seqData.tmp > ${MY_NEXUS_BASENAME}_seqData2.tmp
					cat START.tmp ${MY_NEXUS_BASENAME}_seqData2.tmp > ${MY_NEXUS_BASENAME}_seqData.txt
					rm ./*.tmp

					##--Split sequence data blocks into separate files with "seqDataBlock" base
					##--followed by an integer number:
					awk '/START/{x="seqDataBlock"++i;next}{print > x;}' ${MY_NEXUS_BASENAME}_seqData.txt
				
					##--Extract taxon labels from first seqDataBlock, "seqDataBlock1":
					sed 's/\ [^\ ].*//g' ./seqDataBlock1 > ${MY_NEXUS_BASENAME}_taxonLabels.txt

					##--Rename seqDataBlock files so that they contain the NEXUS basename, and
					##--while you are looping through them, also remove the taxon labels from
					##--each file:
					MY_SEQDATABLOCK_FILES=./seqDataBlock*
					(
						for j in $MY_SEQDATABLOCK_FILES; do
							MY_SEQDATA_BASENAME=$(echo $j | sed 's/\.\///g')
							sed 's/^[A-Za-z0-9\_]*\(\ .*\)/\1/g' $j | sed 's/\ //g' > ${MY_NEXUS_BASENAME}_${MY_SEQDATA_BASENAME}.txt
						done
					)				
				
				##--Paste all ${MY_NEXUS_BASENAME}_${MY_SEQDATA_BASENAME}.txt files together
				##--so that you are left with all the sequence data (from all seqData blocks)
				##--for each individual on a single line, as is required in non-interleaved
				##--phylip format. The resulting file has no taxon labels or spaces.
				MY_EDITED_SEQDATABLOCK_FILES=./${MY_NEXUS_BASENAME}_seqDataBlock*.txt
				paste $MY_EDITED_SEQDATABLOCK_FILES | sed $'s/\t//g' > ${MY_NEXUS_BASENAME}_concatenatedSeqs.txt


				##--Remove all the original "seqDataBlock*" and "seqData.txt" files from 
				##--the working directory:
				rm ./seqDataBlock*
				rm ./*_seqDataBlock*
				rm ./*_seqData.txt
				
			elif [[ "$MY_NEXUS_INTERL_STATUS" = "no" ]] || [[ "$MY_NEXUS_INTERL_STATUS" = "No" ]] || [[ "$MY_NEXUS_INTERL_STATUS" = "NO" ]]; then

			## SECTION B. NON-INTERLEAVED FILES:
				cat $i > ${MY_NEXUS_BASENAME}_concatenatedSeqs.txt
					
			fi
					
		done
	)

}

##--Don't forget to run the function!
seqDataSplitter


echo "INFO      | $(date) | STEP #3: USING TAXON LABEL AND CONCATENATED SEQUENCE FILES GENERATED DURING PREVIOUS STEP TO CREATE ONE "
echo "INFO      | $(date) |          FILE WITH MISSING DATA COUNTS AND PROPORTIONS FOR EACH INDIVIDUAL. "
###### OK, All we need to do at this point is count the number of a) characters and b) missing 
##--character states and gaps (we treat gaps as missing! see Felsenstein's Inferring
##--Phylogenies, https://groups.google.com/forum/#!topic/beast-users/ixrGUA1p4OM/, and other
##--resources for more info) for each individual (=each line, which is a 1:1 match to the 
##--taxonLabels.txt file) in the "...concatenatedSeqs.txt" files; and paste the taxonLabels 
##--and reg/missing character counts together (as "label'\t'count") in a final output file, 
##--with one such file per concatenatedSeqs/NEXUS file. We do this using my summarizeSeqData
##--script, placed in a function, as follows:

################################## summarizeSeqData.sh ###################################

	summarizeSeqData () {

	MY_CONCATSEQ_FILES=./*_concatenatedSeqs.txt
	(
		for k in $MY_CONCATSEQ_FILES; do
			MY_CONCATSEQ_BASENAME=$(echo $k | sed 's/\.\///g; s/\_concatenatedSeqs\.txt//g')

			##### SAVE VALUES IN ENVIRONMENTAL VARIABLES:
			TOTAL_NUM_CHAR="$(awk -F '[A-Za-z\?\-]' '{print NF-1}' $k | sort -rn | head -n1)"
			## echo $TOTAL_NUM_CHAR
			## echo $TOTAL_NUM_CHAR > CharacterTotal.txt
			NUM_MISSING_CHAR="$(awk -F '\?' '{print NF-1}' $k)"
			## echo $NUM_MISSING_CHAR
			NUM_GAP_CHAR="$(awk -F '\-' '{print NF-1}' $k)"
			## echo $NUM_GAP_CHAR
			
			sed 's/\-//g; s/\?//g' $k > ${MY_CONCATSEQ_BASENAME}_regCharFILE.txt
			NUM_REG_CHAR="$(awk -F '[A-Za-z]' '{print NF-1}' ${MY_CONCATSEQ_BASENAME}_regCharFILE.txt)"
			
			##### SAVE VALUES INTO FILES:
			awk -F '\?' '{print NF-1}' $k > ${MY_CONCATSEQ_BASENAME}_missingChar.txt
			awk -F '\-' '{print NF-1}' $k > ${MY_CONCATSEQ_BASENAME}_gapChar.txt
			awk -F '[A-Za-z]' '{print NF-1}' ${MY_CONCATSEQ_BASENAME}_regCharFILE.txt | sed 's/\-1/0/g' > ${MY_CONCATSEQ_BASENAME}_regCharCOUNTS.txt
			awk -F '[A-Za-z\?\-]' '{print NF-1}' $k > ${MY_CONCATSEQ_BASENAME}_totalChar.tmp
			awk -v total="$TOTAL_NUM_CHAR" '{$2=total; print}' ${MY_CONCATSEQ_BASENAME}_totalChar.tmp > ${MY_CONCATSEQ_BASENAME}_totalChar.txt

				##### FIX TAXON LABELS THEN PASTE FILES TOGETHER:
				sed -E $'s/\ +/\t/g' ${MY_CONCATSEQ_BASENAME}_taxonLabels.txt | sed $'s/\t//g' > ${MY_CONCATSEQ_BASENAME}_taxa.txt
				rm ${MY_CONCATSEQ_BASENAME}_taxonLabels.txt
				
				paste ./${MY_CONCATSEQ_BASENAME}_taxa.txt ./${MY_CONCATSEQ_BASENAME}_regCharCOUNTS.txt ./${MY_CONCATSEQ_BASENAME}_missingChar.txt ./${MY_CONCATSEQ_BASENAME}_gapChar.txt ${MY_CONCATSEQ_BASENAME}_totalChar.txt > ${MY_CONCATSEQ_BASENAME}_BIGSummary.txt
				awk -v total="$TOTAL_NUM_CHAR" '{ print $1, 100*$2/total, 100*($3+$4)/total }' ${MY_CONCATSEQ_BASENAME}_BIGSummary.txt > ${MY_CONCATSEQ_BASENAME}_propSummary.txt
				awk -v total="$TOTAL_NUM_CHAR" '{ print $1, 100*$2/total }' ${MY_CONCATSEQ_BASENAME}_BIGSummary.txt > ${MY_CONCATSEQ_BASENAME}_propData.txt
				awk -v total="$TOTAL_NUM_CHAR" '{ print $1, 100*($3+$4)/total }' ${MY_CONCATSEQ_BASENAME}_BIGSummary.txt > ${MY_CONCATSEQ_BASENAME}_propMissing.txt		## This file outputs the correct estimate of proportion missing data, because it (correctly) treats gaps as missing, contributing nothing to phylogenetic estimation etc (see refs/URL listed above).
				
		done
	)
}

##--Don't forget to run the function!
summarizeSeqData



echo "INFO      | $(date) | STEP #4: PRE-PROCESSING MRBAYES CONSENSUS TREE INPUT FILE, IF PRESENT: SPLIT .con.tre FILE, EXTRACT "
echo "INFO      | $(date) |          TERMINAL BRANCH LENGTHS & THEIR CONFIDENCE INTERVALS, AND CREATE BRANCH LENGTH SUMMARY TABLE. "
##--This section helps accomodate consensus tree files with extension "*.con.tre" output 
##--directly from MrBayes into MissingDataFX analysis. Although functions exist for reading 
##--in BEAST trees and extracting a table of BEAST node label data in the 'ips' R package,
##--the only suitable R functions that I've found tailored to MrBayes trees allow reading in
##--the tree; however, there is no function that summarizes or reads in MrBayes tree statistics.
##--So, this section is a necessary evil. It first checks for a MrBayes tree file, and if 
##--one is found, it runs my new mbTreeStatMiner script (placed in a function), which
##--extracts and links tip label and branch length information from a .con.tre file in 
##--current working directory, and outputs a table of the data in R-readable format. 
shopt -s nullglob
if [[ -n $(find . -name "*.con.tre" -type f) ]]; then

	echo "INFO      | $(date) |          Found '.con.tre' file. Reading in MrBayes consensus tree from current working directory... "
	##--Starting from file output directly from MrBayes, with extension "*.con.tre", present
	##--in the current working directory...
	MY_MRBAYES_CONSENSUS=./*.con.tre

	echo "INFO      | $(date) |          Running mbTreeStatMiner script... "

################################### mbTreeStatMiner.sh ###################################

	mbTreeStatMiner () {
		(
			for i in $MY_MRBAYES_CONSENSUS; do											## *IMPORTANT* NOTE: Here, we use a for loop, but we still assume that there is only a *single* MrBayes tree file in pwd (for now; check later to see if you can get to work on multiple files).

				###### COLLECT DATA ON THE INPUT TREE FILE AND SPLIT INTO SEPARATE FILES FOR PROCESSING.
				MY_CONSENSUS_BASENAME="$(ls ${i} | \
				sed 's/\.\///g; s/\.con\.tre$//g')"										## Get basename of MrBayes tree file, for later use.
				cp $i "$(ls $i | sed 's/\.con\.tre/\.tree/g')"							## We rename the MrBayes consensus trees to have the extension (.tree) that we want to use downstream... e.g. in MissingData.
				MY_MB_DOTTREE_FILE=$MY_CONSENSUS_BASENAME.tree
				
				grep "tree\ " $i > ./${MY_CONSENSUS_BASENAME}_treeDataBlock1.txt		## Move tree block to separate file ending in "_treeDataBlock1.txt."

				MY_NTAX="$(grep -n "ntax" $i | awk -F"=" '{print $NF}' | sed 's/\;//g')"
				NLINES_TO_END_TAXLABS_TRANS="$(calc $(sed -n '/^[\ ]*tree\ /{=; q;}' $i)-1)"
				NLINES_TO_END_FIRST_TAXLABS="$(calc $(sed -n '/^end\;/{=; q;}' $i)-2)"
				NLINES_TOTAL=$(wc -l $i | sed 's/\ \.\/[A-Za-z\.0-9\_]*//')		
				
				sed -n 6,"$NLINES_TO_END_FIRST_TAXLABS"p $i > ./taxLabels.tmp
				sed $'s/\t//g' ./taxLabels.tmp > ./${MY_CONSENSUS_BASENAME}_taxLabels.txt	## Save taxon label names in order, delete rest of each line, send to file ending in "_taxLabels.txt."


				###### NOW MODIFY INITIAL treeDataBlock FILE FROM ABOVE TO GET TERMINAL BRANCH LENGTH DATA (ONE SET PER TIP TAXON) INTO USEFUL FORMAT.
				##--Section A. Perform series of gross text changes for reformatting:
					sed 's/(//g' ./${MY_CONSENSUS_BASENAME}_treeDataBlock1.txt | \
					perl -pe 's/\,([0-9]*\[\&)/\n$1/g' | \
					perl -pe 's/\ ([0-9]*\[\&)/\n$1/g' | \
					sed 1,2d | \
					sed 's/)\[\&.*//g' | \
					sed 's/\[\&.*\(\[\&length_mean.*length_median.*\)\}.*/\1/g' | \
					perl -pe 's/\[\&/\t/g; s/\,/\t/g' | \
					sed 's/length\_mean\=//g; s/length\_median\=//g; s/length\_95\%HPD\={//g' > ./${MY_CONSENSUS_BASENAME}_tipBranchData2.txt

				###### FINAL TEXT PROCESSING. Clean file by sorting by taxon label number, so that rows in the tip branch data file match those in the taxLabels.txt file.
				## Also make a final file combining the tip label and tip branch data, with a header indicating column contents.
				sort -n ./${MY_CONSENSUS_BASENAME}_tipBranchData2.txt > ./${MY_CONSENSUS_BASENAME}_tipBranchData3.txt
				perl -pe 's/^[0-9]*\t//g' ./${MY_CONSENSUS_BASENAME}_tipBranchData3.txt > ./${MY_CONSENSUS_BASENAME}_tipBranchData4.txt
				paste ./${MY_CONSENSUS_BASENAME}_taxLabels.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData4.txt > ./${MY_CONSENSUS_BASENAME}_tipBranchData5.txt
				
				##--Make the header:
				echo "tip.label term_brL term_brL_median term_brL_lower_95HPD term_brL_upper_95HPD" > ./header.txt
				cat ./header.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData5.txt > ./${MY_CONSENSUS_BASENAME}_tipBranchData6.txt
				perl -pe 's/\ /\t/g' ./${MY_CONSENSUS_BASENAME}_tipBranchData6.txt | \
				perl -pe 's/^/\"/g; s/\t/\"\t\"/g; s/$/\"/g' | perl -pe 's/^\"\"/\"/g' | \
				perl -pe 's/^\"$//g' > ./${MY_CONSENSUS_BASENAME}_term_brL_Summary.txt

				##--Remove all temporary and intermediary files from the working directory:
				rm ./taxLabels.tmp
				rm ./${MY_CONSENSUS_BASENAME}_treeDataBlock1.txt
				rm ./${MY_CONSENSUS_BASENAME}_taxLabels.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData3.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData4.txt
				rm ./header.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData2.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData5.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData6.txt

			done
		)
}

##--Don't forget to run the function!
mbTreeStatMiner


elif [[ -n $(find . -name "*.tree" -type f) ]]; then
echo "INFO      | $(date) |          Encountered one or more '.tree' files. Testing them... "


    echo "INFO      | $(date) |          No '.con.tre' tree file in current working directory. Checking for '.tree' file... "
    ##--If there is no '.con.tre' treefile from MrBayes in pwd, then 1) there should be a '.tree' 
    ##--file but 2) that file could be of MrBayes or BEAST format. In the case that the file is 
    ##--from BEAST, there would be no problem with continuing with present environment. However, 
    ##--if the .tree file is from MrBayes, then there will be errors since the '*_term_brL_Summary.txt'
    ##--file was only created above for MrBayes .con.tre files. To fix this and allow users to supply
    ##--MrBayes treefiles with the .tree extension, we need to test which program the .tree file in
    ##--pwd is from. If there is no tree file, then echo message and quit. We do this as follows:

		MY_DOTTREE_FILE_VAR="$(find . -name '*.tree' -type f)"
		MY_MRBAYES_TREEFILE_TEST="$(grep -h 'prob\=' $MY_DOTTREE_FILE_VAR | head | wc -l)"

		if [[ "$MY_MRBAYES_TREEFILE_TEST" -eq "0" ]]; then 
			echo "INFO      | $(date) |          Your tree file looks like it may be from BEAST. Assuming BEAST tree(s) available hereafter... "

		elif [[ "$MY_MRBAYES_TREEFILE_TEST" -eq "1" ]]; then 
			echo "INFO      | $(date) |          Your tree file looks like it is from MrBayes. Assuming MrBayes tree(s) available hereafter... "
			##--Now we run a second version of mbTreeStatMiner.sh modified for MrBayes trees having the 
			##--'.tree' file extension:

			echo "INFO      | $(date) |          Running mbTreeStatMiner script modified for '.tree' files... "

			############################### dottree_mbTreeStatMiner.sh ###############################

			dottree_mbTreeStatMiner () {
				(
					for i in $MY_DOTTREE_FILE_VAR; do
		
						###### COLLECT DATA ON THE INPUT TREE FILE AND SPLIT INTO SEPARATE FILES FOR PROCESSING.
						MY_CONSENSUS_BASENAME="$(ls ${i} | sed 's/\.\///g; s/\.tree$//g')"		## Get basename of MrBayes tree file, for later use.
						grep "tree\ " $i > ./${MY_CONSENSUS_BASENAME}_treeDataBlock1.txt		## Move tree block to separate file ending in "_treeDataBlock1.txt."

						MY_NTAX="$(grep -n "ntax" $i | awk -F"=" '{print $NF}' | sed 's/\;//g')"
						NLINES_TO_END_TAXLABS_TRANS="$(calc $(sed -n '/^[\ ]*tree\ /{=; q;}' $i)-1)"
						NLINES_TO_END_FIRST_TAXLABS="$(calc $(sed -n '/^end\;/{=; q;}' $i)-2)"
						NLINES_TOTAL=$(wc -l $i | sed 's/\ \.\/[A-Za-z\.0-9\_]*//')		
				
						sed -n 6,"$NLINES_TO_END_FIRST_TAXLABS"p $i > ./taxLabels.tmp
						sed $'s/\t//g' ./taxLabels.tmp > ./${MY_CONSENSUS_BASENAME}_taxLabels.txt	## Save taxon label names in order, delete rest of each line, send to file ending in "_taxLabels.txt."


						###### NOW MODIFY INITIAL treeDataBlock FILE FROM ABOVE TO GET TERMINAL BRANCH LENGTH DATA (ONE SET PER TIP TAXON) INTO USEFUL FORMAT.
						##--Section A. Perform series of gross text changes for reformatting:
							sed 's/(//g' ./${MY_CONSENSUS_BASENAME}_treeDataBlock1.txt | \
							perl -pe 's/\,([0-9]*\[\&)/\n$1/g' | \
							perl -pe 's/\ ([0-9]*\[\&)/\n$1/g' | \
							sed 1,2d | \
							sed 's/)\[\&.*//g' | \
							sed 's/\[\&.*\(\[\&length_mean.*length_median.*\)\}.*/\1/g' | \
							perl -pe 's/\[\&/\t/g; s/\,/\t/g' | \
							sed 's/length\_mean\=//g; s/length\_median\=//g; s/length\_95\%HPD\={//g' > ./${MY_CONSENSUS_BASENAME}_tipBranchData2.txt

						###### FINAL TEXT PROCESSING. Clean file by sorting by taxon label number, so that rows in the tip branch data file match those in the taxLabels.txt file.
						## Also make a final file combining the tip label and tip branch data, with a header indicating column contents.
						sort -n ./${MY_CONSENSUS_BASENAME}_tipBranchData2.txt > ./${MY_CONSENSUS_BASENAME}_tipBranchData3.txt
						perl -pe 's/^[0-9]*\t//g' ./${MY_CONSENSUS_BASENAME}_tipBranchData3.txt > ./${MY_CONSENSUS_BASENAME}_tipBranchData4.txt
						paste ./${MY_CONSENSUS_BASENAME}_taxLabels.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData4.txt > ./${MY_CONSENSUS_BASENAME}_tipBranchData5.txt
				
						##--Make the header:
						echo "tip.label term_brL term_brL_median term_brL_lower_95HPD term_brL_upper_95HPD" > ./header.txt
						cat ./header.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData5.txt > ./${MY_CONSENSUS_BASENAME}_tipBranchData6.txt
						perl -pe 's/\ /\t/g' ./${MY_CONSENSUS_BASENAME}_tipBranchData6.txt | \
						perl -pe 's/^/\"/g; s/\t/\"\t\"/g; s/$/\"/g' | perl -pe 's/^\"\"/\"/g' | \
						perl -pe 's/^\"$//g' > ./${MY_CONSENSUS_BASENAME}_term_brL_Summary.txt

						##--Remove all temporary and intermediary files from the working directory:
						rm ./taxLabels.tmp
						rm ./${MY_CONSENSUS_BASENAME}_treeDataBlock1.txt
						rm ./${MY_CONSENSUS_BASENAME}_taxLabels.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData3.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData4.txt
						rm ./header.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData2.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData5.txt ./${MY_CONSENSUS_BASENAME}_tipBranchData6.txt

					done
				)
			}

			##--Don't forget to run the function!
			dottree_mbTreeStatMiner

		fi

elif [[ ! -n $(find . -name "*.tree" -type f) ]]; then

		echo "WARNING!  | $(date) |          No tree file present in current working directory. Quitting..."
		exit

fi


echo "INFO      | $(date) | STEP #5: MAKE R SCRIPT THAT A) EXTRACTS PARAMETER ESTIMATES FROM BEAST TREES OR READS IN MRBAYES DATA"
echo "INFO      | $(date) |          TABLE (STEP #3) IN WORKING DIR THEN B) TESTS FOR IMPACT OF MISSING DATA ON PHYLO SUPPORT AND BRANCH LENGTHS. "
###### Before getting started making the script, we first set some environmental variables
##--that FIX issues with echoing shell text containing dollar signs to R:
	MY_TIP_LABEL_VAR=$(echo "\$tip.label")		## Make tip.label variable with '$tip.label' text for Rscript...
	MY_NODE_VAR=$(echo "\$node")				## Make tip.label variable with '$node' text for Rscript...
	MY_POSTERIOR_VAR=$(echo "\$posterior")		## Same as above but for '$posterior'...
	MY_LENGTH_VAR=$(echo "\$length")			## Same as above but for '$length'...
	MY_HEIGHT_VAR=$(echo "\$height")			## Same as above but for '$height'...
	MY_PROP_DATA_VAR=$(echo "\$propData")			
	MY_PROP_MISSING_VAR=$(echo "\$propMissing")			
	MY_VAR_NAMES_VAR=$(echo "\$var_names")
	MY_PVALUE_VAR=$(echo "\$p.value")
	MY_TERM_BRL_VAR=$(echo "\$term_brL")
	
	##--Also prepare some variables specific to MrBayes branch length Summary file (if read):
	MY_MEDIAN_BRL_VAR=$(echo "\$term_brL_median")		## No var is specified here for mean values of terminal branch lengths from MrBayes, because these are already in a column named term_brL, which matches $MY_TERM_BRL_VAR above.
	MY_LOWER_BRL_HPD_VAR=$(echo "\$term_brL_lower_95HPD")
	MY_UPPER_BRL_HPD_VAR=$(echo "\$term_brL_upper_95HPD")

	
############ MAKE R SCRIPT
echo "
#!/usr/bin/env Rscript

################################# missingDataFXTester.R ##################################

############ I. SETUP
setwd('$MY_PATH')
#
##--Load needed packages, R code, or other prelim stuff. Install packages if not present.
packages <- c('tools', 'ape', 'ips', 'phytools', 'phangorn', 'ggfortify', 'ggplot2', 'gplots')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
}

library(tools)
library(ape)		## This is an ips dependency, so load it first.
library(ips)		## ips is the main package we use in this script.
library(phytools)
library(phangorn)	## Only used for 'Ancestor' function, but still useful.
library(ggfortify)	## Load these graphics packages, just in case we need them.
library(ggplot2)
library(gplots)

source('../R/read.mrbayes.FIX.r')	## Source modified version of ips 'read.mrbayes' routine
									## that fixes a bug in the code.

############ II. READ AND PREP TREE AND DATA FILE LISTS
##--Create vector lists of filenames for NEXUS, .tree, .drop, and proportion summary (i.e.
##--'*_propSummary.txt') files:
treeFileNames <- list.files(getwd(), pattern='tree$')
treeFileNames
NEXUSFileNames <- list.files(getwd(), pattern='nex$|NEX$')
NEXUSFileNames
dropFileNames <- list.files(getwd(), pattern='drop$')
dropFileNames
propSumFileNames <- list.files(getwd(), pattern='propSummary.txt$')
propSumFileNames

##--A couple of preliminary file checks:
if (length(treeFileNames)=='0') { print('WARNING! No tree file. Quitting...', quote=FALSE); quit('yes') }				## If if there are no tree or proportion summary files, quit session.
if (length(propSumFileNames)=='0') { print('WARNING! No data summary files. Quitting...', quote=FALSE); quit('yes') }	## 


############ IIIA. IF DROP FILE IS PRESENT, READ IN CHARACTER (PROPORTION) DATA AND DROP FILES
############ AND CONDUCT ANALYSES FOR THIS CASE
##--If a .drop file exists (i.e. if dropFileNames is not empty), then detect that and run 
##--a whole set of analysis code tailored for just when this condition is met. Use the 
##--basename of the .drop file as the same basename for the other file types when reading
##--in the files and assigning them to objects; this way we can ensure that no errors occur
##--due to different sorting/ordering of file groups (actual basenames may vary, in case of
##--multiple NEXUS files, for example) across file lists.
if (!length(dropFileNames)=='0') {
date()
print('DROP FILE DETECTED...CONDUCTING ANALYSES FOR THE DROP FILE CASE. ', quote=FALSE)
	for(i in 1:length(dropFileNames)){
		print('Reading in the .drop basename and contents...', quote=FALSE)
		basename <- file_path_sans_ext(dropFileNames[i])	## Get basename from the .drop file.
#		
		##--Loop through the .drop filenames and read each file in and assign it to an object:
		drop.in <- as.matrix(read.table(file=dropFileNames[i], header=FALSE, sep='	'))
		assign(paste('drop',i,basename,sep='_'), drop.in)
#
		print('Reading in the proportion summary data... ', quote=FALSE)
		##--Loop through the proportion summary filenames to read each file in and assign it to
		##--an object:
		prop.in <- as.matrix(read.table(file=propSumFileNames[i], header=FALSE, sep=' '))
		assign(paste('prop',i,basename,sep='_'), prop.in)
#
		print('Reading in tree file(s)... ', quote=FALSE)
		##--Loop through the .tree filenames to read each tree in and assign it to an appropriately 
		##--named object. However, we will do this using different ips functions depending whether the
		##--tree file is from BEAST or MrBayes ('read.beast' or 'read.mrbayes', respectively). So, we
		##--first need to test this. We'll grep patterns that are exclusive to the format of BEAST vs.
		##--MrBayes consensus NEXUS tree files, then use conditionals to guide how we fill 'tree.in'
		##--var below.
		beast_tree_post_patterns <- sapply(list.files(getwd(), pattern=treeFileNames[i]), FUN=function(x){grep('posterior', readLines(x))})
		mrbayes_tree_post_patterns <- sapply(list.files(getwd(), pattern=treeFileNames[i]), FUN=function(x){grep('prob', readLines(x))})
#
		if(is.numeric(beast_tree_post_patterns) == TRUE) { print('Your tree file looks like it is from BEAST.', quote=FALSE) 
			tree.in <- read.beast(file=treeFileNames[i], digits = NULL)
			treename <- paste('tree',i,basename,sep='_')
			assign(treename, tree.in)
			treename
			tree.in
#
#
		print('Reading in the tree annotation data... ', quote=FALSE)
		treedata.in <- read.beast.table(file=treeFileNames[i], digits = 4)
		treedata_name <- paste('treedata',i,basename,sep='_')
		assign(treedata_name, treedata.in)
		treedata_name
		treedata.in
#
#
print('PRINTING UNPRUNED AND PRUNED TREE GRAPHICS. ', quote=FALSE)
		print('Producing PDFs plotting the original tree with node and tip label numbers (for reference)... ', quote=FALSE)
		pdf(file=sprintf('%s.orig.lab.pdf',treename))
		plot(get(treename), no.margin=T, label.offset=0.1, cex=0.5)
		nodelabels(cex=0.5)
		tiplabels(cex=0.5)
		dev.off()
#
#		pdf(file=sprintf('%s.orig.lab.ladder.pdf',treename))
#		plot(ladderize(get(treename), right=FALSE), no.margin=T, label.offset=0.1, cex=0.5)
#		nodelabels(cex=0.5)
#		tiplabels(cex=0.5)
#		dev.off()
#
#
		print('Pruning tree(s) using the drop file(s)... ', quote=FALSE)
		as.data.frame(get(paste('drop',i,basename,sep='_')))
		pruned_tree <- drop.tip(get(paste('tree',i,basename,sep='_')),c(paste(as.data.frame(get(paste('drop',i,basename,sep='_')))[,1])))
#
		graphics.off()
		quartz()
		print('Printing pruned tree graphics... ', quote=FALSE)
		print('Producing PDFs of simple and ladderized pruned tree plots... ', quote=FALSE)
		pdf(file=sprintf('%s.pruned.pdf',treename))
		plot(pruned_tree, no.margin=F, label.offset=0.1, cex=0.5)
		axisPhylo()
		dev.off()
#
		pdf(file=sprintf('%s.pruned.ladder.pdf',treename))
		plot(ladderize(pruned_tree, right=FALSE), no.margin=F, label.offset=0.1, cex=0.5)
		axisPhylo()
		dev.off()
#
		graphics.off()
		quartz()
		print('Producing PDFs of pruned tree plots showing node and tip label numbers (for reference)... ', quote=FALSE)
		pdf(file=sprintf('%s.pruned.lab.pdf',treename))
		plot(pruned_tree, no.margin=T, label.offset=0.1, cex=0.5)
		nodelabels(cex=0.5)
		tiplabels(cex=0.5)
		dev.off()
#
		pdf(file=sprintf('%s.pruned.lab.ladder.pdf',treename))
		plot(ladderize(pruned_tree, right=FALSE), no.margin=T, label.offset=0.1, cex=0.5)
		nodelabels(cex=0.5)
		tiplabels(cex=0.5)
		dev.off()
#
#		
print('GETTING NODE NUMBERS FOR EACH TIP TAXON USING ORIGINAL, UNPRUNED TREE. ', quote=FALSE)
#
			##--Create function for picking the terminal node number for a given tip taxon. I made this
			##--function, 'terminalNode', by modifying the version of 'fastMRCA' function listed on
			##--this Liam Revell blog post (URL: http://blog.phytools.org/2012/05/function-to-efficiently-return-mrca-of.html/; 
			##--which is *NOT* identical to that of the current phytools distro, Feb 2017). The 'Ancestors'
			##--function is from phangorn, so phytools and phangorn should be loaded at start of script. 
			#
			##--NOTE: If desired, you may want to use the 'fastMRCA' function of phytools to get the
			##--number of the corresponding tip label (=shown in yellow along tips in labeled-tree PDF
			##--plots above, i.e. from calling nodelabels()). fastMRCA just matches the taxon name it
			##--is given to the corresponding tip.label.
			terminalNode <- function(tree, species) {
			    x <- match(species, tree$MY_TIP_LABEL_VAR)
			    y <- match(species, tree$MY_TIP_LABEL_VAR)
			    a <- Ancestors(tree, x)
			    b <- Ancestors(tree, y)
			    z <- a%in%b
			    return(a[min(which(z))])
			}
#
		##--Now that we have the 'terminalNode' function, we can use it to loop through all taxon
		##--tip.labels and save the corresponding terminal node numbers into 'terminal_node_numbers'. 
		terminal_node_numbers <- c()
		for(i in 1:length(get(treename)$MY_TIP_LABEL_VAR)){
			terminal_node_numbers[i] <- terminalNode(get(treename), get(treename)$MY_TIP_LABEL_VAR[i])
		}
#
print('COLLECTING TREE DATA FOR TERMINAL NODES. ', quote=FALSE)
		##--Next, use for loop to 1) match terminal_node_numbers to node column in treedata data frame,
		##--in order to get row numbers in treedata dataframe corresponding to each taxon; 2) get the
		##--a) posterior probability, b) terminal branch length, and c) estimated node height for the
		##--corresponding terminal node/tip; and 3) save all of this data into a list, with appropriate
		##--names for list elements.
		terminal_node_labels <- c()
		terminal_node_row_numbers <- c()
		terminal_node_posteriors <- c()
		terminal_node_lengths <- c()
		terminal_node_heights <- c()
		for(i in 1:length(get(treename)$MY_TIP_LABEL_VAR)){
			terminal_node_labels[i] <- c(get(treename)$MY_TIP_LABEL_VAR[i])
			terminal_node_row_numbers[i] <- match(terminal_node_numbers[i], as.data.frame(get(treedata_name))$MY_NODE_VAR)
			terminal_node_posteriors[i] <- as.data.frame(get(treedata_name))[c(terminal_node_row_numbers), ][i,]$MY_POSTERIOR_VAR
			terminal_node_lengths[i] <- as.data.frame(get(treedata_name))[c(terminal_node_row_numbers), ][i,]$MY_LENGTH_VAR
			terminal_node_heights[i] <- as.data.frame(get(treedata_name))[c(terminal_node_row_numbers), ][i,]$MY_HEIGHT_VAR
		}
#		
		terminal_node_data <- list(terminal_node_labels, terminal_node_numbers, terminal_node_posteriors, terminal_node_lengths, terminal_node_heights)
		names(terminal_node_data) <- c('tip.label', 'node', 'posterior', 'pre_brL', 'term_brL')
#			
##--Save R workspace to file in wd:
save.image(file='MissingDataFX_R_Workspace.RData')
print('GENERATING FINAL DATA FRAME BY MERGING TREE DATA WITH PROPORTION SUMMARY DATA. ', quote=FALSE)
		##--Now, we want to merge the 'prop.in' and 'terminal_node_data' objects/lists. First,
		##--let's give the proportion data (stored as specific to this itertation/set of files
		##--in the prop.in object) column names. Then, let's merge prop.in and the terminal
		##--node data, so that we can link all the following variables back to the tip labels:
		##--a) propData, b) propMissing, c) node, d) posterior, e) pre_brL, and f) term_brL.
		colnames(prop.in) <- c('tip.label', 'propData', 'propMissing')
#		
		##--Ignore next six lines for now, which reflect something I wanted to try (they do no harm, but I'm not really using them yet).
		terminal_node_data_labels <- c()
		terminal_node_data_row_numbers <- c()
		for(i in 1:length(get(treename)$MY_TIP_LABEL_VAR)){
			terminal_node_data_labels[i] <- c(terminal_node_data$MY_TIP_LABEL_VAR[i])
			terminal_node_data_row_numbers[i] <- match(terminal_node_data_labels[i], prop.in[i])
		}
#
		is.data.frame(merge(as.data.frame(prop.in), terminal_node_data, by='tip.label'))		## Merge idea was inspired by code I'd seen in phylogenetics packages, and this post (https://docs.tibco.com/pub/enterprise-runtime-for-R/1.5.0_may_2013/TERR_1.5.0_LanguageRef/base/merge.html).
		## Should produce the following output: 
		## [1] TRUE
		ALL_pruned_data <- merge(as.data.frame(prop.in), terminal_node_data, by='tip.label')
#
#
print('CONDUCTING EXPLORATORY AND CORRELATIONAL ANALYSES. ', quote=FALSE)
		##--We are interested in propData, propMissing, and term_brL as predictors of posterior; and
		##--we are also interested in propData and propMissing as predictors of term_brL. The node
		##--column is not needed and the other possible relationships between variables in the final
		##--ALL_pruned_data df are not of interest.
		#
		##--Let's loop through the four variables of interest above and conduct analyses along two 
		##--lines: 1) if Shapiro-Wilk normality tests suggest the data are normal, we can conduct 
		##--regression and correlation analyses using GLMs, but 2) if Shapiro-Wilk normality tests
		##--show the data are non-normal we should avoid regular correlation and instead test for
		##--correlation using the nonparametric Spearman's rank coefficient method.
		var_names <- c('propData', 'propMissing', 'posterior', 'term_brL')
		SW_pvalues <- c()
		count <- 0
		order <- c(2,3,5,7)
		for(i in 1:length(var_names)){
			SW_pvalues[i] <- shapiro.test(log(as.numeric(ALL_pruned_data[,order[i]])))$MY_PVALUE_VAR
			if (shapiro.test(log(as.numeric(ALL_pruned_data[,order[i]])))$MY_PVALUE_VAR < 0.05) {
				count <- count+1
			}
		}
		SW_pvalues
		print('Writing Shapiro-Wilk normality test results to file for four variables of interest in var_names. ', quote=FALSE)
		write.table(as.data.frame(as.matrix(cbind(var_names, SW_pvalues))), file='Shapiro-Wilk_pvalues.txt', quote=FALSE, sep='\t')
#
##--NOTE: The following section is a complex, potentially problematic set of code that could
##--cause this script to fail. CHECK AND DEBUG thoroughly.
		if ( count == 0 ) {
			## PARAMETRIC/REGULAR CORRELATION TESTS
			print('Testing for correlations between propData, propMissing, posterior, and term_brL using standard correlation tests and GLMs... ', quote=FALSE)
#			
		file.create('corr_test_output.txt')
		zz <- file('corr_test_output.txt', 'wt')	## Idea for this kind of sink (that works) came from https://stat.ethz.ch/pipermail/r-help/2008-November/178925.html/.
		sink(zz, append=TRUE)
			print('#####   GLM MODEL RESULTS   ##### ', quote=FALSE)
			print('##--Proportion data, proportion missing data, and terminal branch length as predictors of posterior support: ', quote=FALSE)
			print('##--GLM: posterior ~ propData', quote=FALSE)
			print(summary(glm(formula = log(ALL_pruned_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)))))
			print('#', quote=FALSE)
			print('##--GLM: posterior ~ propMissing', quote=FALSE)
			print(summary(glm(formula = log(ALL_pruned_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)))))
			print('#', quote=FALSE)
			print('##--GLM: posterior ~ term_brL', quote=FALSE)
			print(summary(glm(formula = log(ALL_pruned_data$MY_POSTERIOR_VAR) ~ log(ALL_pruned_data$MY_TERM_BRL_VAR))))
			print('#', quote=FALSE)
			print('##--Proportion data and missing data as predictors of terminal branch lengths: ', quote=FALSE)
			print('##--GLM: term_brL ~ propData', quote=FALSE)
			print(summary(glm(formula = log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)))))
			print('#', quote=FALSE)
			print('##--GLM: term_brL ~ propMissing', quote=FALSE)
			print(summary(glm(formula = log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)))))
			print('#', quote=FALSE)
			print('#', quote=FALSE)
			print('#####   PEARSON CORRELATION RESULTS   ##### ', quote=FALSE)
			print('##--Pearson correlation: propData, posterior', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)), log(ALL_pruned_data$MY_POSTERIOR_VAR), method='pearson'))
			print('#', quote=FALSE)
			print('##--Pearson correlation: propMissing, posterior', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)), log(ALL_pruned_data$MY_POSTERIOR_VAR), method='pearson'))
			print('#', quote=FALSE)
			print('##--Pearson correlation: term_brL, posterior', quote=FALSE)
			print(cor.test(log(ALL_pruned_data$MY_TERM_BRL_VAR), log(ALL_pruned_data$MY_POSTERIOR_VAR), method='pearson'))
			print('#', quote=FALSE)
			print('##--Pearson correlation: propData, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)), log(ALL_pruned_data$MY_TERM_BRL_VAR), method='pearson'))
			print('#', quote=FALSE)
			print('##--Pearson correlation: propMissing, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)), log(ALL_pruned_data$MY_TERM_BRL_VAR), method='pearson'))
		sink()
		close(zz)	
#
		print('PRINTING SCATTERPLOTS WITH GLM BEST-FIT LINES. ', quote=FALSE)
			graphics.off()
			pdf('GLM_scatterplots.pdf')
			plot(log(ALL_pruned_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)))
			abline(glm(formula = log(ALL_pruned_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR))))
			plot(log(ALL_pruned_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)))
			abline(glm(formula = log(ALL_pruned_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR))))
			plot(log(ALL_pruned_data$MY_POSTERIOR_VAR) ~ log(ALL_pruned_data$MY_TERM_BRL_VAR))
			abline(glm(formula = log(ALL_pruned_data$MY_POSTERIOR_VAR) ~ log(ALL_pruned_data$MY_TERM_BRL_VAR)))
			plot(log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)))
			abline(glm(formula = log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR))))
			plot(log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)))
			abline(glm(formula = log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR))))
			dev.off()
#
#			
		} else {
#
		if ( count > 0 ) {
			## NONPARAMETRIC CORRELATION TESTS
			print('Data are non-normal. Testing for correlations between propData, propMissing, posterior, and term_brL using standard correlation tests and GLMs... ')
#			
		file.create('nonparam_corr_test_output.txt')
		yy <- file('nonparam_corr_test_output.txt', 'wt')
		sink(yy, append=TRUE)
			print('#####   NONPARAMETRIC SPEARMAN CORRELATION RESULTS   ##### ', quote=FALSE)
			print('##--Spearman correlation: propData, posterior', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)), log(ALL_pruned_data$MY_POSTERIOR_VAR), method='spearman'))
			print('#', quote=FALSE)
			print('##--Spearman correlation: propMissing, posterior', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)), log(ALL_pruned_data$MY_POSTERIOR_VAR), method='spearman'))
			print('#', quote=FALSE)
			print('##--Spearman correlation: term_brL, posterior', quote=FALSE)
			print(cor.test(log(ALL_pruned_data$MY_TERM_BRL_VAR), log(ALL_pruned_data$MY_POSTERIOR_VAR), method='spearman'))
			print('#', quote=FALSE)
			print('##--Spearman correlation: propData, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)), log(ALL_pruned_data$MY_TERM_BRL_VAR), method='spearman'))
			print('#', quote=FALSE)
			print('##--Spearman correlation: propMissing, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)), log(ALL_pruned_data$MY_TERM_BRL_VAR), method='spearman'))
		sink()	
		close(yy)			
#
		print('PRINTING BASIC SCATTER PLOTS OF RELATIONSHIPS BETWEEN PROP DATA, PROP MISSING, AND TERMINAL BRANCH LENGTHS. ', quote=FALSE)
			graphics.off()
			pdf('basic_scatterplots.pdf')
			plot(log(as.numeric(ALL_pruned_data$MY_TERM_BRL_VAR)) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)))
			plot(log(as.numeric(ALL_pruned_data$MY_TERM_BRL_VAR)) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)))
			dev.off()
#
#
			}
		}			
rm(count)		
rm(order)		
#
##--Write final data table to file: 
write.table(get(ls(pattern='^ALL')), file=sprintf('%s.ALLData.txt',basename), sep='\t')
#
#
		} else {
#		
		if(is.numeric(mrbayes_tree_post_patterns) == TRUE) { print('Your tree file looks like it is from MrBayes.', quote=FALSE) 
			tree.in <- read.mrbayes.FIX(file=treeFileNames[i], digits = NULL)
			treename <- paste('tree',i,basename,sep='_')
			assign(treename, tree.in)
			treename
			tree.in		
				##--Save R workspace to file in wd before program crashes when trying to read data:
#				save.image(file='MissingDataFX_R_Workspace.RData')
#
#
		##--For the MrBayes tree file case, we need to read in the terminal branch length
		##--summary file ('*_term_brL_Summary.txt') created using mbTreeStatMiner.sh under
		##--overall analysis STEP #4 above.
		print('Reading in the tree annotation data... ', quote=FALSE)
		treedata.in <- read.table(list.files(getwd(), pattern='_term_brL_Summary.txt$'), header=TRUE, sep='\t')
		treedata_name <- paste('treedata',i,basename,sep='_')
		assign(treedata_name, treedata.in)
		treedata_name
		treedata.in
#
#
print('PRINTING UNPRUNED AND PRUNED TREE GRAPHICS. ', quote=FALSE)
		print('Producing PDFs plotting the original tree with node and tip label numbers (for reference)... ', quote=FALSE)
		pdf(file=sprintf('%s.orig.lab.pdf',treename))
		plot(get(treename), no.margin=T, label.offset=0.01, cex=0.5)
		nodelabels(cex=0.5)
		tiplabels(cex=0.5)
		dev.off()
#
#		pdf(file=sprintf('%s.orig.lab.ladder.pdf',treename))
#		plot(ladderize(get(treename), right=FALSE), no.margin=T, label.offset=0.01, cex=0.5)
#		nodelabels(cex=0.5)
#		tiplabels(cex=0.5)
#		dev.off()
#
#
		print('Pruning tree(s) using the drop file(s)... ', quote=FALSE)
		as.data.frame(get(paste('drop',i,basename,sep='_')))
		pruned_tree <- drop.tip(get(paste('tree',i,basename,sep='_')),c(paste(as.data.frame(get(paste('drop',i,basename,sep='_')))[,1])))
#
		graphics.off()
		quartz()
		print('Printing pruned tree graphics... ', quote=FALSE)
		print('Producing PDFs of simple and ladderized pruned tree plots... ', quote=FALSE)
		pdf(file=sprintf('%s.pruned.pdf',treename))
		plot(pruned_tree, no.margin=F, label.offset=0.01, cex=0.5)
		axisPhylo()
		dev.off()
#
		pdf(file=sprintf('%s.pruned.ladder.pdf',treename))
		plot(ladderize(pruned_tree, right=FALSE), no.margin=F, label.offset=0.01, cex=0.5)
		axisPhylo()
		dev.off()
#
		graphics.off()
		quartz()
		print('Producing PDFs of pruned tree plots showing node and tip label numbers (for reference)... ', quote=FALSE)
		pdf(file=sprintf('%s.pruned.lab.pdf',treename))
		plot(pruned_tree, no.margin=T, label.offset=0.01, cex=0.5)
		nodelabels(cex=0.5)
		tiplabels(cex=0.5)
		dev.off()
#
		pdf(file=sprintf('%s.pruned.lab.ladder.pdf',treename))
		plot(ladderize(pruned_tree, right=FALSE), no.margin=T, label.offset=0.01, cex=0.5)
		nodelabels(cex=0.5)
		tiplabels(cex=0.5)
		dev.off()
#
#		
##--Save R workspace to file in wd:
save.image(file='MissingDataFX_R_Workspace.RData')
print('GENERATING FINAL DATA FRAME BY MERGING TREE DATA WITH PROPORTION SUMMARY DATA. ', quote=FALSE)
		##--Now, we want to merge the 'prop.in' and 'terminal_node_data' objects/lists. First,
		##--let's give the proportion data (stored as specific to this itertation/set of files
		##--in the prop.in object) column names. Then, let's merge prop.in and the terminal
		##--node data, so that we can link all the tree_node_data variables back to the tip
		##--labels.
		colnames(prop.in) <- c('tip.label', 'propData', 'propMissing')
#		
		is.data.frame(merge(as.data.frame(prop.in), treedata.in, by='tip.label'))
		## Should produce the following output: 
		## [1] TRUE
		ALL_pruned_data <- merge(as.data.frame(prop.in), treedata.in, by='tip.label')
#
#
print('CONDUCTING EXPLORATORY AND CORRELATIONAL ANALYSES. ', quote=FALSE)
		##--We are interested in propData and propMissing as predictors of term_brL. The node
		##--column is not needed and the other possible relationships between variables in the final
		##--ALL_pruned_data df are not of interest. Unfortunately, this script cannot yet extract
		##--posteriors for MrBayes nodes and link them to tips/labels, so no tests can be conducted
		##--yet for impacts of prop data/missing data on posterior support using MrBayes trees.
		#
		##--Let's loop through the three variables of interest above and conduct analyses along two 
		##--lines: 1) if Shapiro-Wilk normality tests suggest the data are normal, we can conduct 
		##--regression and correlation analyses using GLMs, but 2) if Shapiro-Wilk normality tests
		##--show the data are non-normal we should avoid regular correlation and instead test for
		##--correlation using the nonparametric Spearman's rank coefficient method.
		var_names <- c('propData', 'propMissing', 'term_brL')
		SW_pvalues <- c()
		count <- 0
		order <- c(2,3,4)
		for(i in 1:length(var_names)){
			SW_pvalues[i] <- shapiro.test(log(as.numeric(ALL_pruned_data[,order[i]])))$MY_PVALUE_VAR
			if (shapiro.test(log(as.numeric(ALL_pruned_data[,order[i]])))$MY_PVALUE_VAR < 0.05) {
				count <- count+1
			}
		}
		SW_pvalues
		print('Writing Shapiro-Wilk normality test results to file for four variables of interest in var_names. ', quote=FALSE)
		write.table(as.data.frame(as.matrix(cbind(var_names, SW_pvalues))), file='Shapiro-Wilk_pvalues.txt', quote=FALSE, sep='\t')
#
##--NOTE: The following section is a complex, potentially problematic set of code that could
##--cause this script to fail. CHECK AND DEBUG thoroughly.
		if ( count == 0 ) {
			## PARAMETRIC/REGULAR CORRELATION TESTS
			print('Testing for correlations between propData, propMissing, and term_brL using standard correlation tests and GLMs... ', quote=FALSE)
#			
		file.create('corr_test_output.txt')
		zz <- file('corr_test_output.txt', 'wt')	## Idea for this kind of sink (that works) came from https://stat.ethz.ch/pipermail/r-help/2008-November/178925.html/.
		sink(zz, append=TRUE)
			print('#####   GLM MODEL RESULTS   ##### ', quote=FALSE)
			print('##--Proportion data and missing data as predictors of terminal branch lengths: ', quote=FALSE)
			print('##--GLM: term_brL ~ propData', quote=FALSE)
			print(summary(glm(formula = log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)))))
			print('#', quote=FALSE)
			print('##--GLM: term_brL ~ propMissing', quote=FALSE)
			print(summary(glm(formula = log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)))))
			print('#', quote=FALSE)
			print('#', quote=FALSE)
			print('#####   PEARSON CORRELATION RESULTS   ##### ', quote=FALSE)
			print('##--Pearson correlation: propData, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)), log(ALL_pruned_data$MY_TERM_BRL_VAR), method='pearson'))
			print('#', quote=FALSE)
			print('##--Pearson correlation: propMissing, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)), log(ALL_pruned_data$MY_TERM_BRL_VAR), method='pearson'))
		sink()
		close(zz)	
#
		print('PRINTING SCATTERPLOTS WITH GLM BEST-FIT LINES. ', quote=FALSE)
			graphics.off()
			pdf('GLM_scatterplots.pdf')
			plot(log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)))
			abline(glm(formula = log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR))))
			plot(log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)))
			abline(glm(formula = log(ALL_pruned_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR))))
			dev.off()
#			
#			
		} else {
#
		if ( count > 0 ) {
			## NONPARAMETRIC CORRELATION TESTS
			print('Data are non-normal. Testing for correlations between propData, propMissing, and term_brL using standard correlation tests and GLMs... ')
#			
		file.create('nonparam_corr_test_output.txt')
		yy <- file('nonparam_corr_test_output.txt', 'wt')
		sink(yy, append=TRUE)
			print('#####   NONPARAMETRIC SPEARMAN CORRELATION RESULTS   ##### ', quote=FALSE)
			print('##--Spearman correlation: propData, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)), log(ALL_pruned_data$MY_TERM_BRL_VAR), method='spearman'))
			print('#', quote=FALSE)
			print('##--Spearman correlation: propMissing, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)), log(ALL_pruned_data$MY_TERM_BRL_VAR), method='spearman'))
		sink()	
		close(yy)				
#
		print('PRINTING BASIC SCATTER PLOTS OF RELATIONSHIPS BETWEEN PROP DATA, PROP MISSING, AND TERMINAL BRANCH LENGTHS. ', quote=FALSE)
			graphics.off()
			pdf('basic_scatterplots.pdf')
			plot(log(as.numeric(ALL_pruned_data$MY_TERM_BRL_VAR)) ~ log(as.numeric(ALL_pruned_data$MY_PROP_DATA_VAR)))
			plot(log(as.numeric(ALL_pruned_data$MY_TERM_BRL_VAR)) ~ log(as.numeric(ALL_pruned_data$MY_PROP_MISSING_VAR)))
			dev.off()
#
#
			}
		}			
rm(count)		
rm(order)		
#
##--Write final data table to file: 
write.table(get(ls(pattern='^ALL')), file=sprintf('%s.ALLData.txt',basename), sep='\t')
#
#
		}
}



	}
}



############ IIIB. IF *NO* DROP FILE IS PRESENT, READ IN DATA AND TREE FILE(S), PLOT TREE(S),
############ AND CONDUCT ANALYSES FOR THE REGULAR CASE.
##--Loop through the .tree filenames to read each tree in and assign it to an appropriately 
##--named object, then conduct analyses with the tree and data extracted from it and save
##--the results to file:
if (length(dropFileNames)=='0') {
date()
print('CONDUCTING ANALYSES FOR THE REGULAR CASE. ', quote=FALSE)
	for(i in 1:length(treeFileNames)){
		print('Reading in basename and data/files...', quote=FALSE)
		basename <- file_path_sans_ext(treeFileNames[i])
#
		print('Reading in the proportion summary data... ', quote=FALSE)
		##--Loop through the proportion summary filenames to read each file in and assign it to an
		##--object:
		prop.in <- as.matrix(read.table(file=propSumFileNames[i], header=FALSE, sep=' '))
		assign(paste('prop',i,basename,sep='_'), prop.in)
#
		print('Reading in tree file(s)... ', quote=FALSE)
		##--Loop through the .tree filenames to read each tree in and assign it to an appropriately 
		##--named object. However, we will do this using different ips functions depending whether the
		##--tree file is from BEAST or MrBayes ('read.beast' or 'read.mrbayes', respectively). So, we
		##--first need to test this. We'll grep patterns that are exclusive to the format of BEAST vs.
		##--MrBayes consensus NEXUS tree files, then use conditionals to guide how we fill 'tree.in'
		##--var below.
		beast_tree_post_patterns <- sapply(list.files(getwd(), pattern=treeFileNames[i]), FUN=function(x){grep('posterior', readLines(x))})
		mrbayes_tree_post_patterns <- sapply(list.files(getwd(), pattern=treeFileNames[i]), FUN=function(x){grep('prob', readLines(x))})
#
		if(is.numeric(beast_tree_post_patterns) == TRUE) { print('Your tree file looks like it is from BEAST.', quote=FALSE) 
			tree.in <- read.beast(file=treeFileNames[i], digits = NULL)
			treename <- paste('tree',i,basename,sep='_')
			assign(treename, tree.in)
			treename
			tree.in
#
		print('Reading in the tree annotation data... ', quote=FALSE)
		treedata.in <- read.beast.table(file=treeFileNames[i], digits = 4)
		treedata_name <- paste('treedata',i,basename,sep='_')
		assign(treedata_name, treedata.in)
		treedata_name
		treedata.in
#
#
print('PRINTING REGULAR TREE GRAPHICS. ')
		print('Producing PDFs of simple tree plot... ', quote=FALSE)
		pdf(file=sprintf('%s.pdf',treename))
		plot(get(treename), no.margin=F, label.offset=0.1, cex=0.5)
		axisPhylo()
		dev.off()
#
		pdf(file=sprintf('%s.ladder.pdf',treename))
		plot(ladderize(get(treename), right=FALSE), no.margin=F, label.offset=0.1, cex=0.5)
		axisPhylo()
		dev.off()
#
		graphics.off()
		quartz()
		print('Producing PDFs of tree plots showing node and tip label numbers (for reference)... ')
		pdf(file=sprintf('%s.lab.pdf',treename))
		plot(get(treename), no.margin=T, label.offset=0.1, cex=0.5)
		nodelabels(cex=0.5)
		tiplabels(cex=0.5)
		dev.off()
#
		pdf(file=sprintf('%s.lab.ladder.pdf',treename))
		plot(ladderize(get(treename), right=FALSE), no.margin=T, label.offset=0.1, cex=0.5)
		nodelabels(cex=0.5)
		tiplabels(cex=0.5)
		dev.off()
#
#
print('GETTING NODE NUMBERS FOR EACH TIP TAXON USING ORIGINAL, UNPRUNED TREE. ', quote=FALSE)
#
			##--Create function for picking the terminal node number for a given tip taxon. I made this
			##--function, 'terminalNode', by modifying the version of 'fastMRCA' function listed on
			##--this Liam Revell blog post (URL: http://blog.phytools.org/2012/05/function-to-efficiently-return-mrca-of.html/; 
			##--which is *NOT* identical to that of the current phytools distro, Feb 2017). The 'Ancestors'
			##--function is from phangorn, so phytools and phangorn should be loaded at start of script. 
			#
			##--NOTE: If desired, you may want to use the 'fastMRCA' function of phytools to get the
			##--number of the corresponding tip label (=shown in yellow along tips in labeled-tree PDF
			##--plots above, i.e. from calling nodelabels()). fastMRCA just matches the taxon name it
			##--is given to the corresponding tip.label.
			terminalNode <- function(tree, species) {
			    x <- match(species, tree$MY_TIP_LABEL_VAR)
			    y <- match(species, tree$MY_TIP_LABEL_VAR)
			    a <- Ancestors(tree, x)
			    b <- Ancestors(tree, y)
			    z <- a%in%b
			    return(a[min(which(z))])
			}
#
		##--Now that we have the 'terminalNode' function, we can use it to loop through all taxon
		##--tip.labels and save the corresponding terminal node numbers into 'terminal_node_numbers'. 
		terminal_node_numbers <- c()
		for(i in 1:length(get(treename)$MY_TIP_LABEL_VAR)){
			terminal_node_numbers[i] <- terminalNode(get(treename), get(treename)$MY_TIP_LABEL_VAR[i])
		}
#
print('COLLECTING TREE DATA FOR TERMINAL NODES. ', quote=FALSE)
		##--Next, use for loop to 1) match terminal_node_numbers to node column in treedata data frame,
		##--in order to get row numbers in treedata dataframe corresponding to each taxon; 2) get the
		##--a) posterior probability, b) terminal branch length, and c) estimated node height for the
		##--corresponding terminal node/tip; and 3) save all of this data into a list, with appropriate
		##--names for list elements.
		terminal_node_labels <- c()
		terminal_node_row_numbers <- c()
		terminal_node_posteriors <- c()
		terminal_node_lengths <- c()
		terminal_node_heights <- c()
		for(i in 1:length(get(treename)$MY_TIP_LABEL_VAR)){
			terminal_node_labels[i] <- c(get(treename)$MY_TIP_LABEL_VAR[i])
			terminal_node_row_numbers[i] <- match(terminal_node_numbers[i], as.data.frame(get(treedata_name))$MY_NODE_VAR)
			terminal_node_posteriors[i] <- as.data.frame(get(treedata_name))[c(terminal_node_row_numbers), ][i,]$MY_POSTERIOR_VAR
			terminal_node_lengths[i] <- as.data.frame(get(treedata_name))[c(terminal_node_row_numbers), ][i,]$MY_LENGTH_VAR
			terminal_node_heights[i] <- as.data.frame(get(treedata_name))[c(terminal_node_row_numbers), ][i,]$MY_HEIGHT_VAR
		}
#		
		terminal_node_data <- list(terminal_node_labels, terminal_node_numbers, terminal_node_posteriors, terminal_node_lengths, terminal_node_heights)
		names(terminal_node_data) <- c('tip.label', 'node', 'posterior', 'pre_brL', 'term_brL')
#			
##--Save R workspace to file in wd:
save.image(file='MissingDataFX_R_Workspace.RData')
print('GENERATING FINAL DATA FRAME BY MERGING TREE DATA WITH PROPORTION SUMMARY DATA. ', quote=FALSE)
		##--Now, we want to merge the 'prop.in' and 'terminal_node_data' objects/lists. First,
		##--let's give the proportion data (stored as specific to this itertation/set of files
		##--in the prop.in object) column names. Then, let's merge prop.in and the terminal
		##--node data, so that we can link all the following variables back to the tip labels:
		##--a) propData, b) propMissing, c) node, d) posterior, e) pre_brL, and f) term_brL.
		colnames(prop.in) <- c('tip.label', 'propData', 'propMissing')
#		
		##--Ignore next six lines for now, which reflect something I wanted to try (they do no harm, but I'm not really using them yet).
		terminal_node_data_labels <- c()
		terminal_node_data_row_numbers <- c()
		for(i in 1:length(get(treename)$MY_TIP_LABEL_VAR)){
			terminal_node_data_labels[i] <- c(terminal_node_data$MY_TIP_LABEL_VAR[i])
			terminal_node_data_row_numbers[i] <- match(terminal_node_data_labels[i], prop.in[i])
		}
#
		is.data.frame(merge(as.data.frame(prop.in), terminal_node_data, by='tip.label'))		## Merge idea was inspired by code I'd seen in phylogenetics packages, and this post (https://docs.tibco.com/pub/enterprise-runtime-for-R/1.5.0_may_2013/TERR_1.5.0_LanguageRef/base/merge.html).
		## Should produce the following output: 
		## [1] TRUE
		ALL_data <- merge(as.data.frame(prop.in), terminal_node_data, by='tip.label')
#		
#		
print('CONDUCTING EXPLORATORY AND CORRELATIONAL ANALYSES. ', quote=FALSE)
		##--We are interested in propData, propMissing, and term_brL as predictors of posterior; and
		##--we are also interested in propData and propMissing as predictors of term_brL. The node
		##--column is not needed and the other possible relationships between variables in the final
		##--ALL_data df are not of interest.
		#
		##--Let's loop through the four variables of interest above and conduct analyses along two 
		##--lines: 1) if Shapiro-Wilk normality tests suggest the data are normal, we can conduct 
		##--regression and correlation analyses using GLMs, but 2) if Shapiro-Wilk normality tests
		##--show the data are non-normal we should avoid regular correlation and instead test for
		##--correlation using the nonparametric Spearman's rank coefficient method.
		var_names <- c('propData', 'propMissing', 'posterior', 'term_brL')
		SW_pvalues <- c()
		count <- 0
		order <- c(2,3,5,7)
		for(i in 1:length(var_names)){
			SW_pvalues[i] <- shapiro.test(log(as.numeric(ALL_data[,order[i]])))$MY_PVALUE_VAR
			if (shapiro.test(log(as.numeric(ALL_data[,order[i]])))$MY_PVALUE_VAR < 0.05) {
				count <- count+1
			}
		}
		SW_pvalues
		print('Writing Shapiro-Wilk normality test results to file for four variables of interest in var_names. ')
		write.table(as.data.frame(as.matrix(cbind(var_names, SW_pvalues))), file='Shapiro-Wilk_pvalues.txt', quote=FALSE, sep='\t')
#
##--NOTE: The following section is a complex, potentially problematic set of code that could
##--cause this script to fail. CHECK AND DEBUG thoroughly.
		if ( count == 0 ) {
			## PARAMETRIC/REGULAR CORRELATION TESTS
			print('Testing for correlations between propData, propMissing, posterior, and term_brL using standard correlation tests and GLMs... ', quote=FALSE)
#			
		file.create('corr_test_output.txt')
		zz <- file('corr_test_output.txt', 'wt')
		sink(zz, append=TRUE)
			print('#####   GLM MODEL RESULTS   ##### ', quote=FALSE)
			print('##--Proportion data, proportion missing data, and terminal branch length as predictors of posterior support: ', quote=FALSE)
			print('##--GLM: posterior ~ propData', quote=FALSE)
			print(summary(glm(formula = log(ALL_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR)))))
			print('#', quote=FALSE)
			print('##--GLM: posterior ~ propMissing', quote=FALSE)
			print(summary(glm(formula = log(ALL_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)))))
			print('#', quote=FALSE)
			print('##--GLM: posterior ~ term_brL', quote=FALSE)
			print(summary(glm(formula = log(ALL_data$MY_POSTERIOR_VAR) ~ log(ALL_data$MY_TERM_BRL_VAR))))
			print('#', quote=FALSE)
			print('##--Proportion data and missing data as predictors of terminal branch lengths: ', quote=FALSE)
			print('##--GLM: term_brL ~ propData', quote=FALSE)
			print(summary(glm(formula = log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR)))))
			print('#', quote=FALSE)
			print('##--GLM: term_brL ~ propMissing', quote=FALSE)
			print(summary(glm(formula = log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)))))
			print('#', quote=FALSE)
			print('#', quote=FALSE)
			print('#####   PEARSON CORRELATION RESULTS   ##### ', quote=FALSE)
			print('##--Pearson correlation: propData, posterior', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_DATA_VAR)), log(ALL_data$MY_POSTERIOR_VAR), method='pearson'))
			print('#', quote=FALSE)
			print('##--Pearson correlation: propMissing, posterior', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)), log(ALL_data$MY_POSTERIOR_VAR), method='pearson'))
			print('#', quote=FALSE)
			print('##--Pearson correlation: term_brL, posterior', quote=FALSE)
			print(cor.test(log(ALL_data$MY_TERM_BRL_VAR), log(ALL_data$MY_POSTERIOR_VAR), method='pearson'))
			print('#', quote=FALSE)
			print('##--Pearson correlation: propData, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_DATA_VAR)), log(ALL_data$MY_TERM_BRL_VAR), method='pearson'))
			print('#', quote=FALSE)
			print('##--Pearson correlation: propMissing, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)), log(ALL_data$MY_TERM_BRL_VAR), method='pearson'))
		sink()
		close(zz)	
#
		print('PRINTING SCATTERPLOTS WITH GLM BEST-FIT LINES. ', quote=FALSE)
			graphics.off()
			pdf('GLM_scatterplots.pdf')
			plot(log(ALL_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR)))
			abline(glm(formula = log(ALL_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR))))
			plot(log(ALL_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)))
			abline(glm(formula = log(ALL_data$MY_POSTERIOR_VAR) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR))))
			plot(log(ALL_data$MY_POSTERIOR_VAR) ~ log(ALL_data$MY_TERM_BRL_VAR))
			abline(glm(formula = log(ALL_data$MY_POSTERIOR_VAR) ~ log(ALL_data$MY_TERM_BRL_VAR)))
			plot(log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR)))
			abline(glm(formula = log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR))))
			plot(log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)))
			abline(glm(formula = log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR))))
			dev.off()
#
#			
		} else {
#
		if ( count > 0 ) {
			## NONPARAMETRIC CORRELATION TESTS
			print('Data are non-normal. Testing for correlations between propData, propMissing, posterior, and term_brL using standard correlation tests and GLMs... ')
#			
		file.create('nonparam_corr_test_output.txt')
		yy <- file('nonparam_corr_test_output.txt', 'wt')
		sink(yy, append=TRUE)
			print('#####   NONPARAMETRIC SPEARMAN CORRELATION RESULTS   ##### ', quote=FALSE)
			print('##--Spearman correlation: propData, posterior', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_DATA_VAR)), log(ALL_data$MY_POSTERIOR_VAR), method='spearman'))
			print('#', quote=FALSE)
			print('##--Spearman correlation: propMissing, posterior', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)), log(ALL_data$MY_POSTERIOR_VAR), method='spearman'))
			print('#', quote=FALSE)
			print('##--Spearman correlation: term_brL, posterior', quote=FALSE)
			print(cor.test(log(ALL_data$MY_TERM_BRL_VAR), log(ALL_data$MY_POSTERIOR_VAR), method='spearman'))
			print('#', quote=FALSE)
			print('##--Spearman correlation: propData, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_DATA_VAR)), log(ALL_data$MY_TERM_BRL_VAR), method='spearman'))
			print('#', quote=FALSE)
			print('##--Spearman correlation: propMissing, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)), log(ALL_data$MY_TERM_BRL_VAR), method='spearman'))
		sink()	
		close(yy)	
#			
		print('PRINTING BASIC SCATTER PLOTS OF RELATIONSHIPS BETWEEN PROP DATA, PROP MISSING, AND TERMINAL BRANCH LENGTHS. ', quote=FALSE)
			graphics.off()
			pdf('basic_scatterplots.pdf')
			plot(log(as.numeric(ALL_data$MY_TERM_BRL_VAR)) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR)))
			plot(log(as.numeric(ALL_data$MY_TERM_BRL_VAR)) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)))
			dev.off()
#
#
			}
		}			
rm(count)		
rm(order)
#
##--Write final data table to file: 
write.table(get(ls(pattern='^ALL')), file=sprintf('%s.ALLData.txt',basename), sep='\t')
#
#
		} else {
#		
		if(is.numeric(mrbayes_tree_post_patterns) == TRUE) { print('Your tree file looks like it is from MrBayes.', quote=FALSE) 
			tree.in <- read.mrbayes.FIX(file=treeFileNames[i], digits = NULL)
			treename <- paste('tree',i,basename,sep='_')
			assign(treename, tree.in)
			treename
			tree.in	
				##--Save R workspace to file in wd before program crashes when trying to read data:
#				save.image(file='MissingDataFX_R_Workspace.RData')
#
#
		##--For the MrBayes tree file case, we need to read in the terminal branch length
		##--summary file ('*_term_brL_Summary.txt') created using mbTreeStatMiner.sh under
		##--overall analysis STEP #4 above.
		print('Reading in the tree annotation data... ', quote=FALSE)
		treedata.in <- read.table(list.files(getwd(), pattern='_term_brL_Summary.txt$'), header=TRUE, sep='\t')
		treedata_name <- paste('treedata',i,basename,sep='_')
		assign(treedata_name, treedata.in)
		treedata_name
		treedata.in
#
#
print('PRINTING REGULAR TREE GRAPHICS. ')
		print('Producing PDFs of simple tree plot... ', quote=FALSE)
		pdf(file=sprintf('%s.pdf',treename))
		plot(get(treename), no.margin=F, label.offset=0.01, cex=0.5)
		axisPhylo()
		dev.off()
#
		pdf(file=sprintf('%s.ladder.pdf',treename))
		plot(ladderize(get(treename), right=FALSE), no.margin=F, label.offset=0.01, cex=0.5)
		axisPhylo()
		dev.off()
#
		graphics.off()
		quartz()
		print('Producing PDFs of tree plots showing node and tip label numbers (for reference)... ')
		pdf(file=sprintf('%s.lab.pdf',treename))
		plot(get(treename), no.margin=T, label.offset=0.01, cex=0.5)
		nodelabels(cex=0.5)
		tiplabels(cex=0.5)
		dev.off()
#
		pdf(file=sprintf('%s.lab.ladder.pdf',treename))
		plot(ladderize(get(treename), right=FALSE), no.margin=T, label.offset=0.01, cex=0.5)
		nodelabels(cex=0.5)
		tiplabels(cex=0.5)
		dev.off()
#
#			
##--Save R workspace to file in wd:
save.image(file='MissingDataFX_R_Workspace.RData')
print('GENERATING FINAL DATA FRAME BY MERGING TREE DATA WITH PROPORTION SUMMARY DATA. ', quote=FALSE)
		##--Now, we want to merge the 'prop.in' and 'terminal_node_data' objects/lists. First,
		##--let's give the proportion data (stored as specific to this itertation/set of files
		##--in the prop.in object) column names. Then, let's merge prop.in and the terminal
		##--node data, so that we can link all the tree_node_data variables back to the tip
		##--labels.
		colnames(prop.in) <- c('tip.label', 'propData', 'propMissing')
#		
		is.data.frame(merge(as.data.frame(prop.in), treedata.in, by='tip.label'))
		## Should produce the following output: 
		## [1] TRUE
		ALL_data <- merge(as.data.frame(prop.in), treedata.in, by='tip.label')
#
#
print('CONDUCTING EXPLORATORY AND CORRELATIONAL ANALYSES. ', quote=FALSE)
		##--We are interested in propData and propMissing as predictors of term_brL. The node
		##--column is not needed and the other possible relationships between variables in the final
		##--ALL_data df are not of interest. Unfortunately, this script cannot yet extract
		##--posteriors for MrBayes nodes and link them to tips/labels, so no tests can be conducted
		##--yet for impacts of prop data/missing data on posterior support using MrBayes trees.
		#
		##--Let's loop through the three variables of interest above and conduct analyses along two 
		##--lines: 1) if Shapiro-Wilk normality tests suggest the data are normal, we can conduct 
		##--regression and correlation analyses using GLMs, but 2) if Shapiro-Wilk normality tests
		##--show the data are non-normal we should avoid regular correlation and instead test for
		##--correlation using the nonparametric Spearman's rank coefficient method.
		var_names <- c('propData', 'propMissing', 'term_brL')
		SW_pvalues <- c()
		count <- 0
		order <- c(2,3,4)
		for(i in 1:length(var_names)){
			SW_pvalues[i] <- shapiro.test(log(as.numeric(ALL_data[,order[i]])))$MY_PVALUE_VAR
			if (shapiro.test(log(as.numeric(ALL_data[,order[i]])))$MY_PVALUE_VAR < 0.05) {
				count <- count+1
			}
		}
		SW_pvalues
		print('Writing Shapiro-Wilk normality test results to file for four variables of interest in var_names. ', quote=FALSE)
		write.table(as.data.frame(as.matrix(cbind(var_names, SW_pvalues))), file='Shapiro-Wilk_pvalues.txt', quote=FALSE, sep='\t')
#
##--NOTE: The following section is a complex, potentially problematic set of code that could
##--cause this script to fail. CHECK AND DEBUG thoroughly.
		if ( count == 0 ) {
			## PARAMETRIC/REGULAR CORRELATION TESTS
			print('Testing for correlations between propData, propMissing, and term_brL using standard correlation tests and GLMs... ', quote=FALSE)
#			
		file.create('corr_test_output.txt')
		zz <- file('corr_test_output.txt', 'wt')	## Idea for this kind of sink (that works) came from https://stat.ethz.ch/pipermail/r-help/2008-November/178925.html/.
		sink(zz, append=TRUE)
			print('#####   GLM MODEL RESULTS   ##### ', quote=FALSE)
			print('##--Proportion data and missing data as predictors of terminal branch lengths: ', quote=FALSE)
			print('##--GLM: term_brL ~ propData', quote=FALSE)
			print(summary(glm(formula = log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR)))))
			print('#', quote=FALSE)
			print('##--GLM: term_brL ~ propMissing', quote=FALSE)
			print(summary(glm(formula = log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)))))
			print('#', quote=FALSE)
			print('#', quote=FALSE)
			print('#####   PEARSON CORRELATION RESULTS   ##### ', quote=FALSE)
			print('##--Pearson correlation: propData, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_DATA_VAR)), log(ALL_data$MY_TERM_BRL_VAR), method='pearson'))
			print('#', quote=FALSE)
			print('##--Pearson correlation: propMissing, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)), log(ALL_data$MY_TERM_BRL_VAR), method='pearson'))
		sink()
		close(zz)	
#			
		print('PRINTING SCATTERPLOTS WITH GLM BEST-FIT LINES. ', quote=FALSE)
			graphics.off()
			pdf('GLM_scatterplots.pdf')
			plot(log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR)))
			abline(glm(formula = log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR))))
			plot(log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)))
			abline(glm(formula = log(ALL_data$MY_TERM_BRL_VAR) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR))))
			dev.off()
#
#			
		} else {
#
		if ( count > 0 ) {
			## NONPARAMETRIC CORRELATION TESTS
			print('Data are non-normal. Testing for correlations between propData, propMissing, and term_brL using standard correlation tests and GLMs... ')
#			
		file.create('nonparam_corr_test_output.txt')
		yy <- file('nonparam_corr_test_output.txt', 'wt')
		sink(yy, append=TRUE)
			print('#####   NONPARAMETRIC SPEARMAN CORRELATION RESULTS   ##### ', quote=FALSE)
			print('##--Spearman correlation: propData, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_DATA_VAR)), log(ALL_data$MY_TERM_BRL_VAR), method='spearman'))
			print('#', quote=FALSE)
			print('##--Spearman correlation: propMissing, term_brL', quote=FALSE)
			print(cor.test(log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)), log(ALL_data$MY_TERM_BRL_VAR), method='spearman'))
		sink()	
		close(yy)	
#			
		print('PRINTING BASIC SCATTER PLOTS OF RELATIONSHIPS BETWEEN PROP DATA, PROP MISSING, AND TERMINAL BRANCH LENGTHS. ', quote=FALSE)
			graphics.off()
			pdf('basic_scatterplots.pdf')
			plot(log(as.numeric(ALL_data$MY_TERM_BRL_VAR)) ~ log(as.numeric(ALL_data$MY_PROP_DATA_VAR)))
			plot(log(as.numeric(ALL_data$MY_TERM_BRL_VAR)) ~ log(as.numeric(ALL_data$MY_PROP_MISSING_VAR)))
			dev.off()
#
#
			}
		}			
rm(count)		
rm(order)		
#
##--Write final data table to file: 
write.table(get(ls(pattern='^ALL')), file=sprintf('%s.ALLData.txt',basename), sep='\t')
#
#
		}
}




	}
}


##--Save R workspace to file in wd:
save.image(file='MissingDataFX_R_Workspace.RData')

######################################### END ############################################
" > missingDataFXTester.r


############ IV. FINAL STEPS:
echo "INFO      | $(date) | STEP #6: RUN THE R SCRIPT (WHICH ALSO SAVES RESULTS TO FILE). "
R CMD BATCH ./missingDataFXTester.R


echo "INFO      | $(date) | STEP #7: CLEANUP: ORGANIZE RESULTS, REMOVE UNNECESSARY FILES. "
###### Make dir and organize NEXUS data into "NEXUS_data" sub-folder:
	if [[ ! -n $(find . -name "NEXUS_data" -type d) ]]; then
		mkdir NEXUS_data
		mv ./*_BIGSummary.txt ./*_propData.txt ./*_propMissing.txt ./*_propSummary.txt ./*_taxa.txt ./NEXUS_data/
	else
		mv ./*_BIGSummary.txt ./*_propData.txt ./*_propMissing.txt ./*_propSummary.txt ./*_taxa.txt ./NEXUS_data/
	fi
	
###### mbTreeStatMiner results sub-folder:
	if [[ ! -n $(find . -name "mb_tree_stats" -type d) ]]; then
		mkdir mb_tree_stats
		mv ./*_brL_Summary.txt ./mb_tree_stats/
	else
		mv ./*_brL_Summary.txt ./mb_tree_stats/
	fi
	
###### final data set sub-folder:
	if [[ ! -n $(find . -name "final_dataset" -type d) ]]; then
		mkdir final_dataset
		mv ./*ALLData.txt ./final_dataset/
	else
		mv ./*ALLData.txt ./final_dataset/
	fi
	
###### Make dir and organize R results into "R_results" sub-folder:
	if [[ ! -n $(find . -name "R_results" -type d) ]]; then
		mkdir R_results
		if [ -s ./MissingDataFX_R_Workspace.RData ]; then
			echo "INFO      | $(date) |          R workspace file confirmed in dir: $MY_PATH "
		elif [ ! -s ./MissingDataFX_R_Workspace.RData ]; then
			echo "WARNING!  | $(date) |          Could not find R workspace file in dir: $MY_PATH "
		fi
		mv ./missingDataFXTester.r ./missingDataFXTester.Rout ./MissingDataFX_R_Workspace.RData ./*.pdf ./R_results/
		mv ./*_test_output.txt ./*_pvalues.txt ./R_results/
	else
		if [ -s ./MissingDataFX_R_Workspace.RData ]; then
			echo "INFO      | $(date) |          R workspace file confirmed in dir: $MY_PATH "
		elif [ ! -s ./MissingDataFX_R_Workspace.RData ]; then
			echo "WARNING!  | $(date) |          Could not find R workspace file in dir: $MY_PATH "
		fi
		mv ./missingDataFXTester.r ./missingDataFXTester.Rout ./MissingDataFX_R_Workspace.RData ./*.pdf ./R_results/
		mv ./*_test_output.txt ./*_pvalues.txt ./R_results/
	fi
	
###### Clean up any temporary files remaining in working directory:
	if [[ -n $(find . -name "*headless.nex" -type f) ]]; then rm ./*headless.nex; fi
	if [[ -n $(find . -name "*_concatenatedSeqs.txt" -type f) ]]; then rm ./*_concatenatedSeqs.txt; fi
	if [[ -n $(find . -name "*_totalChar.txt" -type f) ]]; then rm ./*_totalChar.txt; fi
	if [[ -n $(find . -name "*_totalChar.tmp" -type f) ]]; then rm ./*_totalChar.tmp; fi
	if [[ -n $(find . -name "*_regCharFILE.txt" -type f) ]]; then rm ./*_regCharFILE.txt; fi
	if [[ -n $(find . -name "*_regCharCOUNTS.txt" -type f) ]]; then rm ./*_regCharCOUNTS.txt; fi
	if [[ -n $(find . -name "*_gapChar.txt" -type f) ]]; then rm ./*_gapChar.txt; fi
	if [[ -n $(find . -name "*_missingChar.txt" -type f) ]]; then rm ./*_missingChar.txt; fi
		

echo "INFO      | $(date) | STEP #8: TEST RUN OUTPUT USING MDFXTester, WHICH WILL RERUN MissingDataFX IF RESULTS ARE INCOMPLETE. "
	chmod u+x ../shell/MDFXTester.sh
	../shell/MDFXTester.sh


echo "INFO      | $(date) | Done analyzing the amount and potential effects of missing data on phylogenetic support and branch "
echo "INFO      | $(date) | lengths using MissingDataFX. "
echo "INFO      | $(date) | Bye.
"
#
#
#
######################################### END ############################################

exit 0
