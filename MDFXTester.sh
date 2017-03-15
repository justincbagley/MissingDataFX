#!/bin/sh

######################################## START ###########################################
echo "INFO      | $(date) | Starting MDFXTester analysis... "
echo "INFO      | $(date) | STEP #1: DO GRAPHICS CHECK ON BASIC SCATTERPLOTS FROM R; IF FAILED, DO fetchTaxonLabels.sh FIX "
echo "INFO      | $(date) |          AND COMPLETELY RERUN MissingDataFX SCRIPT. "

if [[ -n $(find . -name "basic_scatterplots.pdf" -type f) ]] || [[ -n $(find ../R_results/ -name "basic_scatterplots.pdf" -type f) ]]; then
	echo "INFO      | $(date) |          Passed graphics check. Moving on... "
else
	echo "WARNING!  | $(date) |          FAILED graphics check. Running fetchTaxonLabels script, then re-running MissingDataFX.sh... "
	chmod u+x ../shell/fetchTaxonLabels.sh
	../shell/fetchTaxonLabels.sh
	#
	rm ./missingDataFXTester.Rout ./MissingDataFX_R_Workspace.RData ./*.pdf ./*_test_output.txt ./*_pvalues.txt
	#
	./MissingDataFX.sh
fi


#
#
#
######################################### END ############################################

exit 0
