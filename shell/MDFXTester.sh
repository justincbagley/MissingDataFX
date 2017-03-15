#!/bin/sh

######################################## START ###########################################

echo "
### RUNNING MDFXTester.sh ###
"
echo "DO GRAPHICS CHECK ON BASIC SCATTERPLOTS FROM R; IF FAILED, DO fetchTaxonLabels.sh FIX "
echo "AND COMPLETELY RERUN MissingDataFX SCRIPT. "

if [[ -n $(find . -name "basic_scatterplots.pdf" -type f) ]] || [[ -n $(find ../R_results/ -name "basic_scatterplots.pdf" -type f) ]]; then
	echo "     Passed graphics check. Moving on... "
else
	echo "     WARNING! FAILED graphics check. Running fetchTaxonLabels script, then re-running MissingDataFX.sh... "
	chmod u+x ../shell/fetchTaxonLabels.sh
	../shell/fetchTaxonLabels.sh

	rm ./missingDataFXTester.Rout ./MissingDataFX_R_Workspace.RData ./*.pdf ./*_test_output.txt ./*_pvalues.txt


echo "
### RE-RUNNING MissingDataFX.sh ###
"

./MissingDataFX.sh

fi


#
#
#
######################################### END ############################################

exit 0
