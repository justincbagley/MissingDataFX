#!/bin/sh

######################################## START ###########################################

echo "
### RUNNING MDFXTester.sh ###
"
echo "DO GRAPHICS CHECK ON BASIC SCATTERPLOTS FROM R; IF FAILED, DO fetchTaxonLabels.sh FIX "
echo "AND COMPLETELY RERUN MissingDataFX SCRIPT. "

if [[ -n $(find . -name "basic_scatterplots.pdf" -type f) ]] || [[ -n $(find ../R_results/ -name "basic_scatterplots.pdf" -type f) ]]; then
	echo "     Passed graphics check. Moving on... 
"
else
	echo "     WARNING! FAILED graphics check. Running fetchTaxonLabels script, then re-running MissingDataFX.sh... "
	chmod u+x ../shell/fetchTaxonLabels.sh
	../shell/fetchTaxonLabels.sh

	if [[ -n $(find . -name "missingDataFXTester.Rout" -type f) ]]; then rm ./missingDataFXTester.Rout; fi
	if [[ -n $(find . -name "MissingDataFX_R_Workspace.RData" -type f) ]]; then rm ./MissingDataFX_R_Workspace.RData; fi
	if [[ -n $(find . -name "*.pdf" -type f) ]]; then rm ./*.pdf; fi
	if [[ -n $(find . -name "*_test_output.txt" -type f) ]]; then rm ./*_test_output.txt; fi 
	if [[ -n $(find . -name "*_pvalues.txt" -type f) ]]; then rm ./*_pvalues.txt; fi


echo "
### RE-RUNNING MissingDataFX.sh!!! ###
"

./MissingDataFX.sh

fi


#
#
#
######################################### END ############################################

exit 0
