#!/bin/bash
#
# Sophie COUTANT
# 28/11/2014
#
#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
# Ce script permet l'automatisation du lancement des rapports CANOES                                                                                        #
#-----------------------------------------------------------------------------------------------------------------------------------------------------------#

# usage
function usage
{
	echo -e "#------------------------------------------------------------------------------------------------------------------------------------------------------------#"
	echo -e "# Ce script permet l'automatisation du lancement des rapports CANOES                                                                                         #"
	echo -e "# 6. Lancement du calling des CNV. Script R CANOES.                                                                                                          #"
	echo -e "#------------------------------------------------------------------------------------------------------------------------------------------------------------#"
    echo -e "\nUSAGE: run_makeReport.sh -o <directory> -log <file> -op <prefix>"
    echo "		 -o <output Folder>"
    echo "		 -log <log file>"
    echo "		 -op <output prefix>"
    echo -e "\nEXAMPLE: ./run_makeReport.sh -o /storage/crihan-msa/RunsPlateforme/GAIIx/111125_HWUSI-EAS1884_00002_FC64F86AAXX/CNV/CANOES -log /storage/crihan-msa/RunsPlateforme/GAIIx/111125_HWUSI-EAS1884_00002_FC64F86AAXX/CANOES-151008-23172-AGBedSplit.log"
}

# get the arguments of the command line
if [ $# -lt 6 ]; then
	usage
	exit
else
	while [ "$1" != "" ]; do
	    case $1 in

		-o | --output )         shift
					if [ "$1" != "" ]; then
						outputFolder=$1
					else
						usage
						exit
					fi
		                        ;;
		-log | --log )         shift
					if [ "$1" != "" ]; then
						logFile=$1
					else
						usage
						exit
					fi
		                        ;;
		-op | --outPrefix )         shift
					if [ "$1" != "" ]; then
						outPrefix=$1
					else
						usage
						exit
					fi
		                        ;;		                        
		*)           		usage
		                        exit
		                        ;;
	    esac
	    shift
	done
fi

log=$(basename $logFile)

date=`grep "TIME: BEGIN CANOES" $logFile`
echo -e "\nDate de lancement de CANOES\n$date"

echo -e "\n-------------------------------------------------"

echo -e "\nL'ensemble des fichiers bam analysés (patients+pool) sont :"
echo -e "`grep -o "[0-9A-Za-z-]\+.sorted.dedup.withRG.real.BQSR.bam$" $logFile`"
listBAM=`grep -o "[0-9A-Za-z-]\+.sorted.dedup.withRG.real.BQSR.bam" $logFile`

echo -e "\n-------------------------------------------------"

echo -e "\nLes patients analysés sont :"
echo -e "`cat $outputFolder/$outPrefix"_SampleName.list"`"


#Erreur lors du lancement ?
echo -e "\n-------------------------------------------------"
Erreur=`grep -B 1 "Erreur" $logFile`
Erreur2=`grep -B 1 "ERROR" $logFile | sed s/\t//`
if [ "$Erreur" != "" ]
then
	#recuperer le dernier code patient ayant passé l'analyse
	patientOK=`echo $Erreur | grep -o "[0-9A-Za-z]\+-[0-9A-Za-z]\+"`
	echo -e "\nATTENTION : Erreur détecté lors de l'analyse!"
	#recuperer la ligne correspondante a ce dernier code patient
	LigneOK=`echo "$listBAM" | grep -n $patientOK | awk -F: '{print $1}'`
	#incrémenter de 1 pour avoir la ligne posant problème
	(( LigneError = $LigneOK + 1 ))
	#retrouver l'identifiant à la ligne posant problème
	patientErr=`echo "$listBAM" | tail -n +$LigneError`
	echo -e "Le fichier $patientErr ne peut pas être analysé par CANOES"
	echo -e "\n-------------------------------------------------"
	exit
elif [ "$Erreur2" != "" ]
then
	echo -e "$Erreur2"
	echo -e "CANOES ne peut probablement pas estimer la Variance."
	echo -e "\n-------------------------------------------------"
	exit
else
	echo -e "\nAucune erreur détéctée, tous les patients ont pu être analysés"
fi

echo -e "\n-------------------------------------------------"
echo -e "\nCNVs détectés sur toute la capture (fichier : "$outPrefix"_CNV-result.csv ) :"
echo -e "`grep "CNVs called" $logFile`"

echo -e "\n-------------------------------------------------"
if [ -f $outputFolder/$outPrefix"_CNV-result-AnnotComplil10k.bed" ]; 
then
	CNVs=`wc -m $outputFolder/$outPrefix"_CNV-result-AnnotComplil10k.bed" | awk '{print $1}'`

	if [ $CNVs -eq 60 ]; #Si il y a 60 caractères, il n'y a que le header -> pas de CNVs
	then
		echo -e "\n0 CNV détecté sur les zones d'interet diagnostique."
	else
		echo -e "\nCNVs détectés sur les zones d'interet diagnostique (fichier : "$outPrefix"_CNV-result-AnnotComplil10k.bed) : "
		echo -e "`cat $outputFolder/$outPrefix"_CNV-result-AnnotComplil10k.bed"`"
	fi
else #Si le fichier n'existe pas, rien n'a été détecté
	echo -e "Aucun CNV détecté sur les zones d'interet diagnostique."
fi
echo -e "\n-------------------------------------------------"


