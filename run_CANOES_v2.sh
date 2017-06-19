#!/bin/bash
#
# Sophie COUTANT
# 28/11/2014
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Ce script permet l'automatisation du lancement de CANOES à partir des bam GATK (BQSR.bam only)															 #
#------------------------------------------------------------------------------------------------------------------------------------------------------------#


# usage
function usage
{
	echo -e "#------------------------------------------------------------------------------------------------------------------------------------------------------------#"
	echo -e "# Ce script permet l'automatisation du lancement de CANOES à partir des bam GATK (BQSR.bam only)                                                             #"
	echo -e "# 1. Test si le dossier output existe, si non, création                                                                                                      #"
	echo -e "# 2. Créé le fichier du taux de GC sur le bed spécifié (si non existant)                                                                                     #"
	echo -e "# 3. Calcul la profondeur moyenne de toutes les cibles du .bed pour chaque .bam (CNV.reads.txt)                                                              #"
	echo -e "# 4. Lancement du calling des CNV. Script R CANOES.                                                                                                          #"
	echo -e "#------------------------------------------------------------------------------------------------------------------------------------------------------------#"
    echo -e "\nUSAGE: run_CANOES.sh -i <file> -pool <file> -o <directory> -op <output prefix> -ref <file> -bed <file> -ann <file> -bt2 <path> -can <path> -gatk <path> -type <keyword>"
    echo "		 -i <input Bam List> contain full path of sample bam *BQSR.bam"
	echo "		 -pool <pool bam List> contain full path of pool bam *BQSR.bam"	
    echo "		 -o <output Folder>"
    echo "		 -op <output prefix>"
    echo "		 -ref <fasta genome reference file>"
    echo "		 -bed <target bed file>"
    echo "		 -ann <annotated bed file (pathway)>"
    echo "		 -bt2 <bedtools2 program path>"
    echo "		 -can <CANOES program path + Scripts perl>"
    echo "		 -gatk <GATK program path>"
    echo "		 -type <Explore or Genotype>"
    echo -e "\nEXAMPLE COLON: /opt/pipeline_NGS/CANOES/run_CANOES_v2.sh -i /storage/crihan-msa/RunsPlateforme/MiSeq/161208_M02807_0139_000000000-AWA2L/BWA-GATK_Run139-COLON-MS32/BAM/CNVsampleBamList.txt -pool /storage/crihan-msa/RunsPlateforme/MiSeq/161208_M02807_0139_000000000-AWA2L/BWA-GATK_Run139-COLON-MS32/BAM/CNVpoolBamList.txt -o /storage/crihan-msa/RunsPlateforme/MiSeq/161208_M02807_0139_000000000-AWA2L/CNV/CANOES -op Run139-COLON-MS32 -ref /storage/crihan-msa/RunsPlateforme/Reference/Homo_sapiens/hg19/human_g1k_v37.fasta -bed /storage/crihan-msa/RunsPlateforme/Reference/Capture/MMR/036540_D_BED_20110915-DiagK_colique-U614_TARGET-CNV_PMS2CL.bed -ann /storage/crihan-msa/RunsPlateforme/Reference/Capture/MMR/DiagCapture-11genes_20131009_sansChr_forCNV.bed -bt2 /opt/BEDTOOLS/ -can /opt/pipeline_NGS/CANOES -gatk /opt/GATK -type Explore"
    echo -e "\nPour faire les listes -i et -pool se mettre dans le dossier contenant les bam voulu et taper l'une des commandes suivante:"
    echo -e " Pour le Pool:"    
    echo -e "\t ls -d1 \$PWD/*BQSR.bam > CNVpoolBamList.txt"
    echo -e " Pour les Patients:"
    echo -e "\t ls -d1 \$PWD/*BQSR.bam > CNVsampleBamList.txt"	
    echo -e "\tPour HYCO ou pour faire une liste par patient, taper:"
    echo -e "\t for f in \`ls -d1 \$PWD/*BQSR.bam\`; do name=\$(echo \`basename \$f\` | awk -F. '{print \$1}' ); echo \$f > \$name\"_SampleBamList.txt\"; done"
}

# get the arguments of the command line
if [ $# -lt 22 ]; then
	usage
	exit
else
	while [ "$1" != "" ]; do
	    case $1 in
		-i | --input )         shift
					if [ "$1" != "" ]; then
						sampleBamList=$1
					else
						usage
						exit
					fi
		                        ;;
		-pool | --pool )         shift
					if [ "$1" != "" ]; then
						poolBamList=$1
					else
						usage
						exit
					fi
		                        ;;
		-o | --output )         shift
					if [ "$1" != "" ]; then
						outputFolder=$1
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
		-bed | --bedFile )         shift
					if [ "$1" != "" ]; then
						bedFile=$1
					else
						usage
						exit
					fi
		                        ;;		  
		-ann | --bedAnnot )         shift
					if [ "$1" != "" ]; then
						bedAnnot=$1
					else
						usage
						exit
					fi
		                        ;;
		-ref | --refFile )         shift
					if [ "$1" != "" ]; then
						refFile=$1
					else
						usage
						exit
					fi
		                        ;;
		-bt2 | --bedtools2Path )         shift
					if [ "$1" != "" ]; then
						bedtools2Path=$1
					else
						usage
						exit
					fi
		                        ;;
   		-can | --canoesPath )         shift
					if [ "$1" != "" ]; then
						canoesPath=$1
					else
						usage
						exit
					fi
		                        ;;	 
   		-gatk | --gatkPath )         shift
					if [ "$1" != "" ]; then
						gatkPath=$1
					else
						usage
						exit
					fi
		                        ;;	
   		-type | --type )         shift
					if [ "$1" != "" ]; then
						type=$1
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

logpath=`pwd`
dat=`date +"%Y-%m-%d_%Hh%Mm%Ss"`
log="$outputFolder/CANOES-"$type"_"$outPrefix"_"$dat".log"

echo -e "\nTIME: BEGIN CANOES ".`date`

#Test if the output directory exists, if no, create it
if [ -d $outputFolder ]; then
 touch $log
 echo -e "\tTIME: BEGIN CANOES ".`date` >> $log
 echo -e "\nLOG FILE: $log "
 echo -e "\n\tOUTPUT FOLDER: $outputFolder (folder already exist)"  >> $log
else
 mkdir -p $outputFolder 
 touch $log
 echo -e "\tTIME: BEGIN CANOES ".`date` >> $log
 echo -e "\nLOG FILE: $log "
 echo -e "\n\tOUTPUT FOLDER : $outputFolder (folder created)" >> $log
fi

#créer fichier taux de GC si besoins
bedFich=$(basename $bedFile) #nom du fichier sans path
bedBaseName=${bedFich%.*} #nom du fichier sans path et sans sa dernière extension
gcFile=${bedFile%.*}".gc.txt"

if [ -f $gcFile ]; then
	 echo -e "\n\tGC File for bed: $gcFile (file already exist)"  >> $log
else
	echo -e "\n\tCOMMAND : java -Xmx4g -jar $GATKPath/GenomeAnalysisTK.jar -T GCContentByInterval -L $bedFile -R $refFile -o $gcFile" >> $log
	java -Xmx4g -jar $gatkPath/GenomeAnalysisTK.jar -T GCContentByInterval -L $bedFile -im OVERLAPPING_ONLY -R $refFile -o $gcFile
	echo -e "\n\tGC File for bed: $gcFile (file created)"  >> $log
fi

#Preparer les fichiers de liste de patients
IndOrder=""
SampleOrder=""
cat $poolBamList $sampleBamList > $outputFolder/$outPrefix"_AllBamPath.list"
> $outputFolder/$outPrefix"_AllName.list"
> $outputFolder/$outPrefix"_SampleName.list"

#Preparer le fichier multicov global
cp $bedFile $outputFolder/$outPrefix"_CNV.reads.txt"

##Traitement des fichier du Pool
##-Creer les fichiers multicov pour le pool s'ils n'existent pas déjà
##-concatener les profondeur de chaque bam du pool au fichier multicov (reads.txt) global
echo -e "\n poolBam :" >> $log
while read poolBam ;
do
	echo -e "$poolBam" >> $log
	poolPath=`dirname $poolBam`
	poolName=`basename $poolBam | awk -F. '{print $1}'`
	IndOrder=$IndOrder" "$poolName
	echo $poolName >> $outputFolder/$outPrefix"_AllName.list"
	#-Creer les fichiers multicov pour le pool s'ils n'existent pas déjà
	if [ -f $poolPath/$poolName"_multicov.txt" ]
	then
		echo -e "\tMulticov File: "$poolPath"/"$poolName"_multicov.txt (file already exist)"
		echo -e "\tMulticov File: "$poolPath"/"$poolName"_multicov.txt (file already exist)"  >> $log
	else
	    echo -e "\tCOMMAND : $bedtools2Path/bedtools multicov -bams $poolBam -bed $bedFile -q 20 > "$poolPath"/"$poolName"_multicov.txt"
		echo -e "\tCOMMAND : $bedtools2Path/bedtools multicov -bams $poolBam -bed $bedFile -q 20 > "$poolPath"/"$poolName"_multicov.txt" >> $log
		$bedtools2Path/bedtools multicov -bams $poolBam -bed $bedFile -q 20 > $poolPath/$poolName"_multicov.txt"
	fi
	#-Ajouter le multicov de chaque bam du pool au multicov global au fur et à mesure
	cp $outputFolder/$outPrefix"_CNV.reads.txt" $outputFolder/$outPrefix"_CNV.reads.txt.tmp"
	pr -mts $outputFolder/$outPrefix"_CNV.reads.txt.tmp" <(cut -f4 $poolPath/$poolName"_multicov.txt") > $outputFolder/$outPrefix"_CNV.reads.txt"
done < $poolBamList

##Traitement des fichier des Sample
##-Creer les fichiers multicov pour les samples s'ils n'existent pas déjà
##-concatener les profondeur de chaque bam de sample au fichier multicov (reads.txt) global
echo -e "\n sampleBam :" >> $log
while read sampleBam ;
do
	echo -e "$sampleBam" >> $log
	samplePath=`dirname $sampleBam`
	sampleName=`basename $sampleBam | awk -F. '{print $1}'`
	IndOrder=$IndOrder" "$sampleName
	echo $sampleName >> $outputFolder/$outPrefix"_AllName.list"
	echo $sampleName >> $outputFolder/$outPrefix"_SampleName.list"
	#-Creer les fichiers multicov pour le pool s'ils n'existent pas déjà
	if [ -f $samplePath/$sampleName"_multicov.txt" ]
	then	
		echo -e "\tMulticov File: "$samplePath"/"$sampleName"_multicov.txt (file already exist)"
		echo -e "\tMulticov File: "$samplePath"/"$sampleName"_multicov.txt (file already exist)"  >> $log
	else
	    echo -e "\tCOMMAND : $bedtools2Path/bedtools multicov -bams $sampleBam -bed $bedFile -q 20 > "$samplePath"/"$sampleName"_multicov.txt"
		echo -e "\tCOMMAND : $bedtools2Path/bedtools multicov -bams $sampleBam -bed $bedFile -q 20 > "$samplePath"/"$sampleName"_multicov.txt" >> $log
		$bedtools2Path/bedtools multicov -bams $sampleBam -bed $bedFile -q 20 > $samplePath/$sampleName"_multicov.txt"
	fi
	#-Ajouter le multicov de chaque bam de sample au multicov global au fur et à mesure
	cp $outputFolder/$outPrefix"_CNV.reads.txt" $outputFolder/$outPrefix"_CNV.reads.txt.tmp"
	pr -mts $outputFolder/$outPrefix"_CNV.reads.txt.tmp" <(cut -f4 $samplePath/$sampleName"_multicov.txt") > $outputFolder/$outPrefix"_CNV.reads.txt"
done < $sampleBamList

#Supression fichier temporaire
rm $outputFolder/$outPrefix"_CNV.reads.txt.tmp"


#Correct chrX & chrY -> rename chr23 & chr24
echo -e "\nTEST for ChrX and ChrY name correction" >> $log
chrX=0;
chrY=0;
chrX=$(echo `cat $outputFolder"/"$outPrefix"_CNV.reads.txt" | cut -f 1 | sort -u | grep "X"`)
chrY=$(echo `cat $outputFolder"/"$outPrefix"_CNV.reads.txt" | cut -f 1 | sort -u | grep "Y"`)
if [ "$chrX" != "X" ];
then
	echo -e "\tNo chrX in file : SKIP correction" >> $log
else
	echo -e "\tchrX in file : CORRECTION" >> $log
	echo -e "\tchrX => chr23" >> $log
	cp $outputFolder"/"$outPrefix"_CNV.reads.txt" $outputFolder"/"$outPrefix"_CNV.reads.txt.tmp"
	cp $gcFile $outputFolder"/"$bedBaseName".gc.txt"
	gcFile=$outputFolder"/"$bedBaseName".gc.txt"
	cp $gcFile $gcFile".tmp"
	sed -e s/"X"/"23"/ $outputFolder"/"$outPrefix"_CNV.reads.txt.tmp" > $outputFolder"/"$outPrefix"_CNV.reads.txt"
	sed -e s/"X"/"23"/ $gcFile".tmp" > $gcFile
fi
if [ "$chrY" != "Y" ];
then
	echo -e "\tNo chrY in file : SKIP correction" >> $log
else
	echo -e "\tchrY in file : CORRECTION" >> $log
	echo -e "\tchrY => chr24" >> $log
	cp $outputFolder"/"$outPrefix"_CNV.reads.txt" $outputFolder"/"$outPrefix"_CNV.reads.txt.tmp"
	cp $gcFile $outputFolder"/"$bedBaseName".gc.txt"
	gcFile=$outputFolder"/"$bedBaseName".gc.txt"
	cp $gcFile $gcFile".tmp"
	sed -e s/"Y"/"24"/ $outputFolder"/"$outPrefix"_CNV.reads.txt.tmp" > $outputFolder"/"$outPrefix"_CNV.reads.txt"
	sed -e s/"Y"/"24"/ $gcFile".tmp" > $gcFile
fi

##-----------------------------------------------------------------------------------##
##Se placer dans l'output folder et lancer R
#~ mv $log $outputFolder/$log
cd $outputFolder

##EXPLORE MODE
if [ $type == "Explore" ];
then

	#Rscript
	#	args[1] #path CANOES.R
	#	args[2] #path gc file 
	#	args[3] #path reads.txt file
	#	args[4] #output file name prefix
	#	args[5] #analysis type "Explore" or "Genotype"
	#   args[6] #path to SampleName.list to analyses
	#   args[7] #path to AllName.list (list of all the sample in the reads.txt file)
	echo -e "\n\tRun script CANOES.R" >> $log
	echo -e "\tCOMMAND : Rscript $canoesPath/run_CANOES_v2.R "$canoesPath"/CANOES.R" $gcFile $outPrefix"_CNV.reads.txt" $outPrefix $type $outPrefix"_SampleName.list" $outPrefix"_AllName.list" >> $log
	echo -e "\tCOMMAND : Rscript $canoesPath/run_CANOES_v2.R "$canoesPath"/CANOES.R" $gcFile $outPrefix"_CNV.reads.txt" $outPrefix $type $outPrefix"_SampleName.list" $outPrefix"_AllName.list"
	Rscript $canoesPath/run_CANOES_v2.R $canoesPath"/CANOES.R" $gcFile $outPrefix"_CNV.reads.txt" $outPrefix $type $outPrefix"_SampleName.list" $outPrefix"_AllName.list"  >> $log

	
	if [ -f $outPrefix"_CNV-result.csv" ];
	then
		CNVnb=$(wc -l $outPrefix"_CNV-result.csv" | cut -d" " -f 1)
	else
		echo -e "\tThe Result file "$outPrefix"_CNV-result.csv isn't generated" >> $log
		#~ touch $outPrefix"_CNV-result.csv"
	fi
fi

##GENOTYPE MODE
if [ $type == "Genotype" ];
then
	#étape de normalisation par suppression des régions dont un ou plusieurs individus présente moins de 0 reads
	#~ echo -e "\n\tRemove regions having mean depth < 0" >> $log
	#~ echo -e "\tCOMMAND : perl $canoesPath/cleanCANOESentries.pl --gc $gcFile --reads CNV.reads.txt --outPrefix CNVClean0" >> $log
	#~ perl $canoesPath/cleanCANOESentries.pl --gc $gcFile --reads $outPrefix"_CNV.reads.txt" --outPrefix $outPrefix"_CNVClean0"
	
	#paste line number in the first column using grep
	while read line; do grep -n -e "$line" $outPrefix"_CNV.reads.txt"; done < $outPrefix"_CNV.reads.txt" > $outPrefix".readsLineNb.txt"
	#generate xcnvs file
	awk -F"\t" '{OFS="\t";print $1,$2,$3}' $outPrefix".readsLineNb.txt" | sed -e "s/:/\t/g" | awk -F"\t" 'BEGIN {print "SAMPLE;CNV;INTERVAL;KB;CHR;MID_BP;TARGETS;NUM_TARG;MLCN;Q_SOME"}''{printf "\"SAMPLE\";\"DUP\";\""$2":"$3"-"$4"\";"$4-$3";"$2";%i;\""$1".."$1"\";"$1";3;99\n", ($4+$3)/2}' > $outPrefix".xcnvs.txt"
	
	echo -e "\n\tRun script CANOES_ExploGeno.R" >> $log
	#Rscript
	#   args[1] #path CANOES.R
	#   args[2] #path gc file
	#   args[3] #path reads.txt file
	#   args[4] #output file name prefix	
	#   args[5] #analysis type "Explore" or "Genotype"
	#   args[6] #path to SampleName.list to analyses
	#   args[7] #path to AllName.list (list of all the sample in the reads.txt file)
	#   args[8] #path to xcnvs CNV file
	
	echo -e "\tCOMMAND : Rscript $canoesPath/run_CANOES_v2.R "$canoesPath"/CANOES.R" $gcFile $outPrefix"_CNV.reads.txt" $outPrefix $type $outPrefix"_SampleName.list" $outPrefix"_AllName.list" $outPrefix".xcnvs.txt" >> $log
	Rscript $canoesPath/run_CANOES_v2.R $canoesPath"/CANOES.R" $gcFile $outPrefix"_CNV.reads.txt" $outPrefix $type $outPrefix"_SampleName.list" $outPrefix"_AllName.list" $outPrefix".xcnvs.txt" >> $log
	
	if [ -f $outPrefix".genotype.csv" ];
	then
		CNVnb=$(wc -l $outPrefix".genotype.csv" | cut -d" " -f 1)
	else
		echo -e "\tThe Result file "$outPrefix".genotype.csv isn't generated" >> $log
		#~ touch $outPrefix"_CNV-result.csv"
	fi
	
	if [ ! -d "FichiersTemp" ]; then
		mkdir FichiersTemp
	fi
	mv $outPrefix".readsLineNb.txt" FichiersTemp/.
	mv $outPrefix".xcnvs.txt" FichiersTemp/.
	#mv $gcFile FichiersTemp/.
	mv $outPrefix"_CNV.reads.txt" FichiersTemp/.
	

	
fi

##Fin R
##-----------------------------------------------------------------------------------##
if [ $type != "Genotype" ];
then

	if [ ! -f $outPrefix"_CNV-result.csv" ]; #Si erreur le fichier $outPrefix"_CNV-result.csv" existe mais est vide nb ligne = 0
	then
		echo -e "\n\tERROR DETECTED : "$outPrefix"_CNV-result.csv doesn't exist" >> $log
		echo -e "\n\tERROR DETECTED : "$outPrefix"_CNV-result.csv doesn't exist"
	elif [ $CNVnb -eq 1 ]; #Si pas de CNV il n'y a que le header dans le fichier $outPrefix"_CNV-result.csv" nb ligne = 1
	then
		echo -e "\n\tNO CNV DETECTED" >> $log
		echo -e "\n\tNO CNV DETECTED"
		
		if [ ! -d "FichiersTemp" ]; then
			mkdir FichiersTemp
		fi
		#mv $gcFile FichiersTemp/.
		mv $outPrefix"_CNV.reads.txt" FichiersTemp/.
		
		
	else	#Si CNV détectés création des fichiers de sortie

		#Correct chrX & chrY -> rename chr23 & chr24
		echo -e "\nTEST for ChrX and ChrY REVERT Name Correction" >> $log
		#~ chrX=0;
		#~ chrY=0;
		#~ chrX=$(echo `cat $outPrefix"_CNV-result.csv" | grep "23:"`)
		#~ chrX=$(echo `cat $outPrefix"_CNV-result.csv" | grep "24:"`)
		if [ "$chrX" != "X" ];
		then
			echo -e "\tNo chr23 in file : SKIP Revert Correction" >> $log
		else
			echo -e "\tchr23 in file : REVERT CORRECTION" >> $log
			echo -e "\tchr23 => chrX" >> $log
			cp $outPrefix"_CNV-result.csv" $outPrefix"_CNV-result.csv.tmp"
			sed -e s/"23:"/"X:"/ $outPrefix"_CNV-result.csv.tmp" > $outPrefix"_CNV-result.csv"
		fi
		if [ "$chrY" != "Y" ];
		then
			echo -e "\tNo chr24 in file : SKIP Revert Correction" >> $log
		else
			echo -e "\tchr24 in file : REVERT CORRECTION" >> $log
			echo -e "\tchr24 => chrY" >> $log
			cp $outPrefix"_CNV-result.csv" $outPrefix"_CNV-result.csv.tmp"
			sed -e s/"24:"/"Y:"/ $outPrefix"_CNV-result.csv.tmp" > $outPrefix"_CNV-result.csv"
		fi

		#extraction des colonnes et changement d'ordre pour transformation en .bed
		while read line; do 
			IFS=';' read -ra ADDR <<< "$line"; 
			echo -e "${ADDR[3]}\t${ADDR[1]}\t${ADDR[2]}\t${ADDR[4]}\t${ADDR[8]}\t${ADDR[9]}\t${ADDR[10]}"; 
		done < $outPrefix"_CNV-result.csv" > $outPrefix"_CNV-result.intervals"

		#passage du format intervals vers bed
		echo -e "\n\tCOMMAND : sed -e 's/:/\t/' $outPrefix'_CNV-result.intervals' | sed -e 's/-/\t/' | sed -e s/\"//g | sed 's/INTERVAL/CHR\tSTART\tEND/' > $outPrefix'_CNV-result.bed'" >> $log
		sed -e "s/:/\t/" $outPrefix"_CNV-result.intervals" | sed -e "s/-/\t/" | sed -e s/\"//g | sed "s/INTERVAL/CHR\tSTART\tEND/" > $outPrefix"_CNV-result.bed"
		rm $outPrefix"_CNV-result.intervals"
		
		if [ ! -d "FichiersTemp" ]; then
			mkdir FichiersTemp
		fi
		#mv $gcFile FichiersTemp/.
		mv $outPrefix"_CNV.reads.txt" FichiersTemp/.
		
		if [ -z "$bedAnnot" ];
		then 
			echo -e "\tBed file for annotation not defined : SKIP" >> $log
		else
			echo -e "\tCOMMAND : tail -n +2 "$outPrefix"_CNV-result.bed > "$outPrefix"_CNV-result-Annot.bed.tmp" >> $log
			tail -n +2 $outPrefix"_CNV-result.bed" > $outPrefix"_CNV-result-Annot.bed.tmp"
			
			#extraction des CNV détectés au sein de notre sélection de gènes
			echo -e "\tExtraction des CNV detectes au sein de notre selection de genes" >> $log
			echo -e "\tCOMMAND : $bedtools2Path/intersectBed -a "$outPrefix"_CNV-result-Annot.bed.tmp -b $bedAnnot -wo > "$outPrefix"_CNV-result-Annot.bed"
			$bedtools2Path/intersectBed -a $outPrefix"_CNV-result-Annot.bed.tmp" -b $bedAnnot -wo > $outPrefix"_CNV-result-Annot.bed"
			
			#compilation des CNV proches pour faciliter la lecture (distance de 10Kb, même individus, même DUP/DEL)
			echo -e "\tCompilation des CNV proches pour faciliter la lecture (distance de 10Kb, meme individus, meme DUP/DEL)" >> $log
			echo -e "\tCOMMAND : perl $canoesPath/compilBedFromCANOES.pl --in outPrefix_CNV-result-Annot.bed --out outPrefix_CNV-result-AnnotComplil10k.bed"
			perl $canoesPath/compilBedFromCANOES.pl --in $outPrefix"_CNV-result-Annot.bed" --out $outPrefix"_CNV-result-AnnotComplil10k.bed"
			
			rm *CNV-result-Annot.bed
			rm *CNV-result.bed
			rm $outPrefix"_CNV-result-Annot.bed.tmp"
		fi
	fi

	echo -e "\n\tTIME: END CANOES ".`date` >> $log
	echo -e "\nTIME: END CANOES ".`date`


	#~ mv $log $logpath/$log
	cd $outputFolder
	logName=$(basename $log)

	echo -e "\nGenerate Report: $outputFolder/Report_"$logName".txt"
	echo -e "\tCOMMAND : $canoesPath/run_makeReport.sh -o $outputFolder -log $logName -op $outPrefix > $outputFolder/Report_"$logName".txt"
	$canoesPath/run_makeReport.sh -o $outputFolder -log $logName -op $outPrefix > $outputFolder"/Report_"$logName".txt"
	

	
fi
	rm $outputFolder/*.list
	rm $outputFolder/*.tmp

	#rm $outputFolder/*.pdf
