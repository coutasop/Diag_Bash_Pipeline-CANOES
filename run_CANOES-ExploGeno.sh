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
	echo -e "# 3. Créé fichier listBAM.txt et sample.list à partir des BQSR.bam présents dans le dossier input                                                            #"
	echo -e "# 4. Calcul la profondeur moyenne de toutes les cibles du .bed pour chaque .bam (CNV.reads.txt)                                                              #"
	echo -e "# 5. Etape de normalisation par suppression des régions dont un ou plusieurs individus présente moins de 20 reads                                            #"
	echo -e "# 6. Lancement du calling des CNV. Script R CANOES.                                                                                                          #"
	echo -e "#------------------------------------------------------------------------------------------------------------------------------------------------------------#"
    echo -e "\nUSAGE: run_CANOES_v2.sh -i <directory> -o <directory> -op <output prefix> -ref <file> -bed <file> -ann <file> -type <Explore or Genotype> -bt2 <path> -can <path> -gatk <path>"
    echo "		 -i <input Folder> contain *BQSR.bam"
    echo "		 -o <output Folder>"
    echo "		 -op <output prefix>"
    echo "		 -ref <fasta genome reference file>"
    echo "		 -bed <target bed file>"
    echo "		 -ann <annotated bed file (pathway)>"
    echo "		 -type <type of analysis \"Explore\" ou \"Genotype\">"
    echo "		 -bt2 <bedtools2 program path>"
    echo "		 -can <CANOES program path + Scripts perl>"
    echo "		 -gatk <GATK program path>"
    echo -e "\nEXAMPLE: ./run_CANOES_v2.sh -i /storage/crihan-msa/RunsPlateforme/GAIIx/111125_HWUSI-EAS1884_00002_FC64F86AAXX/BWA-GATK/BAM -o /storage/crihan-msa/RunsPlateforme/GAIIx/111125_HWUSI-EAS1884_00002_FC64F86AAXX/CNV/CANOES -op Run1Diag -ref /storage/crihan-msa/RunsPlateforme/Reference/Homo_sapiens/hg19/human_g1k_v37.fasta -bed /storage/crihan-msa/RunsPlateforme/Reference/Capture/MMR/036540_D_BED_20110915-DiagK_colique-U614_TARGET.bed -ann /storage/crihan-msa/RunsPlateforme/Reference/Capture/MMR/DiagCapture-11genes_20131009_sansChr.bed -bt2 /opt/bedtools2/bin -can /opt/CANOE -gatk /opt/GATK"
}

# get the arguments of the command line
if [ $# -lt 20 ]; then
	usage
	exit
else
	while [ "$1" != "" ]; do
	    case $1 in
		-i | --input )         shift
					if [ "$1" != "" ]; then
						inputFolder=$1
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
		-type | --type )         shift
					if [ "$1" != "" ]; then
						type=$1
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
		*)           		usage
		                        exit
		                        ;;
	    esac
	    shift
	done
fi

echo -e "\tTIME: BEGIN CANOES ".`date`

#Test if the output directory exists, if no, create it
if [ -d $outputFolder ]; then
 echo -e "\n\tOUTPUT FOLDER: $outputFolder (folder already exist)" 
else
 mkdir -p $outputFolder 
 echo -e "\n\tOUTPUT FOLDER : $outputFolder (folder created)"
fi

#créer fichier taux de GC si besoins
bedFich=$(basename $bedFile) #nom du fichier sans path
bedBaseName=${bedFich%.*} #nom du fichier sans path et sans sa dernière extension
gcFile=${bedFile%.*}".gc.txt"

if [ -f $gcFile ]; then
	 echo -e "\n\tGC File for bed: $gcFile (file already exist)" 
else
	echo -e "\n\tCOMMAND : java -Xmx4g -jar $GATKPath/GenomeAnalysisTK.jar -T GCContentByInterval -L $bedFile -R $refFile -o $gcFile"
	java -Xmx4g -jar $gatkPath/GenomeAnalysisTK.jar -T GCContentByInterval -L $bedFile -im OVERLAPPING_ONLY -R $refFile -o $gcFile
	echo -e "\n\tGC File for bed: $gcFile (file created)" 
fi

cd $inputFolder
echo -e "\n\tCreate listBAM.txt"
echo -e "\tCOMMAND : ls *BQSR.bam > $outputFolder/listBAM.txt"
echo -e "\n\tCOMMAND : ls -l $inputFolder/*BQSR.bam | cut -d " " -f 10 > $outputFolder/listBAM.txt"
ls *BQSR.bam > $outputFolder/listBAM.txt

echo -e "\n\tCreate sample.list"
echo -e "\tCOMMAND : cat $outputFolder/listBAM.txt | cut -d "." -f 1 > $outputFolder/sample.list"
cat $outputFolder/listBAM.txt | cut -d "." -f 1 > $outputFolder/sample.list

echo -e "\n\tCreate $outPrefix_CNV.reads.txt"
echo -e "\tCOMMAND : $bedtools2Path/bedtools multicov -bams `cat $outputFolder/listBAM.txt` -bed $bedFile -q 20 > $outputFolder/"$outPrefix"_CNV.reads.txt"
$bedtools2Path/bedtools multicov -bams `cat $outputFolder/listBAM.txt` -bed $bedFile -q 20 > $outputFolder"/"$outPrefix"_CNV.reads.txt"


#Correct chrX & chrY -> rename chr23 & chr24
echo -e "\nTEST for ChrX and ChrY name correction";
chrX=0;
chrY=0;
chrX=$(echo `cat $outputFolder"/"$outPrefix"_CNV.reads.txt" | cut -f 1 | sort -u | grep "X"`)
chrY=$(echo `cat $outputFolder"/"$outPrefix"_CNV.reads.txt" | cut -f 1 | sort -u | grep "Y"`)
if [ "$chrX" != "X" ];
then
	echo -e "\tNo chrX in file : SKIP correction"
else
	echo -e "\tchrX in file : CORRECTION"
	echo -e "\tchrX => chr23";
	cp $outputFolder"/"$outPrefix"_CNV.reads.txt" $outputFolder"/"$outPrefix"_CNV.reads.txt.tmp"
	cp $outputFolder"/"$outPrefix"_CNV.reads.txt" $outputFolder"/"$outPrefix"_CNV.noChrXYCorrection.reads.txt"
	cp $gcFile $outputFolder"/"$bedBaseName".gc.txt"
	gcFile=$outputFolder"/"$bedBaseName".gc.txt"
	cp $gcFile $gcFile".tmp"
	sed -e s/"X"/"23"/ $outputFolder"/"$outPrefix"_CNV.reads.txt.tmp" > $outputFolder"/"$outPrefix"_CNV.reads.txt"
	sed -e s/"X"/"23"/ $gcFile".tmp" > $gcFile
fi
if [ "$chrY" != "Y" ];
then
	echo -e "\tNo chrY in file : SKIP correction"
else
	echo -e "\tchrY in file : CORRECTION"
	echo -e "\tchrY => chr24";
	cp $outputFolder"/"$outPrefix"_CNV.reads.txt" $outputFolder"/"$outPrefix"_CNV.reads.txt.tmp"
	cp $gcFile $outputFolder"/"$bedBaseName".gc.txt"
	gcFile=$outputFolder"/"$bedBaseName".gc.txt"
	cp $gcFile $gcFile".tmp"
	sed -e s/"Y"/"24"/ $outputFolder"/"$outPrefix"_CNV.reads.txt.tmp" > $outputFolder"/"$outPrefix"_CNV.reads.txt"
	sed -e s/"Y"/"24"/ $gcFile".tmp" > $gcFile
fi

##EXPLORE MODE
if [ $type == "Explore" ];
then
	#Rscript
	#	args[1] #path CANOES.R
	#	args[2] #path gc file 
	#	args[3] #path reads.txt file
	#	args[4] #output file name prefix
	#	args[5] #analysis type "Explore" or "Genotype" or "Both"
	echo -e "\n\tRun script CANOES_ExploGeno.R"

	echo -e "\tCOMMAND : Rscript $canoesPath/run_CANOES-ExploGeno.R "$canoesPath"/CANOES.R" $gcFile $outPrefix"_CNV.reads.txt" $outPrefix $type""
	Rscript $canoesPath/run_CANOES-ExploGeno.R $canoesPath"/CANOES.R" $gcFile $outPrefix"_CNV.reads.txt" $outPrefix $type
	
	CNVnb=$(wc -l $outPrefix"_CNV-result.csv" | cut -d" " -f 1)
fi

##GENOTYPE MODE
if [ $type == "Genotype" ];
then
	#étape de normalisation par suppression des régions dont un ou plusieurs individus présente moins de 0 reads
	cd $outputFolder
	echo -e "\n\tRemove regions having mean depth < 0"
	echo -e "\tCOMMAND : perl $canoesPath/cleanCANOESentries.pl --gc $gcFile --reads CNV.reads.txt --outPrefix CNVClean0"
	perl $canoesPath/cleanCANOESentries.pl --gc $gcFile --reads $outPrefix"_CNV.reads.txt" --outPrefix $outPrefix"_CNVClean0"
	
	#paste line number in the first column using grep
	while read line; do grep -n -e "$line" $outPrefix"_CNVClean0.reads.txt"; done < $outPrefix"_CNVClean0.reads.txt" > $outPrefix"_CNVClean0.readsLineNb.txt"
	#generate xcnvs file
	awk -F"\t" '{OFS="\t";print $1,$2,$3}' $outPrefix"_CNVClean0.readsLineNb.txt" | sed -e "s/:/\t/g" | awk -F"\t" 'BEGIN {print "SAMPLE;CNV;INTERVAL;KB;CHR;MID_BP;TARGETS;NUM_TARG;MLCN;Q_SOME"}''{printf "\"SAMPLE\";\"DUP\";\""$2":"$3"-"$4"\";"$4-$3";"$2";%i;\""$1".."$1"\";"$1";3;99\n", ($4+$3)/2}' > $outPrefix"_CNVClean0.xcnvs.txt"
	
	echo -e "\n\tRun script CANOES_ExploGeno.R"
	#Rscript
	#    args[1] #path CANOES.R
	#    args[2] #path gc file
	#    args[3] #path reads.txt file
	#    args[4] #output file name prefix	
	#	 args[5] #analysis type "Explore" or "Genotype" or "Both"
	#    args[6] #path to sample.list to analyses
	#    args[7] #path to xcnvs CNV file
	
	echo -e "\tCOMMAND : Rscript $canoesPath/run_CANOES-ExploGeno.R "$canoesPath"/CANOES.R" $outPrefix"_CNVClean0.gc.txt" $outPrefix"_CNVClean0.reads.txt" $outPrefix $type $outputFolder/sample.list $outPrefix"_CNVClean0.xcnvs.txt"
	Rscript $canoesPath/run_CANOES-ExploGeno.R $canoesPath"/CANOES.R" $outPrefix"_CNVClean0.gc.txt" $outPrefix"_CNVClean0.reads.txt" $outPrefix $type $outputFolder/sample.list $outPrefix"_CNVClean0.xcnvs.txt"
	
	CNVnb=$(wc -l $outPrefix".genotype.csv" | cut -d" " -f 1)
	
	cp $outPrefix".genotype.csv" $outPrefix"_CNV-result.csv"
	
fi

#Si CNV détectés
if [ $CNVnb -eq 1 ];
then
	echo -e "\n\tNO CNV DETECTED"
else

	#Correct chrX & chrY -> rename chr23 & chr24
	echo -e "\nTEST for ChrX and ChrY REVERT Name Correction";
	#~ chrX=0;
	#~ chrY=0;
	#~ chrX=$(echo `cat $outPrefix"_CNV-result.csv" | grep "23:"`)
	#~ chrX=$(echo `cat $outPrefix"_CNV-result.csv" | grep "24:"`)
	if [ "$chrX" != "X" ];
	then
		echo -e "\tNo chr23 in file : SKIP Revert Correction"
	else
		echo -e "\tchr23 in file : REVERT CORRECTION"
		echo -e "\tchr23 => chrX";
		cp $outPrefix"_CNV-result.csv" $outPrefix"_CNV-result.csv.tmp"
		sed -e s/"23:"/"X:"/ $outPrefix"_CNV-result.csv.tmp" > $outPrefix"_CNV-result.csv"
	fi
	if [ "$chrY" != "Y" ];
	then
		echo -e "\tNo chr24 in file : SKIP Revert Correction"
	else
		echo -e "\tchr24 in file : REVERT CORRECTION"
		echo -e "\tchr24 => chrY";
		cp $outPrefix"_CNV-result.csv" $outPrefix"_CNV-result.csv.tmp"
		sed -e s/"24:"/"Y:"/ $outPrefix"_CNV-result.csv.tmp" > $outPrefix"_CNV-result.csv"
	fi


	##EXPLORE MODE
	if [ $type == "Explore" ];
	then
		#Explore
		#extraction des colonnes et changement d'ordre pour transformation en .bed
		while read line; do 
			IFS=';' read -ra ADDR <<< "$line"; 
			echo -e "${ADDR[3]}\t${ADDR[1]}\t${ADDR[2]}\t${ADDR[4]}\t${ADDR[8]}\t${ADDR[9]}\t${ADDR[10]}"; 
		done < $outPrefix"_CNV-result.csv" > $outPrefix"_CNV-result.intervals"

		#passage du format intervals vers bed
		echo -e "\n\tCOMMAND : sed -e 's/:/\t/' $outPrefix'_CNV-result.intervals' | sed -e 's/-/\t/' | sed -e s/\"//g | sed 's/INTERVAL/CHR\tSTART\tEND/' > $outPrefix'_CNV-result.bed'"
		sed -e "s/:/\t/" $outPrefix"_CNV-result.intervals" | sed -e "s/-/\t/" | sed -e s/\"//g | sed "s/INTERVAL/CHR\tSTART\tEND/" > $outPrefix"_CNV-result.bed"
		rm $outPrefix"_CNV-result.intervals"

		if [ -z "$bedAnnot" ];
		then 
			echo -e "\tBed file for annotation not defined : SKIP"
		else
			tail -n +2 $outPrefix"_CNV-result.bed" > $outPrefix"_CNV-result-Annot.bed.tmp"
			
			#extraction des CNV détectés au sein de notre sélection de gènes
			echo -e "\tExtraction des CNV detectes au sein de notre selection de genes"
			#~ echo -e "\tCOMMAND : intersectBed -a outPrefix_CNV-result.bed -b $bedAnnot -wo > outPrefix_CNV-result-Annot.bed"
			intersectBed -a $outPrefix"_CNV-result-Annot.bed.tmp" -b $bedAnnot -wo > $outPrefix"_CNV-result-Annot.bed"
			
			#compilation des CNV proches pour faciliter la lecture (distance de 10Kb, même individus, même DUP/DEL)
			echo -e "\tCompilation des CNV proches pour faciliter la lecture (distance de 10Kb, meme individus, meme DUP/DEL)"
			#~ echo -e "\tCOMMAND : perl $canoesPath/compilBedFromCANOES.pl --in outPrefix_CNV-result-Annot.bed --out outPrefix_CNV-result-AnnotComplil10k.bed"
			perl $canoesPath/compilBedFromCANOES.pl --in $outPrefix"_CNV-result-Annot.bed" --out $outPrefix"_CNV-result-AnnotComplil10k.bed"
			
			rm $outPrefix"_CNV-result-Annot.bed.tmp"
		fi
	fi	
	
	##GENOTYPE MODE
	if [ $type == "Genotype" ];
	then
		
		cp $outPrefix"_CNV-result.csv" $outPrefix".genotype.csv"
		awk -F";" '($5>=90&&$3<=70) {OFS="\t";print $1,"DEL",$0'} $outPrefix"_CNV-result.csv" > $outPrefix".genotype.potentialDel.csv"
		awk -F";" '($3>=90&&$5<=70) {OFS="\t";print $1,"DUP",$0'} $outPrefix"_CNV-result.csv" > $outPrefix".genotype.potentialDup.csv"
		
		#extraction des colonnes et changement d'ordre pour transformation en .bed
		paste <(cat $outPrefix".genotype.potentialDel.csv" | cut -f 1 | cut -d "." -f 2) <(cat $outPrefix".genotype.potentialDel.csv") | sed -e s/\"//g | sed -e s/:/\\t/ | sed -e s/-/\\t/ | sed -e s/_/\\t/ | sed -e s/\\./\\t/g | sed -e s/\;/\\t/g | awk -v OFS="\t" '{print $1,$2,$3,$4,$7,$9,$11,$12,$13,$14}' > $outPrefix".genotype.potentialDel.bed"
		paste <(cat $outPrefix".genotype.potentialDup.csv" | cut -f 1 | cut -d "." -f 2) <(cat $outPrefix".genotype.potentialDup.csv") | sed -e s/\"//g | sed -e s/:/\\t/ | sed -e s/-/\\t/ | sed -e s/_/\\t/ | sed -e s/\\./\\t/g | sed -e s/\;/\\t/g | awk -v OFS="\t" '{print $1,$2,$3,$4,$7,$9,$11,$12,$13,$14}' > $outPrefix".genotype.potentialDup.bed"
		#Annotation
		echo -e "chr\tstart\tend\tsample\ttype\ttarget\tNQDel\tSQDel\tNQDup\tSQDup\tchr\tstart\tend\tNM\tGene\tgeneSize\ttargetSize" > $outPrefix".genotype.potentialDel.annot.bed"
		echo -e "chr\tstart\tend\tsample\ttype\ttarget\tNQDel\tSQDel\tNQDup\tSQDup\tchr\tstart\tend\tNM\tGene\tgeneSize\ttargetSize" > $outPrefix".genotype.potentialDup.annot.bed"	
		intersectBed -a $outPrefix".genotype.potentialDel.bed" -b $bedAnnot -wo >> $outPrefix".genotype.potentialDel.annot.bed"
		intersectBed -a $outPrefix".genotype.potentialDup.bed" -b $bedAnnot -wo >> $outPrefix".genotype.potentialDup.annot.bed"	
		
	fi
	
fi

echo -e "\tTIME: END CANOES ".`date`
