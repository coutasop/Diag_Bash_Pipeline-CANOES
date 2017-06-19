#launch CANOES
args<-commandArgs(TRUE)
print(args)

if((args[5] == "Explore")){
	print("Run \"Explore\" mode") 
}else if ((args[5] == "Genotype")){
	#print("Missing Analysis type : please precise \"Explore\" or \"Genotype\"") 
	print("Run \"Genotype\" mode") 
}else{
	stop("Missing Analysis type : please precise \"Explore\" or \"Genotype\"") 
}

##Explore Mode
if (args[5] == "Explore"){
	print("CANOES launched in Explore mode")
	#	args[1] #path CANOES.R
	#	args[2] #path gc file 
	#	args[3] #path reads.txt file
	#	args[4] #output file name prefix
	#	args[5] #analysis type "Explore" or "Genotype"
	#   args[6] #path to SampleName.list to analyses
	#   args[7] #path to AllName.list (list of all the sample in the reads.txt file)

	source(args[1])
	
	gc <- read.table(args[2])$V2
	canoes.reads <- read.table(args[3])
	sample.names <- unlist(read.table(args[7], stringsAsFactors=FALSE))
	names(canoes.reads) <- c("chromosome", "start", "end", sample.names)
	analyse<-unlist(read.table(args[6],stringsAsFactors=FALSE))
	target <- seq(1, nrow(canoes.reads))
	canoes.reads <- cbind(target, gc, canoes.reads)
	xcnv.list <- vector('list', length(sample.names))

	for (i in 1:length(analyse)){
		xcnv.list[[i]] <- CallCNVs(analyse[i], canoes.reads)
	}

	xcnvs <- do.call('rbind', xcnv.list)
	write.csv2(xcnvs,file=paste(args[4],"_CNV-result.csv",sep = ""))

	if (nrow(xcnvs)!=0){
		pdf(paste(args[4],"_CNVplots.pdf",sep = ""))
		for (i in 1:nrow(xcnvs)){
		  PlotCNV(canoes.reads, xcnvs[i, "SAMPLE"], xcnvs[i, "TARGETS"])
		}
		dev.off()
	}
}

##Genotyping Mode
if (args[5] == "Genotype"){

	print("CANOES launched in Genotype mode")
	#   args[1] #path CANOES.R
	#   args[2] #path gc file
	#   args[3] #path reads.txt file
	#   args[4] #output file name prefix	
	#   args[5] #analysis type "Explore" or "Genotype"
	#   args[6] #path to SampleName.list to analyses
	#   args[7] #path to AllName.list (list of all the sample in the reads.txt file)
	#   args[8] #path to xcnvs CNV file


	source(args[1])

	#gc <- read.table("CNVClean20.gc.txt")$V2
	gc <- read.table(args[2])$V2
	#canoes.reads <- read.table("CNVClean20.reads.txt")
	canoes.reads <- read.table(args[3])
	sample.names <- unlist(read.table(args[7], stringsAsFactors=FALSE))
	print(sample.names)
	names(canoes.reads) <- c("chromosome", "start", "end", sample.names)
	analyse<-unlist(read.table(args[6],stringsAsFactors=FALSE))
	target <- seq(1, nrow(canoes.reads))
	canoes.reads <- cbind(target, gc, canoes.reads)
	#myCNV<-data.frame(SAMPLE=factor("F_SAMPLE"),CNV=factor("DUP"),INTERVAL=factor("3:47960150-47960363"),KB=0.213,CHR=3,MID_BP=47960256,TARGETS=factor("41172..41172"),NUM_TARG=1,MLCN=3,Q_SOME=47)
	x<-read.table(file=args[8],sep=";",header=T)
	geno.list <- list()
	for (i in 1:length(analyse)){
		geno.list[[analyse[i]]] <- GenotypeCNVs(x,analyse[i],canoes.reads,Tnum=2)
	}
	genos <- do.call('rbind',geno.list)

	write.csv2(genos,file=paste0(args[4],".genotype.csv")) 						


}

##COMMANDES OLIVIER
#1#ne garder que les gene ABeta
#intersectBed -a /storage/crihan-msa/DATA/ALZ_201508_CNG/READCOUNT/CNG201508.extend.clean20.reads.txt -b /opt/REFERENCE/ABetaNetwork_2015.bed -wa | awk -F"\t" '{print $1"\t"$2"\t"$3}' | uniq > tmp
#while read line; do grep -n -e "$line" /storage/crihan-msa/DATA/ALZ_201508_CNG/READCOUNT/CNG201508.extend.clean20.reads.txt; done < tmp > tmp2
#2#generer le fichier xcnvs.txt pour les à partir des reads netoyés.
#awk -F"\t" '{OFS="\t";print $1,$2,$3}' tmp2 | sed -e "s/:/\t/g" | awk -F"\t" 'BEGIN {print "SAMPLE;CNV;INTERVAL;KB;CHR;MID_BP;TARGETS;NUM_TARG;MLCN;Q_SOME"}''{print "\"SAMPLE\";\"DUP\";\""$2":"$3"-"$4"\";"$4-$3";"($3+$4)/2";"$2";\""$1".."$1"\";1;3;99"}' > /storage/crihan-msa/DATA/ALZ_201508_CNG/READCOUNT/ABeta2015.target.csv
