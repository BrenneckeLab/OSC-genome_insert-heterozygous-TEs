#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=20g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short
#SBATCH --time=1:00:00

hostname
set -u

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"

TIME=$(date "+%s")

###################################################################################################
#setup-phase

#create path variables
locTMP=${TMPdir}fill_gaps/

#create directories
mkdir $locTMP

#load tools
source ${SCRIPTdir}tools

#calculate THREADS
THREADS=$(($SLURM_CPUS_PER_TASK * 2))
echo $THREADS
###################################################################################################
#setup-phase

locTMP=${TMPdir}map_TEs_to_reads/
mkdir -p ${locTMP}
echo $locTMP

###################################################################################################
#identify and characterize TEs in reads

#---------------------------------------------------------------------------------------------
#create minimap2 index
if [[ ! -f ${locTMP}reads.mmi ]]; then
  minimap2 -x map-ont -t $THREADS -d ${locTMP}reads.mmi $INFILE
  minimap2 -x map-ont -t $THREADS -d ${locTMP}TEs.mmi $TEconsensus
fi

#---------------------------------------------------------------------------------------------
#mapping of nanopore reads against the TE index

minimap2 -x map-ont -t $THREADS ${locTMP}TEs.mmi $INFILE |
  sort -k1,1 -k3,3n >${locTMP}mapped-reads_to_TEs.paf

#---------------------------------------------------------------------------------------------
#generate bed-file for TEs including flanking regions in reads
awk -v OFS="\t" -v TEgap=$TEgap -v blockSIZE=$blockSIZE -v locTMP=$locTMP '

function classify_type(nTE, US, DS)
{
  if(nTE == 1 ){
    if(US > 2000 && DS > 2000){
      return "single-distant"
    }else{
      return "single-close"
    }
  }else{
    if(US > 2000 && DS > 2000){
      return "multiple-distant"
    }else{
      return "multiple-close"
    }
  }
}

{
  if(NR==1){
    #read in all required features of the first read
    oldID=$1
    oldLENGTH=$2
    oldSTART=$3
    oldEND=$4
    usEND=0
    TEcontent=$6"-"$11
    BLOCKsize=$11
    TEcounter=1
  }else{
    if($1 == oldID){
      #code if still the same read
      if($3<oldEND+TEgap){
        #code if TE fragments are very close to each other
        
        #add TE id to the TEcontent 
        TEcontent=TEcontent":!:"$6"-"$11

        #add gap and new alignment to blocksize
        BLOCKsize+=($3-oldEND)+$11

        #set oldEND to the end of the new TE
        oldEND=$4
        TEcounter+=1
      }else{
        #print out old TE block and initiate a new block as too distant to previous

        TYPE=classify_type(TEcounter, oldSTART-usEND, $3-oldEND)
        split(oldID,splitNAME,/:/)
        print oldID,oldSTART,oldEND,splitNAME[1]"::"oldSTART"::TYPE="TYPE"::SIZE="BLOCKsize"::usDIST="oldSTART-usEND"::dsDIST="$3-oldEND,0,"+"
        print splitNAME[1],oldSTART,oldEND,TEcontent > locTMP "TEcontent.txt"
        usEND=oldEND
        oldID=$1
        oldLENGTH=$2
        oldSTART=$3
        oldEND=$4
        TEcontent=$6"-"$11
        BLOCKsize=$11
        TEcounter=1
      }
    }else{
      #code if new read is reached
      TYPE=classify_type(TEcounter, oldSTART-usEND, oldLENGTH-oldEND)
      split(oldID,splitNAME,/:/)
      print oldID,oldSTART,oldEND,splitNAME[1]"::"oldSTART"::TYPE="TYPE"::SIZE="BLOCKsize"::usDIST="oldSTART-usEND"::dsDIST="oldLENGTH-oldEND,0,"+"
      print splitNAME[1],oldSTART,oldEND,TEcontent > locTMP "TEcontent.txt"
      usEND=0
      oldID=$1
      oldLENGTH=$2
      oldSTART=$3
      oldEND=$4
      TEcontent=$6"-"$11
      BLOCKsize=$11
      TEcounter=1
    }
  }
}' ${locTMP}mapped-reads_to_TEs.paf  >  ${locTMP}TEs_in_reads.bed

###################################################################################################
#identify TEs already contained in the genome

#generate flanking sequences
awk -v OFS="\t" '{
  n=split($4, splitNAME,/::|=/)
  if(NR==1){
    for(i=1; i<=n; i++){
      if(splitNAME[i] == "usDIST"){usCOL=i+1}
      if(splitNAME[i] == "dsDIST"){dsCOL=i+1}
      if(splitNAME[i] == "BLOCKsize"){BLOCKsizeCOL=i+1}
    }
  }
  n=split($1, splitREAD,/::|=/)
  if(NR==1){
    for(i=1; i<=n; i++){
      if(splitREAD[i] == "length"){lengthCOL=i+1}
    }
  }

  if(splitNAME[usCOL]>2000){US=$2-2000;}else US=$2-splitNAME[usCOL]
  if(splitNAME[dsCOL]>2000){DS=$3+2000;}else DS=$3+splitNAME[dsCOL]

  #only print if enough space on read
  if( US<$2 && splitREAD[lengthCOL]>DS ){
    print $1,US,$2,$4":@:US",0,"+"
    print $1,$3,DS,$4":@:DS",0,"+"
  }
}' ${locTMP}TEs_in_reads.bed >${locTMP}flanking-regions.bed

#---------------------------------------------------------------------------------------------
#map flanking regions to the genome

#convert reads to single line fasta
rm -rf ${locTMP}reads.fa.fa*
seqkit fq2fa --line-width 0 $INFILE >${locTMP}reads.fa

#extract flanking sequences
bedtools getfasta -name -fi ${locTMP}reads.fa -bed ${locTMP}flanking-regions.bed -fo ${locTMP}flanking-regions.fa

#map and output as paf for later filtering
minimap2 -x map-ont -t $THREADS ${TMPdir}ref.mmi ${locTMP}flanking-regions.fa > ${locTMP}flanking_regions.all.paf

#map and convert to bam-file for visualization
minimap2 -ax map-ont -t $THREADS ${TMPdir}ref.mmi ${locTMP}flanking-regions.fa |
  samtools view -bS -F0x900- |
  samtools sort -o ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/flanking_regions.all.bam -@ $THREADS

samtools index ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/flanking_regions.all.bam

#add bam-file to the hub
#remove entries from hub
while [[ -f ${ASSEMBLYhub}wait.txt ]]; do
  sleep 10s
done

#block hub for other processes
touch ${ASSEMBLYhub}wait.txt
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" '{
if($1 !~ "track TE_flanking") print
}' ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt >${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.tmp
mv ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.tmp ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt

#add track to hub
printf "track TE_flanking\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 10000000\nvisibility squish\nbigDataUrl map-TEs/flanking_regions.all.bam\nshortLabel TE-flanking-regions\nlongLabel TE flanking regions of max 2 kb\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt

rm -rf ${ASSEMBLYhub}wait.txt

###################################################################################################
#filter TEs that are already present in the genome

awk -v OFS="\t" -v locTMP=$locTMP '
{
  split($1,truncNAME,/:@:/)
  n=split(truncNAME[1],splitNAME,/::|=/)
  
  #determine array position of TEsize
  if(NR==1){
    for(i=1;i<=n; i++){ if(splitNAME[i]=="SIZE") SIZEcol=i+1}
  }

  if( $4-$3 >$2*0.9 && $12>0){
    X[truncNAME[1]][truncNAME[2]]["START"]=$8
    X[truncNAME[1]][truncNAME[2]]["STOP"]=$9
    X[truncNAME[1]][truncNAME[2]]["LINE"]=$0
    X[truncNAME[1]]["SIZE"]=splitNAME[SIZEcol]
    X[truncNAME[1]]["N"]+=1
  }else{
    if($12==0){
      print "badALN",$0 > locTMP "TE.bad.txt"
    }else{
      print "shortALN",$0 > locTMP "TE.bad.txt"
    }
  }
}
END{
  for(TE in X){
    #test if both upstream and downstream parts are mapped exactly once each
    if("US" in X[TE] &&  "DS" in X[TE] && X[TE]["N"]==2){
      #*maybe test for different contigs
      #test if TE present or not

      #first test for orientation and gather distanze between end
      if(X[TE]["US"]["START"] < X[TE]["DS"]["START"] ){
        DIST=X[TE]["DS"]["START"]-X[TE]["US"]["STOP"]
      }else{
          DIST=X[TE]["US"]["START"]-X[TE]["DS"]["STOP"]
      }

      #classify according to size
      if(DIST < 50 ) {
        print "TEabsent",TE > locTMP "TE.absent.txt"
      }else{
        if(DIST > X[TE]["SIZE"]*0.8 && DIST < X[TE]["SIZE"]*1.2){
          print "TEpresent", TE > locTMP "TE.present.txt"
        }else{
          if(DIST < X[TE]["SIZE"]*0.8){
            print "GAPshort", TE,DIST > locTMP "TE.bad.txt"
          }else{
            print "GAPlong", TE,DIST > locTMP "TE.bad.txt"
          }
        }
      }
    }else{
      print "not2",TE > locTMP "TE.bad.txt"
    }
  }
}' ${locTMP}flanking_regions.all.paf


#add bam-file to the hub
#remove entries from hub
while [[ -f ${ASSEMBLYhub}wait.txt ]]; do
  sleep 10s
done

touch ${ASSEMBLYhub}wait.txt
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" '{
    if($1 !~ "track TE_classified-flanking") print
}' ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt >${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.tmp
mv ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.tmp ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
  rm -rf ${ASSEMBLYhub}wait.txt

#loop through the categories and filter alignments
for TYPE in absent bad present; do
  echo $TYPE

    sort -k2,2 ${locTMP}TE.${TYPE}.txt | uniq -f 1 | tr -s " " >${locTMP}TE.${TYPE}.sort.txt

  samtools view -h ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/flanking_regions.all.bam |
    awk -v OFS="\t" -v INFILE=${locTMP}TE.${TYPE}.sort.txt '
    BEGIN{
      while((getline LINE < INFILE) > 0) {
        split(LINE,splitLINE,/ |\t/)
        X[splitLINE[2]]=splitLINE[1]
      }
    }
    {
      if($1~"^@"){
        print
      }else{
        split($1,splitNAME,/:@:/)
        if(splitNAME[1] in X){
          split(splitNAME[1],Y,/:/)
          $1=X[splitNAME[1]]":_:"$1
          print $0
        }
      }
    }' | samtools view -bS - >${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/flanking_regions.${TYPE}.bam

  samtools index ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/flanking_regions.${TYPE}.bam

  #add bam-file to the hub
  #remove entries from hub
  while [[ -f ${ASSEMBLYhub}wait.txt ]]; do
    sleep 10s
  done

  #block hub for other processes
  touch ${ASSEMBLYhub}wait.txt

  #add track to hub
  printf "track TE_classified-flanking_${TYPE}\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 10000000\nvisibility squish\nbigDataUrl map-TEs/flanking_regions.${TYPE}.bam\nshortLabel TE-flanking-${TYPE}\nlongLabel TE flanking regions of max 2 kb classified as ${TYPE}\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt

  rm -rf ${ASSEMBLYhub}wait.txt

done

###################################################################################################
#create locations to be replaced
CUTOFF=4


bedtools bamtobed -i ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/flanking_regions.absent.bam >${locTMP}flanking_regions.absent.bed
bedtools genomecov -bg -g ${TMPdir}chrom.size -i ${locTMP}flanking_regions.absent.bed |
  bedtools merge -d 20 -c 4,4 -o max,mean |
  awk -v OFS="\t" -v CUTOFF=$CUTOFF '
  { 
    if($4>CUTOFF ){ 
      print $1,$2,$3,"MAX="$4"_AVG="$5,0,"+"
      
    }
  }' |
  LC_COLLATE=C sort -k1,1 -k2,2n >${locTMP}missingTEs.bed

bedToBigBed ${locTMP}missingTEs.bed ${TMPdir}chrom.size ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/missingTEs.final.bb

bedtools bamtobed -i ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/flanking_regions.absent.bam >${locTMP}flanking_regions.absent.bed
bedtools genomecov -bga -g ${TMPdir}chrom.size -i ${locTMP}flanking_regions.absent.bed | 
  awk -v OFS="\t" -v CUTOFF=$CUTOFF '
  BEGIN{
    MAXval=0
    MAXpos=""
  }
  {
    if($4!=0){
      if(MAXval==0 || MAXval<$4){
        MAXval=$4
        MAXpos=$0
      }
    }else{
      if(MAXval>CUTOFF){
        print MAXpos
        MAXval=0
        MAXpos=""
      }
    }
  }'|
  LC_COLLATE=C sort -k1,1 -k2,2n >${locTMP}missingTEs.max.bed
  
bedToBigBed ${locTMP}missingTEs.max.bed ${TMPdir}chrom.size ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/missingTEs.final.max.bb


#add bam-file to the hub
#remove entries from hub
while [[ -f ${ASSEMBLYhub}wait.txt ]]; do
  sleep 10s
done

#block hub for other processes
touch ${ASSEMBLYhub}wait.txt

awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" '{
  if($1 !~ "track missingTEs.final") print
}' ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt >${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.tmp
mv ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.tmp ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt

#add track to hub
printf "track missingTEs.final\ntype bigBed\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 10000000\nvisibility pack\nbigDataUrl map-TEs/missingTEs.final.bb\nshortLabel missingTEs\nlongLabel regions identified with missing TEs\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
printf "track missingTEs.final.max\ntype bigBed\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 10000000\nvisibility pack\nbigDataUrl map-TEs/missingTEs.final.max.bb\nshortLabel missingTEs.max\nlongLabel regions identified with missing TEs - max position\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt

rm -rf ${ASSEMBLYhub}wait.txt

exit

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################


###################################################################################################
PROCESSED_TIME=$(echo -e $(date "+%s") "$TIME" | mawk '{ print ($1-$2)/60 }')
echo "identification and classification of TEs in reads and genome - processing_time=" "${PROCESSED_TIME}" >>"${OPENdir}time-log.txt"

exit
