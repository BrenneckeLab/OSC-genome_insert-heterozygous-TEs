#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --mem=40g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short
#SBATCH --time=4:00:00

hostname
set -ux

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

#map to genome with minimap
#?tried ngmlr but it is way to slow to be feasable
minimap2 -ax map-ont --MD -t $THREADS $assemblyFASTA $INFILE > ${locTMP}mapped.sam

#sort and output as bam-file
samtools sort -@ $THREADS -O bam -T $locTMP -o ${locTMP}mapped.sort.bam ${locTMP}mapped.sam

#analyze mappings with sniffles
sniffles -n -1 -l 100 --tmp_file $locTMP/sniffles.tmp -m ${locTMP}mapped.sort.bam -v ${locTMP}output.raw.vcf 
#-d 10000 -l 100 

#filter problematic SV types
grep -v SVTYPE=BND ${locTMP}output.raw.vcf | LC_COLLATE=C sort -k1,1 -k2,2n >${locTMP}output.vcf 

#phase the SVs
awk -v OFS="\t" '
BEGIN{
  PHASE="1"
  SWITCH="N"
}
{
  if($1~"^#"){
    print 
  }else{
    if(SWITCH == "N"){ 
      CHR=$1
      n=split($8,splitTAG,/;|=/)
      for(i=1; i<=n; i++){ 
        if(splitTAG[i]~"^RNAMES$"){
          m=split(splitTAG[i+1], splitREADS,/,/)
          for(j=1;j<=m;j++){
            X[splitREADS[j]]=splitREADS[j]
            nREADS+=1
          }
          print splitTAG[i+1] > TMPdir "readIDs.txt"
        }
      }
      SWITCH="Y"
    }else{
      n=split($8,splitTAG,/;|=/)
      for(i=1; i<=n; i++){ 
        if(splitTAG[i]~"^RNAMES$"){
          m=split(splitTAG[i+1], splitREADS,/,/)
          for(j=1;j<=m;j++){
            if(splitREADS[j] in X){
              nOVERLAP+=1
            }else{
              nNOT+=1
            }
          }
          delete X
          for(j=1;j<=m;j++){
            X[splitREADS[j]]=splitREADS[j]
          }
        }
      }
      if(nOVERLAP*100/(nOVERLAP+nNOT)<25){
        PHASE=PHASE+1
      }
      print $0,"PHASEall="PHASE, "nOVERLAP:nNOT:Fraction="nOVERLAP":"nNOT":"nOVERLAP*100/(nOVERLAP+nNOT)<25
      nOVERLAP=0
      nNOT=0
    }
  }
}' ${locTMP}output.vcf  > ${locTMP}phased.vcf

#process phase-blocks into genome browser tracks
awk -v OFS="\t" '{
  if($1!~"^#"){
    print
  }
}' ${locTMP}phased.vcf |
awk -v OFS="\t" -v INFILE=${TMPdir}name-conversion.txt '
BEGIN{
    while((getline LINE < INFILE) > 0) {
      split(LINE,splitLINE,/ |\t/)
      NAME[splitLINE[1]]=splitLINE[2]
    }
}
{
  if(NR==1){
    PHASE=$11
    START=$2
    CHR=$1
    LAST=$2
  }else{
    if(CHR==$1){
      if(PHASE==$11){
        LAST=$2
      }else{
        print NAME[CHR],START,LAST,PHASE,0,"+"
        PHASE=$11
        START=$2
        LAST=$2
        CHR=$1
      }
    }else{
      print NAME[CHR],START,LAST,PHASE,0,"+"
      PHASE=$11
      START=$2
      CHR=$1
      LAST=$2
    }
  }
}' | LC_COLLATE=C sort -k1,1 -k2,2n > ${locTMP}phase_blocks_all.bed
bedToBigBed ${locTMP}phase_blocks_all.bed ${TMPdir}chrom.size ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/phase-blocks.allSV.bb


###################################################################################################
#filter insertions for TE content

#create fasta file from VCF 
mawk '{
  if($0~"INS" && length($5)>50){
    print ">"$3"\n"$5
  }
}' ${locTMP}phased.vcf > ${locTMP}insertions.fa

nLINES=$(grep ">" ${locTMP}insertions.fa | wc -l )
echo "number of insertions >50nt = $nLINES\n" | tee ${OPENdir}stats.txt

#run repeatmasker over ther insertions
RepeatMasker -dir ${locTMP}repeatmasker -gff -s -nolow -no_is  -xsmall  -e  ncbi  -lib  $TEconsensus -pa  $THREADS ${locTMP}insertions.fa

rmsk2bed < ${locTMP}repeatmasker/insertions.fa.out | cut -f 1-6 | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools merge > ${locTMP}repeatmasker.out.insertions.bed
awk -v OFS="\t" -v INFILE=${locTMP}repeatmasker.out.insertions.bed -v TMP=$locTMP  '
  BEGIN{
    while((getline LINE < INFILE) > 0) {
      if(LINE !~ "^#"){
        #print LINE
        split(LINE,splitLINE,/ |\t/)
        X[splitLINE[1]]+=splitLINE[3]-splitLINE[2]
      }
    }
    print "ID","TElength","InsertionLength" > TMP "TEhist.INS.txt"
  }
  {
    if($0~"INS" && $1!~"^#" && length($5)>50 ){
      if($3 in X){
        print $3, X[$3], length($5) > TMP "TEhist.INS.txt"
      }else{
        print $3, 0, length($5) > TMP "TEhist.INS.txt"
      }
    }
  }' < ${locTMP}phased.vcf 

Rscript ${SCRIPTdir}plot_TEinsertion_hist.R INPUT=${locTMP}TEhist.INS.txt OPENdir=$OPENdir EXT=INS



#----------------------------------------------------------------
#deletion evaluation 

#create fasta file from VCF 
mawk '{
  if($0~"DEL" && length($4)>50){
    print ">"$3"\n"$4
  }
}' ${locTMP}phased.vcf > ${locTMP}deletions.fa

nLINES=$(grep ">" ${locTMP}insertions.fa | wc -l )
echo "number of deletions >50nt = $nLINES\n" | tee ${OPENdir}stats.txt

#run repeatmasker over ther insertions
RepeatMasker -dir ${locTMP}repeatmasker -gff -s -nolow -no_is  -xsmall  -e  ncbi  -lib  $TEconsensus -pa  $THREADS ${locTMP}deletions.fa

rmsk2bed < ${locTMP}repeatmasker/deletions.fa.out | cut -f 1-6 > ${locTMP}repeatmasker.out.deletions.bed
awk -v OFS="\t" -v INFILE=${locTMP}repeatmasker.out.deletions.bed -v TMP=$locTMP  '
  BEGIN{
    while((getline LINE < INFILE) > 0) {
      if(LINE !~ "^#"){
        #print LINE
        split(LINE,splitLINE,/ |\t/)
        X[splitLINE[1]]=splitLINE[3]-splitLINE[2]
      }
    }
    print "ID","TElength","InsertionLength" > TMP "TEhist.DEL.txt"
  }
  {
    if($0~"DEL" && $1!~"^#" && length($4)>50 ){
      if($3 in X){
        print $3, X[$3], length($4) > TMP "TEhist.DEL.txt"
      }else{
        print $3, 0, length($4) > TMP "TEhist.DEL.txt"
      }
    }
  }' < ${locTMP}phased.vcf 

Rscript ${SCRIPTdir}plot_TEinsertion_hist.R INPUT=${locTMP}TEhist.DEL.txt OPENdir=$OPENdir EXT=DEL


###################################################################################################
#output TE  containing insertions only for insertion script
awk -v OFS="\t" -v INFILE=${locTMP}TEhist.INS.txt -v INFILE2=${TMPdir}name-conversion.txt -v TMP=$locTMP  '
  BEGIN{
    while((getline LINE < INFILE2) > 0) {
      split(LINE,splitLINE,/ |\t/)
      NAME[splitLINE[1]]=splitLINE[2]
    }
    
    while((getline LINE < INFILE) > 0) {
      if(LINE !~ "ID"){
        split(LINE,splitLINE,/ |\t/)
        if(splitLINE[2]*100/ splitLINE[3] > 50){
          X[splitLINE[1]]=splitLINE[2]*100/ splitLINE[3]
        }
      }
    }
  }
  {
    if($0~"INS" && $1!~"^#" && length($5)>50 ){
      if($3 in X){
        print 
        n=split($8,splitTAG,/;|=/)
        for(i=1; i<=n; i++){ 
          if(splitTAG[i]~"^END$"){
            ENDcorr=splitTAG[i+1]
          }
        }
        print NAME[$1],$2-1,ENDcorr,"INS_"$3,0,"+" > TMP "insertions.bed"
      }
    }
    if($0~"DEL" && $1!~"^#" && length($4)>50 ){
      print
      if($3 in X){
        print 
        n=split($8,splitTAG,/;|=/)
      }
    }
  }' ${locTMP}phased.vcf |
  #try to phase SVs containing missing TEs
  awk -v OFS="\t" '
  BEGIN{
    PHASE="1"
    SWITCH="N"
  }
  {
    if($1~"^#"){
      print 
    }else{
      if(SWITCH == "N" && $0~"DEL"){
        x=b
      }else{
        if(SWITCH=="N"){
          CHR=$1
          n=split($8,splitTAG,/;|=/)
          for(i=1; i<=n; i++){ 
            if(splitTAG[i]~"^RNAMES$"){
              m=split(splitTAG[i+1], splitREADS,/,/)
              for(j=1;j<=m;j++){
                X[splitREADS[j]]=splitREADS[j]
                nREADS+=1
              }
            }
          }
          SWITCH="Y"
        }else{
          if($0~"DEL"){
            PHASE=PHASE+1
            nOVERLAP=0
            nNOT=0
          }else{
            n=split($8,splitTAG,/;|=/)
            for(i=1; i<=n; i++){ 
              if(splitTAG[i]~"^RNAMES$"){
                m=split(splitTAG[i+1], splitREADS,/,/)
                for(j=1;j<=m;j++){
                  if(splitREADS[j] in X){
                    nOVERLAP+=1
                  }else{
                    nNOT+=1
                  }
                }
                delete X
                for(j=1;j<=m;j++){
                  X[splitREADS[j]]=splitREADS[j]
                }
              }
            }
            if(nOVERLAP*100/(nOVERLAP+nNOT)<25){
              PHASE=PHASE+1
            }
          
            print $0,"PHASEteIns="PHASE, "nOVERLAP:nNOT:Fraction="nOVERLAP":"nNOT":"nOVERLAP*100/(nOVERLAP+nNOT)<25
            nOVERLAP=0
            nNOT=0
          }
        }
      }
    }
  }' > ${locTMP}insertions_to_process.raw.txt

LC_COLLATE=C sort -k1,1 -k2,2n ${locTMP}insertions.bed > ${locTMP}insertions.sort.bed
bedToBigBed ${locTMP}insertions.sort.bed ${TMPdir}chrom.size ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/missingTEs.bb
 
#process TE based phase-blocks into genome browser tracks
awk -v OFS="\t" -v INFILE=${TMPdir}name-conversion.txt '
BEGIN{
    while((getline LINE < INFILE) > 0) {
      split(LINE,splitLINE,/ |\t/)
      NAME[splitLINE[1]]=splitLINE[2]
    }
}
{
  if(NR==1){
    PHASE=$13
    START=$2
    CHR=$1
    LAST=$2
  }else{
    if(CHR==$1){
      if(PHASE==$13){
        LAST=$2
      }else{
        print NAME[CHR],START,LAST,PHASE,0,"+"
        PHASE=$13
        START=$2
        CHR=$1
        LAST=$2
      }
    }else{
      print NAME[CHR],START,LAST,PHASE,0,"+"
      PHASE=$13
      START=$2
      CHR=$1
      LAST=$2
    }
  }
}' ${locTMP}insertions_to_process.raw.txt | LC_COLLATE=C sort -k1,1 -k2,2n > ${locTMP}phase_blocks_TE.bed
bedToBigBed ${locTMP}phase_blocks_TE.bed ${TMPdir}chrom.size ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/phase-blocks.TEinsertions.bb

###################################################################################################
#merge TE insertions based on distance and phasing 
cat ${locTMP}insertions_to_process.raw.txt | 
awk -v OFS="\t" '{
  if(NR==1){
    CHR=$1
    START=$2
    PHASEall=$11
    PHASEte=$13
    LINE=$0
  }else{
    if( ($1==CHR && $2 >START+2000 ) || CHR != $1 ){
      print LINE
      CHR=$1
      START=$2
      PHASEall=$11
      PHASEte=$13
      LINE=$0
    }else{ 
      if($11==PHASEall && $13==PHASEte){
        #merge reads togeter
        
        #extract reads from last line 
        m=split(LINE,splitLINE,/ |\t/)
        n=split(splitLINE[8],splitTAG,/;/)
        for(i=1; i<=n; i++){ 
          if(splitTAG[i]~"^RNAMES$"){
            m=split(splitTAG[i+1], splitREADS,/,|=/)
            for(j=2;j<=m;j++){
              X[splitREADS[j]]+=1
            }
          }
          if(splitTAG[i]~"^END="){
            split(splitTAG[i],splitSTOP,/=/)
            STOPlast=splitSTOP[2]
          }
        }

        #extract reads from current line 
        n=split($8,splitTAG,/;/)
        for(i=1; i<=n; i++){ 
          if(splitTAG[i]~"^RNAMES$"){
            m=split(splitTAG[i+1], splitREADS,/,|=/)
            for(j=2;j<=m;j++){
              X[splitREADS[j]]+=1
            }
          }
          if(splitTAG[i]~"^END="){
            split(splitTAG[i],splitSTOP,/=/)
            STOPcurr=splitSTOP[2]
          }
        }

        #create new readID tag
        for(i in X){
          if(READidTAG==""){
            READidTAG="RNAMES="i
          }else{
            READidTAG=READidTAG","i
          }
        }

        #replace TAG in splitLINE 
        splitLINE[8]=READidTAG
        READidTAG=""
        delete X 

        #find STOP of current SV

        #merge back splitTAG
        for(i=1; i<=n;i++){
          if(splitTAG[i]~"^END="){
            splitTAG[i]="END="STOPcurr
          }
          if(TAG==""){
            TAG=splitTAG[i]
          }else{
            TAG=TAG";"splitTAG[i]
          }
        } 
        TAG=TAG";merged"

        #create new merged LINE variable 
        LINE=splitLINE[1] "\t" splitLINE[2] "\t" splitLINE[3]":"$3 "\t" splitLINE[4] "\t" splitLINE[5]":"$5 "\t" splitLINE[6] "\t" splitLINE[7] "\t" TAG  "\t" splitLINE[9] "\t" splitLINE[10] "\t" $11 "\t" $12
        #move the START to the last SV to allow merging of another SV
        START=$2

        TAG="" 
        
      }else{
        print LINE
        CHR=$1
        START=$2
        PHASEall=$11
        PHASEte=$13
        LINE=$0
      }
    }
  }
}
END{
  print LINE
}'  >  ${TMPdir}insertions_to_process.txt
 
awk -v OFS="\t" -v INFILE=${TMPdir}name-conversion.txt '
BEGIN{
    while((getline LINE < INFILE) > 0) {
      split(LINE,splitLINE,/ |\t/)
      NAME[splitLINE[1]]=splitLINE[2]
    }
}
{
  n=split($8,splitTAG,/;/)
  for(i=1; i<=n; i++){ 
    if(splitTAG[i]~"^END="){
      split(splitTAG[i],splitSTOP,/=/)
      STOP=splitSTOP[2]
    }
  }
  $1=NAME[$1]
  print $1,$2-1,STOP,$3,0,"+"
}' ${TMPdir}insertions_to_process.txt | LC_COLLATE=C sort -k1,1 -k2,2n > ${locTMP}insertions.merged.bed
bedToBigBed ${locTMP}insertions.merged.bed ${TMPdir}chrom.size ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/missingTEs.merged.bb


###################################################################################################
#add the SVs and the phase-blocks to the hub

#remove entries from hub
while [[ -f ${ASSEMBLYhub}wait.txt ]]; do
  sleep 10s
done
#block hub for other processes
touch ${ASSEMBLYhub}wait.txt

awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" '{
    if($0 !~ "group map-TEs") print
}' ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt >${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.tmp
mv ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.tmp ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
  rm -rf ${ASSEMBLYhub}wait.txt

#add track to hub
printf "track missingTEs\ntype bigBed\nshowNames on\nmaxWindowToDraw 10000000\nvisibility pack\nbigDataUrl map-TEs/missingTEs.bb\nshortLabel missingTEs\nlongLabel TEs (>50\% of insertion) not present in the genome\n group map-TEs \n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
printf "track missingTEs_merged\ntype bigBed\nshowNames on\nmaxWindowToDraw 10000000\nvisibility pack\nbigDataUrl map-TEs/missingTEs.merged.bb\nshortLabel missingTEs_merged\nlongLabel TEs (>50\% of insertion) not present in the genome - merged by phase and distance\n group map-TEs \n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
printf "track phaseBlocks-all\ntype bigBed\nshowNames on\nmaxWindowToDraw 10000000\nvisibility pack\nbigDataUrl map-TEs/phase-blocks.allSV.bb\nshortLabel phaseBlocks-allSV\nlongLabel phase blocks determined by all SVs reported by sniffles\n group map-TEs \n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
printf "track phaseBlocks-TE\ntype bigBed\nshowNames on\nmaxWindowToDraw 10000000\nvisibility pack\nbigDataUrl map-TEs/phase-blocks.TEinsertions.bb\nshortLabel phaseBlocks-TEins\nlongLabel phase blocks determined by reads supporting TE insertions\n group map-TEs \n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
rm -rf ${ASSEMBLYhub}wait.txt

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
