#!/usr/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#@ #SBATCH --time=24:00:00
#@ #SBATCH --qos=medium
#SBATCH --mem=10g

hostname

set -u
###################################################################################################
#extract variables
VARI=$1
splitVARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$splitVARI"
echo $1 | tr ',' '\n'

TIME=$(date "+%s")
source ${SCRIPTdir}tools
LC_ALL=C

printf "start of pipeline\n\n" >${OPENdir}time-log.txt

###################################################################################################
#setup-phase

TMPdirRAW=$TMPdir

#load tools
source ${SCRIPTdir}tools

#define thread-variable
THREADS=$(($SLURM_CPUS_PER_TASK * 2))

###################################################################################################
#preset hub

#extract name used for assembly
ASSEMBLYversion=$(cat ${ASSEMBLYhub}hub.txt | head -n 1 | awk '{ print $NF}')
VARI="${VARI},ASSEMBLYversion=${ASSEMBLYversion}"

#---------------------------------------------------------------------------------------------
#prepare hub for visualization of TE related tracks
if [[ ! -d ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/ || $FORCE == Y || $FORCEhub == Y ]]; then

  #do not continue until file is unblocked by other process
  while [[ -f ${ASSEMBLYhub}wait.txt ]]; do
    sleep 10s
  done

  #block trackDb from other processes
  touch ${ASSEMBLYhub}wait.txt

  #create folder for files
  rm -rf ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/
  mkdir ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/

  #remove entry from all files
  cp ${ASSEMBLYhub}groups.txt ${ASSEMBLYhub}groups.txt.backup
  awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="map-TEs" '
  {
    if( $0 !~ NAME ) print
  }' ${ASSEMBLYhub}groups.txt >${TMPdir}groups.tmp
  mv ${TMPdir}groups.tmp ${ASSEMBLYhub}groups.txt

  cp ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt.backup
  awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="map-TEs" '
  {
    if( $0 !~ NAME ) print
  }' ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt >${TMPdir}trackDb.tmp
  mv ${TMPdir}trackDb.tmp ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt

  #add AP-hub into groups.txt
  printf "name map-TEs\nlabel map-TEs\npriority 1\n defaultIsClosed 0\n\n" >>${ASSEMBLYhub}groups.txt

  #unblock hub for other tools
  rm -rf ${ASSEMBLYhub}wait.txt
fi

if [[ ! -s ${TMPdir}assembly.fa ]]; then
  seqkit fx2tab ${assemblyFASTA} |
    awk -v OFS="\t" -v TMP=$TMPdir '{
      old=$1
      gsub(":","_",old)
      split($1,splitNAME,/:/)
      $1=splitNAME[1]
      print $1,old > TMP "name-conversion.txt"
      print ">"$1"\n"$NF
    }'  > ${TMPdir}assembly.fa
fi

  seqkit fx2tab --length -n ${assemblyFASTA} | 
  awk -v OFS="\t" '{print $1,$NF}' | tr ':' '_' |
    LC_COLLATE=C sort -k1,1 >${TMPdir}chrom.size

assemblyFASTA=${TMPdir}assembly.fa
VARI=$VARI,assemblyFASTA=${TMPdir}assembly.fa

###################################################################################################
#prepare input nanopore reads

if [[ ! -s ${TMPdir}reads.fa || $FORCE == Y ]]; then
  #copy over reads used for scaffolding
  seqkit seq --min-len 10000 $nanoporeFASTQ | seqkit fq2fa --line-width 0 >${TMPdir}reads.fa

  COMMAND="${SCRIPTdir}NanoPlot.sh"
  VARI="${VARI},STEP=filtered,INFILE=${TMPdir}reads.fa"
  sbatch $COMMAND ${VARI}
fi


###################################################################################################
###################################################################################################
###################################################################################################
#map reads against the genome using minimap2
INFILE=${TMPdir}reads.fa

#?@ #! remove head to allow processing of full files
#?@ head -n 40000 ${TMPdir}reads.fastq >${TMPdir}combined_red.fastq
#?@ INFILE=${TMPdir}combined_red.fastq
###################################################################################################
###################################################################################################
###################################################################################################

#---------------------------------------------------------------------------------------------
#create minimap2 index
if [[ ! -s ${TMPdir}ref.mmi ]]; then
  minimap2 -x map-ont -t $THREADS -d ${TMPdir}ref.mmi ${assemblyFASTA}
fi

###################################################################################################
#map TEs to reads and generate flanking regions
COMMAND="${SCRIPTdir}identify_and_classify_TEs.sh"
VARI="${VARI},INFILE=$INFILE"

if [[ ! -s ${TMPdir}insertions_to_process.txt ]]; then
  if [[ $COMPUTING == C ]]; then
    sbatch --wait $COMMAND ${VARI}
  else
    $COMMAND ${VARI}
  fi
fi

###################################################################################################
#insert TEs into the genome

nMissingTEs=$(cat ${TMPdir}insertions_to_process.txt | wc -l)
COMMAND="${SCRIPTdir}insertTEs.sh"
VARI="${VARI},INFILE=$INFILE"


if [[ $COMPUTING == C || $multiCLUSTER == Y ]]; then
  rm -rf ${TMPdir}insert_TEs/*
  #nMissingTEs=20
  sbatch --array=1-$nMissingTEs --wait $COMMAND ${VARI}
else
  if [[ $multiLOCAL == Y ]]; then
    nMissingTEs=10
    rm -rf ${TMPdir}insert_TEs/*
    for i in $(seq 1 $nMissingTEs); do
      $COMMAND ${VARI},SLURM_ARRAY_TASK_ID=$i | tee ${OPENdir}insertTEs.log.${i}.txt 
    done
  else
    nMissingTEs=1
    $COMMAND ${VARI},SLURM_ARRAY_TASK_ID=230 | tee ${OPENdir}insertTEs.log.txt 
  fi
fi


#assemble logs into big table
rm -rf  ${TMPdir}all.log
rm -rf ${TMPdir}headers.txt
for i in $(seq 1 $nMissingTEs); do
  sort -k1,1n ${TMPdir}insert_TEs/${i}/log.txt  | sed 's/=/\t/g;s/:-:/\t/g' | cut -f 3 | tr '\n' '\t' >>${TMPdir}all.log 
  printf "\n" >>${TMPdir}all.log
done

sed 's/=/\t/g' ${TMPdir}insert_TEs/1/log.txt | cut -f 1 |sort -k1,1n  | uniq | sed 's/:-:/\t/g' | cut -f 2 | tr '\n' '\t' > ${OPENdir}TEinsertions.log 
printf "\n" >> ${OPENdir}TEinsertions.log
cat ${TMPdir}all.log  >> ${OPENdir}TEinsertions.log


tail -n +2 ${OPENdir}TEinsertions.log | cut -f 1  | sort | uniq -c | tee ${OPENdir}result-summary.txt


###################################################################################################

printf '##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##ALT=<ID=INS,Description="Insertion">
##FORMAT=<ID=STATUS,Number=1,Type=String,Description="PredictionStatus">
' > ${TMPdir}toInsert.vcf

seqkit fx2tab --name --length $assemblyFASTA |
  awk -v OFS="\t" '{print $1,$2}' |
  LC_COLLATE=C sort -k1,1  > ${TMPdir}chrom.size.new

awk '{
  print "##contig=<ID="$1",length="$2">"
}' ${TMPdir}chrom.size.new >> ${TMPdir}toInsert.vcf

printf '#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT HapTEs
' | tr ' ' '\t' >> ${TMPdir}toInsert.vcf

rm -rf ${TMPdir}toInsert.tmp
for i in $(seq 1 $nMissingTEs); do  
  if [[ -s ${TMPdir}insert_TEs/${i}/out.vcf ]]; then
    cat ${TMPdir}insert_TEs/${i}/out.vcf | tr ' ' '\t' >> ${TMPdir}toInsert.tmp
  fi
done

sort -k1,1 -k2,2n ${TMPdir}toInsert.tmp >> ${TMPdir}toInsert.vcf
bgzip --force ${TMPdir}toInsert.vcf
bcftools index ${TMPdir}toInsert.vcf.gz
bcftools consensus  --fasta-ref ${assemblyFASTA} --chain ${OPENdir}chain.txt --output ${OPENdir}${assemblyVERSION}.TEfilled.fa ${TMPdir}toInsert.vcf.gz

###################################################################################################
#run annotation script

#start the annotation pipeline for the plished assembly
COMMAND="annotate_assembly.sh -A ${OPENdir}${assemblyVERSION}.TEfilled.fa -N ${assemblyVERSION}_TEinserted -r ${INFILE} -DHF"
echo $COMMAND
eval $COMMAND

#################################################################################################
#collect flank-sequences for TEs that cannot be inserted

rm -rf ${TMPdir}not-inserted-flanks.fa
for i in $(seq 1 $nMissingTEs); do
  if [[ ! -s ${TMPdir}insert_TEs/$i/out.vcf ]]; then
    echo $i
    seqkit fx2tab ${TMPdir}insert_TEs/$i/flanks.fa | 
      awk -v OFS="\t" -v ID=$i '{print ">"ID"_"$1"\n"$NF}' >> ${TMPdir}not-inserted-flanks.fa
  fi
done

#################################################################################################
#wait for hub to be ready

SWITCH=N
ASSEMBLYdir=${OPENdir}annotate-assembly_${assemblyVERSION}_TEinserted/UCSC/

while [[ $SWITCH == N ]]; do
  if [[ -d ${ASSEMBLYdir} ]]; then
    sleep 60s
    SWITCH=Y
  else
    sleep 60s
  fi
done

###################################################################################################
#mark inserted TEs

awk -v OFS="\t" -v INFILE=${TMPdir}name-conversion.txt -v RESULTfile=${TMPdir}xy.tmp '
BEGIN{
    while((getline LINE < INFILE) > 0) {
      split(LINE,splitLINE,/ |\t/)
      NAME[splitLINE[2]]=splitLINE[1]
    }
    while((getline LINE < RESULTfile) > 0) {
      split(LINE,splitLINE,/ |\t/)
      INSERTION[splitLINE[3]]=splitLINE[1]
    }
}
{
  if($4 in INSERTION){
    $1=NAME[$1]
    $2=$2-1
    print
  }
}' ${TMPdir}map_TEs_to_reads/insertions.merged.bed | LC_COLLATE=C sort -k1,1 -k2,2n > ${TMPdir}insertions.old.bed

awk -v OFS="\t" '{
  if($0~"chain"){
    CHR=$3 
    LENGTH=$12
    SWITCH="Y"
    POS=0
  }else{
    if( $0 !=""){
      if(NF==3){
        if(SWITCH==Y){
          print CHR,$1,$1+$3,"TEinserted",0,"+"
          POS=$1+$3
          SWITCH="N"
        }else{
          print CHR,POS+$1,POS+$1+$3,"TEinserted",0,"+"
          POS=POS+$1+$3
        }
      }
    }
  }
}' ${OPENdir}chain.txt | LC_COLLATE=C sort -k1,1 -k2,2n > ${TMPdir}insertedTEs.bed

bedToBigBed ${TMPdir}insertedTEs.bed ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/${assemblyVERSION}_TEinserted.chrom.sizes ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/insertedTEs.bb

#remove entries from hub
while [[ -f ${ASSEMBLYdir}wait.txt ]]; do
  sleep 10s
done
#block hub for other processes
touch ${ASSEMBLYdir}wait.txt

awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" '{
    if($0 !~ "track insertedTEs") print
}' ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/trackDb.txt >${ASSEMBLYdir}${assemblyVERSION}_TEinserted/trackDb.tmp
mv ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/trackDb.tmp ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/trackDb.txt
  rm -rf ${ASSEMBLYhub}wait.txt

#add track to hub
printf "track insertedTEs\ntype bigBed\nshowNames on\nmaxWindowToDraw 10000000\nvisibility pack\nbigDataUrl insertedTEs.bb\nshortLabel insertedTEs\nlongLabel TEs inserted into the assembly via my pipeline \n\n" >>${ASSEMBLYdir}${assemblyVERSION}_TEinserted/trackDb.txt

rm -rf ${ASSEMBLYdir}wait.txt

###################################################################################################
#mark TEs not inserted
minimap2 -ax map-ont --secondary=no -t $THREADS ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/${assemblyVERSION}_TEinserted.fa ${TMPdir}not-inserted-flanks.fa | 
  samtools view -bS | bedtools bamtobed -i -  | LC_COLLATE=C sort -k1,1 -k2,2n > ${TMPdir}not-inserted.bed

 
bedToBigBed ${TMPdir}not-inserted.bed ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/${assemblyVERSION}_TEinserted.chrom.sizes ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/NOTinsertedTEs.bb

#remove entries from hub
while [[ -f ${ASSEMBLYdir}wait.txt ]]; do
  sleep 10s
done
#block hub for other processes
touch ${ASSEMBLYdir}wait.txt

awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" '{
    if($0 !~ "track NOTinsertedTEs") print
}' ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/trackDb.txt >${ASSEMBLYdir}${assemblyVERSION}_TEinserted/trackDb.tmp
mv ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/trackDb.tmp ${ASSEMBLYdir}${assemblyVERSION}_TEinserted/trackDb.txt
  rm -rf ${ASSEMBLYhub}wait.txt

#add track to hub
printf "track NOTinsertedTEs\ntype bigBed\nshowNames on\nmaxWindowToDraw 10000000\nvisibility pack\nbigDataUrl NOTinsertedTEs.bb\nshortLabel NOTinsertedTEs\nlongLabel flanking regions of TEs my pipeline could not insert \n\n" >>${ASSEMBLYdir}${assemblyVERSION}_TEinserted/trackDb.txt

rm -rf ${ASSEMBLYdir}wait.txt

###################################################################################################
###################################################################################################

#@ #---------------------------------------------------------------------------------------------
#@ #map and generate raw tracks
#@ COMMAND="${SCRIPTdir}map_raw_to_genome.sh"
#@ VARI="${VARI},INFILE=$INFILE"

#@ if [[ $COMPUTING == C ]]; then
#@   sbatch $COMMAND ${VARI}
#@ else
#@   $COMMAND ${VARI}
#@ fi
