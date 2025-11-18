#!/usr/bin/bash

#SBATCH --cpus-per-task=14
#SBATCH -e "%x.e.%j-%a.txt"
#SBATCH -o "%x.o.%j-%a.txt"
#SBATCH --qos=medium
#SBATCH --time=24:00:00
#SBATCH --mem=30g

TMPDIR=${SCRATCHDIR}

set -u
HOST=$(hostname)
echo $HOST

###################################################################################################
#extract variables if computing on PIWI

#check if host is piwi and process variables accordingly
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"

if [[ $VARI == *"COMPUTING=C"* ]]; then
  echo ${JOB_ID+x} >>"${TMPdir}jobIDs.txt"
fi

TIME=$(date "+%s")

###################################################################################################
#setup-phase

locTMP=${TMPdir}map_reads/
mkdir -p ${locTMP}

###################################################################################################
#map reads to assembled genome using minimap2

CORES=$(echo $(($SLURM_CPUS_PER_TASK - 1)))
echo $CORES cores

#---------------------------------------------------------------------------------------------
#mapping with default settings
${SINGULARITYdir}minimap2.simg minimap2 -ax map-ont -t $CORES ${TMPdir}ref.mmi $INFILE |
samtools view -bS - |
samtools sort -o ${locTMP}mapped-reads_default.bam -@ $CORES

cp ${locTMP}mapped-reads_default.bam ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_default.bam
samtools index ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_default.bam

#---------------------------------------------------------------------------------------------
#filtering out secondary and hybrid alignments

samtools view -h -F0x900 ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_default.bam |
  samtools view -bS > ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_no2nd_noChim.bam

samtools index ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_no2nd_noChim.bam
#---------------------------------------------------------------------------------------------
#mapping suppressing 2nd mappings
${SINGULARITYdir}minimap2.simg minimap2 --secondary=no -ax map-ont -t $CORES ${TMPdir}ref.mmi $INFILE |
samtools view -bS - |
samtools sort -o ${locTMP}mapped-reads_no2nd.bam -@ $CORES

cp ${locTMP}mapped-reads_no2nd.bam ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_no2nd.bam
samtools index ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_no2nd.bam

#do not continue until file is unblocked by other process
while [[ -f ${ASSEMBLYhub}wait.txt ]]; do
  sleep 10s
done

#block hub for other processes
touch ${ASSEMBLYhub}wait.txt

awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" '{
if($1 !~ "track map-TEs_RAW" && $1!~ "track map-TEs_mod") print 
}' ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt >${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.tmp
mv ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.tmp ${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt

printf "track map-TEs_RAW\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 10000000\nvisibility squish\nbigDataUrl map-TEs//mapped-reads_default.bam\nshortLabel RAWmapping\nlongLabel mapping of raw reads \ngroup map-TEs\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
printf "track map-TEs_RAW_no2ndChim\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 10000000\nvisibility squish\nbigDataUrl map-TEs//mapped-reads_no2nd_noChim.bam\nshortLabel no2ndChim_RAWmapping\nlongLabel mapping of raw reads filtering for 2nd or chimeric mappings\ngroup map-TEs\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
printf "track map-TEs_RAW_no2nd\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 10000000\nvisibility squish\nbigDataUrl map-TEs//mapped-reads_no2nd.bam\nshortLabel no2nd_RAWmapping\nlongLabel mapping of raw reads not allowing 2nd alignments\ngroup map-TEs\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt

rm -rf ${ASSEMBLYhub}wait.txt

###################################################################################################
#mapping with changed settings
${SINGULARITYdir}minimap2.simg minimap2 -ax map-ont -r 50 --no-long-join -K 2G -t $CORES ${TMPdir}ref.mmi $INFILE |
samtools view -bS - |
samtools sort -o ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_mod.bam -@ $CORES

samtools index ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_mod.bam

#---------------------------------------------------------------------------------------------
#filtering out secondary and hybrid alignments

samtools view -h -F0x900 ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_mod.bam |
  samtools view -bS > ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_mod_no2nd_noChim.bam

samtools index ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_mod_no2nd_noChim.bam

#---------------------------------------------------------------------------------------------
#do not continue until file is unblocked by other process
while [[ -f ${ASSEMBLYhub}wait.txt ]]; do
  sleep 10s
done

#block hub for other processes
touch ${ASSEMBLYhub}wait.txt

printf "track map-TEs_RAW_mod\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 10000000\nvisibility squish\nbigDataUrl map-TEs//mapped-reads_mod.bam\nshortLabel modRAWmapping\nlongLabel mapping of raw reads noLongJoing\ngroup map-TEs\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
printf "track map-TEs_RAW_mod_no2ndChim\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 10000000\nvisibility squish\nbigDataUrl map-TEs//mapped-reads_mod_no2nd_noChim.bam\nshortLabel modNo2ndChim_RAWmapping\nlongLabel mapping of raw reads filtering for 2nd or chimeric mappings - noLongJoin\ngroup map-TEs\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt

rm -rf ${ASSEMBLYhub}wait.txt

###################################################################################################
#generate coverage tracks

MEM=$(scontrol show job $SLURM_JOBID | awk '{ if($0~"TRES"){split($1,X,/,|=/); if( X[5]~"G"){sub("G","",X[5]); if(X[5]>20) print X[5]-10;else print 5; }else{print 20 }}}')

bedtools bamtobed -i ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/mapped-reads_mod_no2nd_noChim.bam |
  LC_COLLATE=C sort -S${MEM}G --parallel=$SLURM_CPUS_PER_TASK -k1,1 -k2,2n >${locTMP}mapped-reads_no2nd-sort.bed

#trimm of starting and last nucleotide for proper fraction calculation with end-tracks later on
awk -v OFS="\t" '
{
  print $1,$2+1,$3-1,$4"_START",$5,$6
}
' ${locTMP}mapped-reads_no2nd-sort.bed >${locTMP}mapped-reads_no2nd-sort_trimmed.bed

bedtools genomecov -bg -i ${locTMP}mapped-reads_no2nd-sort_trimmed.bed -g ${TMPdir}chrom.size >${locTMP}full-coverage.bg
bedGraphToBigWig ${locTMP}full-coverage.bg ${TMPdir}chrom.size ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/full-coverage.bw

#block hub for other processes
touch ${ASSEMBLYhub}wait.txt
printf "track map-TEs_RAW_fullCov\ntype bigWig\nautoScale on\nvisibility full\nalwaysZero on\nbigDataUrl map-TEs/full-coverage.bw\nshortLabel FULLcoverage\nlongLabel coverage taking the full read-length into consideration\n\ngroup map-TEs\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
rm -rf ${ASSEMBLYhub}wait.txt

#---------------------------------------------------------------------------------------------
#generate end-coverage
for SITE in start end; do
  #restrict mapping to first and last nucleotide
  awk -v OFS="\t" -v SITE=$SITE '
  {
    if(SITE=="start"){ 
      print $1,$2,$2+1,$4"_START",$5,$6
    }else{
      print $1,$3-1,$3,$4"_END",$5,$6
    }
  }
  ' ${locTMP}mapped-reads_no2nd-sort.bed >${locTMP}mapped-ends.bed

  bedtools genomecov -bg -i ${locTMP}mapped-ends.bed -g ${TMPdir}chrom.size >${locTMP}end-coverage_${SITE}.bg
  bedGraphToBigWig ${locTMP}end-coverage_${SITE}.bg ${TMPdir}chrom.size ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/end-coverage_${SITE}.bw

  touch ${ASSEMBLYhub}wait.txt
  printf "track map-TEs_RAW_endCov_${SITE}\ntype bigWig\nautoScale on\nvisibility full\nalwaysZero on\nbigDataUrl map-TEs/end-coverage_${SITE}.bw\nshortLabel ENDcoverage_${SITE}\nlongLabel coverage taking only alignment ${SITE} ends into consideration\n\ngroup map-TEs\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
  rm -rf ${ASSEMBLYhub}wait.txt

  #---------------------------------------------------------------------------------------------
  #generate END read
  bedtools unionbedg -i ${locTMP}full-coverage.bg ${locTMP}end-coverage_${SITE}.bg |
    awk -v OFS="\t" '
    {
      print $1,$2,$3,$5*100/($4+$5)
    }
    ' >${locTMP}coverage-fraction.bg

  bedGraphToBigWig ${locTMP}coverage-fraction.bg ${TMPdir}chrom.size ${ASSEMBLYhub}${ASSEMBLYversion}/map-TEs/coverage-fraction_${SITE}.bw

  touch ${ASSEMBLYhub}wait.txt
  printf "track map-TEs_RAW_endFrac_${SITE}\ntype bigWig\nautoScale on\nvisibility full\nalwaysZero on\nbigDataUrl map-TEs/coverage-fraction_${SITE}.bw\nshortLabel ENDfract_${SITE}\nlongLabel fraction of reads with ${SITE} ends \n\ngroup map-TEs\n\n" >>${ASSEMBLYhub}${ASSEMBLYversion}/trackDb.txt
  rm -rf ${ASSEMBLYhub}wait.txt
done
exit
###################################################################################################
PROCESSED_TIME=$(echo -e $(date "+%s") "$TIME" | mawk '{ print ($1-$2)/60 }')
echo "mapping of raw-reads to the genome - processing_time=" "${PROCESSED_TIME}" >>"${OPENdir}time-log.txt"

exit

#need to revisit the conversion to bam -- potentially it is better to rename the reads to indicate right away which portion of the read is contained in the mapping... then it would be back map each mapping portion to the TEs and filter TE only mappings.
#	maybe I can also use the mapping number for the fragment, but I don't know how varied splittings for the individual reads can be
#---------------------------------------------------------------------------------------------
#generate bed-files

betools bamToBed -i ${HUBdir}dm6/mapped-reads_no2nd.bam >${locTMP}mapped_raw.bed

genomeCoverageBed -split -ibam ${HUBdir}dm6/mapped-reads_no2nd.bam -strand + -5 -bg -g $chromLENGTH &>${locTMP}error.txt | head
