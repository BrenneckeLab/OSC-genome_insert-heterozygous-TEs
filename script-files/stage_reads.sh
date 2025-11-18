#!/usr/bin/bash 

#SBATCH --cpus-per-task=4
#SBATCH -e "%x.e.%j-%a.txt"
#SBATCH -o "%x.o.%j-%a.txt"   
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --mem=2g

TMPDIR=${SCRATCHDIR}

set -u
HOST=$( hostname )
echo $HOST 
  
###################################################################################################
#extract variables if computing on PIWI


#check if host is piwi and process variables accordingly
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g' )
eval "$VARI"
echo $1 | tr ',' '\n'

if [[ $VARI == *"COMPUTING=C"* ]]
then
  echo ${JOB_ID+x} >> "${TMPdir}jobIDs.txt"
fi 

TIME=$(date "+%s")


###################################################################################################
#setup-phase
RUNs=$( echo $subRUNs | tr '~' '\t'  ) 
currRUN=$( echo $RUNs | awk -v N=$SLURM_ARRAY_TASK_ID '{ print $N }' )

locTMP=${TMPdir}prepare_reads/${currRUN}_individual/
rm -rf $locTMP
mkdir -p ${locTMP}

###################################################################################################
#filter reads
FILES=$(ls -d ${RAW}${currRUN}/reads/${BCversion}/untrimmed/* | tr '\n' '\t' )
rm -rf ${locTMP}orig-count.txt 

run_code () {
  SINGULARITYdir=$1
  FILE=$2
  locTMP=$3
  jobID=$4
  
  gunzip -c $FILE  > ${locTMP}${jobID}_in.fastq
  ${SINGULARITYdir}filtlong.img --min_length 10000  ${locTMP}${jobID}_in.fastq > ${locTMP}${jobID}_out.fastq
}


ml parallel/20171122-foss-2017a
export -f run_code
parallel -j $SLURM_CPUS_PER_TASK run_code $SINGULARITYdir {} $locTMP {#}  ::: $FILES


CORES=$( echo $(( $SLURM_CPUS_PER_TASK)))

cat ${locTMP}*_in.fastq > ${locTMP}input_reads.fastq
touch ${locTMP}selected_reads.fastq
cat ${locTMP}*_out.fastq > ${locTMP}selected_reads.fastq

origCOUNT=$(cat ${locTMP}input_reads.fastq | wc -l )
selCOUNT=$(cat ${locTMP}selected_reads.fastq | wc -l )

awk -v selCOUNT=$selCOUNT -v currRUN=${currRUN} -v  OFS="\t" '{ X+=$1} END{ print currRUN,selCOUNT, selCOUNT * 100 / X }' <(echo ${origCOUNT}) >> ${OPENdir}prepare_reads/selected_reads.txt


 

#---------------------------------------------------------------------------------------------------------
#cleanup 

PROCESSED_TIME=$(echo -e $(date "+%s") "$TIME" | mawk '{ print ($1-$2)/60 }' )
echo "stage reads - processing_time="	"${PROCESSED_TIME}" >> "${OPENdir}time-log.txt"
