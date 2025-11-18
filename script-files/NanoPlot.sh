#!/usr/bin/bash

#SBATCH --cpus-per-task=5
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"   
#SBATCH --qos=short
#SBATCH --time=2:00:00
#SBATCH --mem=4g

hostname
set -u

###################################################################################################
#extract variables 
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g' )
eval "$VARI"

TIME=$(date "+%s")
  printf "start of processing \n\n"> "${OPENdir}time-log.txt"

###################################################################################################
#setup-phase

source ${SCRIPTdir}tools

###################################################################################################


if [[ $STEP == filtered ]] || [[ $STEP == raw ]]
then
  mkdir -p ${OPENdir}prepare_reads/NanoPlot-${STEP}
  NanoPlot -t $SLURM_CPUS_PER_TASK --loglength --readtype 1D -o ${OPENdir}prepare_reads/NanoPlot-${STEP} --format png --fastq $INFILE  
fi


PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | awk '{ print ($1-$2)/60 }' )
printf "NanoPlot - processing-time=	${PROCESSED_TIME} \n\n">> "${OPENdir}time-log.txt"
TIME="$(date "+%s")"
