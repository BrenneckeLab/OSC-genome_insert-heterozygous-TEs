###############################################################################################
set -u

# Argument = -i input -c chunksize -D blast-database -v
usage() {
  cat <<EOF
  usage: $0 options
  
  ###############################################################################
  This tool calculates read counts on a template of choice. It reportes the 
  counts for both, the sense and the antisense strand.
  
  usage: [PATH]/selectUTRs [options] -F 
  
  OPTIONS:
      -h  Show this message
      -A  assembly fasta file
      -V  genome assembly version
      -N  assembly name
            [default is name of fasta-file used as assembly-fasta]
      -n  nanopore reads fastq
          any fastq file from the previous steps would be sufficient
          I mostly use the one from polishing
      -m  minimum detected TE fragment size [ default=500 ]
      -g  maximum gap between TE mappings for merging them [ default=100 ]
      -b  size of flanking regions to be mapped [ default=1000 ]
      -a  path to the hub-folder containing the hub.txt file
      -T  location of the TE consensus file
      -O  supply output folder
            [default is a new subfolder in the directory containing assembly fasta]
      -t  supply path for base-temporary directory for non-permanent storage
      -x  compute all steps in local mode - if not set some steps will be skipped for faster debug-runs
      -M  set flag for local processing of the first 20 insertions (use only if multiple cores available)
      -R  set flag for cluster processing of insertions even if main run local (use only if multiple cores available)
      -C  set flag for local processing (use only if multiple cores available)
      -D  sed debug mode - does not trigger git commit
      -F  force delete everything before restarting the script
      -H  force re-generation of the hub
EOF
}

assemblyFASTA=
assemblyVERSION=
assemblyNAME=
nanoporeFASTQ=
ASSEMBLYhub=
minTE=500
TEgap=500
blockSIZE=1000
TEconsensus="TE_annot_new_WO_DUST.fa"
OPENdir=
TMPdir="/scratch/brennecke/handler/nanopore/map-TEs/"
multiLOCAL=N
multiCLUSTER=N
COMPUTING=C
DEBUG=N
ALL=N
FORCE=
FORCEhub=

while getopts ÒhA:V:N:n:a:m:g:b:T:O:t:xMRCDFH,Ó OPTION; do
  case $OPTION in
  h)
    usage
    exit 1
    ;;
  A)
    assemblyFASTA=$OPTARG
    ;;
  V)
    assemblyVERSION=$OPTARG
    ;;
  N)
    assemblyNAME=$OPTARG
    ;;
  n)
    nanoporeFASTQ=$OPTARG
    ;;
  a)
    ASSEMBLYhub=$OPTARG
    ;;
  m)
    minTE=$OPTARG
    ;;
  g)
    TEgap=$OPTARG
    ;;
  b)
    blockSIZE=$OPTARG
    ;;
  T)
    TEconsensus=$OPTARG
    ;;
  O)
    OPENdir=$OPTARG
    ;;
  t)
    TMPdir=$OPTARG
    ;;
  x)
    ALL=Y
    ;;
  M)
    multiLOCAL=Y
    ;;
  R)
    multiCLUSTER=Y
    ;;
  C)
    COMPUTING=L
    ;;
  D)
    DEBUG=Y
    ;;
  F)
    FORCE=Y
    ;;
  H)
    FORCEhub=Y
    ;;
  ?)
    usage
    exit
    ;;
  esac
done

###################################################################################################
#test supplied variables or set default values

#test assemblyFASTA
if [[ ! -f $assemblyFASTA ]]; then
  printf "
  please supply valid path to assembly fasta file
  the supplied file $assemblyFASTA does not exist
  \n"
  exit
fi

if [[ -z $assemblyVERSION ]]; then
  #usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "       Please provide assembly version in option V!\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
  exit
fi

#determine assemblyNAME from assemblyFASTA if no value supplied
if [[ -z $assemblyNAME ]]; then
  assemblyNAME=$(basename $assemblyFASTA)
fi

#test nanoporeFASTQ
if [[ ! -s $nanoporeFASTQ ]]; then
  printf "
  please supply valid path to nanopore read fasta file
  the supplied file $nanoporeFASTQ does not exist
  \n"
  exit
fi

#test Hub file
if [[ ! -d $ASSEMBLYhub ]]; then
  printf "
  please supply valid path to a hub.txt file for the current genome
  the supplied file $ASSEMBLYhub does not exist
  \n"
  exit
fi

#test TE consensus file
if [[ ! -f $TEconsensus ]]; then
  printf "
  please supply valid path to TE consensus fasta file
  the supplied file $TEconsensus does not exist
  \n"
  exit
fi

#set minimal TE fragment to default size
if [[ -z $minTE ]]; then
  minTE=500
fi

#determine OPENdir from assemblyFASTA if no value supplied
if [[ -z $OPENdir ]]; then
  OPENdir=$(dirname $assemblyFASTA)
  OPENdir=${OPENdir}/
fi


#@!@#
#separator for printing of settings to log
###################################################################################################
#variable-setup

OPENdir=${OPENdir}TEintegration_${assemblyVERSION}_$assemblyNAME/
SINGULARITYdir=${TMPdir}singu/
TMPdir=${TMPdir}TEintegration_${assemblyVERSION}_$assemblyNAME/

#wipe directories completely if clean restart requested
if [[ $FORCE == Y ]]; then
  rm -rf $OPENdir
  rm -rf $TMPdir
fi

#create required directories
mkdir -p $OPENdir
mkdir -p $TMPdir
mkdir -p $SINGULARITYdir

###################################################################################################
#preset scripts

#determine script-location$
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done

SCRIPTdirRAW="$(cd -P "$(dirname "$SOURCE")" && pwd)"
SCRIPTdirRAW="${SCRIPTdirRAW}/"
SCRIPTdir="${SCRIPTdirRAW}script-files/"
UTILITYdir="${SCRIPTdirRAW}utility-files/"

#move scripts to TMP-directory
rm -rf ${TMPdir}script-files
cp -r ${SCRIPTdir} ${TMPdir}script-files
SCRIPTdir=${TMPdir}script-files/

mv ${TMPdir}script-files/main.sh ${TMPdir}script-files/mTE_${assemblyNAME}_${minTE}

###############################################################################################

TEconsensus=${UTILITYdir}${TEconsensus}

###############################################################################################
#download all singularity images required

cd ${SINGULARITYdir}
wget -O MARVEL.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2025_Handler_OSC-genome/Apptainer/MARVEL.app
wget -O minimap2.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2025_Handler_OSC-genome/Apptainer/minimap2.app
wget -O R.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2025_Handler_OSC-genome/Apptainer/R.app
wget -O racon.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2025_Handler_OSC-genome/Apptainer/racon.app
#wget -O flye.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2025_Handler_OSC-genome/Apptainer/flye.app
wget -O wtdbg2.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2025_Handler_OSC-genome/Apptainer/wtdbg2.app
wget -O sniffles.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2025_Handler_OSC-genome/Apptainer/sniffles.app
wget -O basicTools.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2025_Handler_OSC-genome/Apptainer/basicTools.app
wget -O repeatmasker.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2025_Handler_OSC-genome/Apptainer/repeatmasker.app
wget -O whatshap.app https://brenneckelab.imba.oeaw.ac.at/Publication_Data/2025_Handler_OSC-genome/Apptainer/whatshap.app

###################################################################################################
#report settings and initiate log-files

cd ${SCRIPTdirRAW}

if [[ $DEBUG != Y ]]; then

  #ask for commit-message
  while true; do
    read -r -p "Plese specify a commit-message: " msg
    case $msg in
    [Nn]) break ;;
    *)
      commitMESSAGE=$msg
      break
      ;;
    esac
  done

  #commit all changes
  git add .
  if [[ -z $commitMESSAGE ]]; then
    git commit -m "automatic commit on submission"
  else
    git commit -m "$commitMESSAGE"
  fi
  #git push --all
fi
commitID=$(git log -1 --pretty=format:"%h")

printf "\n\n############################################################################\n\n" >>"${OPENdir}log.txt"
printf "CommitID= ${commitID}\n\n" >>"${OPENdir}log.txt"
echo ${OPENdir}

###################################################################################################
mkdir -p ${OPENdir}settings/
cat ${SCRIPTdirRAW}*.sh |
  awk -v RS="#@!@#" '{if (NR==1) print }' >${OPENdir}settings/SETTINGS_${commitID}.log

TIME="$(date "+%s")"

#report settings to log
printf "
assemblyFASTA=$assemblyFASTA
assemblyVERSION=${assemblyVERSION}
assemblyNAME=$assemblyNAME
nanoporeFASTQ=$nanoporeFASTQ
TEconsensus=$TEconsensus
minTE=$minTE
OPENdir=$OPENdir
TMPdir=$TMPdir
">>${OPENdir}log.txt

###################################################################################################
#submit main-run script
LOG=${OPENdir}/LOGs/
if [[ $COMPUTING == C ]]; then
  rm -rf ${LOG}/*
fi
mkdir -p ${LOG}
cd $LOG

COMMAND="${SCRIPTdir}/mTE_${assemblyNAME}_${minTE}"
VARI="OPENdir=${OPENdir},assemblyVERSION=${assemblyVERSION},TMPdir=${TMPdir},SINGULARITYdir=${SINGULARITYdir},COMPUTING=${COMPUTING},SCRIPTdir=${SCRIPTdir},UTILITYdir=${UTILITYdir},LOG=${LOG},FORCEhub=${FORCEhub},FORCE=${FORCE},assemblyFASTA=${assemblyFASTA},assemblyNAME=${assemblyNAME},nanoporeFASTQ=${nanoporeFASTQ},ASSEMBLYhub=${ASSEMBLYhub},TEconsensus=${TEconsensus},minTE=${minTE},TEgap=${TEgap},blockSIZE=${blockSIZE},multiLOCAL=${multiLOCAL},multiCLUSTER=${multiCLUSTER},ALL=${ALL}"

if [[ $COMPUTING == C ]]; then
  sbatch "$COMMAND" $VARI
else
  if [[ -z ${SLURM_JOB_ID+x} ]]; then
    srun --cpus-per-task=10 --mem-per-cpu=5g --qos=short $COMMAND ${VARI}
  else
    $COMMAND ${VARI}
  fi
fi

