#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --mem=15g
#SBATCH -e "%x.e.%j.%A-%a.txt"
#SBATCH -o "%x.o.%j.%A-%a.txt"
#SBATCH --qos=rapid
#SBATCH --time=0:10:00

hostname
set -u

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"

TIME=$(date "+%s")

###################################################################################################
#setup-phase

#load tools
source ${SCRIPTdir}tools

#calculate THREADS
THREADS=$(($SLURM_CPUS_PER_TASK * 2))
echo $THREADS
###################################################################################################
#setup-phase

locTMP=${TMPdir}insert_TEs//$SLURM_ARRAY_TASK_ID/
mkdir -p ${locTMP}
echo $locTMP

#extract sniffles record
VCFline=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${TMPdir}insertions_to_process.txt)


###################################################################################################
#prepare stuff
rm -rf ${locTMP}*.log

#initiate log for current variant
nLOG=1
VARIANTid=$(echo $VCFline | tr ' ' '\t' | cut -f 3 )
echo ${nLOG}:-:JOBid=$SLURM_ARRAY_TASK_ID > ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))
echo ${nLOG}:-:VARIANTid=$VARIANTid >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))
CHR=$(echo $VCFline |  tr ' ' '\t' | cut -f 1 )
echo ${nLOG}:-:CHR=$CHR >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))
START=$(echo $VCFline |  tr ' ' '\t' | cut -f 2 )
echo ${nLOG}:-:START=$START >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))
#get read IDs from the SVC line
ENDpos=$( echo $VCFline | awk -v OFS="\t" -v TMPdir=$locTMP '{
  n=split($8,splitTAG,/;|=/)
  for(i=1; i<=n; i++){ 
    if(splitTAG[i]~"^RNAMES$"){
      gsub(",","\n",splitTAG[i+1])
      print splitTAG[i+1] > TMPdir "readIDs.txt"
    }
    if(splitTAG[i]~"^END$"){
      #hacky way of getting the last end position - better woruld be to fix the position in the TE identification scripty
      STOP=splitTAG[i+1]
    }
  }
  print STOP
  n=split($3, splitID,/:/)
  split($5, splitSEQ,/:/)
  for(i=1;i<=n; i++) {
    print ">"splitID[i]"\n"splitSEQ[i] > TMPdir "insertion.fa"
  }
}'  )

echo ${nLOG}:-:END=$ENDpos >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))

nINSERTION=$(seqkit fx2tab ${locTMP}insertion.fa | wc -l)
echo ${nLOG}:-:nINSERTION=$nINSERTION >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))
INSlength=$(echo $VCFline | awk '{ print length($5)}')
echo ${nLOG}:-:INSERTION_length=$INSlength >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))
VARseq=$(echo $VCFline | tr ' ' '\t' | cut -f 5 )
echo ${nLOG}:-:INSERTION_sequence=$VARseq >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))


###################################################################################################
#preparations

#create flank coordinates
FLANKshift=50

echo $CHR $START $ENDpos | 
  awk -v OFS="\t" -v FLANKshift=$FLANKshift -v FLANKlength=10000 -v TMP=$locTMP '
  {
    print $1,$2-FLANKlength-FLANKshift,$2-FLANKshift,"usFLANK",0,"+"
    print $1,$3+FLANKshift,$3+FLANKshift+FLANKlength,"dsFLANK",0,"+"
    print $1,$2-FLANKshift,$3+FLANKshift,"REPLACE",0,"+" > TMP "to_replace.bed"
  }' > ${locTMP}flanks.bed


#extract flank sequence
bedtools getfasta  -fo ${locTMP}flanks.fa -fi $assemblyFASTA -bed ${locTMP}flanks.bed -name

#determine length of flanks and if they might be too short
flankTooShort=$(seqkit fx2tab --name --length ${locTMP}flanks.fa | 
  awk -v OFS="\t" '{
    if($NF<2000){
      if(X=="")X=$1"_"$NF; else X=X":"$1"_"$NF
    }
  }
  END{
    if(X==""){
      print "OK"
    }else{
      print "flankShort:"X
    }
  }'
  )

echo ${nLOG}:-:FLANK_length=$flankTooShort >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))

if [[ $flankTooShort != OK ]]; then
    echo 0:-:STATUS=FLANKS_nOK >> ${locTMP}log.txt
    exit
fi

#-----------------------e---------------------------------------------------------------------------
#extract reads


#report variant supporting read number and IDs
nREADS=$(cat ${locTMP}readIDs.txt | wc -l)
readIDs=$(cat ${locTMP}readIDs.txt | tr '\n' ',')
echo ${nLOG}:-:READs_number=$nREADS >> ${locTMP}log.txt 
nLOG=$(( nLOG + 1 ))
echo ${nLOG}:-:READs_IDs=$readIDs >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))

#extract reads from input fasta 
if [[ $ALL == Y || ! -s  ${locTMP}reads.fa ]]; then
  seqkit grep --line-width 0 -r -f ${locTMP}readIDs.txt $INFILE > ${locTMP}reads.fa
fi

#determine number of variant supporting reads found in fasta
nreadsFOUND=$(seqkit fx2tab ${locTMP}reads.fa | wc -l)
echo ${nLOG}:-:READs_found=$nreadsFOUND >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))


###################################################################################################
#assemble reads to contig

if [[ $ALL == Y || ! -s ${locTMP}out-corr.fasta ]]; then
  {
    #-k 10 -p 0 --no-read-clip --no-chainning-clip -e 1 -R -A -l 500 -AS 1 -L 1000 --aln-dovetail  -
    wtdbg2 -x ont -g 30k -i ${locTMP}/reads.fa -L 5000 -R -S 1 -A -l 1024 -t $THREADS -fo ${locTMP}wtdbg2-out
    wtpoa -t $THREADS -i ${locTMP}wtdbg2-out.ctg.lay.gz -fo ${locTMP}dbg.raw.fa

    minimap2 -2 -Q -x map-ont --secondary=no -t $THREADS ${locTMP}dbg.raw.fa ${locTMP}/reads.fa >${locTMP}reads_mapped_wtcbg-contig.paf 2>/dev/null
    racon -m 8 -x -6 -g -8 -w 500 -t $THREADS ${locTMP}/reads.fa ${locTMP}reads_mapped_wtcbg-contig.paf ${locTMP}dbg.raw.fa >${locTMP}out-corr.fasta
  } 2>&1 | tee ${locTMP}assembly.log
fi
#determine the number of assembled contigs
nCONTIG=$(seqkit fx2tab ${locTMP}out-corr.fasta | wc -l )
echo ${nLOG}:-:CONTIG_number=$nCONTIG >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))
if [[ $nCONTIG -lt 1 ]]; then
  echo 0:-:STATUS=noContigAssembled >> ${locTMP}log.txt
  exit
fi

###################################################################################################
#determine position of insert within the contigs and generate flanking sequences
##test if insert only contained in one contig

minimap2 -x map-ont --secondary=no -t $THREADS ${locTMP}out-corr.fasta ${locTMP}insertion.fa > ${locTMP}insertion.mapped.paf

if [[ $nCONTIG -gt 1 ]]; then
  nCONTIGwithINSERT=$(cut -f 6 ${locTMP}insertion.mapped.paf | sort | uniq | wc -l)
  if [[ $nCONTIGwithINSERT -eq 1 ]]; then 
    corrCTG=$(cut -f 6 ${locTMP}insertion.mapped.paf | sort | uniq )
    mv ${locTMP}out-corr.fasta ${locTMP}out-corr.multi.fasta
    seqkit grep -p $corrCTG -r ${locTMP}out-corr.multi.fasta > ${locTMP}out-corr.fasta
    minimap2 -x map-ont --secondary=no -t $THREADS ${locTMP}out-corr.fasta ${locTMP}insertion.fa > ${locTMP}insertion.mapped.paf
  elif [[ $nCONTIGwithINSERT -eq 0 ]]; then
    echo 0:-:STATUS=InsertionNotOnContig >> ${locTMP}log.txt
    exit
  else
    echo 0:-:STATUS=tooManyContigs >> ${locTMP}log.txt
    exit
  fi
fi

for INS in $(echo $VARIANTid | tr ':' '\t' ); do
  awk -v OFS="\t" -v INS=$INS 'BEGIN{X=0}{if($1==INS) X+=1}END{print X}' ${locTMP}insertion.mapped.paf >> ${locTMP}nINSERTIOonCONTIG.log
  TMPvar=$(awk -v OFS="\t" -v INS=$INS '{if($1==INS) X[NR]=$5}END{for(i in X) {if(Y=="") {Y=X[i];} else Y=Y"~"X[i] ;} print Y}' ${locTMP}insertion.mapped.paf )
  echo $TMPvar >> ${locTMP}INSERTIOonCONTIG.STRAND.log
  TMPvar=$(awk -v OFS="\t" -v INS=$INS '{if($1==INS) {X[NR]=($4-$3)*100/$2}}END{for(i in X) {if(Y=="") {Y=X[i];} else Y=Y"~"X[i] ;} printf "%3.0f\n",  Y}' ${locTMP}insertion.mapped.paf ) 
  echo $TMPvar
  echo $TMPvar >> ${locTMP}INSERTIOonCONTIG.FRACTIONinsertion.log
  if [[ $TMPvar -lt 80 || $TMPvar -gt 120 ]]; then ERRORvar=Y; fi
  TMPvar=$(awk -v OFS="\t" -v INS=$INS '{if($1==INS) {X[NR]=($9-$8)*100/$2}}END{for(i in X) {if(Y=="") {Y=X[i];} else Y=Y"~"X[i] ;} printf "%3.0f\n",  Y}' ${locTMP}insertion.mapped.paf ) 
  echo $TMPvar >> ${locTMP}INSERTIOonCONTIG.FRACTIONcontig.log
  if [[ $TMPvar -lt 80 || $TMPvar -gt 120 ]]; then ERRORvar=Y; fi
done


nINSERTIOonCONTIG=$(cat ${locTMP}nINSERTIOonCONTIG.log | tr '\n' ':' | sed 's/:$/\n/')
echo ${nLOG}:-:number_INSERTIONS_on_contig=$nINSERTIOonCONTIG >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))

if [[ $nINSERTIOonCONTIG == *"0"* ]]; then
  echo 0:-:STATUS=InsertionNotOnContig >> ${locTMP}log.txt
  exit
fi

INSERTIOonCONTIG_STRAND=$(cat ${locTMP}INSERTIOonCONTIG.STRAND.log | tr '\n' ':' | sed 's/:$/\n/')
echo ${nLOG}:-:strands_INSERTIONS_on_contig=$INSERTIOonCONTIG_STRAND >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))
nSTRANDtypes=$(cat ${locTMP}INSERTIOonCONTIG.STRAND.log | tr '~' '\n'  | sort | uniq | wc -l )
STRANDtype=$(cat ${locTMP}INSERTIOonCONTIG.STRAND.log | tr '~' '\n'  | sort | uniq )

INSERTIOonCONTIG_FRACTIONinsertion=$(cat ${locTMP}INSERTIOonCONTIG.FRACTIONinsertion.log | tr '\n' ':' | sed 's/:$/\n/')
echo ${nLOG}:-:fractionInsertion_INSERTIONS_on_contig=$INSERTIOonCONTIG_FRACTIONinsertion >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))

INSERTIOonCONTIG_FRACTIONcontig=$( cat ${locTMP}INSERTIOonCONTIG.FRACTIONcontig.log |  tr '\n' ':' | sed 's/:$/\n/')
echo ${nLOG}:-:fractionContig_INSERTIONS_on_contig=$INSERTIOonCONTIG_FRACTIONcontig >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))

if [[ $INSERTIOonCONTIG_STRAND == *"~" || ${ERRORvar+x} == Y || $nSTRANDtypes -gt 1 ]]; then
    echo 0:-:STATUS=InsertionOnContigProblem >> ${locTMP}log.txt
  exit
fi


#reverse complement contig if both insertions are on - strand
if [[ $STRANDtype == "-" ]]; then
  mv ${locTMP}out-corr.fasta ${locTMP}out-corr.orig.fasta
  seqkit seq --reverse --complement --line-width 0 ${locTMP}out-corr.orig.fasta > ${locTMP}out-corr.fasta
  mv ${locTMP}insertion.mapped.paf ${locTMP}insertion.mapped.orig.paf
  minimap2 -x map-ont --secondary=no -t $THREADS ${locTMP}out-corr.fasta ${locTMP}insertion.fa > ${locTMP}insertion.mapped.paf
fi


#merge paf intervals to create insertion-coordinate
INSERTIONcoordONcontig=$(sort -k8,8n ${locTMP}insertion.mapped.paf | 
  awk -v OFS="\t" -v TMP=${locTMP} -v FLANKshift=$FLANKshift '{
    if(NR==1){
      CTG=$6
      START=$8
      STOP=$9
      STRAND=$5
      ctgLENGTH=$7
    } else{
      STOP=$9
      if($6!=CTG || $8<START || $5!=STRAND || $9<STOP) print "HELP" > TMP "error.log"
    }
  }
  END{
    print CTG":"START-FLANKshift"-"STOP+FLANKshift"_"STRAND"_"ctgLENGTH

    if(START<30000){ LEFT=START;}else{LEFT=30000}
    if(ctgLENGTH-STOP<30000){ ctgLENGTH-STOP;}else{RIGHT=30000}


  }'
)

echo ${nLOG}:-:INSERTIONcoordONcontig=$INSERTIONcoordONcontig >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))



###################################################################################################
#determine flank locations on the contig

#option -c gives base accuracy alignment
minimap2 -x splice -c -G 10k --splice-flank=no --secondary=no -t $THREADS ${locTMP}out-corr.fasta  ${locTMP}flanks.fa  > ${locTMP}flanks.mapped.paf

#performed for each flank individually to make reporting easier
    
for FLANK in us ds; do
  echo $FLANK
  #extract current flank sequence and map it to the genome
  awk -v OFS="\t" -v FLANK=$FLANK '{
    ID=FLANK"FLANK"
    if($1==ID){
      print
    }
  }' ${locTMP}flanks.mapped.paf > ${locTMP}flanks.mapped.${FLANK}.paf

  #determine how often the flank maps in the genome
  nFLANKmappings=$(cat ${locTMP}flanks.mapped.${FLANK}.paf | wc -l)
  echo ${nLOG}:-:nFLANKmappings_${FLANK}=$nFLANKmappings >> ${locTMP}log.txt
  nLOG=$(( nLOG + 1 ))

  #determine if ends of flanks mapped propperly
  rm -rf ${locTMP}multiple_flank_mappings.${FLANK}.txt
  FLANKendDIST=$( cat ${locTMP}flanks.mapped.${FLANK}.paf |
    awk -v FLANK=$FLANK -v TMPdir=$locTMP  -v OFS="\t" '{ 
      if(FLANK == "us"){ X=$2-$4}
      if(FLANK == "ds"){ X=$3}
      print X
    }' | 
    awk -v FLANK=$FLANK -v TMPdir=$locTMP  '{
      if(NR==1) {
        X=$1
      }else{
         X=X"-"$1
      } 
    }
    END{
      print X
      n=split(X,splitX,/-/)
      if(n>1){
        print X > TMPdir "multiple_flank_mappings."FLANK".txt"
      }
    }' )
  echo ${nLOG}:-:FLANK_unmappedEnd_${FLANK}=$FLANKendDIST >> ${locTMP}log.txt
  nLOG=$(( nLOG + 1 ))

  #!unclear if this is required with genomic flanks mapped to the contig
  if [[ -s ${locTMP}multiple_flank_mappings.${FLANK}.txt  ]]; then
    mv ${locTMP}flanks.mapped.${FLANK}.paf ${locTMP}flanks.mapped.${FLANK}.multi.paf

    awk -v OFS="\t" -v FLANK=$FLANK '{
      if($5=="+"){
        if(FLANK == "us"){
          X[$0]=$2-$4
        }else{
          X[$0]=$3
        }
      }
    }
    END{
      OLD=100000000000000
      for(i in X){
        if(X[i] < OLD){
          OLD=X[i]
          CURRid=i
        }
      }
      print CURRid
    }' ${locTMP}flanks.mapped.${FLANK}.multi.paf  >${locTMP}flanks.mapped.${FLANK}.paf
  fi

  #determine how often the flank maps in the genome after reduction
  nFLANKmappings=$(cat ${locTMP}flanks.mapped.${FLANK}.paf | wc -l)
  echo ${nLOG}:-:nFLANKmappings_${FLANK}-afer-reduction=$nFLANKmappings >> ${locTMP}log.txt
  nLOG=$(( nLOG + 1 ))

  if [[ ! $nFLANKmappings -eq 1 ]]; then
    echo 0:-:STATUS=FLANKS_nOK >> ${locTMP}log.txt
    exit
  fi

  FLANKendDIST2=$( cat ${locTMP}flanks.mapped.${FLANK}.paf |
    awk -v FLANK=$FLANK -v TMPdir=$locTMP  '{ 
      if(FLANK == "us"){ X=$2-$4}
      if(FLANK == "ds"){ X=$3}
      print X
    }' | awk '{if(NR==1) X=$1; else X=X"-"$1} END{print X}' )
  echo ${nLOG}:-:FLANK_unmappedEnd-afterReduction_${FLANK}=$FLANKendDIST2 >> ${locTMP}log.txt
  nLOG=$(( nLOG + 1 ))

  if [[  $FLANKendDIST2 -le 30  && $FLANKendDIST2 -gt -20 ]]; then
    echo ${nLOG}:-:FLANK_dist-STATUS_${FLANK}=OK >> ${locTMP}log.txt
    nLOG=$(( nLOG + 1 ))
  else 
    echo ${nLOG}:-:FLANK_dist-STATUS_${FLANK}=${FLANK}FLANK_tooDistant >> ${locTMP}log.txt
    nLOG=$(( nLOG + 1 ))
    SWITCH=Y
  fi

  FLANKstrand=$(cat ${locTMP}flanks.mapped.${FLANK}.paf |
    awk '{print $5}')
  echo ${nLOG}:-:FLANK_MapStrand_${FLANK}=$FLANKstrand >> ${locTMP}log.txt
  nLOG=$(( nLOG + 1 ))

done

SINGLEflankSTATUS=""
if [[ -n ${SWITCH+x} ]]; then
  SINGLEflankSTATUS=$(
    awk -v OFS="\t" '{
      if($0~"INSERTION_length="){ split($1, splitTAG, /=/); INSERTIONlength=splitTAG[2]}
      if($0~"FLANK_dist-STATUS_us="){ split($1, splitTAG, /=/); FLANKdistSTATUS_us=splitTAG[2]}
      if($0~"FLANK_dist-STATUS_ds="){ split($1, splitTAG, /=/); FLANKdistSTATUS_ds=splitTAG[2]}
      if($0~"FLANK_unmappedEnd-afterReduction_us="){ split($1, splitTAG, /=/); FLANKunmappedEndafterReduction_us=splitTAG[2]}
      if($0~"FLANK_unmappedEnd-afterReduction_ds="){ split($1, splitTAG, /=/); FLANKunmappedEndafterReduction_ds=splitTAG[2]}
    }
    END{
      if(INSERTIONlength > 1000){
        if( FLANKdistSTATUS_us == "OK") {
          if (FLANKunmappedEndafterReduction_us < 25 && FLANKunmappedEndafterReduction_ds < 10000 ){
            print "OK-singleFLANK"
          }
        }else{
          if(FLANKdistSTATUS_ds == "OK" ) {
            if (FLANKunmappedEndafterReduction_ds < 25 && FLANKunmappedEndafterReduction_us < 10000 ){
              print "OK-singleFLANK"
            }

          }
        }
      }
    }' ${locTMP}log.txt
  )
  echo ${SINGLEflankSTATUS}
  if [[ ${SINGLEflankSTATUS} != "OK-singleFLANK" ]]; then
    echo 0:-:STATUS=FLANKS_nOK >> ${locTMP}log.txt
    exit
  fi

fi


#reset SWITCHT variable
SWITCH=""

#determine validity of flank position on the contig

#preset the corrected coordinate variable     
INSERTIONcoordONcontigCORR=$INSERTIONcoordONcontig
FLANKendDISTsum=0


for FLANK in us ds; do
  echo $FLANK

  FLANKdist=$(
    awk -v OFS="\t" -v INSERTIONcoordONcontig=${INSERTIONcoordONcontig} -v FLANK=$FLANK -v FLANKshift=$FLANKshift '{
      split(INSERTIONcoordONcontig, splitCOORD,/:|-|_/)
      if($6==splitCOORD[1]) {
        if(FLANK == "us"){
          print  splitCOORD[2] - $9
        }
        if(FLANK == "ds"){
          print $8 - splitCOORD[3]  
        }
      }else{
        print "CTGerror"
      }
    }' ${locTMP}flanks.mapped.${FLANK}.paf
  )

  echo ${nLOG}:-:FLANK_dist-INS-onCTG_${FLANK}=$FLANKdist >> ${locTMP}log.txt
  nLOG=$(( nLOG + 1 ))


  if [[  $FLANKdist -le 2500  && $FLANKdist -gt -100  && $FLANKdist != "CTGerror" ]]; then
    #extend insertion region on the contig
    INSERTIONcoordONcontigCORR=$(
      echo $INSERTIONcoordONcontigCORR | 
      awk -v OFS="\t" -v FLANK=$FLANK -v FLANKendDIST2=$FLANKdist '{
        if(FLANK=="us"){ 
          split($1,splitCOORD,/:|-/)
          print splitCOORD[1]":"splitCOORD[2]-FLANKendDIST2"-"splitCOORD[3]
        }else{
          split($1,splitCOORD,/-|_/)
          print splitCOORD[1]"-"FLANKendDIST2+splitCOORD[2]"_"splitCOORD[3]"_"splitCOORD[4]
        }
      }' 
      )

    FLANKendDISTsum=$(( $FLANKendDISTsum + $FLANKdist ))

    if [[ $FLANKdist -lt 100 ]]; then
      echo ${nLOG}:-:FLANK_dist-STATUS_${FLANK}=OK >> ${locTMP}log.txt
      nLOG=$(( nLOG + 1 ))
    else
      echo ${nLOG}:-:FLANK_dist-STATUS_${FLANK}=EXTENDED >> ${locTMP}log.txt
      nLOG=$(( nLOG + 1 ))
    fi
  elif [[ $FLANKdist == CTGerror ]]; then 
    echo ${nLOG}:-:FLANK_dist-STATUS_${FLANK}=${FLANK}wrongContig >> ${locTMP}log.txt
    nLOG=$(( nLOG + 1 ))
    SWITCH=Y
  else 
    echo ${nLOG}:-:FLANK_dist-STATUS_${FLANK}=${FLANK}FLANK_tooDistant >> ${locTMP}log.txt
    nLOG=$(( nLOG + 1 ))
    SWITCH=Y
  fi

  if [[ $FLANKdist == CHRerror ]]; then
    SWITCH=Y
  fi
done

echo ${nLOG}:-:INSERTIONcoordONcontigCorrected=$INSERTIONcoordONcontigCORR >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))
echo ${nLOG}:-:FLANK_dist-SUM=$FLANKendDISTsum >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))


if [[ ${SWITCH} == Y ]]; then
  echo 0:-:STATUS=FLANKdist-INS_nOK >> ${locTMP}log.txt
  exit
fi

#finish the sequence to insert into the sequence
#create bed file with corrected insert coordinates on contig
echo $INSERTIONcoordONcontigCORR | 
  awk -v OFS="\t" '{
    split($1, splitCOORD,/:|-|_/)
    print splitCOORD[1],splitCOORD[2],splitCOORD[3],"SequenceToInsert",0,splitCOORD[1]splitCOORD[4] 
  }' >${locTMP}SequenceToInsert.bed

bedtools getfasta -tab -fi ${locTMP}out-corr.fasta -s -bed ${locTMP}SequenceToInsert.bed >${locTMP}SequenceToInsert.tab
lengthSequenceToInsert=$( awk '{print length($2)}' ${locTMP}SequenceToInsert.tab )

REPLACEsize=$( 
  echo $VARseq |
  awk -v OFS="\t" -v START=$START -v ENDpos=$ENDpos -v CHR=$CHR -v FLANK=$FLANK -v FLANKshift=$FLANKshift -v FLANKendDISTsum=$FLANKendDISTsum '
  {
    #calculate the length the replacement sequence should
    gsub(/:/,//,$1)
    VARlength=length($1)
    REPLACElength=VARlength+(ENDpos-START)+FLANKendDISTsum+2*FLANKshift
    print REPLACElength
  }' 
)

echo ${nLOG}:-:TheoreticalSizeToInsert=$REPLACEsize >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))
echo ${nLOG}:-:ActualSizeToInsert=$lengthSequenceToInsert >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))


REPLACEfract=$(( $lengthSequenceToInsert * 100 / $REPLACEsize ))
echo ${nLOG}:-:lengthToInsert-vs-lengthActual=$REPLACEfract >> ${locTMP}log.txt
nLOG=$(( nLOG + 1 ))

if [[ ($REPLACEsize -gt 500 && ($REPLACEfract -gt 130 || $REPLACEfract -lt 90 )) || ($REPLACEsize -le 500 && ($REPLACEfract -gt 2000 || $REPLACEfract -lt 60 )) ]]; then
  echo 0:-:STATUS=ReplaceSizeNotMatching >> ${locTMP}log.txt
  exit
fi

#if all is ok give it good status
if [[ ${SINGLEflankSTATUS} == "OK-singleFLANK" ]]; then
  STATUS=OK-singleFLANK 
else
  STATUS=OK 
fi

echo 0:-:STATUS=$STATUS >> ${locTMP}log.txt



###################################################################################################
#create VCFline
#create sequence to replace
bedtools getfasta -tab -fi $assemblyFASTA -s -bed ${locTMP}to_replace.bed >${locTMP}SequenceToReplace.tab

awk -v TMP=$locTMP -v VARIANTid=$VARIANTid -v STATUS=$STATUS '
BEGIN{
  INFILE=TMP "SequenceToInsert.tab"
  while((getline LINE < INFILE) > 0) {
    split(LINE,splitLINE,/ |\t/)
    SeqInsert=splitLINE[2]
  }
  INFILE=TMP "SequenceToReplace.tab"
  while((getline LINE < INFILE ) > 0) {
    split(LINE,splitLINE,/ |\t/)
    SeqReplace=splitLINE[2]
  }
}
{
  print $1,$2+1,VARIANTid,SeqReplace,SeqInsert,0,"PASS","END="$3,"STATUS",STATUS
}' ${locTMP}to_replace.bed > ${locTMP}out.vcf


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

###################################################################################################
PROCESSED_TIME=$(echo -e $(date "+%s") "$TIME" | mawk '{ print ($1-$2)/60 }')
echo "integration of TEs into the genome - processing_time=" "${PROCESSED_TIME}" >>"${OPENdir}time-log.txt"

exit
