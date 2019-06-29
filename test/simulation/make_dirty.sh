#!/bin/bash

# requires:
# (1) base simulation folder (ex. contamination_simulation/)
# (2) contamination folder (ex. brucei_vs_giardia_10_pct/)
# (3) base dataset id
# (4) contaminator dataset id
# (5) percentage

set -e

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ]
  then
    echo "Wrong arguments supplied"
    exit 1
fi

BASE_FOLDER=$(readlink -f $1)
CONT_FOLDER=$(readlink -f $2)
OUR_ID=$3
FOREIGN_ID=$4
PERCENTAGE=$5

DIRTY_FOLDER=$CONT_FOLDER/$OUR_ID
mkdir -p $DIRTY_FOLDER

DATASETS_FOLDER=$BASE_FOLDER/datasets
mkdir -p $DATASETS_FOLDER
OUR_LEFT_READS=${DATASETS_FOLDER}/${OUR_ID}_1.fastq
OUR_RIGHT_READS=${DATASETS_FOLDER}/${OUR_ID}_2.fastq

FOREIGN_LEFT_READS=${DATASETS_FOLDER}/${FOREIGN_ID}_1.fastq
FOREIGN_RIGHT_READS=${DATASETS_FOLDER}/${FOREIGN_ID}_2.fastq

## Random reads generation
SEED=100

echo 'Generating random reads...'
count=$(awk '{s++}END{print s/4}' ${OUR_LEFT_READS})
count=$(python -c "print(int(round( ($count/(1.0-$PERCENTAGE/100.0)) * ($PERCENTAGE/100.0) )))")
echo "${count} reads to slice..."

echo 'Making left samples'
seqtk sample -s$SEED ${FOREIGN_LEFT_READS} $count > ${DIRTY_FOLDER}/left_foreign.fastq

echo 'Making right samples'
seqtk sample -s$SEED ${FOREIGN_RIGHT_READS} $count > ${DIRTY_FOLDER}/right_foreign.fastq

## Mixing reads
echo 'Mixing...'
cat ${OUR_LEFT_READS} ${DIRTY_FOLDER}/left_foreign.fastq > ${DIRTY_FOLDER}/left_dirty.fastq
cat ${OUR_RIGHT_READS} ${DIRTY_FOLDER}/right_foreign.fastq > ${DIRTY_FOLDER}/right_dirty.fastq

## Assembling dirty transcriptome
echo 'Assembling'
TRINITY_FOLDER=${DIRTY_FOLDER}/${OUR_ID}_dirty_trinity

Trinity --seqType fq --CPU 30 --max_memory 60G \
        --left ${DIRTY_FOLDER}/left_dirty.fastq \
        --right ${DIRTY_FOLDER}/right_dirty.fastq \
        --output ${TRINITY_FOLDER}

sed "s#TRINITY#${OUR_ID}#" ${TRINITY_FOLDER}/Trinity.fasta > ${DIRTY_FOLDER}/dirty_transcriptome.fasta
rm -r $TRINITY_FOLDER

## Map foreign reads to dirty transcriptome
echo 'Dirty reads mapping'
INDEX_FOLDER=${DIRTY_FOLDER}/dirty_bowtie_index
mkdir $INDEX_FOLDER

echo 'Building index...'
bowtie2-build --threads 30 ${DIRTY_FOLDER}/dirty_transcriptome.fasta ${INDEX_FOLDER}/index

echo 'Running mapping...'
bowtie2 --very-sensitive -p 32 \
                         -x ${INDEX_FOLDER}/index \
                         -1 ${DIRTY_FOLDER}/left_foreign.fastq \
                         -2 ${DIRTY_FOLDER}/right_foreign.fastq \
                         -S ${DIRTY_FOLDER}/foreign_alignment.sam
rm -r $INDEX_FOLDER

## Generate bad contigs list
echo 'Calculating depth'
pileup.sh in=${DIRTY_FOLDER}/foreign_alignment.sam \
          ref=${DIRTY_FOLDER}/dirty_transcriptome.fasta \
          out=${DIRTY_FOLDER}/foreign_alignment.pileup \
          overwrite=true

echo 'Generating bad contigs list'
tail --lines=+2 ${DIRTY_FOLDER}/foreign_alignment.pileup | \
  awk -F $'\t'  'BEGIN {OFS=","} { if ($5 > 50) print $1 }' > \
  ${DIRTY_FOLDER}/contaminated_list.txt

echo 'Copying result to output folder'
DECROSS_FOLDER=$CONT_FOLDER/for_decross
mkdir -p $DECROSS_FOLDER

cp ${DIRTY_FOLDER}/left_dirty.fastq ${DECROSS_FOLDER}/${OUR_ID}_1.fastq
cp ${DIRTY_FOLDER}/right_dirty.fastq ${DECROSS_FOLDER}/${OUR_ID}_2.fastq
cp ${DIRTY_FOLDER}/dirty_transcriptome.fasta ${DECROSS_FOLDER}/${OUR_ID}.fasta
