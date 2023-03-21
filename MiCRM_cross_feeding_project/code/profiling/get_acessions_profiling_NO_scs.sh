#!/bin/bash

#conda activate project

# Get sample name
FILENAME=$(basename -- "$0")
SAMPLE="${FILENAME%.*}"

# Defining PATHS
SAMPLE_CSV="../../results/output/"
mkdir ../../results/output/$SAMPLE
MODEL_OUTPUT="../../results/output/$SAMPLE/"

SIZE=$(cat $SAMPLE_CSV${SAMPLE}.csv | wc -l)

echo "Processing sample $SAMPLE with a community of $SIZE species...."


# *** maybe can delete this line below *** #
#cut -d ',' -f 1,3 $SAMPLE_CSV/${SAMPLE}.csv > acessions_${SAMPLE}.txt

OLDIFS="$IFS"
IFS=$'\n'


echo "Carving started..."

for i in $(cat $SAMPLE_CSV${SAMPLE}.csv)
    do
        ACESSION=$(echo $i | cut -d ',' -f 1)

        GENUS=$(echo $i | cut -d ',' -f 2 | awk -F " " '{print $1}')
        SPECIES=$(echo $i | cut -d ',' -f 2 | awk -F " " '{print $2}')

        organism_name=$(echo "${GENUS}_${SPECIES}")

        carve ../genomes_folder/$SAMPLE/$ACESSION.faa.gz -o $MODEL_OUTPUT/${organism_name}.xml

	wait
	
    done

IFS="$OLDIFS"

FINAL_SIZE=$( ls $MODEl_OUTPUT | wc -l )

echo "Sucessfully generate metabolic models for ${FINAL_SIZE} out of ${SIZE} species \
in the sample ${SAMPLE}!!!"

rm -r ../genomes_folder/$SAMPLE
rm $SAMPLE_CSV${SAMPLE}.csv


# echo "Running SMETANA..."

# smetana -v ../../results/output/$SAMPLE/*.xml -o ../../results/output/${SAMPLE}

# smetana -v -d ../../results/output/$SAMPLE/*.xml -o ../../results/output/${SAMPLE}

# wait

# rm -r $MODEl_OUTPUT

# echo "CarveMe and SMETANA Completed for ${SAMPLE}!!!"

# #================================================================#
# echo "Running SMETANA minimising the minimal media composition..."

# smetana --molweight -v ../../results/output/$SAMPLE/*.xml -o ../../results/output/${SAMPLE} --debug

# smetana --molweight -v -d ../../results/output/$SAMPLE/*.xml -o ../../results/output/${SAMPLE}

# wait

# rm -r $MODEl_OUTPUT

# echo "CarveMe and SMETANA Completed for ${SAMPLE}!!!"

#================================================================#
echo "Running SMETANA adding flavor settings to BIGG"

smetana --molweight --flavor bigg -v ../../results/output/$SAMPLE/*.xml -o ../../results/output/${SAMPLE} --debug

smetana --molweight --no-coupling --flavor bigg -v -d ../../results/output/$SAMPLE/*.xml -o ../../results/output/${SAMPLE}_no_scs

wait

rm -r ../../results/output/$SAMPLE/

echo "CarveMe and SMETANA Completed for ${SAMPLE}!!!"


