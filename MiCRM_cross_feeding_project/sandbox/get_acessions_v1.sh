#!/bin/bash

conda activate project


SAMPLE=$(cut -d ',' -f 2 sampled_acessions.csv | uniq)
SIZE=$(cat sampled_acessions.csv | wc -l)

OLDIFS="$IFS"
IFS=$'\n'

cut -d ',' -f 1,3 sampled_acessions.csv > acessions_${SAMPLE}.txt

mkdir output/$SAMPLE

echo "Carving..."

for i in $(cat acessions_${SAMPLE}.txt)
    do
        ACESSION=$(echo $i | cut -d ',' -f 1)

        GENUS=$(echo $i | cut -d ',' -f 2 | awk -F " " '{print $1}')
        SPECIES=$(echo $i | cut -d ',' -f 2 | awk -F " " '{print $2}')

        organism_name=$(echo "${GENUS}_${SPECIES}")

        carve --refseq ${ACESSION} -o output/$SAMPLE/${organism_name}.xml
	wait
	echo "removing ${ACESSION} to trash!!!"
	mv ${ACESSION}.faa* trash/

    done

IFS="$OLDIFS"

FINAL_SIZE=$(wc -l output/$SAMPLE/)

echo "Sucessfully generate metabolic models for ${FINAL_SIZE} out of ${SIZE} species\
in the sample ${SAMPLE}!!!"

rm trash/*faa*
rm acessions_${SAMPLE}.txt


echo "Running SMETANA..."

smetana -v output/$SAMPLE/*.xml -o output/${SAMPLE}

smetana -v -d output/$SAMPLE/*.xml -o output/${SAMPLE}

wait

rm -r output/$SAMPLE/

echo "CarveMe and SMETANA Completed for ${SAMPLE}!!!"
