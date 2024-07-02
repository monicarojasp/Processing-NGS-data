#!/bin/bash                                                                                                                                          
# align reads to the genome

cd ../logs/

array=( )

echo -e 'File\t# reads\t% uniquely mapped\t# of splices\t% mismatch rate per base\t% deletion rate per base\t% insertion rate per base\t% reads mapped to multiple loci\t% reads unmapped' > mapping_to_human_summary.txt

echo "Array size: ${#array[*]}"
echo "Array items:"

for item in ${array[*]}
do

    printf "   %s\n" $item

    # write log                                                                                                                                       
    cd ../mapping/
    item="$item"Log
    OUTPUT1=$(cat $item.final.out | grep "Number of input reads" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT2=$(cat $item.final.out | grep "Uniquely mapped reads %" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT3=$(cat $item.final.out | grep "Number of splices: Total" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT4=$(cat $item.final.out | grep "Mismatch rate per base, %" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT5=$(cat $item.final.out | grep "Deletion rate per base" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT6=$(cat $item.final.out | grep "Insertion rate per base" | sed 's/.*|//' | awk '{$1=$1}{ print }')

    OUTPUT7=$(cat $item.final.out | grep "% of reads mapped to multiple loci" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT8=$(cat $item.final.out | grep "% of reads mapped to too many loci" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')

    OUTPUT9=$(cat $item.final.out | grep "% of reads unmapped: too many mismatches" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT10=$(cat $item.final.out | grep "% of reads unmapped: too short" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT11=$(cat $item.final.out | grep "% of reads unmapped: other" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')

    OUTPUT12=$(echo $OUTPUT7 + $OUTPUT8 | bc)
    OUTPUT13=$(echo $OUTPUT9 + $OUTPUT10 + $OUTPUT11 | bc)

    echo -e ''$item'\t'${OUTPUT1}'\t'${OUTPUT2}'\t'${OUTPUT3}'\t'${OUTPUT4}'\t'${OUTPUT5}'\t'${OUTPUT6}'\t'${OUTPUT12}'%\t'${OUTPUT13}'%\t' >> ../logs/mapping_to_human_summary.txt

    cd ../nohmrRNA/

done
