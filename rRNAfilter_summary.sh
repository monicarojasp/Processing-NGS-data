# !/bin/bash
# generate summary of rRNA filtering
              
array=()

cd ../logs/

echo -e 'file\t# read pairs\t# concordant once\t% concordant once\t# concordant more than once\t% concordant more than once\t# of aligned discordant\t% of aligned disdorant\t# of mates condordant once\t% of mates concordant once\t# of mates concordant more than once\t% of mates concordant more than once\toverall % alignment' > rRNA_summary.txt

echo "Array size: ${#array[*]}"
echo "Array items:"

for item in ${array[*]}
do
        printf "   %s\n" $item

	item="$item"_rRNAfilter
        OUTPUT1=$(cat rRNAfilter/$item.log | grep -E 'reads' | sed 's/ reads.*//')
        OUTPUT2=$(cat rRNAfilter/$item.log | grep -E 'concordantly exactly 1 time' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT3=$(cat rRNAfilter/$item.log | grep -E 'concordantly exactly 1 time' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT4=$(cat rRNAfilter/$item.log | grep -E 'aligned concordantly >1 times' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT5=$(cat rRNAfilter/$item.log | grep -E 'aligned concordantly >1 times' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT6=$(cat rRNAfilter/$item.log | grep -E 'aligned discordantly 1 time' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT7=$(cat rRNAfilter/$item.log | grep -E 'aligned discordantly 1 time' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT8=$(cat rRNAfilter/$item.log | grep -E 'aligned exactly 1 time' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT9=$(cat rRNAfilter/$item.log | grep -E 'aligned exactly 1 time' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT10=$(cat rRNAfilter/$item.log | grep -E 'aligned >1 times' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT11=$(cat rRNAfilter/$item.log | grep -E 'aligned >1 times' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT12=$(cat rRNAfilter/$item.log | grep -E 'overall alignment rate' | sed 's/%.*//')

        echo -e ''$item'\t'${OUTPUT1}'\t'${OUTPUT2}'\t'${OUTPUT3}'\t'${OUTPUT4}'\t'${OUTPUT5}'\t'${OUTPUT6}'\t'${OUTPUT7}'\t'${OUTPUT8}'\t'${OUTPUT9}'\t'${OUTPUT10}'\t'${OUTPUT11}'\t'${OUTPUT12}'' >> rRNA_summary.txt
         
done
