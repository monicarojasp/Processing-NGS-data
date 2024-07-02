# !/bin/bash
# This commands use bowtie2 to filter rRNA by mapping to index of human rRNA

cd ../raw

array=()


echo "Array size: ${#array[*]}"
echo "Array items:"

for item in ${array[*]}
do
        printf "   %s\n" $item

        #bowtie2                                                                                                     
        cmd='/usr/local/bin/bowtie2 -p 50 -x /vol01/genome/rRNA/bowtie2_index/hmrRNA --un-conc-gz ../nohmrRNA/'$item'_nohmrRNA.fastq.gz -1 ../raw/'$item'_R1_001.fastq.gz -2 ../raw/'$item'_R2_001.fastq.gz -S '$item'.sam 1>>../logs/rRNAfilter/'$item'_rRNAfilter.log 2>&1'
        
	echo $cmd
        echo $cmd > ../logs/rRNAfilter/$item\_rRNAfilter.log
        eval $cmd

done
