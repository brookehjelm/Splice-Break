#!/bin/bash
#Splice-Break v1.0.0
#January 25, 2019
#Written by Michelle Webb, michelgw@usc.edu
#University of Southern California

inputDir=$1
outputDir=$2
refDir=$3
filterFails=0
time=`date +%d-%m-%Y-%H-%M-%S`
hostname=`hostname`
echo "Starting $0 at ${time} on ${hostname}"

if [[ -d ${inputDir} && -d ${outputDir} && -d ${refDir} ]] ; then
        echo "MapSplice Dir: ${inputDir}"
        echo "Output Dir: ${outputDir}"
        echo "Reference Dir: ${refDir}"
else
        echo "### ERROR: Check input parameters"
        exit
fi

#Run Samtools mPileup
echo "Running Samtools"
samtools mpileup -d 1000000 ${inputDir}/alignments.bam > ${outputDir}/pileup.txt
if [ $? -gt 0 ]; then
    echo "#ERROR: Please check that samtools is in your PATH"
    exit
fi

#Run CountBases
echo "Running CountBases.py"
python ${refDir}/CountBases.py ${outputDir}/pileup.txt > ${outputDir}/BaseCounts.pileup.txt
if [ $? -gt 0 ]; then
    echo "#ERROR: Please check your Python version or pileup output"
    exit
fi

#Add coverage across bases
echo "Adding per base coverage"
awk 'NR>1{print $2+$3+$4+$5}' ${outputDir}/BaseCounts.pileup.txt > ${outputDir}/init_coverage.txt
sed -i '1i COVERAGE' ${outputDir}/init_coverage.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

#Prints an index for the coverage values
awk -F'\t' -v OFS='\t' '
  NR == 1 {print "NUM",$0; next}
  {print (NR-1), $0}
' ${outputDir}/init_coverage.txt > ${outputDir}/combined.txt

#Joins the reference position and actual coverage
join <( sed '1d' ${outputDir}/combined.txt)  <(sed '1d' ${refDir}/reference_and_primer_positions.txt) | sed $"1 i\rCRS_NC_012920.1\tCoverage\tPrimer_Position" | tr ' ' '\t' | sed -e "s///" > ${outputDir}/Coverage.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

#Calculates benchmark coverage and determines appropriate value
left_benchmark=`awk 'NR>1{print $1"\t"$2}' ${outputDir}/Coverage.txt | grep -w 957 -A249 | cut -f2 | awk '{total += $1} END {print total/NR}'`
right_benchmark=`sed '1d' ${outputDir}/Coverage.txt | awk '{print $1"\t"$2}' | grep -w 15357 -A249 | cut -f2 | awk '{total += $1} END {print total/NR}'`
left_div_right=`echo "scale=3;${left_benchmark} / ${right_benchmark}" | bc -l`
avg_benchmark=`echo "scale=3;(${left_benchmark} + ${right_benchmark}) / 2" | bc -l`
junction_length=`cat ${inputDir}/junctions.txt | wc -l`
if (( $(echo "${left_div_right}  <= 1.4" | bc -l) )) ; then
    echo "Using average_benchmark"
    eval $(echo printf '"${avg_benchmark}%.s\n"' {1..${junction_length}}) | sed '1i Benchmark_Coverage' > ${outputDir}/avgBenchmark.txt
else
    echo "Using left_benchmark"
    eval $(echo printf '"${left_benchmark}%.s\n"' {1..${junction_length}}) | sed '1i Benchmark_Coverage' > ${outputDir}/avgBenchmark.txt
fi

#Assembling Junctions File
cat ${inputDir}/junctions.txt | cut -d$'\t' -f1-3,5,11 | awk -F ',' '{print $1 "\t" $2}' | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2"-"$3"\t"$3 - $2 - 1}' > ${outputDir}/modifiedJunc.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

#Annotating
echo "Annotating"
awk 'NR==FNR{a[$1]=$2;next}{if($1 in a){print $0,a[$1]; delete a[$1]} else print $0,"novel_deletion"}' ${refDir}/deletion_frequency_types.txt <(cut -f7 ${outputDir}/modifiedJunc.txt) > ${outputDir}/breakpoints.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

#Assembling
paste ${outputDir}/modifiedJunc.txt <(awk '{print $2}' ${outputDir}/breakpoints.txt) <( sed '1d' ${outputDir}/avgBenchmark.txt) | column -t -s $'\t' > ${outputDir}/intermediate.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

paste ${outputDir}/intermediate.txt <(awk '{printf "%.4f\n", ($4/$10)*100}' ${outputDir}/intermediate.txt) > ${outputDir}/intermediate_junctions.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

awk '{print $1"\t"$7"\t"$2"\t"$3"\t"$8"\t"$4"\t"$10"\t"$11"\t"$9"\t"$5"\t"$6}' ${outputDir}/intermediate_junctions.txt > ${outputDir}/intermediate_junctions2.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

#Filter1: Remove junction calls with read overhang < 20bp
echo "Applying Filter 1"
awk '(NR>1) && ( ($10>=20) ) && ( ($11>=20) )' ${outputDir}/intermediate_junctions2.txt > ${outputDir}/large_deletions.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

#Filter2: 5' Breakpoint greater than 356, 3' Breakpoint less than 15926
echo "Applying Filter 2"
awk '(NR>1) && ( ($3>=356) && ($4 <= 15926) )' ${outputDir}/large_deletions.txt | sed $"1iReference_Genome\tMapSplice_Breakpoint\t5'_Break\t3'_Break\tDeletion_Size_bp\tDeletion_Reads\tBenchmark_Coverage\tDeletion_Read_%\tAnnotation\tLeft_Overhang\tRight_Overhang" | column -t -s $'\t' > ${outputDir}/Large_Deletions_NC_012920.1_356-15926.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

#Filter2: 5' Breakpoint greater than 106, 3' Breakpoint less than 16176"
awk '(NR>1) && ( ($3>=106) && ($4 <= 16176) )' ${outputDir}/large_deletions.txt | sed $"1iReference_Genome\tMapSplice_Breakpoint\t5'_Break\t3'_Break\tDeletion_Size_bp\tDeletion_Reads\tBenchmark_Coverage\tDeletion_Read_%\tAnnotation\tLeft_Overhang\tRight_Overhang" | column -t -s $'\t' > ${outputDir}/Large_Deletions_NC_012920.1_106-16176.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

sed $"1iReference_Genome\tMapSplice_Breakpoint\t5'_Break\t3'_Break\tDeletion_Size_bp\tDeletion_Reads\tBenchmark_Coverage\tDeletion_Read_%\tAnnotation\tLeft_Overhang\tRight_Overhang" ${outputDir}/large_deletions.txt | column -t -s $'\t' > ${outputDir}/Large_Deletions_NC_012920.1_No-Position-Filter.txt
if [ $? -gt 0 ]; then
    ((filterFails++))
fi

#Delete Intermediary Files
rm ${outputDir}/{init_coverage.txt,avgBenchmark.txt,large_deletions.txt,breakpoints.txt,modifiedJunc.txt,intermediate.txt,intermediate_junctions.txt,intermediate_junctions2.txt,combined.txt,pileup.txt,BaseCounts.pileup.txt}

#Final error check
if [[ ${filterFails} -eq 0 ]]; then
    echo "Junction Filter Steps Success"
else
    echo "#ERROR: ${filterFails} steps failed. Please check your input files."
fi

time=`date +%d-%m-%Y-%H-%M-%S`
echo "Ending $0 at ${time}"
