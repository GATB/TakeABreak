#!/bin/sh
#************

function resume_parameters {
	
	echo "***Parameters values***\n" 
	echo "Output prefix:" $prefix_output
	
	if [ "$graph_file" == "" ]; then
		echo "De Bruijn Graph construction "	
		echo "\t k="$k "  # option -k"
		echo "\t min_coverage="$nks 
	fi
	
	echo "discoInv"
	echo "\t reverse_tolerance="$reverse_tolerance   
	echo "\t LCT="$loc_cmpx
	echo "\t LCS="$lcs_restriction_percentage 
	
	
}

function help {
echo "TakeBreak.sh,"
echo "Usage: ./TakeBreak.sh options"
echo  "\t In/out Options (-r or -g mandatory):"
echo  "\t \t -r STRING. List of reads separated by comma ',' without space. Note that reads may be in fasta or fastq format, gzipped or not. Example: -r data/toy_example_reads.fasta,data/toy_example_with_inv_reads.fasta"
echo  "\t \t \t Incompatible with -g option"
echo  "\t \t -g STRING. Name of the already existing graph file (.h5)." 
echo  "\t \t \t Incompatible with -r option"
echo  "\t\t -p STRING. All out files will start with this prefix. Default: \"TakeABreak_Expe-date\""
echo  "\t De bruijn Graph Options:"
echo  "\t\t -k INT. Set the length of used kmers. Incompatible with -g option. Default=31."
echo  "\t\t -S INT. Set the Solidity of used kmers (minimal number of occurrences a k-mer should have not to be treated as sequencing error). Incompatible with -g option. Default=3."
echo  "\t Inversion detection options:"
echo  "\t\t -c INT. LCT (local complexity threshold): Defaults: 100."
echo  "\t\t -m INT: max_sim: max similarity percentage: Inversions with a and b' (or u and v') whose longuest common subsequence size is bigger than k*(this value)/100 are discarded. Defaults: 80"
echo  "\t\t -r INT: (optimization parameter lower=longer, higher=false negatives) max repeated size suffix of u and v': Defaults: 8"
echo "Any further question: read the readme file or contact us: claire.lemaitre@inria.fr"
}

k=-1
nks=-1

loc_cmpx_LCT=100
lcs_restriction_percentage=80
reverse_tolerance=8

graph_file=""
read_files=""
prefix_output="TakeABreak_Expe"-`date +"%Y-%m-%d-%H:%M"`


#######################################################################
#################### GET OPTIONS                #######################
#######################################################################
while getopts "hr:g:k:c:s:p:r:S:" opt; do
case $opt in

h)
help
exit
;;

r)
echo "use read set: $OPTARG"
read_files=$OPTARG
;;

g)
echo "consider the graph $OPTARG already constructed"
graph_file=$OPTARG
;;

k)
echo "use k=$OPTARG"
k=$OPTARG
;;

c)
echo "use LCT=$OPTARG (local complexity threshold)"
loc_cmpx_LCT=$OPTARG
;;

s)
echo "use LCS=$OPTARG (Longuest Common Subsequence percentage of the read threshold)"
lcs_restriction_percentage=$OPTARG
;;

r)
echo "use reverse_tolerance=$OPTARG" 
reverse_tolerance=$OPTARG
;;

p)
echo "use prefix=$OPTARG for output files" 
prefix_output=$OPTARG
;;


S)
echo "use solide kmer threshold=$OPTARG" 
nks=$OPTARG
;;

\?)
echo "Invalid option: -$OPTARG"
exit 1
;;

:)
echo "Option -$OPTARG requires an argument." 
exit 1
;;
esac
done
#######################################################################
#################### END GET OPTIONS            #######################
#######################################################################

if [ "$graph_file" == "" ]&&[ "$read_files" == "" ]; then
	echo "-g or -r is mandatory"
	help
	exit
fi

if [ "$graph_file" != "" ]&&[ "$read_files" != "" ]; then
	echo "-g and -r options are not compatible, exit"
	help
	exit
fi	

if [ "$k" != -1 ]&&[ "$graph_file" != "" ]; then
	echo "-k and -g option are not compatible, exit"
	help
	exit
fi

if [ "$nks" != -1 ]&&[ "$graph_file" != "" ]; then
	echo "-S and -g option are not compatible, exit"
	help
	exit
fi

if [ "$k" == -1 ]; then
	k=31
fi

if [ "$nks" == -1 ]; then
	nks=3
fi



log=$prefix_output.log

param_resume=$(resume_parameters)   # or result=`myfunc`


echo "\033[1m" $param_resume "\033[0m"
echo $param_resume > $log

## Generating DBG
if [ "$graph_file" == "" ]; then
	graph_file=$prefix_output
echo "\033[31m ./bin/dbgh5 -in $read_files -nks $nks -kmer-size $k -out $graph_file >> $log\033[0m"
T="$(date +%s)"
 ./bin/dbgh5 -in $read_files -nks $nks -kmer-size $k -out $graph_file  >> $log
T="$(($(date +%s)-T))"
echo "\033[32m generating the graph took ${T} seconds\033[0m"
echo "generating the graph took ${T} seconds" >> $log
# -nks : solidity threshold	
fi
	
	
outfile=$prefix_output".r2sv"
fastafile=$prefix_output".fasta"
## Running BreakFinder
echo "\033[31m ./bin/TakeABreak -i $graph_file -o $prefix_output -r $reverse_tolerance -m $lcs_restriction_percentage -c $loc_cmpx_LCT\033[0m"
T="$(date +%s)"
./bin/TakeABreak -i $graph_file -o $outfile -r $reverse_tolerance -m $lcs_restriction_percentage -c $loc_cmpx_LCT  >> $log
tail -n 2 $log
T="$(($(date +%s)-T))"
echo "\033[32m finding inversions took ${T} seconds\033[0m"
echo "finding inversions took ${T} seconds" >> $log

echo "\033[31m sort -u $outfile > temp"
sort -u $outfile > temp
echo " awk -F \";\" '{print \">inv_\" NR \" \"\$1\"\\\n\"\$2\"\\\n\"\$3\"\\\n\"\$4\"\\\n\"\$5\"\\\n\"\$6\"\\\n\"\$7\"\\\n\"\$8}' temp > $fastafile"
awk -F ";" '{print ">inv_" NR " "$1"\n"$2"\n"$3"\n"$4"\n"$5"\n"$6"\n"$7"\n"$8}' temp > $fastafile

echo " rm -f temp $outfile \033[0m"
rm -f temp
rm -f $outfile


echo "\033[1m\t results in " $fastafile "\033[0m"

