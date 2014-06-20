#!/bin/sh
#************

function resume_parameters {
	
	printf '%s\n\n' "***Parameters values***"
	printf '%s\n' "Output prefix: "$prefix_output
	
	if [ "$graph_file" == "" ]; then
		printf '%s\n' "De Bruijn Graph construction "	
		printf '\t%s\n' "k="$k
		printf '\t%s\n' "min_coverage="$nks
	fi
	
	printf '%s\n' "Inversion detection"
	printf '\t%s\n' "reverse_tolerance="$reverse_tolerance
	printf '\t%s\n' "LCT="$loc_cmpx_LCT
	printf '\t%s\n' "LCS="$lcs_restriction_percentage
    printf '%s\n\n' "***********************"
	
	
}

function help {
printf '%s\n' "TakeBreak.sh,"
printf '%s\n' "Usage: ./TakeBreak.sh options"
printf '\t%s\n'  "In/out Options (-r or -g mandatory):"
printf '\t\t%s\n'  "-r STRING. List of reads separated by comma ',' without space. Note that reads may be in fasta or fastq format, gzipped or not. Example: -r data/toy_example_reads.fasta,data/toy_example_with_inv_reads.fasta"
printf '\t\t\t%s\n'  "Incompatible with -g option"
printf '\t\t%s\n'  "-g STRING. Name of the already existing graph file (.h5)."
printf '\t\t\t%s\n'  "Incompatible with -r option"
printf '\t\t%s\n'  "-p STRING. All out files will start with this prefix. Default: \"TakeABreak_Expe-date\""
printf '\t%s\n'  "De bruijn Graph Options:"
printf '\t\t%s\n'  "-k INT. Set the length of used kmers. Incompatible with -g option. Default=31."
printf '\t\t%s\n'  "-S INT. Set the Solidity of used kmers (minimal number of occurrences a k-mer should have not to be treated as sequencing error). Incompatible with -g option. Default=3."
printf '\t%s\n'  "Inversion detection options:"
printf '\t\t%s\n'  "-c INT. LCT (local complexity threshold): Defaults: 100."
printf '\t\t%s\n'  "-m INT: max_sim: max similarity percentage: Inversions with a and b' (or u and v') whose longuest common subsequence size is bigger than k*(this value)/100 are discarded. Defaults: 80"
printf '\t\t%s\n'  "-r INT: (optimization parameter lower=longer, higher=false negatives) max repeated size suffix of u and v': Defaults: 8"
printf '\t\t%s\n'  "-a INT: number of cores to be used for computation : Defaults: 0, ie. all available cores will be used"
printf '%s\n' "Any further question: read the readme file or contact us: claire.lemaitre@inria.fr"
}

k=-1
nks=-1

loc_cmpx_LCT=100
lcs_restriction_percentage=80
reverse_tolerance=8
nb_cores=0

graph_file=""
read_files=""
prefix_output="TakeABreak_Expe"-`date +"%Y-%m-%d-%H:%M"`


#######################################################################
#################### GET OPTIONS                #######################
#######################################################################
while getopts "hr:g:k:c:s:p:r:a:S:" opt; do
case $opt in

h)
help
exit
;;

r)
printf '%s\n' "use read set: $OPTARG"
read_files=$OPTARG
;;

g)
printf '%s\n' "consider the graph $OPTARG already constructed"
graph_file=$OPTARG
;;

k)
printf '%s\n' "use k=$OPTARG"
k=$OPTARG
;;

c)
printf '%s\n' "use LCT=$OPTARG (local complexity threshold)"
loc_cmpx_LCT=$OPTARG
;;

s)
printf '%s\n' "use LCS=$OPTARG (Longuest Common Subsequence percentage of the read threshold)"
lcs_restriction_percentage=$OPTARG
;;

r)
printf '%s\n' "use reverse_tolerance=$OPTARG" 
reverse_tolerance=$OPTARG
;;

a)
printf '%s\n' "use nb_cores=$OPTARG"
nb_cores=$OPTARG
;;

p)
printf '%s\n' "use prefix=$OPTARG for output files" 
prefix_output=$OPTARG
;;


S)
printf '%s\n' "use solide kmer threshold=$OPTARG" 
nks=$OPTARG
;;

\?)
printf '%s\n' "Invalid option: -$OPTARG"
exit 1
;;

:)
printf '%s\n' "Option -$OPTARG requires an argument." 
exit 1
;;
esac
done
#######################################################################
#################### END GET OPTIONS            #######################
#######################################################################

if [ "$graph_file" == "" ]&&[ "$read_files" == "" ]; then
	printf '%s\n' "-g or -r is mandatory"
	help
	exit
fi

if [ "$graph_file" != "" ]&&[ "$read_files" != "" ]; then
	printf '%s\n' "-g and -r options are not compatible, exit"
	help
	exit
fi	

if [ "$k" != -1 ]&&[ "$graph_file" != "" ]; then
	printf '%s\n' "-k and -g option are not compatible, exit"
	help
	exit
fi

if [ "$nks" != -1 ]&&[ "$graph_file" != "" ]; then
	printf '%s\n' "-S and -g option are not compatible, exit"
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


printf '%s\n' "$param_resume"
printf '%s\n' "$param_resume" > $log

## Generating DBG
if [ "$graph_file" == "" ]; then
	graph_file=$prefix_output
printf '%s\n' "./bin/dbgh5 -in $read_files -nks $nks -kmer-size $k -out $graph_file >> $log"
T="$(date +%s)"
 ./bin/dbgh5 -in $read_files -nks $nks -kmer-size $k -out $graph_file  >> $log
T="$(($(date +%s)-T))"
printf '%s\n' "generating the graph took ${T} seconds"
printf '%s\n' "generating the graph took ${T} seconds" >> $log
# -nks : solidity threshold	
fi
	
	
outfile=$prefix_output".out"
fastafile=$prefix_output".fasta"
## Running BreakFinder
printf '%s\n' "./bin/TakeABreak -i $graph_file -o $outfile -r $reverse_tolerance -m $lcs_restriction_percentage -c $loc_cmpx_LCT"
T="$(date +%s)"
./bin/TakeABreak -i $graph_file -o $outfile -r $reverse_tolerance -m $lcs_restriction_percentage -c $loc_cmpx_LCT -a $nb_cores >> $log
tail -n 2 $log
T="$(($(date +%s)-T))"
printf '%s\n' "finding inversion motifs took ${T} seconds"
printf '%s\n' "finding inversion motifs took ${T} seconds" >> $log

#printf '%s\n' "sort -u $outfile > temp"
sort -u $outfile > temp
#printf '%s\n' " awk -F \";\" '{print \">inv_\" NR \" \"\$1\"\\\n\"\$2\"\\\n\"\$3\"\\\n\"\$4\"\\\n\"\$5\"\\\n\"\$6\"\\\n\"\$7\"\\\n\"\$8}' temp > $fastafile"
awk -F ";" '{print ">inv_" NR " "$1"\n"$2"\n"$3"\n"$4"\n"$5"\n"$6"\n"$7"\n"$8}' temp > $fastafile

nbInv=$(cat temp | wc -l)

#printf '%s\n' "rm -f temp $outfile"
rm -f temp
rm -f $outfile



printf '%s\n' $nbInv" inversions were found"
printf '%s\n' "results in "$fastafile

printf '%s\n' $nbInv" inversions were found" >> $log
printf '%s\n' "results in "$fastafile >> $log

