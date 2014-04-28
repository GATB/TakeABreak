#!/bin/bash
#==============================================================================
#                    TakeABreak scalability tests
#==============================================================================
#
# Usage:
#     On the IGRIDA frontend (igrida-oar-frontend), launch:
#     	   oarsub -S ./job_scalability_igrida.sh
#==============================================================================
#
#------------------------------------------------------------------------------
# Job parameters
#------------------------------------------------------------------------------
#OAR -n TakeABreak
#OAR -l {cluster='lambda'}/nodes=1,walltime=4:00:00
#OAR -O /temp_dd/igrida-fs1/cdeltel/run.%jobid%.out
#OAR -E /temp_dd/igrida-fs1/cdeltel/run.%jobid%.out

set -xv

# >>>>>>>>>>>>>>>>>>>>>> edit me >>>>>>>>>>>>>>>>>>>>>>

tool=TakeABreak
#exp_code=block_nbits10
#exp_code=block_nbits14
exp_code=validation_before_commit
data_dir=/temp_dd/igrida-fs1/cdeltel/bioinfo/Misc

time_stamp="`date +"%Y-%m-%d-%H:%M"`"

scratch_dir=/temp_dd/igrida-fs1/cdeltel/scratchdir_${tool}_${exp_code}_${time_stamp}

fasta_list="aphid_400000.fa aphid_662451seq.fa"
#fasta_list="aphid_400000.fa"
graph_generation=1
compilation=1

nbCoresMax="`cat /proc/cpuinfo|grep processor|wc -l`"

graph_list="aphid_400000.h5 aphid_662451seq.h5"
#graph_list="aphid_400000.h5"

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
tools_dir=/udd/cdeltel/bioinfo/anr-gatb/git-gatb-tools/gatb-tools/tools/$tool

#-- Compilation
cd $tools_dir/build

if [ $compilation -eq 1 ] ; then
	make clean
	make -j 8
fi

#-- scratch_dir preparation 
mkdir -p $scratch_dir
cd $scratch_dir

ln -sf $data_dir/aphid_662451seq.fa .
ln -sf $data_dir/aphid_400000.fa .

cp $tools_dir/build/TakeABreak               . || { echo "Error 1"; exit; }
cp $tools_dir/build/ext/gatb-core/bin/dbgh5  . || { echo "Error 2"; exit; }

#-- Graph generation
if [ $graph_generation -eq 1 ] ; then
	for fasta_file in $fasta_list; do
		./dbgh5 -in $fasta_file
		rm -f cfp
	done
fi

#-- TakeABreak
nbCores=1

while [ $nbCores -le $nbCoresMax ]; do  
	
	for graph_file in $graph_list; do
		
		# ========================= run =========================
		base_name="`basename $graph_file .h5`"
		prefix_output="TakeABreak_Expe_${base_name}_${nbCores}cores"
		
		log=$prefix_output.log
		outfile=$prefix_output".r2sv"
		fastafile=$prefix_output".fasta"
		
		T="$(date +%s)"
		#./bin/TakeABreak -i $graph_file -o $outfile -r $reverse_tolerance -m $lcs_restriction_percentage -c $loc_cmpx_LCT  >> $log
		./TakeABreak -a $nbCores -i $graph_file -o $outfile >> $log
		tail -n 2 $log
		time_elapsed="$(($(date +%s)-T))"
					
		echo "finding inversions took ${time_elapsed} seconds" >> $log

		# ==================== post-processing ==================								
		sort -u $outfile > temp
		awk -F ";" '{ print ">inv_" NR " "$1"\n"$2"\n"$3"\n"$4"\n"$5"\n"$6"\n"$7"\n"$8 }' temp > $fastafile
		
		# ===================== scalability =====================
		[ $nbCores -eq 1 ] && { eval "time_sequential_${base_name}=$time_elapsed"; }
		eval "speedup=\`echo \$time_sequential_${base_name} / $time_elapsed | bc -l\`"
		echo "$base_name $nbCores $time_elapsed $speedup" >> ${base_name}-scalability.txt
											
	done
 
	(( nbCores = nbCores+1 ))
done

#-- Scalability plots

gnuplot << EOF
set terminal postscript
set output "scalability.ps"
set grid
plot "aphid_400000-scalability.txt" u 2:4 w lp, "aphid_662451seq-scalability.txt" u 2:4 w lp, x
EOF

#egrep 'inversions.*found' *log







