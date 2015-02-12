#!/bin/sh
#************

## Claire Lemaitre
## 13/01/2015

## Several tests to make sure we obtain the same results

## dir for the results obtained with the current version
outdir="newtests"

## dir with the good results, obtained with a clean version
gold="gold"

## file with the evaluation of each test PASS/FAILED :
stats="test-results.txt"
>$stats

if [ -d  $outdir ]; then
  rm -rf $outdir
fi

mkdir $outdir

min=1
max=6
## The user can make only the test i iff only one argument
if (($# == 1)); then
    min=$1
    max=$1
fi

## The user can make only the test from I to j iff only two arguments
if (($# == 2)); then
    min=$1
    max=$2
fi


## Test1 : toy example
## -------------------
num=1
if (( $min<=$num )) && (($max>=$num)); then
echo "Test$num : toy example..."
test1="test$num"

../build/TakeABreak -in ../data/toy_example_reads.fasta,../data/toy_example_with_inv_reads.fasta -out $outdir/$test1 > $outdir/$test1.out

../build/ext/gatb-core/bin/dbginfo -in $outdir/$test1.h5 > $outdir/$test1.h5.info

## evaluation

dif=$(diff $outdir/$test1.fasta $gold/$test1.fasta | wc -l | grep -o "[0-9]\+")
if [ "$dif" == "0" ]; then
    echo "Test$num : PASS" >> $stats
else
    echo "Test$num : FAILED result files differ" >> $stats
fi

fi


## Test2 : toy example from the graph file
## ---------------------------------------
num=2
if (( $min<=$num )) && (($max>=$num)); then
echo "Test$num : toy example from graph input..."
test2="test$num"

../build/TakeABreak -graph $outdir/$test1.h5 -out $outdir/$test2 > $outdir/$test2.out

## evaluation
dif=$(diff $outdir/$test2.fasta $gold/$test2.fasta | wc -l | grep -o "[0-9]\+")
if [ "$dif" == "0" ]; then
    echo "Test$num : PASS" >> $stats
else
    echo "Test$num : FAILED result files differ" >> $stats
fi
fi
## remove heavy file
rm -f $outdir/$test1.h5



## Test3 : on a large dataset = ch22
## ---------------------------------------
num=3
if (( $min<=$num )) && (($max>=$num)); then
echo "Test$num : larger dataset ch22..."
ch22dir="/Users/clemaitr/DATA/NGSdata/TakeABreak/ch22_new"
test3="test$num"

../build/TakeABreak -in $ch22dir/humch22c_reads.fasta,$ch22dir/humch22c.inv_reads.fasta -out $outdir/$test3 > $outdir/$test3.out

../build/ext/gatb-core/bin/dbginfo -in $outdir/$test3.h5 > $outdir/$test3.h5.info

## evaluation
dif=$(diff $outdir/$test3.fasta $gold/$test3.fasta | wc -l | grep -o "[0-9]\+")

nb1=$(grep "nb_solid_kmers" $outdir/$test3.out | grep -o "[0-9]\+")
nb2=$(grep "nb_branching_nodes" $outdir/$test3.out | grep -o "[0-9]\+")

nbG1=$(grep "nb_solid_kmers" $gold/$test3.out | grep -o "[0-9]\+")
nbG2=$(grep "nb_branching_nodes" $gold/$test3.out | grep -o "[0-9]\+")


if [ "$dif" == "0" ] && [ "$nb1" == "$nbG1" ] && [ "$nb2" == "$nbG2" ]; then
echo "Test$num : PASS" >> $stats
else
echo "Test$num : FAILED different solutions or different nb of solid and branching ($nb1 vs $nbG1 solid and $nb2 vs $nbG2 branching)" >> $stats
fi

fi


## Test4 : on ch22 with sensitive parameters
## ---------------------------------------
num=4
if (( $min<=$num )) && (($max>=$num)); then
echo "Test$num : ch22 with sensitive parameters..."
ch22dir="/Users/clemaitr/DATA/NGSdata/TakeABreak/ch22_new"
test4="test$num"

../build/TakeABreak -graph $outdir/$test3.h5 -out $outdir/$test4 -max-sim 95 -repeat 15 -lct 1000 > $outdir/$test4.out

## evaluation
dif=$(diff $outdir/$test4.fasta $gold/$test4.fasta | wc -l | grep -o "[0-9]\+")
if [ "$dif" == "0" ]; then
    echo "Test$num : PASS" >> $stats
else
    echo "Test$num : FAILED result files differ" >> $stats
fi
fi

## remove heavy file
rm -f $outdir/$test3.h5



## Test5 : on large reads, testing k=127
## ---------------------------------------
num=5
if (( $min<=$num )) && (($max>=$num)); then
echo "Test$num : large reads k=127..."
largeKfq="/Users/clemaitr/DATA/NGSdata/LbFV_raw_data/Varaldi_LbFV_TTCAGC_L001_R1_001.fastq"
test5="test$num"

../build/TakeABreak -in $largeKfq -out $outdir/$test5 -kmer-size 127 -abundance-min 500 -repeat 110 -max-sim 99 > $outdir/$test5.out

## evaluation :

dif=$(diff $outdir/$test5.fasta $gold/$test5.fasta | wc -l | grep -o "[0-9]\+")

nb1=$(grep "nb_solid_kmers" $outdir/$test5.out | grep -o "[0-9]\+")
nb2=$(grep "nb_branching_nodes" $outdir/$test5.out | grep -o "[0-9]\+")

nbG1=$(grep "nb_solid_kmers" $gold/$test5.out | grep -o "[0-9]\+")
nbG2=$(grep "nb_branching_nodes" $gold/$test5.out | grep -o "[0-9]\+")


if [ "$dif" == "0" ] && [ "$nb1" == "$nbG1" ] && [ "$nb2" == "$nbG2" ]; then
    echo "Test$num : PASS" >> $stats
else
    echo "Test$num : FAILED different solutions or different nb of solid and branching ($nb1 vs $nbG1 solid and $nb2 vs $nbG2 branching)" >> $stats
fi

fi
## remove heavy file
rm -f $outdir/$test5.h5



## Test6 : on large reads, testing k=200 (warning need to recompile)
## ---------------------------------------
num=6
if (( $min<=$num )) && (($max>=$num)); then
echo "Test$num : large reads k=200..."
largeKfq="/Users/clemaitr/DATA/NGSdata/LbFV_raw_data/Varaldi_LbFV_TTCAGC_L001_R1_001.fastq"
test6="test$num"

../build/TakeABreak -in $largeKfq -out $outdir/$test6 -kmer-size 200 -abundance-min 300 > $outdir/$test6.out

## evaluation :

dif=$(diff $outdir/$test6.fasta $gold/$test6.fasta | wc -l)

nb1=$(grep "nb_solid_kmers" $outdir/$test6.out | grep -o "[0-9]\+")
nb2=$(grep "nb_branching_nodes" $outdir/$test6.out | grep -o "[0-9]\+")

nbG1=$(grep "nb_solid_kmers" $gold/$test6.out | grep -o "[0-9]\+")
nbG2=$(grep "nb_branching_nodes" $gold/$test6.out | grep -o "[0-9]\+")


if [ "$dif" == "0" ] && [ "$nb1" == "$nbG1" ] && [ "$nb2" == "$nbG2" ]; then
    echo "Test$num : PASS" >> $stats
else
    echo "Test$num : FAILED different solutions or nb of solid and branching ($nb1 vs $nbG1 solid and $nb2 vs $nbG2 branching)" >> $stats
fi

fi
## remove heavy file
rm -f $outdir/$test6.h5



## AT THE END
echo "*****************************"
echo "Test results : "
cat $stats
echo "*****************************"
