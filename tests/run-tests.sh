#!/bin/sh
#************

## Claire Lemaitre
## 13/01/2015

## Several tests to make sure we obtain the same results

## dir with the results obtained with the current version
outdir="newtests"

## dir with the good results, obtained with a clean version
gold="gold"

## file with test PASS/FAILED :
stats="test-results.txt"

if [ -d  $outdir ]; then
  rm -rf $outdir
fi

mkdir $outdir

## Test1 : toy example
## -------------------
echo "Test1 : toy example..."
test1="test1"

../build/TakeABreak -in ../data/toy_example_reads.fasta,../data/toy_example_with_inv_reads.fasta -out $outdir/$test1 > $outdir/$test1.out

../build/ext/gatb-core/bin/dbginfo -in $outdir/$test1.h5 > $outdir/$test1.h5.info

## evaluation

dif=$(diff $outdir/$test1.fasta $gold/$test1.fasta | wc -l | grep -o "[0-9]\+")
if [ "$dif" == "0" ]; then
    echo "Test1 : PASS" > $stats
else
    echo "Test1 : FAILED result files differ" > $stats
fi




## Test2 : toy example from the graph file
## ---------------------------------------
echo "Test2 : toy example from graph input..."
test2="test2"

../build/TakeABreak -graph $outdir/$test1.h5 -out $outdir/$test2 > $outdir/$test2.out

## evaluation
dif=$(diff $outdir/$test2.fasta $gold/$test2.fasta | wc -l | grep -o "[0-9]\+")
if [ "$dif" == "0" ]; then
    echo "Test2 : PASS" >> $stats
else
    echo "Test2 : FAILED result files differ" >> $stats
fi

## remove heavy file
rm -f $outdir/$test1.h5


## Test3 : on a large dataset = ch22
## ---------------------------------------
echo "Test3 : larger dataset ch22..."
ch22dir="/Users/clemaitr/DATA/NGSdata/TakeABreak/ch22_new"
test3="test3"

../build/TakeABreak -in $ch22dir/humch22c_reads.fasta,$ch22dir/humch22c.inv_reads.fasta -out $outdir/$test3 > $outdir/$test3.out

../build/ext/gatb-core/bin/dbginfo -in $outdir/$test3.h5 > $outdir/$test3.h5.info

## evaluation
dif=$(diff $outdir/$test3.fasta $gold/$test3.fasta | wc -l | grep -o "[0-9]\+")
if [ "$dif" == "0" ]; then
    echo "Test3 : PASS" >> $stats
else
    echo "Test3 : FAILED result files differ" >> $stats
fi


## Test4 : on ch22 with sensitive parameters
## ---------------------------------------
echo "Test4 : ch22 with sensitive parameters..."
ch22dir="/Users/clemaitr/DATA/NGSdata/TakeABreak/ch22_new"
test4="test4"

../build/TakeABreak -graph $outdir/$test3.h5 -out $outdir/$test4 -max-sim 95 -repeat 15 -lct 1000 > $outdir/$test4.out

## evaluation
dif=$(diff $outdir/$test4.fasta $gold/$test4.fasta | wc -l | grep -o "[0-9]\+")
if [ "$dif" == "0" ]; then
    echo "Test4 : PASS" >> $stats
else
    echo "Test4 : FAILED result files differ" >> $stats
fi

## remove heavy file
rm -f $outdir/$test3.h5



## Test5 : on large reads, testing k=127
## ---------------------------------------
echo "Test5 : large reads k=127..."
largeKfq="/Users/clemaitr/DATA/NGSdata/LbFV_raw_data/Varaldi_LbFV_TTCAGC_L001_R1_001.fastq"
test5="test5"

../build/TakeABreak -in $largeKfq -out $outdir/$test5 -kmer-size 127 -abundance 500 -repeat 110 -max-sim 99 > $outdir/$test5.out

## evaluation :

dif=$(diff $outdir/$test5.fasta $gold/$test5.fasta | wc -l | grep -o "[0-9]\+")

nb1=$(grep "nb_solid_kmers" $outdir/$test5.out | grep -o "[0-9]\+")
nb2=$(grep "nb_branching_nodes" $outdir/$test5.out | grep -o "[0-9]\+")

nbG1=$(grep "nb_solid_kmers" $gold/$test5.out | grep -o "[0-9]\+")
nbG2=$(grep "nb_branching_nodes" $gold/$test5.out | grep -o "[0-9]\+")


if [ "$dif" == "0" ] && [ "$nb1" == "$nbG1" ] && [ "$nb2" == "$nbG2" ]; then
    echo "Test5 : PASS" >> $stats
else
    echo "Test5 : FAILED solutions or different nb of solid and branching ($nb1 vs $nbG1 solid and $nb2 vs $nbG2 branching)" >> $stats
fi

## remove heavy file
rm -f $outdir/$test5.h5



## Test6 : on large reads, testing k=200 (warning need to recompile)
## ---------------------------------------
echo "Test6 : large reads k=200..."
largeKfq="/Users/clemaitr/DATA/NGSdata/LbFV_raw_data/Varaldi_LbFV_TTCAGC_L001_R1_001.fastq"
test6="test6"

../build/TakeABreak -in $largeKfq -out $outdir/$test6 -kmer-size 200 -abundance 300 > $outdir/$test6.out

## evaluation :

dif=$(diff $outdir/$test6.fasta $gold/$test6.fasta | wc -l)

nb1=$(grep "nb_solid_kmers" $outdir/$test6.out | grep -o "[0-9]\+")
nb2=$(grep "nb_branching_nodes" $outdir/$test6.out | grep -o "[0-9]\+")

nbG1=$(grep "nb_solid_kmers" $gold/$test6.out | grep -o "[0-9]\+")
nbG2=$(grep "nb_branching_nodes" $gold/$test6.out | grep -o "[0-9]\+")


if [ "$dif" == "0" ] && [ "$nb1" == "$nbG1" ] && [ "$nb2" == "$nbG2" ]; then
    echo "Test6 : PASS" >> $stats
else
    echo "Test6 : FAILED different solutions or nb of solid and branching ($nb1 vs $nbG1 solid and $nb2 vs $nbG2 branching)" >> $stats
fi

## remove heavy file
rm -f $outdir/$test6.h5



## AT THE END
echo "*****************************"
echo "Test results : "
cat $stats
echo "*****************************"
