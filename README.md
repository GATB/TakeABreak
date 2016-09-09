# TakeABreak 
					  
| **Linux** | **Mac OSX** |
|-----------|-------------|
[![Build Status](https://ci.inria.fr/gatb-core/view/TakeABreak/job/tool-takeabreak-build-debian7-64bits-gcc-4.7/badge/icon)](https://ci.inria.fr/gatb-core/view/TakeABreak/job/tool-takeabreak-build-debian7-64bits-gcc-4.7/) | [![Build Status](https://ci.inria.fr/gatb-core/view/TakeABreak/job/tool-takeabreak-build-macos-10.9.5-gcc-4.2.1/badge/icon)](https://ci.inria.fr/gatb-core/view/TakeABreak/job/tool-takeabreak-build-macos-10.9.5-gcc-4.2.1/)

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

# What is TakeABreak?

TakeABreak detects inversion breakpoints without a reference genome by looking for fixed size topological patterns in the de Bruijn graph. TakeABreak requires small memory configuration and processes data quickly. For example,  Illumina reads simulated at 2x40x coverage from human chromosome 22 can be processed  in  2 hours and 1GB of memory.

Lemaitre C., Ciortuz L. and Peterlongo P. (2014) [Mapping-Free and Assembly-Free Discovery of Inversion Breakpoints from Raw NGS Reads](http://link.springer.com/chapter/10.1007%2F978-3-319-07953-0_10).  AlCoB 2014 (Tarragona), Lecture Notes in Computer Science vol. 8542, pp. 119--130.

# Installation instructions

## Requirements

CMake 2.6+; see http://www.cmake.org/cmake/resources/software.html

c++ compiler; compilation was tested with gcc and g++ version>=4.5 (Linux) and clang version>=4.1 (Mac OSX).

## Getting the latest source code with git

    # get a local copy of source code
    git clone --recursive https://github.com/GATB/TakeABreak.git
    
    # compile the code an run a simple test on your computer
    cd TakeABreak
    sh INSTALL
    # the binary file is located in directory build/bin/
    ./build/bin/MindTheGap -help
    # run a simple example
    ./build/bin/TakeABreak -in tests/data/toy_example_reads.fasta,tests/data/toy_example_with_inv_reads.fasta

Note: when updating your local repository with `git pull`, if you see that thirdparty/gatb-core has changed, you have to run also : `git submodule update`. 

## Installing a stable release

Retrieve a binary archive file from one of the official TakeABreak releases (see "Releases" tab on the Github web page, or https://colibread.inria.fr/software/takeabreak/ for releases older than 1.1.2); file name is `TakeABreak-vX.Y.Z-bin-Linux.tar.gz` (for Linux) or `TakeABreak-vX.Y.Z-bin-Darwin.tar.gz` (for MacOs).

    tar -zxf TakeABreak-vX.Y.Z-bin-Darwin.tar.gz
    cd TakeABreak-vX.Y.Z-bin-Darwin
    chmod u+x bin/TakeABreak

    # run a simple example
    ./bin/TakeABreak -in data/toy_example_reads.fasta,data/toy_example_with_inv_reads.fasta


In case the software does not run appropriately on your system, you should consider to install it from its source code. Retrieve the source archive file `TakeABreak-vX.Y.Z-Source.tar.gz`.

    tar -zxf TakeABreak-vX.Y.Z-Source.tar.gz
    cd TakeABreak-vX.Y.Z-Source
    sh INSTALL
    # the binary file is located in directory build/bin/
    ./build/bin/TakeABreak -help
    # run a simple example
    ./build/bin/TakeABreak -in tests/data/toy_example_reads.fasta,tests/data/toy_example_with_inv_reads.fasta


# USER MANUAL	 
								
## Description

TakeABreak detects inversion breakpoints directly from raw NGS reads, without the need of any reference genome and without de novo assembling the genomes. Its implementation has a very limited memory footprint (less than 6GB for analyzing a full human NGS dataset) and acceptable runtime.
	
		
## Usage
TakeABreak now comes as a single executable, combining the de bruijn graph creation and the inversion breakpoint detection algorithm.

1. **Basic command lines**

        TakeABreak (-in <reads.fq> | -graph <graph.h5>) [-out filePrefix] [options]
        #To get help and see all options:
        TakeABreak -help
	
2. **Input data**

    If one or several read sets are provided (option `-in`) TakeABreak pipelines the de Bruijn graph creation with the inversion breakpoint detection phase. In this case, the user can provide the de Bruijn graph creation options plus the breakpoint detection options. 

    If only a graph file is provided (option `-graph`) TakeABreak only computes the breakpoint detections based on this graph. Only the breakpoint detection options can thus be provided.

    NOTE: options `-in` and `-graph` are mutually exclusive, and one of these is mandatory.
	
    If the input is composed of several read files, they can be provided as a list of file paths separated by a comma or as a "file of file" (fof), that is a text file containing on each line the path to each read file. Read file format can be fasta, fastq or gzipped. 
		
3. **de Bruijn graph creation options**

    In addition to input read set(s), the de Bruijn graph creation uses two main parameters, `-kmer-size` and `-abundance-min`: 

    * `-kmer-size`: the k-mer size [default '31']. By default, the largest kmer-size allowed is 128. To use k>128, you will need to re-compile TakeABreak with the two following commands in the build directory: `cmake -DKSIZE_LIST="32 64 96 256" ..` and then `make`. To go back to default, replace 256 by 128. Note that increasing the range between two consecutive kmer-sizes in the list can have an impact on the size of the output h5 files (but none on the results).
    * `-abundance-min`: the minimal abundance threshold, k-mers having less than this number of occurrences are discarded from the graph [default 'auto', ie. automatically inferred from the dataset]. If several datasets are given, this parameter can be a list of thresholds, one for each dataset if `solidity-kind` is set to 'one' or 'all' (see below and section "Details on dealing with several input read files").
    * `-abundance-max`: the maximal abundance threshold, k-mers having more than this number of occurrences are discarded from the graph [default '2147483647' ie. no limit]
    * `-solidity-kind`: the way to consider a solid kmer with several input datasets (sum, one or all) [default 'one']. Details : with 'sum', a kmer is solid if the sum of its abundances in all input datasets respects the abundance-min and abundance-max conditions; with 'one' (resp. 'all') a kmer is solid if its abundance in at least one dataset (resp. in all datasets) is in the interval [min-max]. 
	
4. **Breakpoint detection options**

    The breakpoint detection algorithm uses the following optional parameters:

    * `-lct`: the local complexity threshold, this limits the search, particularly in complex parts of the graph [default '100']. Warning : depending on the graph, increasing this parameter could lead to very long runtime.
    * `-max-sim`: the max similarity percentage, inversions with a and b' (or u and v') whose longest common subsequence size is larger than k*(this value)/100 are discarded [default '80']. Warning : increasing this parameter may lead to numerous false positives (genomic approximate repeats).
    * `-repeat`: the maximal repeat size at the breakpoint (ie. the longest common suffix size of a and reverse complement of b) [default '8']. To be fully effective, the -max-sim parameter should be fixed accordingly. (this option was formerly called reverse tolerance)
			
5. **Computational resources options**
    Additional options are related to computational runtime and memory:

    * `-nb-cores` : number of cores to be used for computation (graph creation and breakpoint detection) [default '0', ie. all available cores will be used].
    * `-max-memory` : max RAM memory for the graph creation (in MBytes)  [default '2000']. Increasing the memory will speed up the graph creation phase.
    * `-max-disk` : max usable disk space for the graph creation (in MBytes)  [default '0', ie. automatically set]. Kmers are counted by writting temporary files on the disk, to speed up the counting you can increase the usable disk space.
		
	
## TakeABreak Output

TakeABreak generates the following output files: 

* a graph file (`.h5`). This is a binary file, to obtain information stored in it, you can use the utility program dbginfo located in your bin directory or in ext/gatb-core/bin/.
* a fasta file (`.fasta`) containing the canonical representations of the detected inversion breakpoints. Each inversion corresponds to 4 entries in the fasta file : the first two correspond to the breakpoint sequences (a-u,v-b) (canonical representation) that should be present in one genome and the last two are the corresponding breakpoint sequences in the other genome (a-revcomp(v),revcomp(u)-b). Additionally, the size of the exact repeat detected at the breakpoints is indicated in each entry header (e.g. |rep_4). 

All the output results are prefixed either by a default name: "TakeABreak_Expe-[date:YY:MM:DD-HH:mm]" or by a user defined prefix (option `-out` of TakeABreak)
	

## Details on dealing with several input read files

Depending on the value of the parameter `-solidity-kind` all input read datasets can be considered as one single dataset (ie. concatenating the files) or as separated datasets. This impacts the step of filtering out kmers with sequencing errors, whether this is done conjointly by summing abundances over all datasets ('sum') or independently for each dataset ('one' or 'all'). This can lead to significant differences in terms of results when numerous individuals are compared or when the sequencing effort is variable between datasets; in these cases, we recommend to treat each dataset independently ('one').

Usually, one dataset corresponds to one individual or one biological condition. If for one individual or condition, one has several read files (for instance, paired read files or more than one sequencing lane), these files should be considered as one in order to increase the read depth (and thus the power of variant detection) per individual. In a typical case with several individuals having each several read files, some files should be concatenated but not all, the solution is to use an arborescent file of file (fof) : in the master fof (given as input to the `-in` option) each line will be the path to another fof representing one individual. The individual fofs will contain paths of files to be concatenated.
	
Example : assume you have 3 individuals (indivA, indivB, indivC) having each two read files (R1 and R2), those files are all located in a directory path/DATA. The command line could be for instance:

	./TakeABreak -in path/DATA/data_3indiv.fof -kmer-size 31 -abundance-min \ 
	3,7,auto -solidity-kind one -out result_3indiv 

The file data_3indiv.fof contains 3 lines such as:

    indivA.fof
    indivB.fof
    indivC.fof

With indivA.fof file, being located in the path/DATA directory and containing the following two lines:

    indivA_readR1.fastq
    indivA_readR2.fastq
		
Note that fof files do not need to be located in the same directory as the files they refer to. If this is not the case, the path (absolute or relative from the location of the fof file) must be added for each file.
	
 
## Full example

This toy example can be run with the provided data.

* Example from raw input reads:

        ./TakeABreak -in tests/data/toy_example_reads.fasta,tests/data/toy_example_with_inv_reads.fasta \
        -out MyFirstTakeABreakExperiment

    NOTE: the input read files are simply separated by comma without spaces.

    This command line first computes the de Bruijn graph (saved in the file MyFirstTakeABreakExperiment.h5) before searching for the 6 artificial inversions contained in the input read sets and finally it outputs the inversion breakpoints in the file MyFirstTakeABreakExperiment.fasta. 
	
* Example from an already created de Bruijn graph 

    (MyFirstTakeABreakExperiment.h5)

        ./TakeABreak -graph MyFirstTakeABreakExperiment.h5 -out MySecondTakeABreakExperiment

    This command line uses the already computed de Bruijn graph (MyFirstTakeABreakExperiment.h5) to detect the 6 artificial inversions and outputs them in MySecondTakeABreakExperiment.fasta. 
		
		
## Utility programs

Either in your `bin/` directory or in `ext/gatb-core/bin/`, you can find additional utility programs:

* `dbginfo`: to get information about a graph stored in a .h5 file
* `dbgh5`: to build a graph from read set(s) and obtain a .h5 file
* `h5dump`: to extract all data stored in a .h5 file
	
# Contact

To contact a developer, request help, or for any feedback on TakeABreak, please use the issue form of github: https://github.com/GATB/TakeABreak/issues

You can see all issues concerning TakeABreak [here](https://github.com/GATB/TakeABreak/issues) and GATB [here](https://www.biostars.org/t/GATB/).

If you do not have any github account, you can also send an email to claire dot lemaitre at inria dot fr
 
