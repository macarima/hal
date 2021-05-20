# HAL_assembly pipeline: *De novo* assembly pipeline for a low-depth sequencing of human genome
HAL_assembly uses wtdbg2, racon, medaka, and RaGOO.

## Installation
Prerequisties: wtdbg2, racon, medaka, RaGOO  
These programs must exist on your path.

No need to installation. Just execute 'HAL_0.1.sh or HAL_0.1 binary' on your linux machine.

## Usage
```
*****************************************************
HAL de novo assembly pipeline

USAGE: ./hal_0.1 [options] <path_to_input_reads> <name_of_genome> <name_of_reference_genome>

DESCRIPTION:
This is a script to assemble draft genome as choromsome level with a reference genome

ARGUMENTS:
path_to_input_reads	Specify file path to sequenced long read fastq file
name_of_genome	Specify the name of genome which you want ot use
name_of_reference_genome Specify the name of the reference genome for the chromosome assembly

OPTIONS:
-d|--depth	Specify the depth of the sequenced long read (Default: D2)
	D1 : below 15x
	D2 : 15 ~ 30x

-c|--cpu	Specify number of threads for processing assembly (Default: 1)
-s|--seq_type	Type of long read sequencer. ONT / PacBio Hifi-CCS or Corrected PacBio reads (Default: ONT)
-g|--gpu	GPU usage for Medaka. Specify number of GPU or X (Default: X)
-m|--medaka_model (Default: r941_prom_high_g303)
-h|--help	Shows this help page.
*****************************************************
```

## Output files
On working directory, you can find the final result with extention of '~.chr.fa'.

## Contact
For questions, bugs and issues, please post on issue page of Github.

## Citation
This page: https://github.com/macarima/hal

## Copyright
HAL_assembly is freely available for academic use. It is provided without warranty of any kind, including but not limited to the warranties of merchantability, fitness for a particular purpose and non-infringement. No liability for the software usage is assumed. Redistribution is allowed. Redistribution of modified version is not allowed. 

For commercial use, please contact me.
