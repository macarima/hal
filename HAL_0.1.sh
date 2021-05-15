#!/bin/bash

#############
# The GNU General Public License version 3.0 (GPLv3)
#
# Copyright (c) 2020 Hui-Su Kim
#
#############
#
# HAL de novo assembly pipeline
#

version=20200607

echo `readlink -f $0`" "$*

echo "version: "${version}

USAGE_short="
*****************************************************
HAL de novo assembly pipeline

USAGE: ./hal_0.1 [options] <path_to_input_reads> <name_of_genome> <name_of_reference_genome>

DESCRIPTION:
This is a script to assemble draft genome with a chromosome-level

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
-m|--medaka_model
-h|--help	Shows this help page.
*****************************************************
"

pipeline=`cd "$( dirname $0)" && pwd`


## default parameter setup

depth="D1"
cpu=1
seq_type="ONT"
gpu="X"
medaka_model="r941_prom_high_g303"

############### HANDLE OPTIONS ###############


while :; do
		case $1 in
				-h|--help)
						echo "$USAGE_short" >&1
						exit 0
		;;
		-d|--depth) OPTARG=$2
						if [ "$OPTARG" == "D1" ] || [ "$OPTARG" == "D2" ]; then
						echo " Depth option $OPTARG is selected." >&1
						assemble_stage=$OPTARG
						else
								echo ":( Wrong syntax for sequencing depth. Exiting!" >&2
						exit 1
						fi
				shift
		;;
		-s|--seq_type) OPTARG=$2
						if [ "$OPTARG" == "ONT" ] || [ "$OPTARG" == "PBCCS" ]; then
						echo " Sequencing platform $OPTARG is selected." >&1
						ecc_stage=$OPTARG
						else
								echo ":( Wrong syntax for sequencing platform. Exiting!" >&2
						exit 1
						fi
				shift
		;;
		-c|--cpu) OPTARG=$2
						re='^[0-9]+$'
						if [[ $OPTARG =~ $re ]]; then
								echo " $OPTARG of cpu will be used for HAL pipeline."
								cpu=$OPTARG
						else
								echo ":( Wrong syntax for cpu parameter value. Using the default value ${cpu}." >&2
						fi
        		shift 
		;;
		-g|--gpu) OPTARG=$2
						if [[ "$OPTARG" == "X" ]]; then
								echo " Medaka error correction with CPU is selected."
								gpu=$OPTARG
						else
								gpu=$OPTARG
								echo " Medaka error correction with GPU is selected. GPU $gpu will be used." >&2
						fi
			shift
		;;
		-m|--medaka_model) OPTARG=$2
						medaka_model=$OPTARG
			shift
		;;
		*)
			break
	esac
	shift
done

############### HANDLE ARGUMENTS ###############

orig_fastq=$1
orig_genome=$2
orig_ref=$3

if [ ! -e "$orig_fastq" ]; then
        echo "$USAGE_short" >&1
        exit 1
fi


if [ "$assemble_stage" == "D1" ] && [ "$ecc_stage" == "ONT" ];then

############### 15x ONT assembly ###############

			wtdbg2 -t ${cpu} -x ont -g 3g -X 100.0 -L 3000 -i ${orig_fastq} -fo ${orig_genome} --edge-min 2 --rescue-low-cov-edges
			wtpoa-cns -t ${cpu} -i ${orig_genome}.ctg.lay.gz -fo ${orig_genome}.ctg.lay.fa

			minimap2 -t ${cpu} -x map-ont ${orig_genome}.ctg.lay.fa ${orig_fastq} > ${orig_genome}.mapping_minimap21x.paf
			racon -t ${cpu} -m 8 -x -6 -g -8 -w 500 ${orig_fastq} ${orig_genome}.mapping_minimap21x.paf ${orig_genome}.ctg.lay.fa > ${orig_genome}.wtdbg2_racon1x.fa

			mkdir medaka_bam
			mini_align -i ${orig_fastq} -r ${orig_genome}.wtdbg2_racon1x.fa -P -m -p medaka_bam/calls_to_draft -t ${cpu}

			if [ "$gpu" == "X" ]; then
					medaka consensus medaka_bam/calls_to_draft.bam medaka_bam/concensus.hdf --batch 80 --threads 40 --model ${medaka_model}
					medaka stitch --jobs ${cpu} medaka_bam/concensus.hdf ${orig_genome}.wtdbg2_racon1x_medaka.fa
			else
					CUDA_VISIBLE_DEVICES=${gpu} medaka consensus medaka_bam/calls_to_draft.bam medaka_bam/concensus.hdf --batch 80 --threads 40 --model ${medaka_model}
					medaka stitch --jobs ${cpu} medaka_bam/concensus.hdf ${orig_genome}.wtdbg2_racon1x_medaka.fa

			fi

			ragoo.py -t ${cpu} -C ${orig_genome}.wtdbg2_racon1x_medaka.fa ${orig_ref}
			cp ragoo_output/ragoo.fasta ${orig_genome}.wtdbg2_racon1x_medaka.chr.fa

			echo "############### Assembly job is done. ###############"

elif [ "$assemble_stage" == "D2" ] && [ "$ecc_stage" == "ONT" ];then

############### 30x ONT assembly ###############

			wtdbg2 -t ${cpu} -x ont -g 3g -X 70.0 -L 5000 -i ${orig_fastq} -fo ${orig_genome} --edge-min 2 --rescue-low-cov-edges
			wtpoa-cns -t ${cpu} -i ${orig_genome}.ctg.lay.gz -fo ${orig_genome}.ctg.lay.fa

			minimap2 -t ${cpu} -x map-ont ${orig_genome}.ctg.lay.fa ${orig_fastq} > ${orig_genome}.mapping_minimap21x.paf
			racon -t ${cpu} -m 8 -x -6 -g -8 -w 500 ${orig_fastq} ${orig_genome}.mapping_minimap21x.paf ${orig_genome}.ctg.lay.fa > ${orig_genome}.wtdbg2_racon1x.fa

			mkdir medaka_bam
			mini_align -i ${orig_fastq} -r ${orig_genome}.wtdbg2_racon1x.fa -P -m -p medaka_bam/calls_to_draft -t ${cpu}

			if [ "$gpu" == "X" ]; then
					medaka consensus medaka_bam/calls_to_draft.bam medaka_bam/concensus.hdf --batch 80 --threads 40 --model ${medaka_model}
					medaka stitch --jobs ${cpu} medaka_bam/concensus.hdf ${orig_genome}.wtdbg2_racon1x_medaka.fa
			else
					CUDA_VISIBLE_DEVICES=${gpu} medaka consensus medaka_bam/calls_to_draft.bam medaka_bam/concensus.hdf --batch 80 --threads 40 --model ${medaka_model}
					medaka stitch --jobs ${cpu} medaka_bam/concensus.hdf ${orig_genome}.wtdbg2_racon1x_medaka.fa

			fi

			ragoo.py -t ${cpu} -C ${orig_genome}.wtdbg2_racon1x_medaka.fa ${orig_ref}
			cp ragoo_output/ragoo.fasta ${orig_genome}.wtdbg2_racon1x_medaka.chr.fa

			echo "############### Assembly job is done. ###############"

elif [ "$assemble_stage" == "D1" ] && [ "$ecc_stage" == "PBCCS" ];then

############### 15x PacBio_CCS assembly ###############

			wtdbg2 -t ${cpu} -x ccs -g 3g -X 100.0 -L 3000 -i ${orig_fastq} -fo ${orig_genome} --edge-min 2 --rescue-low-cov-edges
			wtpoa-cns -t ${cpu} -i ${orig_genome}.ctg.lay.gz -fo ${orig_genome}.ctg.lay.fa

			minimap2 -t ${cpu} -x asm20  ${orig_genome}.ctg.lay.fa ${orig_fastq} > ${orig_genome}.mapping_minimap21x.paf
			racon -t ${cpu} -m 8 -x -6 -g -8 -w 500 ${orig_fastq} ${orig_genome}.mapping_minimap21x.paf ${orig_genome}.ctg.lay.fa > ${orig_genome}.wtdbg2_racon1x.fa

			ragoo.py -t ${cpu} -C ${orig_genome}.wtdbg2_racon1x.fa ${orig_ref}
			cp ragoo_output/ragoo.fasta ${orig_genome}.wtdbg2_racon1x.chr.fa

			echo "############### Assembly job is done. ###############"

else

############### 30x PacBio_CCS assembly ###############

			wtdbg2 -t ${cpu} -x ccs -g 3g -X 70.0 -L 5000 -i ${orig_fastq} -fo ${orig_genome} --edge-min 2 --rescue-low-cov-edges
			wtpoa-cns -t ${cpu} -i ${orig_genome}.ctg.lay.gz -fo ${orig_genome}.ctg.lay.fa

			minimap2 -t ${cpu} -x asm20  ${orig_genome}.ctg.lay.fa ${orig_fastq} > ${orig_genome}.mapping_minimap21x.paf
			racon -t ${cpu} -m 8 -x -6 -g -8 -w 500 ${orig_fastq} ${orig_genome}.mapping_minimap21x.paf ${orig_genome}.ctg.lay.fa > ${orig_genome}.wtdbg2_racon1x.fa

			ragoo.py -t ${cpu} -C ${orig_genome}.wtdbg2_racon1x.fa ${orig_ref}
			cp ragoo_output/ragoo.fasta ${orig_genome}.wtdbg2_racon1x.chr.fa

			echo "############### Assembly job is done. ###############"

fi

