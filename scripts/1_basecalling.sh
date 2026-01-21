#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --mail-user=aanakamo@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --exclude=phoenix-09
#SBATCH --mem=400gb
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=24
#SBATCH --output=slurm-%j.out             # Standard output and error log
#SBATCH --error=slurm-%j.err              # Standard output and error log
#SBATCH --time=3-00:00:00

set -o pipefail
set -e
set -u
set -o xtrace

POD5DIR=$1      ## directory of POD5 files
OUTDIR=$2       ## directory to output files to
MODEL=$3        ## fast,hac,sup

dorado basecaller ${MODEL} ${POD5DIR} \
        --device cuda:all --kit-name SQK-NBD114-24 | dorado demux \
        --emit-fastq --no-classify -o ${OUTDIR}/demux
