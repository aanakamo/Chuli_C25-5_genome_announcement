#!/bin/bash
#SBATCH --job-name=funannotate
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200gb
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=aanakamo@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=slurm-%j.out             # Standard output and error log
#SBATCH --error=slurm-%j.err              # Standard output and error log

### Annotate C. huliohia isolate C25-2 genome with funannotate

genome_name=${1}    ## ie. C25-5
working_dir=${2}    ## ie. /private/groups/corbettlab/aanakamo/ROD_HAWAII/R491_GENOME_ANNOUNCEMENT/000_FINAL_ASSEMBLY/3_annotate/C25-5_funannotate_predict

# Run antiSMASH (optional) using the website: https://fungismash.secondarymetabolites.org/#!/start

## run phobius
cd ${working_dir}
phobius_out=Ceratocystis_huliohia_C25-5.proteins.phobius.out
# phobius.pl -short predict_results/Ceratocystis_huliohia_C25-5.proteins.fa > ${phobius_out}

## run interpro scan
cd ${working_dir}
mkdir -p ${genome_name}_interpro
# interproscan.sh -i predict_results/Ceratocystis_huliohia_C25-5.proteins.fa -d ${genome_name}_interpro -cpu 32 

## Annotate Genome
cd ..
phobius_out=${working_dir}/Ceratocystis_huliohia_C25-5.proteins.phobius.out
iprscan_out=${working_dir}/${genome_name}_interpro/Ceratocystis_huliohia_C25-5.proteins.fa.xml
antismash_out=${working_dir}/Ceratocystis_huliohia_C25-5.antismash.gbk
sbt_file=/private/groups/corbettlab/aanakamo/ROD_HAWAII/R491_GENOME_ANNOUNCEMENT/7_funannotate/template.sbt
source activate /private/home/aanakamo/miniforge3/envs/funannotate
funannotate annotate -i ${genome_name}_funannotate_predict --cpus 8 --sbt ${sbt_file} \
            --phobius ${phobius_out} --iprscan ${iprscan_out} --antismash ${antismash_out}
conda deactivate
