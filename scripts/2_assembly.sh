#!/bin/bash
#SBATCH --job-name=hifiasm_assembly
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

### Assemble C. huliohia isolate C25-2 genome from raw Nanopore reads

fastq=${1}             ## ie. R491.merged.rmdup.fastq.gz
genome_name=${2}       ## ie. C25-5
working_dir=${3}       

cd ${working_dir}
mkdir -p 1_hifiasm; cd 1_hifiasm
source activate /private/home/aanakamo/miniforge3/envs/hifiasm
threads=${SLURM_CPUS_PER_TASK}
assembly=${genome_name}.bp.p_ctg.fa
echo "*** Assembling ${genome_name} based on nanopore reads: ${fastq}"
hifiasm -o ${genome_name} -t ${threads} ${fastq} --ont
echo "*** Extracting fa from gfa for ${genome_name}"
awk '/^S/{print ">"$2;print $3}' ${genome_name}.bp.p_ctg.gfa > ${assembly}
conda deactivate

## clean/remove duplicate contigs
echo "*** Cleaning/removing duplicate contigs"
source activate /private/home/aanakamo/miniforge3/envs/funannotate
cleaned_fa=${genome_name}.bp.p_ctg.funannotate_clean.fa
funannotate clean -i ${assembly} -o ${cleaned_fa} > funannotate_clean.log
conda deactivate

## correct misassemblies around telomeres
## 1. Generate coverage per base tsv
echo "*** Correcting misassemblies at telomeres"
sorted_bam=${genome_name}.bp.p_ctg.funannotate_clean.sorted.bam
primary_bam=${genome_name}.bp.p_ctg.funannotate_clean.sorted.primary.bam
cov_tsv=${genome_name}.bp.p_ctg.funannotate_clean.sorted.primary.COV.tsv
minimap2 -ax map-ont ${cleaned_fa} ${fastq} | \
samtools view -b | \
samtools sort -o ${sorted_bam}
samtools view -b -F 2308 ${sorted_bam} > ${primary_bam}
samtools index ${primary_bam}
samtools depth -aa ${primary_bam} > ${cov_tsv}
## 2. Fix misassemblies
python ~/CD_Lab/rod_hawaii/scripts/telomere_curation.py \
        --fasta ${cleaned_fa} --coverage ${cov_tsv} > telomere_curation.log
## 3. apply the breaks and trimming
source activate /private/home/aanakamo/miniforge3/envs/bedtools
corrected_fa=${genome_name}.bp.p_ctg.funannotate_clean.telo_cur.fa
bedtools getfasta \
  -fi ${cleaned_fa} \
  -bed keep_regions.bed \
  -fo ${corrected_fa}
conda deactivate

## blast against available C4194 C. huli mito genome (MT331822.1) and remove mito contigs
source activate /private/home/aanakamo/miniforge3/envs/blast
makeblastdb -in MT331822.1_Ceratocystis_huliohia_C4194.fasta -parse_seqids -dbtype nucl -out MT331822.1
blastn -query ${corrected_fa} -db MT331822.1 -out C25-5_vs_MT331822.1.BLAST.out -task blastn \
            -evalue 1e-5 -perc_identity 80 -qcov_hsp_perc 50 \
            -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
conda deactivate
cat C25-5_vs_MT331822.1.BLAST.out | awk '{ print $1; }' | sort | uniq > mito_contigs.txt
no_mito_fa=${genome_name}.bp.p_ctg.funannotate_clean.telo_cur.no_mito.fa
source activate /private/home/aanakamo/miniforge3/envs/seqkit
seqkit grep -v -n -f mito_contigs.txt ${corrected_fa} > ${no_mito_fa}
conda deactivate

## sort and rename contigs
echo "*** Sorting and renaming contigs"
source activate /private/home/aanakamo/miniforge3/envs/funannotate
sorted_fa=${working_dir}/1_hifiasm/${genome_name}.bp.p_ctg.funannotate_clean_sort.telo_cur.fa
funannotate sort -i ${no_mito_fa} -o ${sorted_fa} -b contig --minlen 0

## repeat-mask before funannotate
echo "*** Masking repeats"
masked_fa=${working_dir}/1_hifiasm/${genome_name}.bp.p_ctg.funannotate_clean_sort_masked.telo_cur.fa
funannotate mask -i ${sorted_fa} -o ${masked_fa} --cpus ${SLURM_CPUS_PER_TASK}
conda deactivate

echo "   Final assembly: ${masked_fa}"
echo -e "   Checking ${genome_name} assembly stats:\n"
faSize ${masked_fa}
fa_size=${genome_name}.bp.p_ctg.faSize
faSize -detailed ${masked_fa} > ${fa_size}
echo -e "\n   Contig names and sizes here: ${fa_size}\n"

## QC
echo "*** Starting QC ***"
masked_fa=${working_dir}/1_hifiasm/${genome_name}.bp.p_ctg.funannotate_clean_sort_masked.telo_cur.fa

## coverage
cd ${working_dir}
mkdir -p 2_QC; cd 2_QC
echo "*** Coverage"
out_sam=${genome_name}.ONT.minimap2.sam
minimap2 -x map-ont -a -t ${SLURM_CPUS_PER_TASK} ${masked_fa} ${fastq} > ${out_sam}
out_bam=${genome_name}.ONT.minimap2.bam
samtools view -h -Sb ${out_sam} | samtools sort - -O 'bam' -T samtools_temp -o ${out_bam}
filt_bam=${genome_name}.ONT.minimap2.filt.bam
samtools view -h -b ${out_bam} | samtools sort - -O 'bam' -T samtools_temp -o ${filt_bam}
samtools index ${filt_bam}
source activate /private/home/aanakamo/miniforge3/envs/mosdepth
mosdepth_out=${genome_name}
mosdepth -t ${SLURM_CPUS_PER_TASK} ${mosdepth_out} ${filt_bam}
conda deactivate

## GC
echo "*** GC content"
source activate /private/home/aanakamo/miniforge3/envs/emboss
gc_out=geecee_${genome_name}
geecee -sequence ${masked_fa} -outfile ${gc_out}
conda deactivate

source activate /private/home/aanakamo/miniforge3/envs/matplotlib
python ~/CD_Lab/rod_hawaii/scripts/plot_gc_cov_len.py ${gc_out} ${mosdepth_out}.mosdepth.summary.txt ${genome_name}
conda deactivate

## tidk
echo "*** Telomere identification"
source activate /private/home/aanakamo/miniforge3/envs/tidk
tidk search --string AACCCT --output ${genome_name} --dir . ${masked_fa}
tidk plot --tsv ${genome_name}_telomeric_repeat_windows.tsv 
conda deactivate

## Quast
echo "*** QUAST"
source activate /private/home/aanakamo/miniforge3/envs/quast
echo "Running quast for ${genome_name}"
quast ${masked_fa} --nanopore ${fastq} -t 32 -o QUAST_${genome_name} --circos --k-mer-stats --glimmer --conserved-genes-finding --rna-finding --est-ref-size 30000000 --fungus 
conda deactivate

## BUSCO
echo "*** BUSCO"
source activate /private/home/aanakamo/miniforge3/envs/busco
busco -i ${masked_fa} -o busco_${genome_name}_sord -l sordariomycetes_odb10 -m genome -c 32 -f --download_path busco_${genome_name}_sord
busco -i ${masked_fa} -o busco_${genome_name}_asco -l ascomycota_odb10 -m genome -c 32 -f --download_path busco_${genome_name}_asco
busco -i ${masked_fa} -o busco_${genome_name}_fungi -l fungi_odb10 -m genome -c 32 -f --download_path busco_${genome_name}_fungi
conda deactivate

## predict genes
echo "*** Initial gene prediction ***"
cd ${working_dir}
mkdir -p 3_annotate; cd 3_annotate
source activate /private/home/aanakamo/miniforge3/envs/funannotate
predicted_out=${genome_name}_funannotate_predict
funannotate predict -i ${masked_fa} -o ${predicted_out} --species "Ceratocystis huliohia" --strain ${genome_name}
## Add species parameters to database:
# funannotate species -s ceratocystis_huliohia_c25-5 -a ${genome_name}_funannotate_predict/predict_results/ceratocystis_huliohia_c25-5.parameters.json
conda deactivate

cd ${working_dir}
gff_file=${working_dir}/3_annotate/${predicted_out}/predict_results/Ceratocystis_huliohia_C25-5.gff3
gc_repeat_out=${working_dir}/1_hifiasm/${genome_name}.bp.p_ctg.funannotate_clean_sort_masked.telo_cur.GC_REPEAT.tsv
python ~/CD_Lab/rod_hawaii/scripts/windowed_gc_repeat.py \
            -g ${masked_fa} -w 50000 -a ${gff_file} -o ${gc_repeat_out}
