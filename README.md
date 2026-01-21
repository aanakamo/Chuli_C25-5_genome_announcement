# Complete genome sequence of *Ceratocystis huliohia*, a fungal pathogen of the native ‘ōhi‘a tree in Hawai‘i

[![BioProject](https://img.shields.io/badge/BioProject-PRJNA1400107-blue)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1400107/)

This repository contains scripts used to produce and analyze the complete *de novo* genome sequence of *Ceratocystis huliohia*, one of two fungal pathogens causing the Rapid ‘Ōhi‘a Death disease of the native ‘ōhi‘a tree in Hawai‘i. This assembly was generated using long-read Nanopore sequencing of *C. huliohia* isolate C25-5, collected on the island of Maui in April 2025.


## Introduction
Hawai‘i’s native forests are threatened by Rapid ‘Ōhi‘a Death (ROD), a disease of the keystone and culturally significant native ‘ōhi‘a tree (*Metrosideros polymorpha*) species. ROD was first characterized in 2014, and is caused by two novel, non-native fungal pathogen species, *Certaocystis lukuohia* and *huliohia*. Genomic resources for these pathogens remain limited, and while *C. lukuohia* has a high quality, long-read reference genome (GCF_044167205.1), no such genome has been available for *C. huliohia*. Here, we present a complete, long-read genome assembly of *C. huliohia* isolate C25-5 collected in April 2025 at Pu‘u Kukui Elementary School on the island of Maui. This resource will serve as a representative long-read genome for the Asian-Australian Clade of *Ceratocystis*, and facilitate further investigation into the origins of *C. huliohia* in Hawai‘i.

**Isolate Information**
| <!-- -->    | <!-- -->    |
|------|-------|
| Isolate | C25-5 |
| Sample date | April 11, 2025 |
| Sample location | Pu‘u Kukui Elementary School, island of Maui |
| Latitude | 20.88 |
| Longitude | -156.51 |

## Methods Overview

### Data Generation
- **Isolate:** C. huliohia C25-5 
- **Growth conditions** 
    - liquid media (2% malt extract, 0.2% yeast extract)
    - circular shaker at 150 RPM
    - 25°C incubation for six days
- **Preparation of Material**
    - collected mycelium by filtration through miracloth (Millipore #475855)
    - lyophilized for two days
    - ground to a powder with beads on a Geno/Grinder 2010 at 1600 RPM for two minutes
- **DNA Extraction**: Wizard HMW DNA Extraction Kit (Promega #A2920)
- **Library Prep:** Native Barcoding Kit V14 (Oxford Nanopore, SQK-NBD114-24)
- **Sequencing**: Nanopore MinION Mk1B, R10 version flow cell
- **Resulting Raw Reads:** 1.57M (N50 3.23kb)

### Data Analysis
Scripts in this repository are linked below for each step, and code is provided to show how they were run. Necessary programs and their versions are also provided.

-  [1_basecalling.sh]()
    - **Basecalling:** Dorado v1.1.1, sup model
-  [2_assembly.sh]()
    - **Assembly:** Hifiasm v0.25.0-r726, ONT mode
    - **Removal of duplicate contigs:** Funannotate v1.8.17
    - **Telomere Identification:** tidk v0.2.65
    - **Correction of misassemblies around telomeres:** [script here]()
    - **Mitochondrial contig identification:** blastn v2.16.0+ against MT331822.1
    - **Repeat masking:** Funannotate v1.8.17
    - **Quality Assessment:** QUAST v5.3.0 and BUSCO v6.0.0
3. [3_annotation.sh]()
    - **Gene annotation:** Funannotate v1.8.17

## Results
Statistics for the final C. huliohia C25-5 genome assembly (GCA_054512535.1) are shown in the table below. Our assembly demonstrates high completeness, with a BUSCO score indicating 99.7% of fungal orthologs present. The majority of bases (97.6%) and coding sequences (98.6%) reside in the eight largest contigs (>1Mb), falling in the expected chromosome number range of ~7-9 for Ceratocystis species. All eight contigs >1Mb have a telomeric repeat present on at least one end, with five of them being chromosome-level and having telomeres on both ends.

**C25-5 Assembly summary**
| <!-- -->    | <!-- -->    |
|------|-------|
| Total number of contigs | 16 |
| Total assembly length | 30,006,571 bp (30.0 Mb) |
| Number of near-chromosome-level contigs (>1Mb) | 8 |
| Total length of >1Mb contigs | 29,299,846 bp (29.3 Mb) |
| Percentage of bases in >1Mb contigs | 97.6% |
| Contig N50 | 4,906,038 bp (4.9 Mb) |
| Contig L50 | 3 |
| Avg. coverage | 57.5× |
| GC content | 48.9% |
| CDS content | 37.2% |
| Repeat content | 4.77% (1,430,314 bp) |
| BUSCO completeness (fungi_odb10) | 99.7% |
| BUSCO completeness (ascomycota_odb10) | 98.2% |
| BUSCO completeness (sordariomycetes_odb10) | 93.5% |

**C25-5 Annotation summary**
| <!-- -->    | <!-- -->    |
|------|-------|
| Number of genes | 7,356 (7,006 mRNA, 350 tRNA) |
| Avg. gene length | 1,710.6 bp |
| Number of CDS | 17,890 |
| Number of CDS in >1Mb contigs | 17,645 (98.6%) |
| Avg. exon length | 557.7 bp |
| Total annotations | 35,473 |
| Secretome annotations | 558 |
| Transmembrane annotations | 1,330 |
| Secondary metabolite biosynthesis gene clusters | 11 |
| Biosynthetic enzymes | 11 |
| smCOGs | 23 |
| CAZYmes | 192 |


## Data Availability
- BioProject:
- BioSample:
- Raw sequencing data: 
- WGS Accession: 
- Nuclear Genome Assembly & Annotation:
- Mitochondrial Genome Assembly & Annotation:
- Analysis code: https://github.com/aanakamo/Chuli_C25-5_genome_announcement (this GitHub repository)

## Citation
If you use this genome assembly, associated data, or analysis code in your research, please cite:

~~~
Nakamoto A, Keith L, Yu Q, Sugiyama L, Wu X, Luiz B, Villalun M, Jacobs J, Corbett-Detig R, 
Cisneros A, Heath H, Shanks C, Okamoto F, Alburo A, Henricson K, Lan E, Moore H, Seligmann W, 
Zybina Y. Complete genome sequence of Ceratocystis huliohia, a fungal pathogen of the native 
‘ōhi‘a tree in Hawai‘i. (preprint). 2026.
~~~

## Acknowledgments
We acknowledge the University of California Santa Cruz Genomics Institute for providing support for this project, including use of the Phoenix computational cluster. We also thank Rion Parsons for supporting our use of the Hummingbird computational cluster provided by the University of California Santa Cruz. We acknowledge the USDA-ARS, Daniel K. Inouye U.S. Pacific Basin Agricultural Research Center in Hilo, HI for use of laboratory space, and thank Jon Suzuki, Katelin Branco, and Andrew Paresa for their support. We thank the University of California Santa Cruz Baskin Engineering Lab Support team for providing laboratory space, equipment, and support. We also thank Christopher Vollmers and James Letchinger each for the use of their computers to perform Nanopore sequencing, and Honey Mekonen for technical support. Funding for this project was provided by the NIH (T32 HG012344) awarded to AN, AC, HH, CS, and FO, and the NSF-GRFP awarded to AN.

## Contact
- Corresponding Author: Anne Nakamoto (aanakamo@ucsc.edu)
- Lab Website: [Corbett-Detig Lab at UC Santa Cruz](https://corbett-lab.github.io/)
