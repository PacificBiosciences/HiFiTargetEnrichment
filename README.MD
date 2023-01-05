# HiFi target enrichment workflow

## About
This is a snakemake workflow for analyzing targeted HiFi sequence datasets. The workflow was designed for Twist gene panels sequenced on Sequel IIe or Revio systems. Learn more in this [App Note](https://www.pacb.com/wp-content/uploads/Application-Brief-HiFi-Target-Enrichment-Best-Practices.pdf). Please note that this workflow is still in development. There may be updates to the input / output files and workflow steps which may affect the behavior of the program.

## Getting started

### Dependencies
This workflow uses `conda` to install dependencies and set up the compute environment. The workflow is demonstrated with `conda v22+` on a `SLURM` cluster and we recommend 16 cores and 48GB of memory on a single machine. Otherwise, use at your own risk.

### Installation

```bash
# create the base conda environment
$ conda create \
    --channel conda-forge \
    --channel bioconda \
    --prefix ./conda_env \
    python=3.9 snakemake mamba pysam lockfile

# activate the base conda environment
$ conda activate ./conda_env

# clone the github repo
$ git clone https://github.com/PacificBiosciences/HiFiTargetEnrichment.git workflow
```

### Quick start
We recommend using this GRCh38 [reference](https://downloads.pacbcloud.com/public/dataset/HiFiTE_Revio/Nov_2022/TwistAllianceDarkGene/human_GRCh38_no_alt_analysis_set.fasta). Refer to the FAQ for more details on input files and formats.

```bash
# create a couple directories
$ mkdir reference annotation cluster_logs

# drop your reference.fasta and reference.fasta.fai into reference directory
# and adjust the [ref][fasta|index] paths in workflow/config.yaml

# Run full workflow starting with demultiplexing for batch <batch_name>
$ sbatch workflow/run_snakemake.sh <batch_name> <biosample_csv> <hifi_reads> <target_bed> [<probe_bed>]
# Use the <target_bed> as the probes bed if you do not have the probes.  This will allow for generation of HS_Metrics

```

### Understanding output

#### BAM and VCF files

Haplotagged aligned BAM files and VCF files are available for each sample in the `whatshap` subdirectory for each sample. These files were produced using Deep Variant and Whatshap. See FAQ for a detailed directory structure. 

If [pbsv](https://github.com/PacificBiosciences/pbsv) is used in the pipeline, a VCF for structral variants can be found in the `pbsv` subdirectory for each sample.

#### Demultiplexed output

The `demux` directory contains unmapped BAM files for each sample. These files are named with the barcode ID. The directory also contains all of the output for [lima](https://github.com/PacificBiosciences/barcoding), the demultiplexing software for PacBio.

#### Target enrichment stats

See the `stats` directory for many useful tables and summary figures. The file `hs_metrics_consolidated_quickview.tsv` is a good place to start to understand run performance. Definitions of each statistics can be found [here](https://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics).

`read_categories.png` gives a graphical summary of read lengths and the percentage of reads that are demultiplexed, PCR duplicates, unmapped, on and off target.


## FAQ

#### 1) What is the format for the biosample_csv file?

This file provides a mapping of barcodes to samples. It is used for sample demultiplexing, to trim adapter sequences, and to rename files with sample ID. The file is comma-delimited where the first column is a pair of barcode IDs (forward and reverse, must match the headers in the barcodes fasta) and the second column in a sample ID (must be unique). The csv must contain a header row, but the header names can be arbitrary.

Here is an example:
```text
Barcode,Biosample
UDI_1_A01_f--UDI_1_A01_r,HG002
UDI_2_B01_f--UDI_2_B01_r,HG003
UDI_3_C01_f--UDI_3_C01_r,HG004
UDI_4_D01_f--UDI_4_D01_r,HG005
```

#### 2) What is the file format for HiFi reads?

Please use hifi_reads.bam, the native format for PacBio HiFi data, which is an unmapped BAM file ("ubam"). Learn more [here](https://ccs.how/faq/reads-bam#how-to-get-hifi-reads).

#### 3) Can I change the barcode file?

By default the pipeline uses the Tru-Seq barcodes from Twist. You can change you barcode file by providing a different fasta file in `workflow/barcodes`. For example, the barcode file "Sequel_96_barcodes_v1.IDT_pad.fasta" was used for this [dataset](https://www.pacb.com/wp-content/uploads/Portik_ASHG_2021.pdf) and in this [preprint](https://www.biorxiv.org/content/10.1101/2022.10.25.513646v1). Adjust the path to the barcode file in `workflow/config.yaml` and make sure you update your `biosample_csv` file.

#### 4) What is the format of the `target.bed` file?

This is a standard [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file that specifies the coordinates of the regions of interest. It is used by [HS Metrics](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/collecthsmetrics.html) to compute target coverage, on-target mapping rate, and other statistics.

Here is an example `targets.bed`:
```text
chr1	155234451	155244627	GBA
chr5	70049638	70078522	SMN2
chr5	70925030	70953942	SMN1
chr22	42126498	42130810	CYP2D6
```

#### 5) What is a `probes.bed` file? What if I don't have one?

The `probes.bed` file specifies the coordinate for the probes use to prepare the target capture library. The file usually contains hundreds or even thousands of probes depending on the size of the panel. You may not have access to the probe design due to it being proprietary. In this case, you can use the `target.bed` file in place of the `probes.bed` file (you will use `target.bed` twice when you call the workflow command).

Here is an example of `probes.bed`:

```text
chr1	48037726	48037846
chr1	48038919	48039039
chr1	48040111	48040231
chr1	48041516	48041636
```

#### 6) What software is used in the pipeline?
```text
python v3.9+
conda v22+
lima v2.5+
pbmarkdups v1+
pbmm2 v1.7+
deep variant v1.4+
whatshap v1.1+
glnexus v1.4.1+
pbsv v2.8+
htslib v1.13+
hsmetrics v1.3.1
pharmCAT
pangu v0.2.1
```
#### 7) What are the steps in the workflow?
1) demultiplex HiFi reads with [lima](https://github.com/pacificbiosciences/barcoding/)
2) Mark PCR duplicate HiFi reads by sample with [pbmarkdups](https://github.com/PacificBiosciences/pbmarkdup/)
3) align HiFi reads to reference with [pbmm2](https://github.com/PacificBiosciences/pbmm2)
4) call small variants with [DeepVariant v1.4.0](https://github.com/google/deepvariant)
5) phase small variants with [WhatsHap](https://github.com/whatshap/whatshap)
6) haplotag aligned BAMs with WhatsHap
7) call SV with pbsv
8) jointly call all variants (excl pbsv) with glnexus
9) [optionally] call PGx star (*) alleles with [PharmCAT](https://pharmcat.org/) and [pangu](https://github.com/PacificBiosciences/pangu)
10) [optionally] annotate output gVCF with dbsnp or other database containing variant IDs
11) Run some QC reports including hsMetrics


#### 8) What is the full directory structure?

```text
.
├── cluster_logs  # slurm stderr/stdout logs
├── reference
│   ├── reference.chr_lengths.txt  # cut -f1,2 reference.fasta.fai > reference.chr_lengths.txt
│   ├── reference.fasta
│   └── reference.fasta.fai
├── batches
│   └── <batch_name>  # batch_id
│       ├── benchmarks/  # cpu time per task
│       ├── demux/  # demultiplexed hifi reads
│       ├── glnexus/  # intermediate cohort vcf files
│       ├── logs/  # per-rule stdout/stderr logs
│       ├── picard/  # interval lists for hsmetrics
│       ├── stats/  # batch-wide collated reports, including HS Metrics summary
│       ├── whatshap_cohort/  # joint-called, phased SNV 
│       ├── merged_gvcf/  # [optional] annotated gvcf with all batch samples included 
│       ├── <sample_id 1>/  # per-sample results, one for each sample
│       :        ...
│       └── <sample_id n>/  # per-sample results, one for each sample
│           ├── coverage/ # read coverage by target beds
│           ├── deepvariant/  # intermediate DV vcf, incl g.vcf per sample
│           ├── hs_metrics/ # picard hsmetrics for this sample
│           ├── markdup/ # unaligned reads with PCR dups marked
│           ├── pangu/  # [optional] HiFi CYP2D6 star calling results
│           ├── pbsv/  # structural variant calls
│           ├── pharmcat/  # [optional] pharmcat results
│           ├── read_metrics/  # per-read information
│           └── whatshap/  # phased small variants; merged haplotagged alignments
│ 
└── workflow  # clone of this repo
         
```

## Detailed run guidance

```bash
# create the base conda environment
$ conda create \
    --channel conda-forge \
    --channel bioconda \
    --prefix ./conda_env \
    python=3.9 snakemake mamba pysam lockfile

# activate the base conda environment
$ conda activate ./conda_env

# clone the github repo
$ git clone https://github.com/PacificBiosciences/HiFiTargetEnrichment.git workflow

# create a couple directories
$ mkdir reference annotation cluster_logs

# drop your reference.fasta and reference.fasta.fai into reference 
# and adjust the [ref][fasta|index] paths in workflow/config.yaml

# Annotation [optional]
# drop your annotation file and index into annotation
# and adjust the [annotate][variants] path in workflow/config.yaml
# ensure that [annotate][gVCF] is set to True in workflow/config.yaml

# PharmCAT [optional]
# Set [pharmcat][run_analysis] to True in workflow/config.yaml

# run the full workflow including demux/markdup/mapping from a single HiFi movie for batch <batch_name>
# Use the <target_bed> as the probes bed if you do not have the probes.  This will allow for generation of HS_Metrics
$ sbatch workflow/run_snakemake.sh <batch_name> <biosample_csv> <hifi_reads> <target_bed> [<probe_bed>]

# run just variant calling and phasing for a set of bams following demux/markdup/mapping on SL
# <hifi_reads> can be a directory of bams or a textfile with one bam path per line (fofn)
$ sbatch workflow/run_snakemake_SLmapped.sh <batch_name> <hifi_reads> <target_bed> [<probe_bed>]
```
