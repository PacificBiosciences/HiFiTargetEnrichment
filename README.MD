# HiFi Targeted probe-capture workflow

## Workflow steps
1) demulltiplex HiFi reads with [lima](https://github.com/pacificbiosciences/barcoding/)
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


## Directory structure within basedir

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

## To run the pipeline

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

# drop your reference.fasta and reference.fasta.fai into reference and adjust the [ref][fasta|index] paths in workflow/config.yaml

# Annotation [optional]
# drop your annotation file and index into annotation and adjust the [annotate][variants] path in workflow/config.yaml
# ensure that [annotate][gVCF] is set to True

# PharmCAT [optional]
# Set [pharmcat][run_analysis] to True

# run the full workflow including demux/markdup/mapping from a single HiFi movie for batch <batch_name>
# Use the <target_bed> as the probes bed if you do not have the probes.  This will allow for generation of HS_Metrics
$ sbatch workflow/run_snakemake.sh <batch_name> <biosample_csv> <hifi_reads> <target_bed> [<probe_bed>]

# run just variant calling and phasing for a set of bams following demux/markdup/mapping on SL
# <hifi_reads> can be a directory of bams or a textfile with one bam path per line (fofn)
$ sbatch workflow/run_snakemake_SLmapped.sh <batch_name> <hifi_reads> <target_bed> [<probe_bed>]
```
