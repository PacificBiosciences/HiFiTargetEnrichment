# HiFi Targeted probe-capture workflow

## Workflow steps
1) demulltiplex HiFi reads with [lima](https://github.com/pacificbiosciences/barcoding/)
2) Mark PCR duplicate HiFi reads by sample with [pbmarkdups](https://github.com/PacificBiosciences/pbmarkdup/)
3) align HiFi reads to reference with [pbmm2](https://github.com/PacificBiosciences/pbmm2)
4) call small variants with [DeepVariant v1.4.0](https://github.com/google/deepvariant)
5) phase small variants with [WhatsHap](https://github.com/whatshap/whatshap)
6) haplotag aligned BAMs with WhatsHap
7) jointly call all variants with glnexus
8) [optionally] Run some QC reports including hsMetrics


## Directory structure within basedir

```text
.
├── cluster_logs  # slurm stderr/stdout logs
├── reference
│   ├── reference.chr_lengths.txt  # cut -f1,2 reference.fasta.fai > reference.chr_lengths.txt
│   ├── reference.fasta
│   └── reference.fasta.fai
├── batches
│   └── <batch_name>  # bpatch_id
│       ├── benchmarks/  # cpu time per task
│       ├── demux/  # demultiplexed hifi reads
│       ├── glnexus/  # intermediate cohort vcf files
│       ├── logs/  # per-rule stdout/stderr logs
│       ├── picard/  # interval lists for hsmetrics
│       ├── stats/  # batch-wide collated reports, including HS Metrics summary
│       ├── whatshap_cohort/  # joint-called, phased SNV 
│       ├── <sample_id 1>/  # per-sample results, one for each sample
│       :        ...
│       └── <sample_id n>/  # per-sample results, one for each sample
│           ├── aligned/  # intermediate
│           ├── counts/ # read counts by target beds
│           ├── coverage/ # read coverage by target beds
│           ├── deepvariant/  # intermediate DV vcf, incl g.vcf per sample
│           ├── hs_metrics/ # picard hsmetrics for this sample
│           ├── markdup/ # unaligned reads with PCR dups marked
│           ├── read_metrics/  # per-read information
│           └── whatshap/  # phased small variants; merged haplotagged alignments
│ 
└── workflow  # clone of this repo
         
```

## To run the pipeline

```bash
# create the base conda environment
$ conda create \
    --channel bioconda \
    --channel conda-forge \
    --prefix ./conda_env \
    python=3 snakemake=6.15.5 mamba lockfile

# activate the base conda environment
$ conda activate ./conda_env

# clone the github repo
$ git clone https://github.com/jrharting/hifi-probe-capture-workflow.git workflow

# create a couple directories
$ mkdir reference cluster_logs

# drop your reference.fasta and reference.fasta.fai into reference and adjust the path in workflow/config.yaml

# run the workflow for batch <batch_name>
$ sbatch workflow/run_snakemake.sh <batch_name> <biosample_csv> <hifi_reads> <target_bed> <probe_bed>
```
