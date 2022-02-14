# HiFi DeepVariant + WhatsHap workflow

## Workflow steps

1) align HiFi reads to reference with [pbmm2](https://github.com/PacificBiosciences/pbmm2)
2) call small variants with [DeepVariant](https://github.com/google/deepvariant), using two-pass method (DeepVariant :arrow_right: WhatsHap `phase` :arrow_right: WhatsHap `haplotag` :arrow_right: DeepVariant)
3) phase small variants with [WhatsHap](https://github.com/whatshap/whatshap)
4) haplotag aligned BAMs with WhatsHap and merge

## Directory structure within basedir

```text
.
├── cluster_logs  # slurm stderr/stdout logs
├── reference
│   ├── reference.chr_lengths.txt  # cut -f1,2 reference.fasta > reference.chr_lengths.txt
│   ├── reference.fasta
│   └── reference.fasta.fai
├── samples
│   └── <sample_id>  # sample_id regex: r'[A-Za-z0-9_-]+'
│       ├── whatshap/  # phased small variants; merged haplotagged alignments
│       ├── logs/  # per-rule stdout/stderr logs
│       ├── aligned/  # intermediate
│       ├── deepvariant/  # intermediate
│       ├── deepvariant_intermediate/  # intermediate
│       └── whatshap_intermediate/  # intermediate
├── smrtcells
│   ├── done  # move folders from smrtcells/ready to smrtcells/done to prevent re-processing
│   └── ready
│       └── <sample_id>  # uBAMs or FASTQs per sample
│                        # filename regex: r'm\d{5}[Ue]?_\d{6}_\d{6}).(ccs|hifi_reads).bam' or r'm\d{5}[Ue]?_\d{6}_\d{6}).fastq.gz'
└── workflow  # clone of this repo
```

## To run the pipeline

```bash
# create the base conda environment
$ conda create \
    --channel bioconda \
    --channel conda-forge \
    --prefix ./conda_env \
    python=3 snakemake mamba lockfile

# activate the base conda environment
$ conda activate ./conda_env

# clone the github repo
$ git clone https://github.com/jrharting/hifi-probe-capture-workflow.git workflow

# create a few directories
$ mkdir -p reference smrtcells/capture/{ready,done} cluster_logs

# drop your reference.fasta and reference.fasta.fai into reference and adjust the path in workflow/config.yaml

# copy or link some data into a smrtcells/ready/<sample_id> directory
$ ln -s /path/to/some/hifi/data/m00000_000000_000000.hifi_reads.bam smrtcells/capture/ready/<sample_id>/m00000_000000_000000.hifi_reads.bam
$ ln -s /path/to/some/hifi/data/m11111_111111_111111.hifi_reads.bam smrtcells/capture/ready/<sample_id>/m11111_111111_111111.hifi_reads.bam

# run the workflow for sample <sample_id>
$ sbatch workflow/run_snakemake.sh <sample_id>
```
