import re

# get config file from run-script args for testing multiple panels
#configfile: "workflow/config.yaml"


## prefix every task command with:
# set -o pipefail  # trace ERR through pipes
# umask 002  # group write permissions
# export TMPDIR={config['tmpdir']}  # configure temp directory
# export SINGULARITY_TMPDIR={config['tmpdir']}  # configure temp directory
shell.prefix(
    f"set -o pipefail; umask 002; export TMPDIR={config['tmpdir']}; export SINGULARITY_TMPDIR={config['tmpdir']}; "
)

batch = config["batch"]
movie = config["ccsReads"]
# Get expected barcodes and sample map from standard barcode_biosample_csv used in SL
# map both directions
# Assumes csv has header (actual names ignored).  
# two columns barcode,sampleId
barcode2sample = dict( row.split(',') 
                       for i,row in enumerate( open( config["biosamples"] ).read().strip().split('\n') )
                       if i > 0 )
sample2barcode = { v:k for k,v in barcode2sample.items() }
ref = config["ref"]["shortname"]
print(f"Processing batch {batch} with reference {ref}.")

targets = [
            f'batches/{batch}/demux/demultiplex.{barcode}.bam'
            for barcode in sample2barcode.values()
           ]


include: "rules/common.smk"
include: "rules/demux.smk"
include: "rules/preprocess.smk"
include: "rules/pbmm2.smk"
include: "rules/deepvariant.smk"
include: "rules/whatshap.smk"
include: "rules/pbsv.smk"
include: "rules/hifiasm.smk"
include: "rules/glnexus.smk"
include: "rules/hsmetrics.smk"

# DV
targets.extend(
    [
        f"batches/{batch}/{sample}/downsampled_{maxreads}/deepvariant/{sample}.{ref}.deepvariant.{suffix}"
        for suffix in [
            "vcf.gz",
            "vcf.gz.tbi",
            "g.vcf.gz",
            "g.vcf.gz.tbi",
            "visual_report.html",
            "vcf.stats.txt",
        ]
        for sample in sample2barcode.keys()
        for maxreads in config['downsample'] + ['all']
    ]
)
# WH
targets.extend(
    [
        f"batches/{batch}/{sample}/downsampled_{maxreads}/whatshap/{sample}.{ref}.deepvariant.{suffix}"
        for suffix in [
            "phased.vcf.gz",
            "phased.vcf.gz.tbi",
            "phased.gtf",
            "phased.tsv",
            "phased.blocklist",
            "haplotagged.bam",
            "haplotagged.bam.bai",
        ]
        for sample in sample2barcode.keys()
        for maxreads in config['downsample'] + ['all']
    ]
)
#SV
targets.extend(
    [
        f"batches/{batch}/{sample}/downsampled_{maxreads}/pbsv/{sample}.{ref}.pbsv.vcf"
        for sample in sample2barcode.keys()
        for maxreads in config['downsample'] + ['all']
    ]
)

# gVCF/cohort
for maxreads in config['downsample'] + ['all']:
    targets.extend(
        [
         f"batches/{batch}/whatshap_cohort/downsampled_{maxreads}/{batch}.{ref}.deepvariant.glnexus.phased.vcf.gz",
         f"batches/{batch}/whatshap_cohort/downsampled_{maxreads}/{batch}.{ref}.deepvariant.glnexus.phased.gtf",
         f"batches/{batch}/whatshap_cohort/downsampled_{maxreads}/{batch}.{ref}.deepvariant.glnexus.phased.tsv",
         f"batches/{batch}/whatshap_cohort/downsampled_{maxreads}/{batch}.{ref}.deepvariant.glnexus.phased.blocklist"
        ]
    )

# HiFiasm
targets.extend(
    [
        f"batches/{batch}/{sample}/downsampled_{maxreads}/hifiasm/{sample}.asm.{ref}.htsbox.vcf.stats.txt"
        for sample in sample2barcode.keys()
        for maxreads in config['downsample'] + ['all']
    ]
)


# QC extras
if config['QC']['runQC']:
    include: "rules/qc.smk"

ruleorder: demux_ubam > demux_fastq
ruleorder: downsample_bam > downsample_fastq
ruleorder: pbmm2_align_ubam > pbmm2_align_fastq
ruleorder: deepvariant_postprocess_variants_round2 > deepvariant_postprocess_variants_round1 > tabix_vcf


rule all:
    input:
        targets,
