import re

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

# checkpoint for samples that fail to yield data in demux
def _get_demuxed_samples( wildcards ):
    '''Some samples may not have yield (ie failed), so update samples after demuxing'''
    demuxdir       = checkpoints.demux_ubam.get( **wildcards ).output.odir
    outputBarcodes = glob_wildcards( f'{demuxdir}/demultiplex.{{barcode}}.bam' ).barcode
    return [ barcode2sample[ bc ] for bc in outputBarcodes if bc != 'unbarcoded' ]

targets = []

include: "rules/common.smk"
include: "rules/demux.smk"
include: "rules/preprocess.smk"
include: "rules/pbmm2.smk"
include: "rules/deepvariant.smk"
include: "rules/whatshap.smk"
include: "rules/pbsv.smk"
if config["run_cohort"]:
    include: "rules/glnexus.smk"
if config[ "probes" ] != "None":
    include: "rules/hsmetrics.smk"
if config[ "pharmcat" ][ "run_analysis" ]:
    include: "rules/pharmcat.smk"
    include: "rules/pangu_cyp2d6.smk"
if config[ "annotate" ][ "gVCF" ]:
    include: "rules/annotate.smk"

# DV targets
targets.append(
    lambda wildcards: \
            [
                f"batches/{batch}/{sample}/deepvariant/{sample}.{ref}.deepvariant.{suffix}"
                for suffix in [
                    "vcf.gz",
                    "vcf.gz.tbi",
                    "g.vcf.gz",
                    "g.vcf.gz.tbi",
                    "visual_report.html",
                    "vcf.stats.txt",
                ]
                for sample in _get_demuxed_samples( wildcards)
            ]
)

# WH targets
targets.append(
    lambda wildcards: \
           [
             f"batches/{batch}/{sample}/whatshap/{sample}.{ref}.deepvariant.{suffix}"
             for suffix in [
                 "phased.vcf.gz",
                 "phased.vcf.gz.tbi",
                 "phased.gtf",
                 "phased.tsv",
                 "phased.blocklist",
                 "haplotagged.bam",
                 "haplotagged.bam.bai",
             ]
             for sample in _get_demuxed_samples( wildcards )
           ]
)

# pbsv targets
targets.append(
    lambda wildcards:
        [
         f"batches/{batch}/{sample}/pbsv/{sample}.{ref}.pbsv.vcf"
         for sample in _get_demuxed_samples( wildcards )
        ]
)

# QC extras
if config['QC']['runQC']:
    include: "rules/qc_cov.smk"
    include: "rules/qc_ext.smk"

ruleorder: deepvariant_postprocess_variants > tabix_vcf

rule all:
    input:
        targets
