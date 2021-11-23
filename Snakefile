import re
from pathlib import Path


configfile: "workflow/config.yaml"
## prefix every task command with:
# set -o pipefail  # trace ERR through pipes
# umask 002  # group write permissions
# export TMPDIR={config['tmpdir']}  # configure temp directory
# export SINGULARITY_TMPDIR={config['tmpdir']}  # configure temp directory
shell.prefix(f"set -o pipefail; umask 002; export TMPDIR={config['tmpdir']}; export SINGULARITY_TMPDIR={config['tmpdir']}; ")

sample = config['sample']
ref = config['ref']['shortname']
print(f"Processing sample {sample} with reference {ref}.")

# scan smrtcells/ready for inputs
# uBAMs have priority over FASTQs if both are available
# samples and files are expected to match one of the following patterns:
ubam_pattern = re.compile(r'smrtcells/ready/(?P<sample>[A-Za-z0-9_-]+)/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6}).(ccs|hifi_reads).bam')
ubam_dict = {}
fastq_pattern = re.compile(r'smrtcells/ready/(?P<sample>[A-Za-z0-9_-]+)/(?P<movie>m\d{5}[Ue]?_\d{6}_\d{6}).fastq.gz')
fastq_dict = {}
for infile in Path('smrtcells/ready').glob('**/*.bam'):
    ubam_match = ubam_pattern.search(str(infile))
    if ubam_match:
        # create a dict-of-dict to link samples to movie context to uBAM filenames
        ubam_dict[ubam_match.group('movie')] = str(infile)
for infile in Path('smrtcells/ready').glob('**/*.fastq.gz'):
    fastq_match = fastq_pattern.search(str(infile))
    if fastq_match:
        # create a dict-of-dict to link samples to movie context to FASTQ filenames
        fastq_dict[fastq_match.group('movie')] = str(infile)

# create a list of movies from the uBAMs and FASTQs
movies = list(set(list(ubam_dict.keys()) + list(fastq_dict.keys())))
# create a list of aBAMs that will be generated
abams = [f"samples/{sample}/aligned/{movie}.{ref}.bam" for movie in movies]
abam_dict = {movie: f"samples/{sample}/aligned/{movie}.{ref}.bam" for movie in movies}

# build a list of targets
targets = []
include: 'rules/common.smk'

# add ubams and fastqs to targets so that md5sums will be written
targets.extend(ubam_dict.values())
targets.extend(fastq_dict.values())

# call small variants with DeepVariant
# phase small variants with WhatsHap
include: 'rules/deepvariant.smk'
include: 'rules/whatshap.smk'
targets.extend([f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.{suffix}"
                for suffix in ['vcf.gz', 'vcf.gz.tbi', 'g.vcf.gz', 'g.vcf.gz.tbi',
                               'visual_report.html', 'vcf.stats.txt']])
targets.extend([f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.{suffix}"
                for suffix in ['phased.vcf.gz', 'phased.vcf.gz.tbi', 'phased.gtf',
                               'phased.tsv', 'phased.blocklist',
                               'haplotagged.bam', 'haplotagged.bam.bai']])


ruleorder: pbmm2_align_ubam > pbmm2_align_fastq


rule all:
    input: targets
