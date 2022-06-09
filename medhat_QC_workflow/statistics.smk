# import Lib
############
import os, glob
from snakemake.utils import min_version

############################

# Snake Version
###############
min_version("6.15.5")
##############

# Config File
#############
cpath = "workflow/medhat_QC_workflow/config.yaml"
if os.path.isfile( cpath ):
    configfile: cpath
else:
    sys.exit("Looks like there is no config.yaml file in " + os.getcwd() + " make sure there is one or at least specify one with the --configfile commandline parameter.")
#############


# GET WORKING DIRECTORY DEFAULT IS CURRENT DIRECTORY
####################################################
script_location = os.getcwd()
data_dir =  config["sample_directory"] if config['sample_directory'] else os.getcwd()
REFERENCES = config["reference"]
#############

# Preparing conda environements.
###############################
STAT_ENV=script_location+"/envs/stat.yaml"
#############


rule all:
    # input: "{myfile}.mosdepth.done",
    # input: "{myfile}_coverage.bed",
    # input: "{myfile}_zero_covered.bed",
    # input: "{myfile}_average_zero_cover.bed",
    # input: "{myfile}_average_zero_min_max_cover.bed", "{myfile}_SNVs_benchamrk.done"
    # input: "{myfile}_average_zero_min_max_cover.bed", "{myfile}_SNVs_benchamrk_values.tsv"
    # input: "{myfile}_coverage_and_benchamrk.tsv"
    input: "{myfile}_coverage_and_benchamrk.tsv", "{myfile}_average_zero_min_max_cover_exon.bed"
    output: touch("{myfile}.stat")

rule Mosdepth:
    input: "{myfile}.bam"
    output: "{myfile}.mosdepth.done"
    params:
        bed_file=config['bed_file'],
        coverage_prefix="{myfile}"
    threads: config['mosdepth_threads']
    log: "{myfile}.mosdepth.log"
    conda: STAT_ENV
    shell:"""
        mosdepth -t {threads}  --by {params.bed_file} {params.coverage_prefix} {input} > {log} 2>&1 && touch "{output}"
    """

rule MergeCoverageWithBam:
    input: "{myfile}.mosdepth.done"
    output: "{myfile}_coverage.bed"
    message: "merge average coverage per gene with gene statistics"
    params:
        bed_file=config['bed_region']
    conda: STAT_ENV
    shell:"""
        tail -n +2 {params.bed_file} | bedtools intersect -wo -a - -b  {wildcards.myfile}.regions.bed.gz -f 1 -r | cut -f1-12,17 > {output}
    """

rule ZeroCoverage:
    input: "{myfile}_coverage.bed"
    output: "{myfile}_zero_covered.bed"
    message: "Calculating zero covered percentage"
    log: "{myfile}_zero_covered.log"
    conda: STAT_ENV
    params:
        covmin = config['covmin'],
        cov0 = config['cov0'],
        covmin_value = config['covmin_value'],
    shell:"""
        bedtools intersect -wo -a {wildcards.myfile}.per-base.bed.gz -b {input} | awk 'BEGIN{{OFS="\\t";}}$4==0 {{print $5,$6,$7,$9, $(NF)}}' | datamash -s -g 1,2,3,4 sum 5 1> {params.cov0}.bed 2>>{log} &&\
        bedtools intersect -wo -a {wildcards.myfile}.per-base.bed.gz -b {input} | awk -v value={params.covmin_value} 'BEGIN{{OFS="\\t";}}$4<value {{print $5,$6,$7,$9, $(NF)}}' | datamash -s -g 1,2,3,4 sum 5 1> {params.covmin}.bed 2>>{log} &&\
        bedtools intersect -wao -a {params.covmin}.bed -b {params.cov0}.bed -r  -f 1 | bedtools sort | cut -f 1-5,10 | sed 's/\s-1/\t0/' 1> {output} 2>> {log} &&\
        rm {params.covmin}.bed {params.cov0}.bed
    """

rule MergeZeroCoverageWithStat:
    """
    Count both bases with zero coverage and bases less than X coverage, X value from params
    """
    input:
        average_cover = "{myfile}_coverage.bed",
        zero_cover = "{myfile}_zero_covered.bed",
    output: "{myfile}_average_zero_cover.bed"
    message: "Running {rule} to merge coverage with statistics genes file"
    log: "{myfile}_average_zero_cover.log"
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wao -a {input.average_cover} -b {input.zero_cover} -f 1 -r |  cut -f1-13,18-19 | sed 's/\t-1/\t0/;s/\t\./\t0/g' 1> {output} 2> {log}
    """

rule MinAndMaxCover:
    input: "{myfile}_average_zero_cover.bed"
    output: "{myfile}_average_zero_min_max_cover.bed"
    message: "Calculating all coverages for genes"
    log: "{myfile}_average_zero_min_max_cover.log"
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wo -a {wildcards.myfile}.per-base.bed.gz -b {input} |  datamash -s -g 5,6,7,9 min 4 max 4 | bedtools intersect -wo -a {input} -b - -f 1 -r |  cut -f1-15,20,21 1>{output} 2>{log}
    """

rule CollectSNVsBenchmark:
    input: lambda wildcards: expand(["{{SNV_bench_dir}}/{gene}__{{myfile}}".format(gene=os.path.basename(i).split('.')[0]) for i in glob.glob(config['genes_dir']+"/*.bed")], SNV_bench_dir=config['output_dir_name'], myfile=wildcards.myfile)
    output: touch("{myfile}_SNVs_benchamrk.done")

rule SNVsBenchmark:
    input:
        bed_file = "{genes_dir}/{{gene}}.bed".format(genes_dir=config['genes_dir']),
        snv = "{myfile}_SNVs/phased_merge_output.vcf.gz"
    output: directory("{SNV_bench_dir}/{{gene}}__{{myfile}}".format(SNV_bench_dir=config['output_dir_name']))
    threads: config['SNVs_bench_threads']
    params:
        GOLD_SNP=config['gold_vcf'],
        ref=config['REF_38'],
    conda: STAT_ENV
    shell:"""
      rtg RTG_MEM=16G vcfeval --baseline {params.GOLD_SNP} --bed-regions {input.bed_file} -c {input.snv} -o {output} -t {params.ref} -T {threads}
    """

rule WriteSNVsBenchmark:
    input:
        Bench_done = "{myfile}_SNVs_benchamrk.done",
        SNVs_benchamrk_dir_list = glob.glob(config['output_dir_name']+"/*")
    output: "{myfile}_SNVs_benchamrk_values.tsv"
    message: "Writing SNVs Benchamrk for {rule}"
    conda: STAT_ENV
    params:
        genes=config['genes_dir']
    shell:"""
        for i in {input.SNVs_benchamrk_dir_list};
            do
            filename=$(basename $i)
            gene_name=$(echo $filename | sed 's/__{wildcards.myfile}//g')
            variants=$(bcftools view -H {wildcards.myfile}_SNVs/phased_merge_output.vcf.gz -R {params.genes}/${{gene_name}}.bed | wc -l)
            snps=$(bcftools view -H -v snps {wildcards.myfile}_SNVs/phased_merge_output.vcf.gz -R {params.genes}/${{gene_name}}.bed | wc -l)
            indels=$(bcftools view -H -v indels {wildcards.myfile}_SNVs/phased_merge_output.vcf.gz -R {params.genes}/${{gene_name}}.bed | wc -l)
            Bench_values=$(tail -n +3  ${{i}}/summary.txt | head -n 1 | awk 'BEGIN{{OFS="\t"}}{{print $(NF-2), $(NF-1), $(NF)}}')
            echo "${{gene_name}}\t${{variants}}\t${{snps}}\t${{indels}}\t${{Bench_values}}";
            done  >> {output}
    """

rule UpdateSNVsStatWithBenchmark:
    input:
        benchmark="{myfile}_SNVs_benchamrk_values.tsv",
        stat="{myfile}_average_zero_min_max_cover.bed"
    output: "{myfile}_coverage_and_benchamrk.tsv"
    message: "merging benchamrk with coverage"
    conda: STAT_ENV
    params:
        maxcov = config['covmin_value']
    shell:"""
        python {data_dir}/merge_coverage_bench.py {input.benchmark} {input.stat} {params.maxcov} > {output}
    """

rule ExoneAverageCoverage:
    input: "{myfile}.bam"
    output: "{myfile}.mosdepth.exon.done"
    params:
        bed_file=config['exon_file'],
        coverage_prefix="{myfile}.exon"
    threads: config['mosdepth_threads']
    log: "{myfile}.mosdepth.exon.log"
    message: "Running {rule} to calculate exon coverage"
    conda: STAT_ENV
    shell:"""
        mosdepth -t {threads}  --by {params.bed_file} {params.coverage_prefix} {input} > {log} 2>&1 && touch "{output}"
    """

rule MergeExonCoverageWithBam:
    input: "{myfile}.mosdepth.exon.done"
    output: "{myfile}_coverage_exon.bed"
    message: "Running {rule} merge average coverage per exon with exon file"
    params:
        bed_file=config['exon_file']
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wo -a {params.bed_file} -b  {wildcards.myfile}.exon.regions.bed.gz -f 1 -r | cut -f 1-7,12 | uniq > {output}
    """

rule ZeroExonCoverage:
    input: "{myfile}_coverage_exon.bed"
    output: "{myfile}_zero_covered_exon.bed"
    message: "Running {rule} calculating zero covered percentage for exon"
    log: "{myfile}_zero_covered_exon.log"
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wo -a {wildcards.myfile}.exon.per-base.bed.gz -b {input} | awk 'BEGIN{{OFS="\\t";}}$4==0 {{print $5,$6,$7,$11, $(NF)}}' | datamash -g 1,2,3,4 sum 5 1> {output} 2> {log}
    """

rule MergeZeroExonCoverageWithStat:
    input:
        average_cover = "{myfile}_coverage_exon.bed",
        zero_cover = "{myfile}_zero_covered_exon.bed",
    output: "{myfile}_average_zero_cover_exon.bed"
    message: "Running {rule} to merge exon coverage with statistics genes file"
    log: "{myfile}_average_zero_cover_exon.log"
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wao -a {input.average_cover} -b {input.zero_cover} -f 1 -r |  cut -f1-8,13 | sed 's/\s-1/\t0/' 1> {output} 2> {log}
    """

rule MinAndMaxExonCover:
    input: "{myfile}_average_zero_cover_exon.bed"
    output: "{myfile}_average_zero_min_max_cover_exon.bed"
    message: "Running {rule} calculating all coverages for exon"
    log: "{myfile}_average_zero_min_max_cover_exon.log"
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wo -a {wildcards.myfile}.exon.per-base.bed.gz -b {input} |  datamash -s -g  5,6,7,11 min 4 max 4 | bedtools intersect -wo -a {input} -b - -f 1 -r |  cut -f1-9,14,15 | sed '1i Chr\\tStart\\tEnd\\tgene_id\\ttranscript_id\\tgene_name\\texon_id\\tAverage_cov\\tZero_cover\\tMin_bp_cover\\tMax_bp_cover' 1>{output} 2>{log}
    """

## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##

onsuccess:
	shell("mkdir -p snake_log && find . -maxdepth 1  \( -name 'snakejob*' -or -name 'slurm*' \) -type f -exec mv -t snake_log {{}} \;")

onerror:
	shell("mkdir -p snake_log && find . -maxdepth 1  \( -name 'snakejob*' -or -name 'slurm*' \) -type f -exec mv -t snake_log {{}}  \;")
