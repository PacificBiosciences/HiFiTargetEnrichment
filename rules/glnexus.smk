def _get_regions(wildcards):
    regions = {}
    with checkpoints.make_regions.get().output[0].open() as bed:
        for line in bed.read().strip().split('\n'):
            chrm,start,stop,name = line.strip().split()[:4]
            regions[name] = f'{chrm}:{start}-{stop}'
        return regions

def _get_region(wildcards):
    return _get_regions(wildcards)[wildcards.region]
    
checkpoint make_regions:
    input:
        config['targets'],
    output:
        f'batches/{batch}/glnexus/regions.bed',
    params:
        distance=config['merge'],
    log:
        f"batches/{batch}/logs/glnexus/make_regions.log"
    conda:
        'envs/samtools.yaml'
    message:
        'Generating process regions for {input}'
    shell:
        '''
        (bedtools merge -i <(sort -k1,1 -k2,2n {input}) -d {params.distance} -c 4 -o collapse -delim - > {output}) 2> {log}
        '''


rule glnexus:
    input:
        gvcf=expand(f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.g.vcf.gz",sample=sample2barcode.keys()),
        tbi=expand(f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.g.vcf.gz.tbi",sample=sample2barcode.keys()),
        bed=f'batches/{batch}/glnexus/regions.bed',
    output:
        bcf=temp(f"batches/{batch}/glnexus/{batch}.{ref}.deepvariant.glnexus.bcf"),
        scratch_dir=temp(directory(f"batches/{batch}/glnexus/{batch}.{ref}.GLnexus.DB/"))
    log: 
        f"batches/{batch}/logs/glnexus/{batch}.{ref}.log"
    benchmark: 
        f"batches/{batch}/benchmarks/glnexus/{batch}.{ref}.tsv"
    container: 
        f"docker://ghcr.io/dnanexus-rnd/glnexus:{config['GLNEXUS_VERSION']}"
    threads: 24
    message: 
        f"Executing {{rule}}: Joint calling variants from {batch}"
    shell:
        """
        (rm -rf {output.scratch_dir} && \
        glnexus_cli --threads {threads} \
                    --dir {output.scratch_dir} \
                    --config DeepVariant_unfiltered \
                    --bed {input.bed} \
                    {input.gvcf} > {output.bcf}) 2> {log}
        """

rule bcftools_bcf2vcf:
    input: 
        f"batches/{batch}/glnexus/{{prefix}}.bcf"
    output: 
        f"batches/{batch}/glnexus/{{prefix}}.vcf.gz"
    log: 
        f"batches/{batch}/logs/bcftools/glnexus/{{prefix}}.log"
    benchmark: 
        f"batches/{batch}/benchmarks/bcftools/glnexus/{{prefix}}.tsv"
    params: "--threads 4 -Oz"
    threads: 4
    conda: 
        "envs/bcftools.yaml"
    message: 
        "Executing {rule}: Converting GLnexus BCF to VCF for {input}."
    shell: 
        '''
        (bcftools view {params} {input} -o {output}) > {log} 2>&1
        '''


rule split_glnexus_vcf:
    input:
        vcf=f"batches/{batch}/glnexus/{batch}.{ref}.deepvariant.glnexus.vcf.gz",
        tbi=f"batches/{batch}/glnexus/{batch}.{ref}.deepvariant.glnexus.vcf.gz.tbi",
    output: 
        f"batches/{batch}/whatshap_cohort/regions/{batch}.{ref}.{{region}}.deepvariant.glnexus.vcf.gz"
    log: 
        f"batches/{batch}/logs/tabix/query/{batch}.{ref}.{{region}}.glnexus.vcf.log"
    benchmark: 
        f"batches/{batch}/benchmarks/tabix/query/{batch}.{ref}.{{region}}.glnexus.vcf.tsv"
    params: 
        vcf=f"batches/{batch}/whatshap_cohort/regions/{batch}.{ref}.{{region}}.deepvariant.glnexus.vcf",
        region=_get_region,
        extra='-h'
    conda: 
        "envs/htslib.yaml"
    message: 
        "Executing {rule}: Extracting {wildcards.region} variants from {input}."
    shell: 
        '''
        (tabix {params.extra} {input.vcf} {params.region} > {params.vcf} && bgzip {params.vcf}) 2> {log}
        '''


rule whatshap_phase_cohort:
    input:
        reference = config['ref']['fasta'],
        vcf = f"batches/{batch}/whatshap_cohort/regions/{batch}.{ref}.{{region}}.deepvariant.glnexus.vcf.gz",
        tbi = f"batches/{batch}/whatshap_cohort/regions/{batch}.{ref}.{{region}}.deepvariant.glnexus.vcf.gz.tbi",
        phaseinput=expand(f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",sample=sample2barcode.keys()),
        phasebai=expand(f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",sample=sample2barcode.keys())
    output: 
        temp(f"batches/{batch}/whatshap_cohort/regions/{batch}.{ref}.{{region}}.deepvariant.glnexus.phased.vcf.gz")
    log: 
        f"batches/{batch}/logs/whatshap/phase/{batch}.{ref}.{{region}}.log"
    benchmark: 
        f"batches/{batch}/benchmarks/whatshap_cohort/phase/{batch}.{ref}.{{region}}.tsv"
    params:
        extra="--indels"
    conda: 
        "envs/whatshap.yaml"
    message: 
        "Executing {rule}: Phasing {input.vcf} using {input.phaseinput} for targets {wildcards.region}."
    shell:
        """
        (whatshap phase {params.extra} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.phaseinput}) > {log} 2>&1
        """

## to include pedigree information in phasing
# extra = """--indels --ped cohorts/{cohort}/{cohort}.ped \
#            --no-genetic-haplotyping --recombination-list cohorts/{cohort}/whatshap/{cohort}.recombination.list"""

def _get_vcf(index=False):
    def f(wildcards):
        ext = '.tbi' if index else '' 
        return expand(f"batches/{batch}/whatshap_cohort/regions/{batch}.{ref}.{{region}}.deepvariant.glnexus.phased.vcf.gz{ext}",
                      region=_get_regions(wildcards).keys())
    return f

rule whatshap_bcftools_concat:
    input:
        calls=_get_vcf(),
        indices=_get_vcf(index=True)
    output: 
        f"batches/{batch}/whatshap_cohort/{batch}.{ref}.deepvariant.glnexus.phased.vcf.gz"
    log: 
        f"batches/{batch}/logs/bcftools/concat/{batch}.{ref}.whatshap.log"
    benchmark: 
        f"batches/{batch}/benchmarks/bcftools/concat/{batch}.{ref}.whatshap.tsv"
    params: 
        "-a -Oz"
    conda: 
        "envs/bcftools.yaml"
    message: 
        "Executing {rule}: Concatenating WhatsHap phased VCFs: {input.calls}"
    shell: 
        '''
        (bcftools concat {params} -o {output} {input.calls}) > {log} 2>&1
        '''


rule whatshap_stats_cohort:
    input:
        vcf = f"batches/{batch}/whatshap_cohort/{batch}.{ref}.deepvariant.glnexus.phased.vcf.gz",
        tbi = f"batches/{batch}/whatshap_cohort/{batch}.{ref}.deepvariant.glnexus.phased.vcf.gz.tbi",
        chr_lengths = config['ref']['chr_lengths']
    output:
        gtf = f"batches/{batch}/whatshap_cohort/{batch}.{ref}.deepvariant.glnexus.phased.gtf",
        tsv = f"batches/{batch}/whatshap_cohort/{batch}.{ref}.deepvariant.glnexus.phased.tsv",
        blocklist = f"batches/{batch}/whatshap_cohort/{batch}.{ref}.deepvariant.glnexus.phased.blocklist"
    log: 
        f"batches/{batch}/logs/whatshap/cohort_stats/{batch}.{ref}.log"
    benchmark: 
        f"batches/{batch}/benchmarks/whatshap/cohort_stats/{batch}.{ref}.tsv"
    conda: 
        "envs/whatshap.yaml"
    message: 
        "Executing {rule}: Calculating phasing stats for {input.vcf}."
    shell:
        """
        (whatshap stats \
            --gtf {output.gtf} \
            --tsv {output.tsv} \
            --block-list {output.blocklist} \
            --chr-lengths {input.chr_lengths} \
            {input.vcf}) > {log} 2>&1
        """
