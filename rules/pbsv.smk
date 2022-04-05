localrules: bcftools_concat_pbsv_vcf


rule pbsv_discover:
    input:
        bam=f"batches/{batch}/{{sample}}/downsampled_{{maxreads}}/aligned/{{sample}}.{ref}.bam",
        bai=f"batches/{batch}/{{sample}}/downsampled_{{maxreads}}/aligned/{{sample}}.{ref}.bam.bai",
        tr_bed = config['ref']['tr_bed'],
    output: 
        f"batches/{batch}/{{sample}}/downsampled_{{maxreads}}/pbsv/svsig/{movie}.{ref}.svsig.gz",
    log: 
        f"batches/{batch}/logs/pbsv/discover/{{sample}}.{{maxreads}}.{ref}.log"
    benchmark: 
        f"batches/{batch}/benchmarks/pbsv/discover/{{sample}}.{{maxreads}}.{ref}.tsv"
    params:
        extra = "--hifi",
        loglevel = "INFO",
        #region = lambda wildcards: wildcards.region,
    conda: 
        "envs/pbsv.yaml"
    message: 
        "Executing {rule}: Discovering structural variant signatures in {wildcards.sample} from {input.bam}."
    shell:
        """
        (pbsv discover {params.extra} \
            --log-level {params.loglevel} \
            --tandem-repeats {input.tr_bed} \
            {input.bam} {output}) > {log} 2>&1
        """


rule pbsv_call:
    input:
        svsigs=f"batches/{batch}/{{sample}}/downsampled_{{maxreads}}/pbsv/svsig/{movie}.{ref}.svsig.gz",
        reference = config['ref']['fasta'],
    output: 
        f"batches/{batch}/{{sample}}/downsampled_{{maxreads}}/pbsv/{{sample}}.{ref}.pbsv.vcf",
    log: 
        f"batches/{batch}/logs/pbsv/call/{{sample}}.{{maxreads}}.{ref}.log",
    benchmark: 
        f"batches/{batch}/benchmarks/pbsv/call/{{sample}}.{{maxreads}}.{ref}.tsv",
    params:
        extra = "--hifi -m 20",
        loglevel = "INFO"
    threads: 8
    conda: 
        "envs/pbsv.yaml"
    message: 
        "Executing {rule}: Calling structural variants in {wildcards.sample} from SVSIG: {input.svsigs}"
    shell:
        """
        (pbsv call {params.extra} \
            --log-level {params.loglevel} \
            --num-threads {threads} \
            {input.reference} {input.svsigs} {output}) > {log} 2>&1
        """
