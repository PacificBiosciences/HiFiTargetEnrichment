rule pbsv_discover:
    input:
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        bai=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
        tr_bed = config['ref']['tr_bed'],
    output: 
        temp( f"batches/{batch}/{{sample}}/pbsv/svsig/{{sample}}.{ref}.svsig.gz" ),
    log: 
        f"batches/{batch}/logs/pbsv/discover/{{sample}}.{ref}.log"
    benchmark: 
        f"batches/{batch}/benchmarks/pbsv/discover/{{sample}}.{ref}.tsv"
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
        svsigs=f"batches/{batch}/{{sample}}/pbsv/svsig/{{sample}}.{ref}.svsig.gz",
        reference = config['ref']['fasta'],
    output: 
        f"batches/{batch}/{{sample}}/pbsv/{{sample}}.{ref}.pbsv.vcf",
    log: 
        f"batches/{batch}/logs/pbsv/call/{{sample}}.{ref}.log",
    benchmark: 
        f"batches/{batch}/benchmarks/pbsv/call/{{sample}}.{ref}.tsv",
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
