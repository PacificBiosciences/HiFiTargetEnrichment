rule pbmm2_align_ubam:
    input:
        reference=config["ref"]["fasta"],
        ref_index=config["ref"]["index"],
        query=lambda wildcards: f'batches/{batch}/{{sample}}/dedup/deduplicated.{sample2barcode[wildcards.sample]}.bam',
    output:
        bam=(f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam"),
        bai=(f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai"),
    log:
        f"batches/{batch}/logs/pbmm2/align/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.pbmm2.bam.tsv"
    params:
        sample="{sample}",
        preset="CCS",
        extra="--sort --unmapped -c 0 -y 70",
        loglevel="INFO",
    threads: 24
    conda:
        "envs/pbmm2.yaml"
    message:
        "Aligning {input.query} to {input.reference} using pbmm2 with --preset {params.preset} {params.extra}."
    shell:
        """
        (pbmm2 align --num-threads {threads} \
            --preset {params.preset} \
            --sample {params.sample} \
            --log-level {params.loglevel} \
            {params.extra} \
            {input.reference} \
            {input.query} \
            {output.bam}) > {log} 2>&1
        """

rule pbmm2_align_fastq:
    input:
        reference=config["ref"]["fasta"],
        ref_index=config["ref"]["index"],
        query=lambda wildcards: f'batches/{batch}/{{sample}}/dedup/deduplicated.{sample2barcode[wildcards.sample]}.fastq',
    output:
        bam=(f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam"),
        bai=(f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai"),
    log:
        f"batches/{batch}/logs/pbmm2/align/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.pbmm2.fastq.tsv"
    params:
        sample="{sample}",
        preset="CCS",
        extra="--sort --unmapped -c 0 -y 70",
        loglevel="INFO",
    threads: 24
    conda:
        "envs/pbmm2.yaml"
    message:
        "Aligning {input.query} to {input.reference} using pbmm2 with --preset {params.preset} {params.extra}."
    shell:
        """
        (pbmm2 align --num-threads {threads} \
            --preset {params.preset} \
            --sample {params.sample} \
            --log-level {params.loglevel} \
            {params.extra} \
            {input.reference} \
            {input.query} \
            {output.bam}) > {log} 2>&1
        """
