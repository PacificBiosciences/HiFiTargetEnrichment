rule dedup_ubam:
    input:
        lambda wildcards: expand( f'batches/{batch}/demux/{{movie}}/demultiplex.{sample2barcode[wildcards.sample]}.bam',
                                   movie=movies ),
    output:
        dedup=f'batches/{batch}/{{sample}}/dedup/deduplicated.{{barcode}}.bam',
        dups=f'batches/{batch}/{{sample}}/dedup/duplicates.{{barcode}}.bam',
        idx=f'batches/{batch}/{{sample}}/dedup/deduplicated.{{barcode}}.bam.bai',
    log:
        f"batches/{batch}/logs/pbmarkdup/dedup/{{sample}}.{{barcode}}.log",
    params:
        options='--cross-library',
        loglevel='INFO',
    threads:
        16
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.{{barcode}}.dedup.bam.tsv"
    conda:
        'envs/pbmarkdup.yaml'
    message:
        "Deduplicating {input} using pbmarkdup with option {params.options}."
    shell:
        '''
        (pbmarkdup --log-level {params.loglevel} \
                  -j {threads} \
                  {params.options} \
                  --dup-file {output.dups} \
                  {input} {output.dedup}) > {log} 2>&1
        samtools index {output.dedup}
        '''

rule dedup_fastq:
    input:
        lambda wildcards: expand( f'batches/{batch}/demux/{{movie}}/demultiplex.{sample2barcode[wildcards.sample]}.fastq',
                                   movie=movies ),
    output:
        dedup=f'batches/{batch}/{{sample}}/dedup/deduplicated.{{barcode}}.fastq',
        dups=f'batches/{batch}/{{sample}}/dedup/duplicates.{{barcode}}.fastq',
        idx=f'batches/{batch}/{{sample}}/dedup/deduplicated.{{barcode}}.fastq.fai',
    log:
        f"batches/{batch}/logs/pbmarkdup/dedup/{{sample}}.{{barcode}}.log",
    params:
        loglevel='INFO',
    threads:
        16
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.{{barcode}}.dedup.fastq..tsv"
    conda:
        'envs/pbmarkdup.yaml'
    message:
        "Deduplicating {input} using pbmarkdup with option {params.options}."
    shell:
        '''
        (pbmarkdup --log-level {params.loglevel} \
                  -j {threads} \
                  --cross-library \
                  --dup-file {output.dups} \
                  {input} {output.dedup}) > {log} 2>&1
        samtools fqidx {output.dedup}
        '''
