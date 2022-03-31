rule markdup_ubam:
    input:
        lambda wildcards: f'batches/{batch}/demux/demultiplex.{sample2barcode[wildcards.sample]}.bam',
    output:
        markdup=f'batches/{batch}/{{sample}}/markdup/markdups.{{barcode}}.bam',
        idx=f'batches/{batch}/{{sample}}/markdup/markdups.{{barcode}}.bam.bai',
    log:
        f"batches/{batch}/logs/pbmarkdup/{{sample}}.{{barcode}}.log",
    params:
        options='--cross-library',
        loglevel='INFO',
    threads:
        16
    benchmark:
        f"batches/{batch}/benchmarks/pbmarkdup/{{sample}}.{{barcode}}.bam.tsv"
    conda:
        'envs/pbmarkdup.yaml'
    message:
        "Marking duplicates in {input} using pbmarkdup with option {params.options}."
    shell:
        '''
        (pbmarkdup --log-level {params.loglevel} \
                   -j {threads} \
                   {params.options} \
                   {input} {output.markdup}
        samtools index {output.markdup}) > {log} 2>&1
        '''

rule markdup_fastq:
    input:
        lambda wildcards: f'batches/{batch}/demux/demultiplex.{sample2barcode[wildcards.sample]}.fastq',
    output:
        markdup=f'batches/{batch}/{{sample}}/markdup/markdups.{{barcode}}.fastq',
        idx=f'batches/{batch}/{{sample}}/markdup/markdups.{{barcode}}.fastq.fai',
    log:
        f"batches/{batch}/logs/pbmarkdup/{{sample}}.{{barcode}}.log",
    params:
        options='--cross-library',
        loglevel='INFO',
    threads:
        16
    benchmark:
        f"batches/{batch}/benchmarks/pbmarkdup/{{sample}}.{{barcode}}.fastq.tsv"
    conda:
        'envs/pbmarkdup.yaml'
    message:
        "Marking duplicates in {input} using pbmarkdup with option {params.options}."
    shell:
        '''
        (pbmarkdup --log-level {params.loglevel} \
                   -j {threads} \
                   {params.options} \
                   {input} {output.markdup}
        samtools fqidx {output.dedup}) > {log} 2>&1
        '''

rule downsample_bam:
    input:
        lambda wildcards: f'batches/{batch}/{{sample}}/markdup/markdups.{sample2barcode[wildcards.sample]}.bam',
    output:
        f"batches/{batch}/{{sample}}/downsampled/markdups.bam",
    log:
        f"batches/{batch}/logs/pbcoretools/bamsieve/downsample.{{sample}}.log",
    benchmark:
        f"batches/{batch}/benchmarks/pbtools/downsample.{{sample}}.bam.tsv"
    params:
        maxreads=int(config['downsample']),
        seed=42,
    threads: 1
    conda:
        "envs/pbtools.yaml"
    message:
        "Downsampling {input} to at most {params.maxreads}."
    shell:
        '''
        (
        pbindex {input}
        bamsieve -s {params.seed} -n {params.maxreads} {input} {output}
        ) > {log} 2>&1
        '''

rule downsample_fastq:
    input:
        lambda wildcards: f'batches/{batch}/{{sample}}/markdup/markdups.{sample2barcode[wildcards.sample]}.fastq',
    output:
        temp(f"batches/{batch}/{{sample}}/downsampled/markdups.fastq"),
    log:
        f"batches/{batch}/logs/pbcoretools/bamsieve/downsample.{{sample}}.log",
    benchmark:
        f"batches/{batch}/benchmarks/seqtk/downsample.{{sample}}.fastq.tsv"
    params:
        maxreads=int(config['downsample']),
        seed=42,
    threads: 1
    conda:
        "envs/seqtk.yaml"
    message:
        "Downsampling {input} to at most {params.maxreads}."
    shell:
        '''
        ( seqtk sample -s {params.seed} {input} {params.maxreads} > {output} ) > {log} 2>&1
        '''
