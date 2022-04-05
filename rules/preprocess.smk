rule markdup_ubam:
    input:
        lambda wildcards: f'batches/{batch}/demux/demultiplex.{sample2barcode[wildcards.sample]}.bam',
    output:
        markdup=f'batches/{batch}/{{sample}}/downsampled_all/markdup/markdups.bam',
        idx=f'batches/{batch}/{{sample}}/downsampled_all/markdup/markdups.bam.bai',
    log:
        f"batches/{batch}/logs/pbmarkdup/{{sample}}.log",
    params:
        options='--cross-library',
        loglevel='INFO',
    threads:
        16
    benchmark:
        f"batches/{batch}/benchmarks/pbmarkdup/{{sample}}.bam.tsv"
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

rule pbindex_bam:
    input:
        f'batches/{batch}/{{sample}}/downsampled_all/markdup/markdups.bam',
    output:
        f'batches/{batch}/{{sample}}/downsampled_all/markdup/markdups.bam.pbi',
    conda:
        "envs/pbtools.yaml"
    shell:
        '''
        pbindex {input}
        '''

rule markdup_fastq:
    input:
        lambda wildcards: f'batches/{batch}/demux/demultiplex.{sample2barcode[wildcards.sample]}.fastq',
    output:
        markdup=f'batches/{batch}/{{sample}}/downsampled_all/markdup/markdups.fastq',
        idx=f'batches/{batch}/{{sample}}/downsampled_all/markdup/markdups.fastq.fai',
    log:
        f"batches/{batch}/logs/pbmarkdup/{{sample}}.log",
    params:
        options='--cross-library',
        loglevel='INFO',
    threads:
        16
    benchmark:
        f"batches/{batch}/benchmarks/pbmarkdup/{{sample}}.fastq.tsv"
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
        bam=f'batches/{batch}/{{sample}}/downsampled_all/markdup/markdups.bam',
        pbi=f'batches/{batch}/{{sample}}/downsampled_all/markdup/markdups.bam.pbi',
    output:
        f"batches/{batch}/{{sample}}/downsampled_{{maxreads}}/markdup/markdups.bam",
    log:
        f"batches/{batch}/logs/pbcoretools/bamsieve/downsample.{{maxreads}}.{{sample}}.log",
    benchmark:
        f"batches/{batch}/benchmarks/pbtools/downsample.{{maxreads}}.{{sample}}.bam.tsv"
    params:
        maxreads='{maxreads}',
        seed=42,
    threads: 1
    conda:
        "envs/pbtools.yaml"
    message:
        "Downsampling {input.bam} to at most {params.maxreads}."
    shell:
        '''
        bamsieve -s {params.seed} -n {params.maxreads} {input.bam} {output} > {log} 2>&1
        '''

rule downsample_fastq:
    input:
        f'batches/{batch}/{{sample}}/downsampled_all/markdup/markdups.fastq',
    output:
        temp(f"batches/{batch}/{{sample}}/downsampled_{{maxreads}}/markdup/markdups.fastq"),
    log:
        f"batches/{batch}/logs/pbcoretools/bamsieve/downsample.{{maxreads}}.{{sample}}.log",
    benchmark:
        f"batches/{batch}/benchmarks/seqtk/downsample.{{maxreads}}.{{sample}}.fastq.tsv"
    params:
        maxreads='{maxreads}',
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
