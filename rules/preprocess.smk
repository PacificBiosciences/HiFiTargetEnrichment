rule markdup_ubam:
    input:
        lambda wildcards: f'batches/{batch}/demux/demultiplex.{sample2barcode[wildcards.sample]}.bam',
    output:
        markdup=temp( f'batches/{batch}/{{sample}}/markdup/markdups.bam' ),
        idx=temp( f'batches/{batch}/{{sample}}/markdup/markdups.bam.bai' ),
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
        f'batches/{batch}/{{sample}}/markdup/markdups.bam',
    output:
        temp( f'batches/{batch}/{{sample}}/markdup/markdups.bam.pbi' ),
    conda:
        "envs/pbtools.yaml"
    shell:
        '''
        pbindex {input}
        '''
