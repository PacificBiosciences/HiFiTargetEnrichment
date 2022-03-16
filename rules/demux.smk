rule demux_ubam:
    input:
        reads=lambda wildcards: ubam_dict[wildcards.movie],
        barcodes=config["barcodes"],
        biosamples=config["biosamples"],
    output:
        expand( f'batches/{batch}/demux/' + '{{movie}}/demultiplex.{barcode}.bam', barcode=barcode2sample.keys() ),
        f'batches/{batch}/demux/{{movie}}/demultiplex.lima.report',
    log:
        f"batches/{batch}/logs/lima/demux/{{movie}}.log",
    params:
        odir=f'batches/{batch}/demux/{{movie}}',
        #preset='hifi-symmetric',
        #filters='--ccs --min-score 40 --min-end-score 40 --min-ref-span 0.75 --same',
        filters='--min-score 40 --min-end-score 40 --min-ref-span 0.75 --same --ignore-missing-adapters',
        loglevel='INFO',
    threads:
        24
    benchmark:
        f"batches/{batch}/benchmarks/{{movie}}.demux.bam.tsv",
    conda:
        'envs/lima.yaml'
    message:
        "Demultiplexing {input.reads} with {input.barcodes} using lima with options {params.filters}."
    shell:
        '''
        (lima {params.filters} \
             --split-named \
             --dump-removed \
             -j {threads} \
             --log-level {params.loglevel} \
             --biosample-csv {input.biosamples} \
             {input.reads} {input.barcodes} {params.odir}/demultiplex.bam) > {log} 2>&1
        '''        

rule demux_fastq:
    input:
        reads=lambda wildcards: fastq_dict[wildcards.movie],
        barcodes=config["barcodes"],
        biosamples=config["biosamples"],
    output:
        expand( f'batches/{batch}/demux' + '{{movie}}/demultiplex.{barcode}.fastq', barcode=barcode2sample.keys() ),
        f'batches/{batch}/demux/{{movie}}/demultiplex.lima.report',
    log:
        f"batches/{batch}/logs/lima/demux/{{movie}}.log",
    params:
        odir=f'batches/{batch}/demux/{{movie}}',
        #preset='hifi-symmetric',
        filters='--min-score 40 --min-end-score 40 --min-ref-span 0.75 --same --ignore-missing-adapters',
        loglevel='INFO',
    threads:
        24
    benchmark:
        f"batches/{batch}/benchmarks/{{movie}}.demux.fastq.tsv",
    conda:
        'envs/lima.yaml'
    message:
        "Demultiplexing {input.reads} with {input.barcodes} using lima with options {params.filters}."
    shell:
        '''
        (lima {params.filters} \
             --split-named \
             --dump-removed \
             -j {threads} \
             --log-level {params.loglevel} \
             --biosample-csv {input.biosamples} \
             {input.reads} {input.barcodes} {params.odir}/demultiplex.fastq) > {log} 2>&1
        '''        
