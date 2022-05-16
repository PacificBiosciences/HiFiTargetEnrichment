rule demux_ubam:
    input:
        reads=movie,
        barcodes=config["barcodes"],
        biosamples=config["biosamples"],
    output:
        temp( expand( f'batches/{batch}/demux/demultiplex.{{barcode}}.bam', barcode=barcode2sample.keys() ) ),
        f'batches/{batch}/demux/demultiplex.lima.report',
    log:
        f"batches/{batch}/logs/lima/demux.log",
    params:
        odir=f'batches/{batch}/demux/',
        #preset='hifi-symmetric',
        #filters='--ccs --min-score 40 --min-end-score 40 --min-ref-span 0.75 --same',
        filters='--min-score 40 --min-end-score 40 --min-ref-span 0.75 --same --ignore-missing-adapters',
        loglevel='INFO',
    threads:
        24
    benchmark:
        f"batches/{batch}/benchmarks/lima/demux.tsv",
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
        reads=movie,
        barcodes=config["barcodes"],
        biosamples=config["biosamples"],
    output:
        temp( expand( f'batches/{batch}/demux/demultiplex.{{barcode}}.fastq', barcode=barcode2sample.keys() ) ),
        f'batches/{batch}/demux/demultiplex.lima.report',
    log:
        f"batches/{batch}/logs/lima/demux.log",
    params:
        odir=f'batches/{batch}/demux/',
        #preset='hifi-symmetric',
        filters='--min-score 40 --min-end-score 40 --min-ref-span 0.75 --same --ignore-missing-adapters',
        loglevel='INFO',
    threads:
        24
    benchmark:
        f"batches/{batch}/benchmarks/lima/demux.tsv",
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
