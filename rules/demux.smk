checkpoint demux_ubam:
    input:
        reads=movie,
        barcodes=config["barcodes"],
        biosamples=config["biosamples"],
    output:
        odir=directory( f'batches/{batch}/demux' ),
        report=f'batches/{batch}/demux/demultiplex.lima.report',
    log:
        f"batches/{batch}/logs/lima/demux.log",
    params:
        odir=f'batches/{batch}/demux/',
        filters='--hifi-preset ASYMMETRIC --min-score 80 --min-qv 20',
        #filters='--min-score 40 --min-end-score 40 --min-ref-span 0.75 --same --ignore-missing-adapters',
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
        (mkdir -p {output.odir}
        lima {params.filters} \
             --split-named \
             --dump-removed \
             -j {threads} \
             --log-level {params.loglevel} \
             --biosample-csv {input.biosamples} \
             {input.reads} {input.barcodes} {params.odir}/demultiplex.bam) > {log} 2>&1
        '''        

rule check_demux_fail:
    input:
        f'batches/{batch}/demux/',
    params:
        samples=_get_demuxed_samples,
    log:
        f"batches/{batch}/logs/lima/demux_fail.log",
    run:
        missing = sample2barcode.keys() - params.samples
        if len( missing ):
            with open( f"batches/{batch}/demux_FAIL.txt", 'w' ) as ofile:
                for sample in sorted( missing ):
                    ofile.write( f'{sample2barcode[sample]},{sample}\n' ) 


targets.append( f'batches/{batch}/logs/lima/demux_fail.log' )
