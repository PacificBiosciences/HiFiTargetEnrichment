
rule create_dict:
    input:
        config['ref']['fasta'],
    output:
        temp(f"batches/{batch}/picard/reference.dict"),
    log:
        f"batches/{batch}/logs/picard/create_dict.log",
    params:
        extra="", 
    resources:
        mem_mb=1024,
    wrapper:
        "v1.3.1/bio/picard/createsequencedictionary"


rule bed_to_interval_list:
    input:
        bed=lambda wildcards: config[wildcards.bedfile],
        dict=f"batches/{batch}/picard/reference.dict",
    output:
        temp( f"batches/{batch}/picard/{{bedfile}}.interval_list" ),
    log:
        f"batches/{batch}/logs/picard/bedtointervallist/{{bedfile}}.log",
    params:
        extra="--SORT true",  # sort output interval list before writing
    resources:
        mem_mb=1024,
    wrapper:
        "v1.3.1/bio/picard/bedtointervallist"

rule picard_collect_hs_metrics:
    input:
        bam=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.haplotagged.bam",
        reference=config["ref"]["fasta"],
        # Baits and targets should be given as interval lists. These can
        # be generated from bed files using picard BedToIntervalList.
        bait_intervals=f"batches/{batch}/picard/probes.interval_list",
        target_intervals=f"batches/{batch}/picard/targets.interval_list",
    output:
        f"batches/{batch}/{{sample}}/hs_metrics/hs_metrics.txt",
    params:
        # Optional extra arguments. Here we reduce sample size
        # to reduce the runtime in our unit test.
        extra=f"--SAMPLE_SIZE {config['picard']['sample_size']} --NEAR_DISTANCE {config['picard']['near_distance']}",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    log:
        f"batches/{batch}/logs/picard/collect_hs_metrics/{{sample}}.log",
    wrapper:
        "v1.3.1/bio/picard/collecthsmetrics"

rule consolidate_hsmetrics:
    input:
        expand( f'batches/{batch}/' + '{sample}/hs_metrics/hs_metrics.txt', sample=sample2barcode.keys() ),
    output:
        f"batches/{batch}/stats/hs_metrics_consolidated.tsv",
    shell:
        '''
        awk '/METRICS/ {{getline; print;
                         getline; split(FILENAME,s,"/"); $67=s[3]; print}}' {input} \
        | awk 'NR==1 || !/^BAIT_SET/' > {output}
        '''

quickviewColumns = [
                    'SAMPLE',
                    'PCT_OFF_BAIT',
                    'MEAN_BAIT_COVERAGE',
                    'MEAN_TARGET_COVERAGE',
                    'FOLD_ENRICHMENT',
                    'ZERO_CVG_TARGETS_PCT',
                    'PCT_EXC_DUPE',
                    'FOLD_80_BASE_PENALTY',
                    'PCT_TARGET_BASES_20X',
                    'PCT_TARGET_BASES_30X',
                    'AT_DROPOUT',
                    'GC_DROPOUT'
                   ] 


rule hsmetrics_quickview:
    input:
        f"batches/{batch}/stats/hs_metrics_consolidated.tsv",
    output:
        f"batches/{batch}/stats/hs_metrics_consolidated_quickview.tsv",
    run:
        import csv
        from operator import itemgetter
        quickgetter = itemgetter(*quickviewColumns)
        with open( input[0], newline='' ) as hsmetrics,\
             open( output[0], 'w', newline='' ) as quickview:
            reader = csv.DictReader( map( lambda s: s.replace('\t',' '), hsmetrics), dialect='unix', delimiter=' ')
            writer = csv.DictWriter( quickview, quickviewColumns, dialect='unix', delimiter='\t', quoting=0 )
            writer.writeheader()
            for row in reader:
                writer.writerow( dict( zip( quickviewColumns, quickgetter(row) ) ) )
                       

targets.extend(
    [
        f"batches/{batch}/{sample}/hs_metrics/hs_metrics.txt"
        for sample in sample2barcode.keys()
    ]
)

targets.extend(
    [
         f"batches/{batch}/stats/hs_metrics_consolidated.tsv"
    ]
)

targets.extend(
    [
         f"batches/{batch}/stats/hs_metrics_consolidated_quickview.tsv"
    ]
)
