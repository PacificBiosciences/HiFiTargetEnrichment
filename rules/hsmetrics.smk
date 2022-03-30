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
        f"batches/{batch}/picard/{{bedfile}}.interval_list",
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
        expand( f"batches/{batch}/{{sample}}/hs_metrics/hs_metrics.txt",sample=sample2barcode.keys() ),
    output:
        f"batches/{batch}/stats/hs_metrics_consolidated.tsv",
    shell:
        '''
        awk '/METRICS/ {{getline; print;
                         getline; split(FILENAME,s,"/"); $67=s[3]; print}}' {input} \
        | awk 'NR==1 || !/^BAIT_SET/' > {output}
        '''

targets.extend(
    [
        f"batches/{batch}/{sample}/hs_metrics/hs_metrics.txt"
        for sample in sample2barcode.keys()
    ]
)

targets.append( f"batches/{batch}/stats/hs_metrics_consolidated.tsv" )
