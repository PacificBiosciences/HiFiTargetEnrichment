def _get_gene_list():
    '''5th col of bed == "paraphase"'''
    rows = map( lambda s: s.split('\t'), open( config["targets"] ).read().strip().split('\n') )
    return '-g ' + ','.join( row[3] for row in rows if row[4] == "paraphase" )

rule call_paraphase:
    input:
        bam=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.haplotagged.bam",
        index=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.haplotagged.bam.bai",
    output:
        json=f"batches/{batch}/{{sample}}/paraphase/{{sample}}.json",
        bam=f"batches/{batch}/{{sample}}/paraphase/{{sample}}_realigned_tagged.bam",
        bai=f"batches/{batch}/{{sample}}/paraphase/{{sample}}_realigned_tagged.bam.bai",
    params:
        odir=f"batches/{batch}/{{sample}}/paraphase/",
        genes=_get_gene_list() if config["paraphase"]["labeled_bed"] else '',
        reference=config["ref"]["fasta"],
    log:
        f"batches/{batch}/logs/paraphase/{{sample}}.log",
    benchmark:
        f"batches/{batch}/benchmarks/paraphase/{{sample}}.tsv"
    threads: 1
    conda:
        "envs/paraphase.yaml"
    message:
        "Calling paraphase SMN1/2 for {wildcards.sample}"
    shell:
        "(paraphase {params.genes} -b {input.bam} -r {params.reference} -o {params.odir}) > {log} 2>&1"

targets.append(
    lambda wildcards: 
        [
            f"batches/{batch}/{sample}/paraphase/{sample}.json"
            for sample in _get_demuxed_samples( wildcards )
        ]
)