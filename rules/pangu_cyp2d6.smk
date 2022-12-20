localrules: fix_pangu_tsv

rule call_cyp2d6:
    input:
        bam=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.haplotagged.bam",
        index=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.haplotagged.bam.bai",
    output:
        json=f"batches/{batch}/{{sample}}/pangu/{{sample}}_report.json",
        tsv=f"batches/{batch}/{{sample}}/pangu/{{sample}}_pharmcat.tsv",
    params:
        prefix=f"batches/{batch}/{{sample}}/pangu/{{sample}}/",
    log:
        f"batches/{batch}/logs/pangu/{{sample}}.cyp2d6.log",
    benchmark:
        f"batches/{batch}/benchmarks/pangu/{{sample}}.call_cyp2d6.tsv"
    threads: 1
    conda:
        "envs/pangu.yaml"
    message:
        "Calling CYP2D6 for {wildcards.sample}"
    shell:
        "(pangu -m capture -p {params.prefix} {input.bam}) > {log} 2>&1"


rule fix_pangu_tsv:
    input: 
        tsv=f"batches/{batch}/{{sample}}/pangu/{{sample}}_pharmcat.tsv",
    output:
        f"batches/{batch}/{{sample}}/pangu/{{sample}}_pharmcat_fix.tsv",
    threads:
        1
    message:
        "Temp fix for pangu output with missing call for {wildcards.sample}"
    shell:
        '''
        awk 'BEGIN {{OFS="\t"}} !($2 ~ /\//) {{$2=$2"/[]"}} 1' {input} > {output}
        '''

targets.append(
    lambda wildcards: \
        [
            f"batches/{batch}/{sample}/pangu/{sample}_{fname}"
            for fname in [ "report.json", "pharmcat.tsv" ]
            for sample in _get_demuxed_samples( wildcards )
        ]
)
