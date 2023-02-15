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
        "(paraphase -b {input.bam} -o {params.odir}) > {log} 2>&1"

targets.append(
    lambda wildcards: \
        [
            f"batches/{batch}/{sample}/paraphase/{sample}.json"
            for sample in _get_demuxed_samples( wildcards )
        ]
)