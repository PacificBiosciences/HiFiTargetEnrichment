def _get_bam_by_type( index=False ):
    ext = '.bai' if index else ''
    def getter( wildcards ):
        if wildcards.source == "asm":
            return f"batches/{batch}/{wildcards.sample}/hifiasm/{wildcards.sample}.asm.{ref}.bam{ext}",
        else:
            return f"batches/{batch}/{wildcards.sample}/aligned/{wildcards.sample}.{ref}.bam{ext}",
    return getter


rule sniffles_discover:
    input:
        bam=_get_bam_by_type(),
        bai=_get_bam_by_type(index=True),
    output:
        f"batches/{batch}/{{sample}}/sniffles/{{sample}}.{ref}.{{source}}.sv.vcf.gz",
    log:
        f"batches/{batch}/logs/sniffles/{{sample}}.{{source}}.{ref}.log"
    benchmark:
        f"batches/{batch}/benchmarks/sniffles/{{sample}}.{{source}}.{ref}.tsv"
    params:
        min_sv_len=config['min_sv_len'],
        sv_threads=config['sv_threads'],
    conda:
        "envs/sniffles.yaml"
    message:
        "Executing {rule}: Discovering structural variant signatures in {wildcards.sample} from {input.bam}."
    shell:
        """
            sniffles \
                --minsvlen {params.min_sv_len} \
                --sample-id {wildcards.sample} \
                -t {params.sv_threads} \
                --input {input.bam} \
                --vcf {output} > {log} 2>&1
        """

targets.append(
    lambda wildcards: \
        [
            f"batches/{batch}/{sample}/sniffles/{sample}.{ref}.{source}.sv.vcf.gz.tbi"
            for sample in _get_demuxed_samples( wildcards )
            for source in [ "asm","mapped"]
       ]
)
