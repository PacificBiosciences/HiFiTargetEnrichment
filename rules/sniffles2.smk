rule sniffles_discover:
    input:
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        bai=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
    output:
        temp( f"batches/{batch}/{{sample}}/sniffles/{{sample}}.{ref}.sv.gz" ),
    log:
        f"batches/{batch}/logs/sniffles/{{sample}}.{ref}.log"
    benchmark:
        f"batches/{batch}/benchmarks/sniffles/{{sample}}.{ref}.tsv"
    params:
        min_sv_len=config['min_sv_len'],
        sv_threads=config['sv_threads'],
    conda:
        "envs/sniffles.yaml"
    message:
        "Executing {rule}: Discovering structural variant signatures in {wildcards.sample} from {input.bam}."
    shell:
        """
            sniffles --minsvlen {params.min_sv_len} --sample-id {wildcards.sample} -t {params.sv_threads} --input {input.bam} --vcf {output} > {log} 2>&1
        """