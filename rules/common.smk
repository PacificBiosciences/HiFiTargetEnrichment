localrules:
    tabix_vcf,

rule tabix_vcf:
    input:
        f"batches/{batch}/{{sample}}/{{prefix}}.vcf.gz",
    output:
        f"batches/{batch}/{{sample}}/{{prefix}}.vcf.gz.tbi",
    log:
        f"batches/{batch}/logs/tabix/index/{{sample}}.{{prefix}}.log",
    params:
        "-p vcf",
    conda:
        "envs/htslib.yaml"
    message:
        "Executing {rule}: Indexing {input}."
    shell:
        "tabix {params} {input} > {log} 2>&1"
