localrules:
    tabix_vcf,
    #bgzip_vcf_common,


#rule bgzip_vcf_common:
#    input:
#        f"batches/{batch}/{{sample}}/{{prefix}}.vcf",
#    output:
#        f"batches/{batch}/{{sample}}/{{prefix}}.vcf.gz",
#    log:
#        f"batches/{batch}/logs/bgzip/{{sample}}.{{prefix}}.log",
#    threads: 2
#    conda:
#        "envs/htslib.yaml"
#    message:
#        "Executing {rule}: Compressing {input}."
#    shell:
#        "(bgzip --threads {threads} {input}) > {log} 2>&1"


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
