localrules: bgzip_vcf, tabix_vcf


rule bgzip_vcf:
    input: f"samples/{sample}/{{prefix}}.vcf"
    output: f"samples/{sample}/{{prefix}}.vcf.gz"
    log: f"samples/{sample}/logs/bgzip/{{prefix}}.log"
    threads: 2
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Compressing {input}."
    shell: "(bgzip --threads {threads} {input}) > {log} 2>&1"


rule tabix_vcf:
    input: f"samples/{sample}/{{prefix}}.vcf.gz"
    output: f"samples/{sample}/{{prefix}}.vcf.gz.tbi"
    log: f"samples/{sample}/logs/tabix/index/{{prefix}}.log"
    params: "-p vcf"
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Indexing {input}."
    shell: "tabix {params} {input} > {log} 2>&1"
