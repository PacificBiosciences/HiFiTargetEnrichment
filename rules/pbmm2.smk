rule pbmm2_align_ubam:
    input:
        reference = config['ref']['fasta'],
        ref_index = config['ref']['index'],
        query = lambda wildcards: ubam_dict[wildcards.movie]
    output:
        bam = temp(f"samples/{sample}/aligned/{{movie}}.{ref}.bam"),
        bai = temp(f"samples/{sample}/aligned/{{movie}}.{ref}.bam.bai")
    log: f"samples/{sample}/logs/pbmm2/align/{{movie}}.{ref}.log"
    params:
        sample = f"{sample}",
        preset = "CCS",
        extra = "--sort --unmapped -c 0 -y 70",
        loglevel = "INFO"
    threads: 24
    conda: "envs/pbmm2.yaml"
    message: "Aligning {input.query} to {input.reference} using pbmm2 with --preset {params.preset} {params.extra}."
    shell:
        """
        (pbmm2 align --num-threads {threads} \
            --preset {params.preset} \
            --sample {params.sample} \
            --log-level {params.loglevel} \
            {params.extra} \
            {input.reference} \
            {input.query} \
            {output.bam}) > {log} 2>&1
        """


rule pbmm2_align_fastq:
    input:
        reference = config['ref']['fasta'],
        ref_index = config['ref']['index'],
        query = lambda wildcards: fastq_dict[wildcards.movie]
    output:
        bam = temp(f"samples/{sample}/aligned/{{movie}}.{ref}.bam"),
        bai = temp(f"samples/{sample}/aligned/{{movie}}.{ref}.bam.bai")
    log: f"samples/{sample}/logs/pbmm2/align/{{movie}}.{ref}.log"
    params:
        sample = f"{sample}",
        preset = "CCS",
        extra = "--sort --unmapped -c 0 -y 70",
        loglevel = "INFO"
    threads: 24
    conda: "envs/pbmm2.yaml"
    message: "Aligning {input.query} to {input.reference} using pbmm2 with --preset {params.preset} {params.extra}."
    shell:
        """
        (pbmm2 align --num-threads {threads} \
            --preset {params.preset} \
            --sample {params.sample} \
            --log-level {params.loglevel} \
            {params.extra} \
            {input.reference} \
            {input.query} \
            {output.bam}) > {log} 2>&1
        """
