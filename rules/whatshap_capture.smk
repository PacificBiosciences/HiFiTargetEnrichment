rule whatshap_phase_round1:
    input:
        reference=config["ref"]["fasta"],
        vcf=f"batches/{batch}/{{sample}}/deepvariant_intermediate/{{sample}}.{ref}.deepvariant.vcf.gz",
        tbi=f"batches/{batch}/{{sample}}/deepvariant_intermediate/{{sample}}.{ref}.deepvariant.vcf.gz.tbi",
        phaseinput=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        phaseinputindex=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
    output:
        (
            f"batches/{batch}/{{sample}}/whatshap_intermediate/{{sample}}.{ref}.deepvariant.phased.vcf.gz"
        ),
    log:
        f"batches/{batch}/logs/whatshap/phase/{{sample}}.{ref}.whatshap_intermediate.log",
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.wh_phase1.tsv"
    conda:
        "envs/whatshap.yaml"
    message:
        "Phasing {input.vcf} using {input.phaseinput}."
    shell:
        """
        (whatshap phase \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.phaseinput}) > {log} 2>&1
        """


rule whatshap_haplotag_round1:
    input:
        reference=config["ref"]["fasta"],
        vcf=f"batches/{batch}/{{sample}}/whatshap_intermediate/{{sample}}.{ref}.deepvariant.phased.vcf.gz",
        tbi=f"batches/{batch}/{{sample}}/whatshap_intermediate/{{sample}}.{ref}.deepvariant.phased.vcf.gz.tbi",
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        bai=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
    output:
        (
            f"batches/{batch}/{{sample}}/whatshap_intermediate/{{sample}}.{ref}.deepvariant.haplotagged.bam"
        ),
    log:
        f"batches/{batch}/logs/whatshap/haplotag/{{sample}}.{ref}.whatshap_intermediate.log",
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.wh_haplotaq1.tsv"
    params:
        "--tag-supplementary",
    conda:
        "envs/whatshap.yaml"
    message:
        "Haplotagging {input.bam} using phase information from {input.vcf}."
    shell:
        """
        (whatshap haplotag {params} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.bam}) > {log} 2>&1
        """


rule samtools_index_bam_haplotag1:
    input:
        f"batches/{batch}/{{sample}}/whatshap_intermediate/{{sample}}.{ref}.deepvariant.haplotagged.bam",
    output:
        (
            f"batches/{batch}/{{sample}}/whatshap_intermediate/{{sample}}.{ref}.deepvariant.haplotagged.bam.bai"
        ),
    log:
        f"batches/{batch}/logs/samtools/index/whatshap_intermediate/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.wh_index_bam_haplotag1.tsv"
    threads: 4
    conda:
        "envs/samtools.yaml"
    message:
        "Indexing {input}."
    shell:
        "(samtools index -@ 3 {input}) > {log} 2>&1"


rule whatshap_phase_round2:
    input:
        reference=config["ref"]["fasta"],
        vcf=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.vcf.gz",
        tbi=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.vcf.gz.tbi",
        phaseinput=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        phaseinputindex=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
    output:
        f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.vcf.gz",
    log:
        f"batches/{batch}/logs/whatshap/phase/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.wh_phase2.tsv"
    params:
        extra="--indels",
    conda:
        "envs/whatshap.yaml"
    message:
        "Phasing {input.vcf} using {input.phaseinput}."
    shell:
        """
        (whatshap phase {params.extra} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} \
            {input.phaseinput}) > {log} 2>&1
        """


rule whatshap_stats:
    input:
        vcf=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.vcf.gz",
        tbi=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.vcf.gz.tbi",
        chr_lengths=config["ref"]["chr_lengths"],
    output:
        gtf=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.gtf",
        tsv=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.tsv",
        blocklist=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.blocklist",
    log:
        f"batches/{batch}/logs/whatshap/stats/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.wh_stats.tsv"
    conda:
        "envs/whatshap.yaml"
    message:
        "Calculating phasing stats for {input.vcf}."
    shell:
        """
        (whatshap stats \
            --gtf {output.gtf} \
            --tsv {output.tsv} \
            --block-list {output.blocklist} \
            --chr-lengths {input.chr_lengths} \
            {input.vcf}) > {log} 2>&1
        """


rule whatshap_haplotag_round2:
    input:
        reference=config["ref"]["fasta"],
        vcf=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.vcf.gz",
        tbi=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.vcf.gz.tbi",
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        bai=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
    output:
        f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.haplotagged.bam",
    log:
        f"batches/{batch}/logs/whatshap/haplotag/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.wh_haplotag2.tsv"
    params:
        "--tag-supplementary",
    conda:
        "envs/whatshap.yaml"
    message:
        "Haplotagging {input.bam} using phase information from {input.vcf}."
    shell:
        """
        (whatshap haplotag {params} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.bam}) > {log} 2>&1
        """


rule samtools_index_bam_haplotag2:
    input:
        f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.haplotagged.bam",
    output:
        f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.haplotagged.bam.bai",
    log:
        f"batches/{batch}/logs/samtools/index/whatshap/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/{{sample}}.wh_index_bam_haplotag2.tsv"
    threads: 4
    conda:
        "envs/samtools.yaml"
    message:
        "Indexing {input}."
    shell:
        "(samtools index -@ 3 {input}) > {log} 2>&1"
