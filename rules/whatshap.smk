rule whatshap_phase:
    input:
        reference=config["ref"]["fasta"],
        vcf=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.vcf.gz",
        tbi=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.vcf.gz.tbi",
        phaseinput=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        phaseinputindex=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
    output:
        vcf=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.vcf.gz",
        vcf_index=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.vcf.gz.tbi",
    log:
        f"batches/{batch}/logs/whatshap/phase/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/whatshap/{{sample}}.phase.tsv"
    params:
        extra="--indels",
    conda:
        "envs/whatshap.yaml"
    message:
        "Phasing {input.vcf} using {input.phaseinput}."
    shell:
        """
        (whatshap phase {params.extra} \
            ---output {output.vcf} \
            --reference {input.reference} \
            {input.vcf} \
            {input.phaseinput} && tabix {output.vcf}) > {log} 2>&1
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
        f"batches/{batch}/benchmarks/whatshap/{{sample}}.stats.tsv"
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


rule whatshap_haplotag:
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
        f"batches/{batch}/benchmarks/whatshap/{{sample}}.haplotag.tsv"
    params:
        "--tag-supplementary --ignore-read-groups",
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


rule samtools_index_bam_haplotag:
    input:
        f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.haplotagged.bam",
    output:
        f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.haplotagged.bam.bai",
    log:
        f"batches/{batch}/logs/samtools/index/whatshap/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/whatshap/{{sample}}.index_bam_haplotag.tsv"
    threads: 4
    conda:
        "envs/samtools.yaml"
    message:
        "Indexing {input}."
    shell:
        "(samtools index -@ 3 {input}) > {log} 2>&1"
