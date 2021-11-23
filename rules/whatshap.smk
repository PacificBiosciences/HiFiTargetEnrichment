rule whatshap_phase_round1:
    input:
        reference = config['ref']['fasta'],
        vcf = f"samples/{sample}/deepvariant_intermediate/{sample}.{ref}.deepvariant.vcf.gz",
        tbi = f"samples/{sample}/deepvariant_intermediate/{sample}.{ref}.deepvariant.vcf.gz.tbi",
        phaseinput = abams,
        phaseinputindex = [f"{x}.bai" for x in abams]
    output: temp(f"samples/{sample}/whatshap_intermediate/{sample}.{ref}.deepvariant.phased.vcf.gz")
    log: f"samples/{sample}/logs/whatshap/phase/{sample}.{ref}.whatshap_intermediate.log"
    conda: "envs/whatshap.yaml"
    message: "Phasing {input.vcf} using {input.phaseinput}."
    shell:
        """
        (whatshap phase \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.phaseinput}) > {log} 2>&1
        """


rule whatshap_haplotag_round1:
    input:
        reference = config['ref']['fasta'],
        vcf = f"samples/{sample}/whatshap_intermediate/{sample}.{ref}.deepvariant.phased.vcf.gz",
        tbi = f"samples/{sample}/whatshap_intermediate/{sample}.{ref}.deepvariant.phased.vcf.gz.tbi",
        bam = lambda wildcards: abam_dict[wildcards.movie],
        bai = lambda wildcards: f"{abam_dict[wildcards.movie]}.bai"
    output: temp(f"samples/{sample}/whatshap_intermediate/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam")
    log: f"samples/{sample}/logs/whatshap/haplotag/{sample}.{ref}.{{movie}}.whatshap_intermediate.log"
    params: "--tag-supplementary"
    conda: "envs/whatshap.yaml"
    message: "Haplotagging {input.bam} using phase information from {input.vcf}."
    shell:
        """
        (whatshap haplotag {params} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.bam}) > {log} 2>&1
        """


rule samtools_index_bam_haplotag1:
    input: f"samples/{sample}/whatshap_intermediate/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam"
    output: temp(f"samples/{sample}/whatshap_intermediate/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam.bai")
    log: f"samples/{sample}/logs/samtools/index/whatshap_intermediate/{sample}.{ref}.{{movie}}.log"
    threads: 4
    conda: "envs/samtools.yaml"
    message: "Indexing {input}."
    shell: "(samtools index -@ 3 {input}) > {log} 2>&1"


rule whatshap_phase_round2:
    input:
        reference = config['ref']['fasta'],
        vcf = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz",
        tbi = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz.tbi",
        phaseinput = abams,
        phaseinputindex = [f"{x}.bai" for x in abams]
    output: temp(f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz")
    log: f"samples/{sample}/logs/whatshap/phase/{sample}.{ref}.log"
    params: extra = "--indels"
    conda: "envs/whatshap.yaml"
    message: "Phasing {input.vcf} using {input.phaseinput}."
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
        vcf = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz",
        tbi = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz.tbi",
        chr_lengths = config['ref']['chr_lengths']
    output:
        gtf = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.gtf",
        tsv = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.tsv",
        blocklist = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.blocklist"
    log: f"samples/{sample}/logs/whatshap/stats/{sample}.{ref}.log"
    conda: "envs/whatshap.yaml"
    message: "Calculating phasing stats for {input.vcf}."
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
        reference = config['ref']['fasta'],
        vcf = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz",
        tbi = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.phased.vcf.gz.tbi",
        bam = lambda wildcards: abam_dict[wildcards.movie],
        bai = lambda wildcards: f"{abam_dict[wildcards.movie]}.bai"
    output: temp(f"samples/{sample}/whatshap/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam")
    log: f"samples/{sample}/logs/whatshap/haplotag/{sample}.{ref}.{{movie}}.log"
    params: "--tag-supplementary"
    conda: "envs/whatshap.yaml"
    message: "Haplotagging {input.bam} using phase information from {input.vcf}."
    shell:
        """
        (whatshap haplotag {params} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.bam}) > {log} 2>&1
        """


rule samtools_index_bam_haplotag2:
    input: f"samples/{sample}/whatshap/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam"
    output: temp(f"samples/{sample}/whatshap/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam.bai")
    log: f"samples/{sample}/logs/samtools/index/whatshap/{sample}.{ref}.{{movie}}.log"
    threads: 4
    conda: "envs/samtools.yaml"
    message: "Indexing {input}."
    shell: "(samtools index -@ 3 {input}) > {log} 2>&1"


rule merge_haplotagged_bams:
    input: expand(f"samples/{sample}/whatshap/{sample}.{ref}.{{movie}}.deepvariant.haplotagged.bam", movie=movies)
    output: f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam"
    log: f"samples/{sample}/logs/samtools/merge/{sample}.{ref}.haplotag.log"
    threads: 8
    conda: "envs/samtools.yaml"
    message: "Merging {input}."
    shell: "(samtools merge -@ 7 {output} {input}) > {log} 2>&1"


rule samtools_index_merged_bam:
    input: f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam"
    output: f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam.bai"
    log: f"samples/{sample}/logs/samtools/index/whatshap/{sample}.{ref}.log"
    threads: 4
    conda: "envs/samtools.yaml"
    message: "Indexing {input}."
    shell: "(samtools index -@ 3 {input}) > {log} 2>&1"
