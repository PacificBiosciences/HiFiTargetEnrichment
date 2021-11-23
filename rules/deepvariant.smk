shards = [f"{x:05}" for x in range(config['N_SHARDS'])]


rule deepvariant_make_examples_round1:
    input:
        bams = abams,
        bais = [f"{x}.bai" for x in abams],
        reference = config['ref']['fasta']
    output:
        tfrecord = temp(f"samples/{sample}/deepvariant_intermediate/examples/examples.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz")
    log: f"samples/{sample}/logs/deepvariant_intermediate/make_examples/{sample}.{ref}.{{shard}}-of-{config['N_SHARDS']:05}.log"
    container: f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params:
        vsc_min_fraction_indels = "0.12",
        shard = lambda wildcards: wildcards.shard,
        reads = ','.join(abams)
    message: "DeepVariant make_examples {wildcards.shard} for {input.bams}."
    shell:
        f"""
        (/opt/deepvariant/bin/make_examples \
            --norealign_reads \
            --vsc_min_fraction_indels {{params.vsc_min_fraction_indels}} \
            --alt_aligned_pileup=diff_channels \
            --add_hp_channel \
            --sort_by_haplotypes \
            --parse_sam_aux_fields \
            --mode calling \
            --ref {{input.reference}} \
            --reads {{params.reads}} \
            --examples samples/{sample}/deepvariant_intermediate/examples/examples.tfrecord@{config['N_SHARDS']}.gz \
            --task {{wildcards.shard}}) > {{log}} 2>&1
        """


rule deepvariant_call_variants_gpu_round1:
    input: expand(f"samples/{sample}/deepvariant_intermediate/examples/examples.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz", shard=shards)
    output: temp(f"samples/{sample}/deepvariant_intermediate/{sample}.{ref}.call_variants_output.tfrecord.gz")
    log: f"samples/{sample}/logs/deepvariant_intermediate/call_variants/{sample}.{ref}.log"
    container: f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params: model = "/opt/models/pacbio/model.ckpt"
    threads: 8
    message: "DeepVariant call_variants for {input}."
    shell:
        f"""
        (/opt/deepvariant/bin/call_variants \
            --outfile {{output}} \
            --examples samples/{sample}/deepvariant_intermediate/examples/examples.tfrecord@{config['N_SHARDS']}.gz \
            --checkpoint {{params.model}}) > {{log}} 2>&1
        """


rule deepvariant_postprocess_variants_round1:
    input:
        tfrecord = f"samples/{sample}/deepvariant_intermediate/{sample}.{ref}.call_variants_output.tfrecord.gz",
        reference = config['ref']['fasta']
    output:
        vcf = temp(f"samples/{sample}/deepvariant_intermediate/{sample}.{ref}.deepvariant.vcf.gz"),
        vcf_index = temp(f"samples/{sample}/deepvariant_intermediate/{sample}.{ref}.deepvariant.vcf.gz.tbi"),
        report = temp(f"samples/{sample}/deepvariant_intermediate/{sample}.{ref}.deepvariant.visual_report.html")
    log: f"samples/{sample}/logs/deepvariant_intermediate/postprocess_variants/{sample}.{ref}.log"
    container: f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    threads: 4
    message: "DeepVariant postprocess_variants for {input.tfrecord}."
    shell:
        """
        (/opt/deepvariant/bin/postprocess_variants \
            --ref {input.reference} \
            --infile {input.tfrecord} \
            --outfile {output.vcf}) > {log} 2>&1
        """


haplotagged_abams = [f"samples/{sample}/whatshap_intermediate/{sample}.{ref}.{movie}.deepvariant.haplotagged.bam" for movie in movies]


rule deepvariant_make_examples_round2:
    input:
        bams = haplotagged_abams,
        bais = [f"{x}.bai" for x in haplotagged_abams],
        reference = config['ref']['fasta']
    output:
        tfrecord = temp(f"samples/{sample}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz"),
        nonvariant_site_tfrecord = temp(f"samples/{sample}/deepvariant/examples/gvcf.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz")
    log: f"samples/{sample}/logs/deepvariant/make_examples/{sample}.{ref}.{{shard}}-of-{config['N_SHARDS']:05}.log"
    container: f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params:
        vsc_min_fraction_indels = "0.12",
        shard = lambda wildcards: wildcards.shard,
        reads = ','.join(haplotagged_abams)
    message: "DeepVariant make_examples {wildcards.shard} for {input.bams}."
    shell:
        f"""
        (/opt/deepvariant/bin/make_examples \
            --norealign_reads \
            --vsc_min_fraction_indels {{params.vsc_min_fraction_indels}} \
            --alt_aligned_pileup=diff_channels \
            --add_hp_channel \
            --sort_by_haplotypes \
            --parse_sam_aux_fields \
            --mode calling \
            --ref {{input.reference}} \
            --reads {{params.reads}} \
            --examples samples/{sample}/deepvariant/examples/examples.tfrecord@{config['N_SHARDS']}.gz \
            --gvcf samples/{sample}/deepvariant/examples/gvcf.tfrecord@{config['N_SHARDS']}.gz \
            --task {{wildcards.shard}}) > {{log}} 2>&1
        """


rule deepvariant_call_variants_gpu_round2:
    input: expand(f"samples/{sample}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz", shard=shards)
    output: temp(f"samples/{sample}/deepvariant/{sample}.{ref}.call_variants_output.tfrecord.gz")
    log: f"samples/{sample}/logs/deepvariant/call_variants/{sample}.{ref}.log"
    container: f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params: model = "/opt/models/pacbio/model.ckpt"
    message: "DeepVariant call_variants for {input}."
    threads: 8
    shell:
        f"""
        (echo "CUDA_VISIBLE_DEVICES=" $CUDA_VISIBLE_DEVICES; \
        /opt/deepvariant/bin/call_variants \
            --outfile {{output}} \
            --examples samples/{sample}/deepvariant/examples/examples.tfrecord@{config['N_SHARDS']}.gz \
            --checkpoint {{params.model}}) > {{log}} 2>&1
        """


rule deepvariant_postprocess_variants_round2:
    input:
        tfrecord = f"samples/{sample}/deepvariant/{sample}.{ref}.call_variants_output.tfrecord.gz",
        nonvariant_site_tfrecord = expand(f"samples/{sample}/deepvariant/examples/gvcf.tfrecord-{{shard:05}}-of-{config['N_SHARDS']:05}.gz",
                                          shard=range(config['N_SHARDS'])),
        reference = config['ref']['fasta']
    output:
        vcf = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz",
        vcf_index = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz.tbi",
        gvcf = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.g.vcf.gz",
        gvcf_index = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.g.vcf.gz.tbi",
        report = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.visual_report.html"
    log: f"samples/{sample}/logs/deepvariant/postprocess_variants/{sample}.{ref}.log"
    container: f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    message: "DeepVariant postprocess_variants for {input.tfrecord}."
    threads: 4
    shell:
        f"""
        (/opt/deepvariant/bin/postprocess_variants \
            --ref {{input.reference}} \
            --infile {{input.tfrecord}} \
            --outfile {{output.vcf}} \
            --nonvariant_site_tfrecord_path samples/{sample}/deepvariant/examples/gvcf.tfrecord@{config['N_SHARDS']}.gz \
            --gvcf_outfile {{output.gvcf}}) > {{log}} 2>&1
        """


rule deepvariant_bcftools_stats:
    input: f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz"
    output: f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.stats.txt"
    log: f"samples/{sample}/logs/bcftools/stats/{sample}.{ref}.deepvariant.vcf.log"
    params: f"--fasta-ref {config['ref']['fasta']} --apply-filters PASS -s {sample}"
    threads: 4
    conda: "envs/bcftools.yaml"
    message: "Calculating VCF statistics for {input}."
    shell: "(bcftools stats --threads 3 {params} {input} > {output}) > {log} 2>&1"
