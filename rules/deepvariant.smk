shards = [f"{x:05}" for x in range(config["N_SHARDS"])]


rule deepvariant_make_examples_round1:
    input:
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        bai=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
        reference=config["ref"]["fasta"],
    output:
        tfrecord=(
            f"batches/{batch}/{{sample}}/deepvariant_intermediate/examples/examples.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz"
        ),
    log:
        f"batches/{batch}/logs/deepvariant_intermediate/make_examples/{{sample}}.{ref}.{{shard}}-of-{config['N_SHARDS']:05}.log",
    benchmark:
        f"batches/{batch}/benchmarks/deepvariant/{{sample}}.{{shard}}.dv_make_examples1.tsv"
    container:
        f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params:
        vsc_min_fraction_indels="0.12",
        pileup_image_width=199,
        shard=lambda wildcards: wildcards.shard,
        examples=f"batches/{batch}/{{sample}}/deepvariant_intermediate/examples/examples.tfrecord@{config['N_SHARDS']}.gz",
    message:
        "DeepVariant make_examples {wildcards.shard} for {input.bam}."
    shell:
        """
        (/opt/deepvariant/bin/make_examples \
            --norealign_reads \
            --vsc_min_fraction_indels {params.vsc_min_fraction_indels} \
            --pileup_image_width {params.pileup_image_width} \
            --alt_aligned_pileup=diff_channels \
            --add_hp_channel \
            --sort_by_haplotypes \
            --parse_sam_aux_fields \
            --mode calling \
            --ref {input.reference} \
            --reads {input.bam} \
            --examples {params.examples} \
            --task {params.shard}) > {log} 2>&1
        """


rule deepvariant_call_variants_gpu_round1:
    input:
        expand(
            f"batches/{batch}/" + "{{sample}}/deepvariant_intermediate/examples/examples.tfrecord-{shard}" + f"-of-{config['N_SHARDS']:05}.gz",
            shard=shards,
        ),
    output:
        (
            f"batches/{batch}/{{sample}}/deepvariant_intermediate/{{sample}}.{ref}.call_variants_output.tfrecord.gz"
        ),
    log:
        f"batches/{batch}/logs/deepvariant_intermediate/call_variants/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/deepvariant/{{sample}}.dv_call_variants1.tsv"
    container:
        f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params:
        model="/opt/models/pacbio/model.ckpt",
        examples=f"batches/{batch}/{{sample}}/deepvariant_intermediate/examples/examples.tfrecord@{config['N_SHARDS']}.gz",
    threads: 8
    message:
        "DeepVariant call_variants for {input}."
    shell:
        """
        (/opt/deepvariant/bin/call_variants \
            --outfile {output} \
            --examples {params.examples} \
            --checkpoint {params.model}) > {log} 2>&1
        """


rule deepvariant_postprocess_variants_round1:
    input:
        tfrecord=f"batches/{batch}/{{sample}}/deepvariant_intermediate/{{sample}}.{ref}.call_variants_output.tfrecord.gz",
        reference=config["ref"]["fasta"],
    output:
        vcf=(
            f"batches/{batch}/{{sample}}/deepvariant_intermediate/{{sample}}.{ref}.deepvariant.vcf.gz"
        ),
        vcf_index=(
            f"batches/{batch}/{{sample}}/deepvariant_intermediate/{{sample}}.{ref}.deepvariant.vcf.gz.tbi"
        ),
        report=(
            f"batches/{batch}/{{sample}}/deepvariant_intermediate/{{sample}}.{ref}.deepvariant.visual_report.html"
        ),
    log:
        f"batches/{batch}/logs/deepvariant_intermediate/postprocess_variants/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/deepvariant/{{sample}}.dv_postprocess1.tsv"
    container:
        f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    threads: 4
    message:
        "DeepVariant postprocess_variants for {input.tfrecord}."
    shell:
        """
        (/opt/deepvariant/bin/postprocess_variants \
            --ref {input.reference} \
            --infile {input.tfrecord} \
            --outfile {output.vcf}) > {log} 2>&1
        """


rule deepvariant_make_examples_round2:
    input:
        bam=f"batches/{batch}/{{sample}}/whatshap_intermediate/{{sample}}.{ref}.deepvariant.haplotagged.bam",
        bai=f"batches/{batch}/{{sample}}/whatshap_intermediate/{{sample}}.{ref}.deepvariant.haplotagged.bam.bai",
        reference=config["ref"]["fasta"],
    output:
        tfrecord=(
            f"batches/{batch}/{{sample}}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz"
        ),
        nonvariant_site_tfrecord=(
            f"batches/{batch}/{{sample}}/deepvariant/examples/gvcf.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz"
        ),
    log:
        f"batches/{batch}/logs/deepvariant/make_examples/{{sample}}.{ref}.{{shard}}-of-{config['N_SHARDS']:05}.log",
    benchmark:
        f"batches/{batch}/benchmarks/deepvariant/{{sample}}.{{shard}}.dv_make_examples2.tsv"
    container:
        f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params:
        vsc_min_fraction_indels="0.12",
        pileup_image_width=199,
        shard=lambda wildcards: wildcards.shard,
        examples=f"batches/{batch}/{{sample}}/deepvariant/examples/examples.tfrecord@{config['N_SHARDS']}.gz",
        gvcf=f"batches/{batch}/{{sample}}/deepvariant/examples/gvcf.tfrecord@{config['N_SHARDS']}.gz",
    message:
        "DeepVariant make_examples {wildcards.shard} for {input.bam}."
    shell:
        """
        (/opt/deepvariant/bin/make_examples \
            --norealign_reads \
            --vsc_min_fraction_indels {params.vsc_min_fraction_indels} \
            --pileup_image_width {params.pileup_image_width} \
            --alt_aligned_pileup=diff_channels \
            --add_hp_channel \
            --sort_by_haplotypes \
            --parse_sam_aux_fields \
            --mode calling \
            --ref {input.reference} \
            --reads {input.bam} \
            --examples {params.examples} \
            --gvcf {params.gvcf} \
            --task {params.shard}) > {log} 2>&1
        """


rule deepvariant_call_variants_gpu_round2:
    input:
        expand(
            f"batches/{batch}/" + "{{sample}}/deepvariant/examples/examples.tfrecord-{shard:05}" + f"-of-{config['N_SHARDS']:05}.gz",
            shard=shards,
        ),
    output:
        (
            f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.call_variants_output.tfrecord.gz"
        ),
    log:
        f"batches/{batch}/logs/deepvariant/call_variants/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/deepvariant/{{sample}}.dv_call_variants2.tsv"
    container:
        f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params:
        model="/opt/models/pacbio/model.ckpt",
        shards=config["N_SHARDS"],
        examples=f"batches/{batch}/{{sample}}/deepvariant/examples/examples.tfrecord@{config['N_SHARDS']}.gz"
    message:
        "DeepVariant call_variants for {input}."
    threads: 8
    shell:
        """
        (echo "CUDA_VISIBLE_DEVICES=" $CUDA_VISIBLE_DEVICES; \
        /opt/deepvariant/bin/call_variants \
            --outfile {output} \
            --examples {params.examples} \
            --checkpoint {params.model}) > {log} 2>&1
        """


rule deepvariant_postprocess_variants_round2:
    input:
        tfrecord=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.call_variants_output.tfrecord.gz",
        nonvariant_site_tfrecord=expand(
            f"batches/{batch}/" + "{{sample}}/deepvariant/examples/gvcf.tfrecord-{shard:05}" f"-of-{config['N_SHARDS']:05}.gz",
            shard=range(config["N_SHARDS"]),
        ),
        reference=config["ref"]["fasta"],
    output:
        vcf=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.vcf.gz",
        vcf_index=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.vcf.gz.tbi",
        gvcf=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.g.vcf.gz",
        gvcf_index=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.g.vcf.gz.tbi",
        report=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.visual_report.html",
    log:
        f"batches/{batch}/logs/deepvariant/postprocess_variants/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/deepvariant/{{sample}}.dv_postprocess2.tsv"
    params:
        tfrecpath=f'batches/{batch}/{{sample}}/deepvariant/examples/gvcf.tfrecord@{config["N_SHARDS"]}.gz',
    container:
        f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    message:
        "DeepVariant postprocess_variants for {input.tfrecord}."
    threads: 4
    shell:
        """
        (/opt/deepvariant/bin/postprocess_variants \
            --ref {input.reference} \
            --infile {input.tfrecord} \
            --outfile {output.vcf} \
            --nonvariant_site_tfrecord_path {params.tfrecpath} \
            --gvcf_outfile {output.gvcf}) > {log} 2>&1
        """


rule deepvariant_bcftools_stats:
    input:
        f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.vcf.gz",
    output:
        f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.vcf.stats.txt",
    log:
        f"batches/{batch}/logs/bcftools/stats/{{sample}}.{ref}.deepvariant.vcf.log",
    benchmark:
        f"batches/{batch}/benchmarks/deepvariant/{{sample}}.dv_bcftools_stats.tsv"
    params:
        f"--fasta-ref {config['ref']['fasta']} --apply-filters PASS -s {{sample}}",
    threads: 4
    conda:
        "envs/bcftools.yaml"
    message:
        "Calculating VCF statistics for {input}."
    shell:
        "(bcftools stats --threads 3 {params} {input} > {output}) > {log} 2>&1"
