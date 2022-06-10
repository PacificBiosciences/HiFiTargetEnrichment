shards = [f"{x:05}" for x in range(config["N_SHARDS"])]


rule deepvariant_make_examples:
    input:
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        bai=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
        reference=config["ref"]["fasta"],
    output:
        tfrecord=(
            f"batches/{batch}/{{sample}}/deepvariant/examples/examples.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz"
        ),
        nonvariant_site_tfrecord=f"batches/{batch}/{{sample}}/deepvariant/examples/gvcf.tfrecord-{{shard}}-of-{config['N_SHARDS']:05}.gz",
    log:
        f"batches/{batch}/logs/deepvariant/make_examples/{{sample}}.{ref}.{{shard}}-of-{config['N_SHARDS']:05}.log",
    benchmark:
        f"batches/{batch}/benchmarks/deepvariant/{{sample}}.{{shard}}.dv_make_examples.tsv"
    container:
        f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params:
        vsc_min_fraction_indels="0.12",
        pileup_image_width=199,
        shard='{shard}',
        examples=f"batches/{batch}/{{sample}}/deepvariant/examples/examples.tfrecord@{config['N_SHARDS']}.gz",
        gvcf=f"batches/{batch}/{{sample}}/deepvariant/examples/gvcf.tfrecord@{config['N_SHARDS']}.gz",
    message:
        "DeepVariant make_examples {wildcards.shard} for {input.bam}."
    shell:
        """
        (/opt/deepvariant/bin/make_examples \
            --add_hp_channel \
            --alt_aligned_pileup=diff_channels \
            --min_mapping_quality=1 \
            --parse_sam_aux_fields \
            --partition_size=25000 \
            --max_reads_per_partition=600 \
            --phase_reads \
            --pileup_image_width {params.pileup_image_width} \
            --norealign_reads \
            --sort_by_haplotypes \
            --track_ref_reads \
            --vsc_min_fraction_indels {params.vsc_min_fraction_indels} \
            --mode calling \
            --ref {input.reference} \
            --reads {input.bam} \
            --examples {params.examples} \
            --gvcf {params.gvcf} \
            --task {params.shard}) > {log} 2>&1
        """

rule deepvariant_call_variants_gpu:
    input:
        expand(
            f"batches/{batch}/" + "{{sample}}/deepvariant/examples/examples.tfrecord-{shard}" + f"-of-{config['N_SHARDS']:05}.gz",
            shard=shards,
        ),
    output:
        (
            f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.call_variants_output.tfrecord.gz"
        ),
    log:
        f"batches/{batch}/logs/deepvariant/call_variants/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/deepvariant/{{sample}}.dv_call_variants.tsv"
    container:
        f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params:
        model="/opt/models/pacbio/model.ckpt",
        examples=f"batches/{batch}/{{sample}}/deepvariant/examples/examples.tfrecord@{config['N_SHARDS']:05}.gz",
    threads: 8
    message:
        "DeepVariant call_variants for {input}."
    shell:
        """
        (echo "CUDA_VISIBLE_DEVICES=" $CUDA_VISIBLE_DEVICES; \
         /opt/deepvariant/bin/call_variants \
            --outfile {output} \
            --examples {params.examples} \
            --checkpoint {params.model}) > {log} 2>&1
        """


rule deepvariant_postprocess_variants:
    input:
        tfrecord=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.call_variants_output.tfrecord.gz",
        nonvariant_site_tfrecord=expand(
            f"batches/{batch}/" + "{{sample}}/deepvariant/examples/gvcf.tfrecord-{shard}" + f"-of-{config['N_SHARDS']:05}.gz",
            shard=shards
        ),
        reference=config["ref"]["fasta"],
    output:
        vcf=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.vcf.gz",
        vcf_index=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.vcf.gz.tbi",
        gvcf=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.g.vcf.gz",
        gvcf_index=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.g.vcf.gz.tbi",
        report=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.visual_report.html",
    params:
        nonvariant_tfrecs=f"batches/{batch}/{{sample}}/deepvariant/examples/gvcf.tfrecord@{config['N_SHARDS']}.gz",
    log:
        f"batches/{batch}/logs/deepvariant/postprocess_variants/{{sample}}.{ref}.log",
    benchmark:
        f"batches/{batch}/benchmarks/deepvariant/{{sample}}.dv_postprocess.tsv"
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
            --outfile {output.vcf} \
            --nonvariant_site_tfrecord_path {params.nonvariant_tfrecs} \
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
