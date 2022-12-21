rule expand_gvcf_over_targets:
    input:
        gvcf=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.g.vcf.gz",
        gvcf_index=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.g.vcf.gz.tbi",
        reference=config["ref"]["fasta"],
    output:
        temp(f"batches/{batch}/{{sample}}/annotate/{{sample}}.{ref}.expanded.g.vcf.gz"),
    params:
        variants=config["annotate"]["variants"],
    log:
        f"batches/{batch}/logs/annotate/bcftools_convert/{{sample}}.{ref}.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "annotation: expanding gvcf for {wildcards.sample}"
    shell:
        '''
        (bcftools convert --gvcf2vcf \
         -f {input.reference} \
         -R {params.variants} \
         -Oz \
         -o {output} \
         {input.gvcf}) > {log} 2>&1
        '''

rule annotate_expanded_gvcf:
    input:
        gvcf=f"batches/{batch}/{{sample}}/annotate/{{sample}}.{ref}.expanded.g.vcf.gz",
        gvcf_index=f"batches/{batch}/{{sample}}/annotate/{{sample}}.{ref}.expanded.g.vcf.gz.tbi",
    output:
        temp(f"batches/{batch}/{{sample}}/annotate/{{sample}}.{ref}.expanded.annotated.g.vcf.gz"),
    params:
        region=f'batches/{batch}/glnexus/regions.bed',
        variants=config["annotate"]["variants"],
    log:
        f"batches/{batch}/logs/annotate/bcftools_annotate/{{sample}}.{ref}.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "annotating gvcf for {wildcards.sample}"
    shell:
        '''
        (bcftools annotate -c ID \
                           -R {params.region} \
                           -a {params.variants} \
                           -Oz -o {output} \
                           {input.gvcf}) > {log} 2>&1
        '''

rule enforce_annotate_all_rows:
    input: 
        vcf=f"batches/{batch}/{{sample}}/annotate/{{sample}}.{ref}.expanded.annotated.g.vcf.gz",
        vcf_index=f"batches/{batch}/{{sample}}/annotate/{{sample}}.{ref}.expanded.annotated.g.vcf.gz.tbi",
    output:
        temp(f"batches/{batch}/{{sample}}/annotate/{{sample}}.{ref}.dv.annotated.merged.g.vcf.gz"),
    params:
        region=f'batches/{batch}/glnexus/regions.bed',
        variants=config["annotate"]["variants"],
    log:
        f"batches/{batch}/logs/annotate/bcftools_merge_indels/{{sample}}.{ref}.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "annotate: enforcing inclusion for annotations for {wildcards.sample}"
    shell:
        '''
        (bcftools merge \
                  -R {params.region} \
                  -Oz -o {output} \
                  {input.vcf} {params.variants}) > {log} 2>&1
        '''

def _get_annotated_gvcfs( wildcards ):
    return [ f"batches/{batch}/{sample}/annotate/{sample}.{ref}.dv.annotated.merged.g.vcf.gz"
             for sample in _get_demuxed_samples( wildcards ) ]

def _get_annotated_gvcfs_idxs( wildcards ):
    return [ f"batches/{batch}/{sample}/annotate/{sample}.{ref}.dv.annotated.merged.g.vcf.gz.tbi"
             for sample in _get_demuxed_samples( wildcards ) ]

rule merge_annotated_gvcfs:
    input:
        gvcfs=_get_annotated_gvcfs,
        idxs=_get_annotated_gvcfs_idxs,
    output:
        f"batches/{batch}/merged_gvcf/allSamples.{ref}.dv.annotated.merged.g.vcf.gz",
    params:
        region=f'batches/{batch}/glnexus/regions.bed',
    log:
        f"batches/{batch}/logs/annotate/bcftools_merge_all.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "annotate: merging annotated vcfs"
    shell:
        '''
        (bcftools merge \
                  -R {params.region} \
                  -Oz -o {output} \
                  {input.gvcfs}) > {log} 2>&1
        '''

targets.append( f"batches/{batch}/merged_gvcf/allSamples.{ref}.dv.annotated.merged.g.vcf.gz.tbi" )
