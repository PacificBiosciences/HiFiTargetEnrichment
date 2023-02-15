ruleorder: merge_analysis_vcf > bgzip_vcf


def _get_vcf_by_type( index=False ):
    ext = '.tbi' if index else ''
    def getter( wildcards ):
        if wildcards.source == "dv":
            return f"batches/{batch}/{wildcards.sample}/whatshap/{wildcards.sample}.{ref}.deepvariant.phased.vcf.gz{ext}"
        else:
            return f"batches/{batch}/{wildcards.sample}/hifiasm/{wildcards.sample}.asm.{ref}.htsbox.vcf.gz{ext}"
    return getter

rule intersect_targets_by_analysis:
    input:
        vcf=_get_vcf_by_type(),
        idx=_get_vcf_by_type(index=True),
    output:
        temp(f"batches/{batch}/{{sample}}/combined_vcf/intersect/{{sample}}.{{source}}.vcf"),
    params:
        targets=config['targets'],
        chr_len=config['ref']['chr_lengths'],
        buffer=10000,
        source='{source}',
    log:
        f"batches/{batch}/logs/combine_vcf/bedtools/{{sample}}.intesect_{{source}}.log",
    conda:
        "envs/samtools.yaml"
    message:
        "Intersecting {wildcards.source} variant targets for {wildcards.sample}"
    shell:
        '''
        (
        bedtools intersect \
            -wa \
            -header \
            -a {input.vcf} \
            -b <( \
                 bedtools slop \
                        -b {params.buffer} \
                        -i {params.targets} \
                        -g {params.chr_len} \
                 | awk -v src={params.source} '$5 ~ src' \
                ) \
            > {output}
        ) > {log} 2>&1
        '''

rule annotate_variant_source:
    input:
        vcf=f"batches/{batch}/{{sample}}/combined_vcf/intersect/{{sample}}.{{source}}.vcf",
    output:
        temp(f"batches/{batch}/{{sample}}/combined_vcf/annotate/{{sample}}.{{source}}.vcf"),
    params:
        old_desc='"Added by +fill-tags expression VS=[12]"',
        new_desc='"Variant Source (1=deepvariant, 2=hifiasm)"',
        vs=lambda wcs: {"dv":"1","asm":"2"}[wcs.source],
    log:
        f"batches/{batch}/logs/combine_vcf/bedtools/{{sample}}.annotate_{{source}}.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "Annotating {wildcards.source} variant targets for {wildcards.sample}"
    shell:
        '''
        (
        bcftools +fill-tags \
            {input.vcf} \
            -- \
            -t \'VS={params.vs}\' \
        | sed 's/{params.old_desc}/{params.new_desc}/' \
        | bcftools view \
            -Oz \
            -o {output} \
        ) > {log} 1>&2  
        '''
    
rule merge_analysis_vcf:
    input:
        vcf=[f"batches/{batch}/{{sample}}/combined_vcf/annotate/{{sample}}.{source}.vcf.gz" 
             for source in ["asm","dv"]],
        idx=[f"batches/{batch}/{{sample}}/combined_vcf/annotate/{{sample}}.{source}.vcf.gz.tbi" 
             for source in ["asm","dv"]],
    output:
        f"batches/{batch}/{{sample}}/combined_vcf/{{sample}}.{ref}.deepvariant_hifiasm.vcf.gz",
    params:
    log:
        f"batches/{batch}/logs/combine_vcf/bcftools/{{sample}}.merge.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "Merging deepvariant and assembly variants for {wildcards.sample}"
    shell:
        '''
        (bcftools merge -Oz {input.vcf} > {output}) > {log} 2>&1  
        '''

targets.append(
    lambda wildcards: \
        [
            f"batches/{batch}/{sample}/combined_vcf/{sample}.{ref}.deepvariant_hifiasm.vcf.gz.tbi"
            for sample in _get_demuxed_samples( wildcards )
       ]
)
