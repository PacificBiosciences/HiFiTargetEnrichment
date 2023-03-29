def _get_vcf_by_type( index=False ):
    '''Default is to take VCF data from DV'''
    ext = '.tbi' if index else ''
    def getter( wildcards ):
        if wildcards.source == "asm":
            return f"batches/{batch}/{wildcards.sample}/hifiasm/{wildcards.sample}.asm.{ref}.htsbox.vcf.gz{ext}"
        else:
            #return DV if the source is not specifically called out
            return f"batches/{batch}/{wildcards.sample}/whatshap/{wildcards.sample}.{ref}.deepvariant.phased.vcf.gz{ext}"
    return getter

    
rule buffer_target_bed:
    input:
        config['targets'],
    output:
        temp(f"batches/{batch}/targets_buffered.bed"),
    params:
        buffer=10000,
        chr_len=config["ref"]["chr_lengths"],
    log:
        f"batches/{batch}/logs/combine_vcf/bedtools/buffer_slop.log",
    conda:
        "envs/samtools.yaml"
    message:
        "Buffering target bed for merging analyses"
    shell:
        '''
        (bedtools slop \
                  -i {input} \
                  -b {params.buffer} \
                  -g {params.chr_len} > {output}) > {log} 2>&1  
        '''

rule compress_buffered_targets:
    input:
        f"batches/{batch}/targets_buffered.bed",
    output:
        bed=temp(f"batches/{batch}/targets_buffered.bed.gz"),
        idx=temp(f"batches/{batch}/targets_buffered.bed.gz.gzi"),
    log:
        f"batches/{batch}/logs/combine_vcf/htslib/compress_index_bed.log",
    conda:
        "envs/htslib.yaml"
    message:
        "Compressing target bed for analysis merging"
    shell:
        '''
        (bgzip -i {input}) > {log} 2>&1  
        '''

rule annotate_vcf_by_method:
    input:
        vcf=_get_vcf_by_type(),
        idx=_get_vcf_by_type(index=True),
        bed=f"batches/{batch}/targets_buffered.bed.gz",
        bidx=f"batches/{batch}/targets_buffered.bed.gz.gzi",
    output:
        temp(f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{{source}}.vcf.gz"),
    params:
        fields='CHROM,FROM,TO,INFO/GENE,INFO/V',
        header=r''' '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Name">\n##INFO=<ID=V,Number=1,Type=Integer,Description="V=1 Assembly V=2 Deepvariant">' ''',
        vs=lambda wcs: {"asm":"1","dv":"2"}[wcs.source],
    log:
        f"batches/{batch}/logs/combine_vcf/bcftools/{{sample}}.annotate_{{source}}.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "Annotating {wildcards.source} variant targets for {wildcards.sample}"
    shell:
        '''
        (bcftools annotate \
                  -a {input.bed} \
                  -c {params.fields} \
                  -h <(echo -e {params.header}) \
                  {input.vcf} | \
         bcftools view \
                  -i "V=={params.vs}" \
                  -Oz \
                  -o {output}) > {log} 2>&1
        '''
    
rule concat_combined_analyses:
    input:
        vcf=[f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.dv.vcf.gz",
             f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.asm.vcf.gz"],
        idx=[f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.dv.vcf.gz.tbi",
             f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.asm.vcf.gz.tbi"],
    output:
        f"batches/{batch}/{{sample}}/combined_vcf/{{sample}}.{ref}.deepvariant_hifiasm.vcf.gz",
    log:
        f"batches/{batch}/logs/combine_vcf/bcftools/{{sample}}.concat.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "Merging variants by different methods for {wildcards.sample}"
    shell:
        '''
        (bcftools concat -a {input.vcf} -Oz -o {output}) > {log} 2>&1
        '''

rule filter_to_snv:
    input:
        vcf=f"batches/{batch}/{{sample}}/combined_vcf/{{sample}}.{ref}.deepvariant_hifiasm.vcf.gz",
        idx=f"batches/{batch}/{{sample}}/combined_vcf/{{sample}}.{ref}.deepvariant_hifiasm.vcf.gz.tbi",
    output:
        f"batches/{batch}/{{sample}}/combined_vcf/{{sample}}.{ref}.deepvariant_hifiasm.filtered.vcf.gz",
    params:
        maxlen=50,
    log:
        f"batches/{batch}/logs/combine_vcf/bcftools/{{sample}}.filter.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "Filtering to SNV for {wildcards.sample}"
    shell:
        '''
        (bcftools view \
                  -e 'strlen(ALT) >= {params.maxlen} || strlen(REF) >= {params.maxlen}' \
                  {input.vcf} \
                  -Oz \
                  -o {output}) > {log} 2>&1
        '''
    
targets.append(
    lambda wildcards: \
        [
            f"batches/{batch}/{sample}/combined_vcf/{sample}.{ref}.deepvariant_hifiasm.filtered.vcf.gz.tbi"
            for sample in _get_demuxed_samples( wildcards )
       ]
)
