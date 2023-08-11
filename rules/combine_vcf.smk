#ruleorder: filter_to_type > tabix_vcf
#ruleorder: concat_combined_analyses > tabix_vcf 
#ruleorder: annotate_vcf_by_method > tabix_vcf


def _get_vcf_by_type_and_source( index=False ):
    '''Default is to take VCF data from DV'''
    ext = '.tbi' if index else ''
    def getter( wildcards ):
        return {
            ("asm",True): f"batches/{batch}/{wildcards.sample}/hifiasm/{wildcards.sample}.asm.{ref}.htsbox.vcf.gz{ext}",
            ("asm",False): f"batches/{batch}/{wildcards.sample}/hifiasm/{wildcards.sample}.asm.{ref}.htsbox.vcf.gz{ext}",
            ("hifi",True): f"batches/{batch}/{wildcards.sample}/whatshap/{wildcards.sample}.{ref}.deepvariant.phased.vcf.gz{ext}",
            ("hifi",False): f"batches/{batch}/{wildcards.sample}/sniffles/{wildcards.sample}.{ref}.{wildcards.source}.sv.vcf.gz{ext}"
        }[ (wildcards.source, "SNV" in wildcards.vartype) ]
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
        bed=(f"batches/{batch}/targets_buffered.bed.gz"),
        idx=(f"batches/{batch}/targets_buffered.bed.gz.gzi"),
#        bed=temp(f"batches/{batch}/targets_buffered.bed.gz"),
#        idx=temp(f"batches/{batch}/targets_buffered.bed.gz.gzi"),
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

#def is_valid_sample(sample):
#    return '/' not in sample 
    
#wildcard_constraints:
 #   tool    = ["deepvariant_hifiasm","sniffles"],
   # vartype = ["SNV","SV"],
  #  source  = ["asm","hifi"],
#    sample  = is_valid_sample,

rule annotate_vcf_by_method:
    input:
        vcf=_get_vcf_by_type_and_source(),
        idx=_get_vcf_by_type_and_source(index=True),
        bed=f"batches/{batch}/targets_buffered.bed.gz",
        bidx=f"batches/{batch}/targets_buffered.bed.gz.gzi",
    output:
        vcf=(f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{{vartype}}.{{source}}.vcf.gz"),
        idx=(f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{{vartype}}.{{source}}.vcf.gz.tbi"),
    params:
        fields=lambda wcs: 'CHROM,FROM,TO,INFO/GENE,INFO/V,-' if 'SNV' in wcs.vartype else 'CHROM,FROM,TO,INFO/GENE,-,INFO/V',
        header=r''' '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Name">\n##INFO=<ID=V,Number=1,Type=Integer,Description="V=1 Assembly V=2 HiFi_Reads">' ''',
        vs=lambda wcs: {"asm":"1","hifi":"2"}[wcs.source],
    log:
        f"batches/{batch}/logs/combine_vcf/bcftools/{{sample}}.{{vartype}}.annotate_{{source}}.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "Annotating {wildcards.vartype} from {wildcards.source} variant targets for {wildcards.sample}"
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
                  -o {output.vcf} \
         && bcftools index -t {output.vcf}) > {log} 2>&1
        '''
    
rule concat_combined_analyses:
    input:
        vcf=[f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{{vartype}}.hifi.vcf.gz",
             f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{{vartype}}.asm.vcf.gz"],
        idx=[f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{{vartype}}.hifi.vcf.gz.tbi",
             f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{{vartype}}.asm.vcf.gz.tbi"],
    output:
        vcf=temp(f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{ref}.{{vartype}}.{{tool}}.vcf.gz"),
        idx=temp(f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{ref}.{{vartype}}.{{tool}}.vcf.gz.tbi"),
    log:
        f"batches/{batch}/logs/combine_vcf/bcftools/{{sample}}.{{vartype}}.{{tool}}.concat.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "Merging {wildcards.vartype} variants by different methods for {wildcards.sample}"
    shell:
        '''
        (bcftools concat -a {input.vcf} -Oz -o {output.vcf} && bcftools index -t {output.vcf}) > {log} 2>&1
        '''

rule filter_to_type:
    input:
        vcf=f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{ref}.{{vartype}}.{{tool}}.vcf.gz",
        idx=f"batches/{batch}/{{sample}}/combined_vcf/intermed/{{sample}}.{ref}.{{vartype}}.{{tool}}.vcf.gz.tbi",
    output:
        vcf=f"batches/{batch}/{{sample}}/combined_vcf/{{vartype}}/{{sample}}.{ref}.{{tool}}.filtered.vcf.gz",
        idx=f"batches/{batch}/{{sample}}/combined_vcf/{{vartype}}/{{sample}}.{ref}.{{tool}}.filtered.vcf.gz.tbi",
    params:
        exclusion=lambda wcs: \
            f'{"-e" if "SNV" in wcs["vartype"] else "-i"} \'strlen(ALT) >= 50 || strlen(REF) >= 50\'',
    log:
        f"batches/{batch}/logs/combine_vcf/bcftools/{{sample}}.{{vartype}}.{{tool}}.filter.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "Filtering to {wildcards.vartype} for {wildcards.sample}"
    shell:
        '''
        (bcftools view \
                  {params.exclusion} \
                  {input.vcf} \
                  -Oz \
                  -o {output.vcf} \
        && bcftools index -t {output.vcf}) > {log} 2>&1
        '''
    
targets.append(
    lambda wildcards: \
        [
            f"batches/{batch}/{sample}/combined_vcf/{vartype}/{sample}.{ref}.{tool}.filtered.vcf.gz.tbi"
            for sample in _get_demuxed_samples( wildcards )
            for vartype,tool in zip(["SNV","SV"],["deepvariant_hifiasm","sniffles_hifiasm"])
       ]
)
