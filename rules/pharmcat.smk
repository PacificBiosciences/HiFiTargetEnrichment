rule expand_gvcf:
    input:
        gvcf=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.g.vcf.gz",
        gvcf_index=f"batches/{batch}/{{sample}}/deepvariant/{{sample}}.{ref}.deepvariant.g.vcf.gz.tbi",
        reference=config["ref"]["fasta"],
    output:
        temp(f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.{ref}.deepvariant.g.expanded.vcf.gz"),
    params:
        regions=config["pharmcat"]["positions"],
    log:
        f"batches/{batch}/logs/pharmcat/bcftools_convert/{{sample}}.{ref}.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "pharmcat: expanding gvcf for {wildcards.sample}"
    shell:
        '''
        (bcftools convert --gvcf2vcf \
         -f {input.reference} \
         -R {params.regions} \
         -Oz \
         -o {output} \
         {input.gvcf}) > {log} 2>&1
        '''

rule enforce_indels:
    input: 
        vcf=f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.{ref}.deepvariant.g.expanded.vcf.gz",
        vcf_index=f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.{ref}.deepvariant.g.expanded.vcf.gz.tbi",
    output:
        temp(f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.{ref}.deepvariant.g.merged.vcf.gz"),
    params:
        regions=config["pharmcat"]["positions"],
    log:
        f"batches/{batch}/logs/pharmcat/bcftools_merge_indels/{{sample}}.{ref}.log",
    conda:
        "envs/bcftools.yaml"
    message:
        "pharmcat: expanding gvcf for {wildcards.sample}"
    shell:
        '''
        (bcftools merge {params.regions} {input.vcf} |\
         bcftools view -s ^PharmCAT -Oz -o {output}) > {log} 2>&1
        '''

rule pharmcat_preprocess_vcf:
    input:
        vcf=f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.{ref}.deepvariant.g.merged.vcf.gz",
        vcf_index=f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.{ref}.deepvariant.g.merged.vcf.gz.tbi",
        reference=config["ref"]["fasta"],
    output:
        f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.preprocessed.vcf",
    log:
        f"batches/{batch}/logs/pharmcat/preprocess_vcf/{{sample}}.{ref}.log"
    container:
        "docker://pgkb/pharmcat"
    params:
        odir=f"batches/{batch}/{{sample}}/pharmcat/",
        regions=config["pharmcat"]["positions"],
    message:
        "pharmcat: preprocess vcf for {wildcards.sample}"
    shell:
        '''
        (/pharmcat/pharmcat_vcf_preprocessor.py \
            -vcf {input.vcf} \
            -refFna {input.reference} \
            -refVcf {params.regions} \
            -o {params.odir} ) > {log} 2>&1
        '''

rule run_pharmcat:
    input:
        vcf=f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.preprocessed.vcf",
        cyp2d6=f"batches/{batch}/{{sample}}/pangu/{{sample}}_pharmcat_fix.tsv",
    output:
        [ f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.preprocessed.{ftype}"
            for ftype in [ 'match.json', 'phenotype.json', 'report.json', 'report.html' ] ],
    params:
        odir=f"batches/{batch}/{{sample}}/pharmcat/",
    log:
        f"batches/{batch}/logs/pharmcat/run_pharmcat/{{sample}}.log"   
    container:
        "docker://pgkb/pharmcat"
    message:
        "pharmcat: running pharmcat for {wildcards.sample}"
    shell:
        '''
        (/pharmcat/pharmcat \
            -vcf {input.vcf} \
            -reporterJson \
            -po {input.cyp2d6} \
            -o {params.odir}) > {log} 2>&1
        '''

targets.append(
    lambda wildcards: \
        [
            f"batches/{batch}/{sample}/pharmcat/{sample}.preprocessed.report.json"
            for sample in _get_demuxed_samples( wildcards )
        ]
)

# Thanks to Binglan Li for help constructing these rules
