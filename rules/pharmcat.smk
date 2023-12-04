rule pharmcat_preprocess_fill_missing:
    input:
        vcf=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.vcf.gz",
        vcf_index=f"batches/{batch}/{{sample}}/whatshap/{{sample}}.{ref}.deepvariant.phased.vcf.gz.tbi",
        reference=config["ref"]["fasta"],
    output:
        temp(f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.preprocessed.vcf.bgz"),
    log:
        f"batches/{batch}/logs/pharmcat/preprocess_vcf/{{sample}}.{ref}.log"
    container:
        "docker://pgkb/pharmcat:2.8.2"
    params:
        odir=f"batches/{batch}/{{sample}}/pharmcat/",
        regions=config["pharmcat"]["positions"],
        basefile="{sample}",
    message:
        "pharmcat: preprocess vcf for {wildcards.sample}"
    shell:
        '''
        (/pharmcat/pharmcat_vcf_preprocessor.py \
            --missing-to-ref \
            -vcf {input.vcf} \
            -refFna {input.reference} \
            -refVcf {params.regions} \
            -bf {params.basefile} \
            -o {params.odir} ) > {log} 2>&1
        '''

rule pharmcat_remove_positions_with_no_coverage:
    input:
        vcf=f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.preprocessed.vcf.bgz",
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        bamindex=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
        genome=config["ref"]["chr_lengths"],
    output:
        vcf=f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.preprocessed.filtered.vcf",
    log:
        f"batches/{batch}/logs/pharmcat/preprocess_vcf/{{sample}}.remove_nocov.log"
    params:
        mincov=config["pharmcat"]["mincov"],
    conda:
        "envs/samtools.yaml"
    message:
        "pharmcat: remove ref-calls with low mean coverage for {wildcards.sample}"
    shell:
        '''
        (bedtools coverage \
                  -sorted \
                  -g {input.genome} \
                  -f 1 \
                  -header \
                  -mean \
                  -a {input.vcf} \
                  -b {input.bam} \
         | ( sed  -u '/^#CHROM/q' ; awk '$11 >= {params.mincov}' | cut -f1-10 ) > {output}) > {log} 2>&1
        '''

rule run_pharmcat:
    input:
        vcf=f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.preprocessed.filtered.vcf",
        cyp2d6=f"batches/{batch}/{{sample}}/pangu/{{sample}}_pharmcat_fix.tsv",
    output:
        [ f"batches/{batch}/{{sample}}/pharmcat/{{sample}}.preprocessed.filtered.{ftype}"
            for ftype in [ 'match.json', 'phenotype.json', 'report.json', 'report.html' ] ],
    params:
        odir=f"batches/{batch}/{{sample}}/pharmcat/",
    log:
        f"batches/{batch}/logs/pharmcat/run_pharmcat/{{sample}}.log"   
    container:
        "docker://pgkb/pharmcat:2.3.0"
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
            f"batches/{batch}/{sample}/pharmcat/{sample}.preprocessed.filtered.report.json"
            for sample in _get_demuxed_samples( wildcards )
        ]
)
