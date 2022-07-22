rule copy_aligned_bams:
    input:
        lambda wcs: sample2reads[ wcs.sample ]
    output:
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
        bai=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam.bai",
    conda:
        "envs/samtools.yaml"
    shell:
        '''
        cp $( readlink -f {input} ) {output.bam}
        samtools index {output.bam}
        '''
