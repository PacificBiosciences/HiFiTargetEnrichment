localrules: 
    merge_read_metrics,
    consolidate_read_data,
    consolidate_coverage_data,
    count_issue_element,
    get_coverage_fraction,


rule count_ontarget:
    input:
        bed=lambda wildcards: config['beds'][wildcards.bedfile],
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
    output:
        by_elem=(f'batches/{batch}/{{sample}}/counts/readcount_{{bedfile}}.bed'),
        by_read=(f'batches/{batch}/{{sample}}/counts/{{bedfile}}_per_read.csv'),
    params:
        elem='{bedfile}',
    conda:
        "envs/samtools.yaml"
    message:
        "Counting intersections of {wildcards.bedfile} with {input.bam}"
    shell:
        '''
        #by element count of reads
        bedtools intersect -a {input.bed} -b <(samtools view -hbF 0x900 {input.bam}) -c > {output.by_elem}
        #by read count of elements
        bedtools intersect -a <(samtools view -hbF 0x900 {input.bam}) -b {input.bed} -bed -c \
        | awk -v elem={params.elem} \
              'BEGIN {{ OFS="," ; print "readname,chr,start,stop,"elem }} \
               {{ print $4,$1,$2,$3,$NF }}' \
        | ( sed -u 1q; sort ) > {output.by_read} 
        '''

rule get_coverage:
    input:
        bed=lambda wildcards: config['beds'][wildcards.bedfile],
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
    output:
        ontarget=(f'batches/{batch}/{{sample}}/coverage/{{bedfile}}_ontarget.bam'),
        bg=f'batches/{batch}/{{sample}}/coverage/{{bedfile}}_coverage.bedgraph',
        cov=f'batches/{batch}/{{sample}}/coverage/{{bedfile}}_base_coverage.csv',
    params:
        sample='{sample}'
    conda:
        "envs/samtools.yaml"
    message:
        "Getting coverage of {input.bed} from {input.bam}"
    shell:
        '''
        bedtools intersect -a {input.bam} -b {input.bed} > {output.ontarget}
        bedtools genomecov -bga -ibam {output.ontarget} > {output.bg}
        bedtools intersect -a {output.bg} -b {input.bed} -wb \
        | awk -v s={params.sample} \
              'BEGIN {{ OFS = ","; print "sample,chr,start,stop,coverage,target" }} \
              {{ print s,$1,$2,$3,$4,$8 }}' \
        > {output.cov} 
        '''

rule get_coverage_fraction:
    input:
        f'batches/{batch}/{{sample}}/coverage/{{bedfile}}_base_coverage.csv',
    output:
        f'batches/{batch}/{{sample}}/coverage/{{bedfile}}_base_coverage_fraction.csv',
    params:
        sample='{sample}',
        lowcov=config["QC"]["lowcov"],
    shell:
        '''
        awk -F, -v s={params.sample} -v low={params.lowcov} \
            'BEGIN {{ OFS=","; print "sample,target,totalBp,lowcovBp,fracLow" }} 
             NR>1 {{ 
                    span = $4 - $3; 
                    split( $6, t, "_" ); target = t[1];
                    total += span; 
                    subset[target] += span; 
                    if( $5 <= low ) {{ 
                                        lowcov += span; 
                                        lowsubset[target] += span; 
                                      }} 
                   }} 
             END {{ 
                    for ( tg in subset ) {{ print s, tg, subset[tg], lowsubset[tg], lowsubset[tg]/subset[tg]; }} 
                    print s,"all",total,lowcov,lowcov/total; 
                 }}' {input} > {output}
        '''

rule summarize_coverage_fraction:
    input:
        expand( f'batches/{batch}/' + '{sample}/coverage/{{bedfile}}_base_coverage_fraction.csv', sample=sample2barcode.keys() )
    output:
        f'batches/{batch}/stats/{{bedfile}}_covered_fraction.csv',
    shell:
        '''
        tail -qn1 {input} \
        | awk -F, \
            'BEGIN {{ OFS=","; print "sample,target,totalBp,lowcovBp,fracLow" }} \
            {{ \
                total += $3; \
                lowcov += $4; \
                print; \
            }} \
            END {{ \
                    print "Total","allTargets",total,lowcov,lowcov/total;
                }}' \
        > {output}
        '''

rule get_read_metrics:
    input:
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
    output:
        stats=(f'batches/{batch}/{{sample}}/read_metrics/reads.csv'),
    conda:
        "envs/samtools.yaml",
    message:
        "Getting samtools target read metrics from {input.bam}"
    shell:
        '''
        samtools view -F 0x900 {input.bam} \
        | awk 'BEGIN {{ OFS="," ; print "readname,length,qual,ismapped,duplicates,softclip" }} \
               {{ if( match( $0, /ds:i:[0-9]/ ) ) {{ split( substr( $0, RSTART, RLENGTH ),d,":"); dup=d[3] }} \
                  else                            {{ dup=0 }} \
               }} ; \
               {{ soft = match( $6, /[0-9]+S/ ) ? \
                         substr( $6, RSTART, RLENGTH-1 ) : \
                         0 \
               }} ; \
               {{ mapped = !and($2,4) }} ; \
               {{ match( $0, /rq:f:[0-9.]+/ ); \
                  split( substr( $0, RSTART, RLENGTH ),q,":") ; \
                  print $1, length($10), q[3], mapped, dup, soft \
               }}' \
        | ( sed -u 1q; sort) > {output.stats}
        '''

rule get_read_target_metrics:
    input:
        bed=config['beds']['targets'],
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
    output:
        target=(f'batches/{batch}/{{sample}}/read_metrics/read_target.csv'),
    conda:
        "envs/samtools.yaml",
    message:
        "Getting target read metrics from {input.bam}"
    shell:
        ''' 
        bedtools intersect -a <(samtools view -hbF 0x900 {input.bam}) -b {input.bed} -bed -loj \
        | awk 'BEGIN {{ print "readname,target" }} \
               {{ print $4","$NF }}' \
        | ( sed -u 1q; sort -t, -u -k1,1 ) > {output.target}
        ''' 

rule get_duplicate_lengths:
    input: 
        lambda wildcards: f"batches/{batch}/{wildcards.sample}/dedup/duplicates.{sample2barcode[wildcards.sample]}.bam",
    output:
        f'batches/{batch}/{{sample}}/read_metrics/duplicate_lengths.csv',
    params:
        sample='{sample}'
    message:
        "Getting duplicate readlengths for {input}"
    conda:
        "envs/samtools.yaml",
    shell:
        '''
        samtools view {input} \
        | awk -v s={params.sample} '{{print s","length($10)",Duplicates"}}' \
        > {output}
        '''

rule merge_duplengths:
    input:
        expand( f'batches/{batch}/{{sample}}/read_metrics/duplicate_lengths.csv', sample=sample2barcode.keys() ),
    output:
        f'batches/{batch}/stats/all_duplicate_lengths.csv'
    shell:
        ''' 
        awk ' BEGIN {{ print "sample,length,category" }}
              {{ print }}' {input} > {output}
        '''

rule endmost_probe:
    input:
        bed=config['beds']['probes'],
        bam=f"batches/{batch}/{{sample}}/aligned/{{sample}}.{ref}.bam",
    output:
        (f'batches/{batch}/{{sample}}/read_metrics/endmost_probe.csv')
    conda:
        "envs/samtools.yaml"
    message:
        "Finding endmost probe for {input.bam}"
    shell:
        '''
        bedtools intersect -a <(samtools view -hbF 0x900 {input.bam}) -b {input.bed} -bed \
        | awk 'function abs(v) {{return v < 0 ? -v : v}} \
               BEGIN {{ OFS=","; print "readname,probe2end,probeStart,probeStop" }} \
              {{ m=1e9; \
                 for (r=2;r<4;r++) {{ \
                                     for (p=7;p<9;p++) {{ \
                                                         if ( abs($r-$p)<m ) {{ m=abs($r-$p) }} \
                                                       }} \
                                   }}; \
                print $4, m, $7, $8 \
              }}' \
        | ( sed -u 1q; sort -t, -k1,1 -k2,2n | sort -t, -u -k1,1 ) \
        > {output}
        '''

rule merge_read_metrics:
    input:
        stats=f'batches/{batch}/{{sample}}/read_metrics/reads.csv',
        target=f'batches/{batch}/{{sample}}/read_metrics/read_target.csv',
        probes=f'batches/{batch}/{{sample}}/counts/probes_per_read.csv', 
        endprb=f'batches/{batch}/{{sample}}/read_metrics/endmost_probe.csv',
        exons=f'batches/{batch}/{{sample}}/counts/exons_per_read.csv',
    output:
        f'batches/{batch}/{{sample}}/read_metrics/read_data.csv',
    params:
        sample='{sample}'
    shell:
        '''
        join -t, -e. -a1 -o auto \
             --header \
            {input.stats} \
            <( \
              join -t, \
                   --header \
                  <( \
                    join -t, \
                         --header \
                        <( \
                          join -t, \
                               --header \
                               -o 1.1,1.2,1.3,1.4,1.5,2.5 \
                               {input.probes} \
                               {input.exons} \
                         ) \
                        {input.target} \
                   ) \
                  {input.endprb} \
              ) \
        | sed -e '1s/^/sample,/' -e '2,$s/^/{params.sample},/' \
        > {output}
        '''

rule consolidate_read_data:
    input: 
        expand(
                f'batches/{batch}/{{sample}}/read_metrics/read_data.csv',sample=sample2barcode.keys()
              ),
    output:
        f'batches/{batch}/stats/all_read_data.csv',
    shell:
        '''
        awk 'NR==1 || FNR>1' {input} > {output}
        '''

rule consolidate_coverage:
    input:
        expand(
                f'batches/{batch}/' + '{sample}/coverage/{{bedfile}}_base_coverage.csv',sample=sample2barcode.keys()
              ),
    output:
        f'batches/{batch}/stats/all_{{bedfile}}_coverage.csv',
    shell:
        '''
        awk 'NR==1 || FNR>1' {input} > {output}
        '''

rule count_issue_elements:
    input:
        expand(
                f'batches/{batch}/' + '{sample}/counts/readcount_{{elem}}.bed', sample=sample2barcode.keys()
              ),
    output:
        f'batches/{batch}/stats/{{status}}/{{elem}}.csv',
    params:
        low= lambda wildcards: 0 if wildcards.status == 'dropped' else config['QC']['lowcov']
    shell:
        '''
        awk -v min={params.low} '$5 <= min' {input} \
        | sort -k1,1 -k2,2n \
        | uniq -c \
        | awk 'BEGIN {{ OFS=","; print "chr,start,stop,target,samplesDropped" }} \
                     {{ print $2,$3,$4,$5,$1 }}' \
        > {output}
        '''

rule gc_content:
    input:
        bed=lambda wildcards: config['beds'][wildcards.bedfile],
        ref=config["ref"]["fasta"],
    output:
        f'batches/{batch}/stats/{{bedfile}}_frac_gc.csv',
    conda:
        "envs/samtools.yaml"
    shell:
        '''
        bedtools nuc -fi {input.ref} -bed {input.bed} \
        | awk 'BEGIN {{ OFS=","; print "chr,start,stop,target,frac_gc" }} \
               NR>1  {{ print $1,$2,$3,$4,$(NF-7) }}' \
        > {output}
        '''

rule plot_read_metrics:
    input:
        csv=f'batches/{batch}/{{sample}}/read_metrics/read_data.csv',
        bed=config['beds']['targets'],
    output:
        [f'batches/{batch}/{{sample}}/read_metrics/{p}' 
          for p in 
            ['mean_base_coverage_by_target.png',
             'readlength_hist.png',
             'readlength_hist_by_target.png',
             'readlength_violins_by_target.png',
             'dedup_count_by_target.png',
             'dedup_rate_by_target.png',
             'readqual_by_target.png',
             'off-target_length_hist.png'
            ]
        ],
    params:
        script=f'{config["scripts"]}/plot_read_data.py',
        odir=f'batches/{batch}/{{sample}}/read_metrics',
    shell:
        '''
        {params.script} {input.csv} {input.bed} {params.odir}        
        '''

rule plot_coverage:
    input:
        f'batches/{batch}/{{sample}}/coverage/targets_base_coverage.csv',
    output:
        f'batches/{batch}/{{sample}}/coverage/coverage_by_target.png',
    params:
        script=f'{config["scripts"]}/plot_coverage.py',
        odir=f'batches/{batch}/{{sample}}/coverage',
    shell:
        '''
        {params.script} {input} {params.odir}
        '''
rule plot_multi_coverage:
    input:
        f'batches/{batch}/stats/all_targets_coverage.csv',
    output:
        f'batches/{batch}/stats/multi_coverage_by_target.png',
    params:
        script=f'{config["scripts"]}/plot_multi_coverage.py',
        odir=f'batches/{batch}/stats',
    shell:
        '''
        {params.script} {input} {params.odir}
        '''

rule plot_multi_read:
    input:
        csv=f'batches/{batch}/stats/all_read_data.csv',
        bed=config['beds']['targets'],
    output:
        [
          f'batches/{batch}/stats/{p}' 
          for p in
            ['mean_base_coverage.png',
             'mean_base_coverage_by_sample.png',
             'dedup_length_by_target.png',
             'dedup_rate_by_target.png',
             'readlength_hist_by_target.png'
            ]
        ],
    params:
        script=f'{config["scripts"]}/plot_multi_reads.py',
        odir=f'batches/{batch}/stats',
    shell:
        '''
        {params.script} {input.csv} {input.bed} {params.odir}
        '''

rule plot_read_categories:
    input:
        readCsv=f'batches/{batch}/stats/all_read_data.csv',
        lima=f'batches/{batch}/demux/{{movie}}/demultiplex.lima.report',
        dups=f'batches/{batch}/stats/all_duplicate_lengths.csv',
    output:
        f'batches/{batch}/stats/read_categories_{{movie}}.png',
    params:
        script=f'{config["scripts"]}/plot_read_cats.py',
    shell:
        '''
        {params.script} {input.readCsv} {input.lima} {input.dups} {output}
        '''
    

# extand targets
targets.extend(
    [
        f'batches/{batch}/{sample}/counts/readcount_{bedfile}.bed'
        for sample in sample2barcode.keys()
        for bedfile in config['beds'].keys()
    ]
)

targets.extend(
    [
        f'batches/{batch}/{sample}/counts/{bedfile}_per_read.csv'
        for sample in sample2barcode.keys()
        for bedfile in config['beds'].keys()
    ]
)

targets.extend(
    [
        f'batches/{batch}/{sample}/read_metrics/reads.csv'
        for sample in sample2barcode.keys()
    ]
)

targets.extend(
    [
        f'batches/{batch}/{sample}/read_metrics/read_data.csv'
        for sample in sample2barcode.keys()
    ]
)

targets.extend(
    [
        f'batches/{batch}/stats/{bedfile}_frac_gc.csv'
        for bedfile in config['beds'].keys()
    ]
)

targets.extend(
    [
        f'batches/{batch}/stats/all_{bedfile}_coverage.csv'
        for bedfile in ['targets','exons']
    ]
)

targets.extend(
    [
        f'batches/{batch}/{sample}/read_metrics/mean_base_coverage_by_target.png'
        for sample in sample2barcode.keys()
    ]
)

targets.extend(
    [
        f'batches/{batch}/{sample}/coverage/coverage_by_target.png'
        for sample in sample2barcode.keys()
    ]
)

targets.extend(
    [
        f'batches/{batch}/stats/{status}/{elem}.csv'
        for status in ['dropped','lowcov']
        for elem in config['beds'].keys()
    ]
)

targets.extend(
    [
        f'batches/{batch}/{sample}/coverage/{bedfile}_base_coverage_fraction.csv'
        for sample in sample2barcode.keys()
        for bedfile in ['targets','exons']
    ]
)

targets.extend(
    [
        f'batches/{batch}/stats/{bedfile}_covered_fraction.csv'
        for bedfile in ['targets','exons']
    ]
)

targets.extend(
    [
        f'batches/{batch}/stats/all_read_data.csv',
        f'batches/{batch}/stats/mean_base_coverage.png',
        f'batches/{batch}/stats/multi_coverage_by_target.png',
        f'batches/{batch}/stats/all_duplicate_lengths.csv'
    ]
)

targets.extend(
    [
        f'batches/{batch}/{sample}/read_metrics/duplicate_lengths.csv'
        for sample in sample2barcode.keys()
    ]
)

targets.extend(
    [
        f'batches/{batch}/stats/read_categories_{movie}.png'
        for movie in movies
    ]
)
