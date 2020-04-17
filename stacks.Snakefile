#!/usr/bin/env python3

import multiprocessing
import pandas
import tempfile


#############
# FUNCTIONS #
#############


###########
# GLOBALS #
###########

filtered_popmap = 'output/000_config/filtered_population_map.txt'
ref = 'data/draft_genome.fasta'

# containers
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
samtools = 'shub://TomHarrop/TomHarrop/align-utils:samtools_1.9'
stacks_container = 'shub://TomHarrop/variant-utils:stacks_2.53'

#########
# SETUP #
#########

all_cpus = multiprocessing.cpu_count()

# read popmap
try:
    popmap = pandas.read_csv(filtered_popmap,
                             delimiter='\t',
                             header=None,
                             names=['individual', 'population'])
    all_indivs = sorted(set(popmap['individual']))
except FileNotFoundError:
    print(f'ERROR. {filtered_popmap} not found.')
    print('       Run process_reads.Snakefile')
    raise FileNotFoundError

#########
# RULES #
#########

subworkflow process_reads:
    snakefile: 'process_reads.Snakefile'

rule target:
    input:
        'output/050_stacks/populations/populations.snps.vcf'

# 8 run populations once to generate stats for filtering
rule populations:
    input:
        'output/050_stacks/catalog.fa.gz',
        'output/050_stacks/catalog.calls',
        map = filtered_popmap
    output:
        'output/050_stacks/populations/populations.snps.vcf'
    params:
        stacks_dir = 'output/050_stacks',
        outdir = 'output/050_stacks/populations'
    log:
        'output/logs/populations.log'
    singularity:
        stacks_container
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-O {params.outdir} '
        '-r 0 '
        '--genepop --vcf '
        '&> {log}'

# 7 run gstacks
rule gstacks:
    input:
        bam_files = expand(
            'output/050_stacks/bamfiles/{individual}.bam',
            individual=all_indivs),
        popmap = filtered_popmap
    output:
        'output/050_stacks/catalog.fa.gz',
        'output/050_stacks/catalog.calls'
    params:
        stacks_i = 'output/050_stacks/bamfiles',
        stacks_o = 'output/050_stacks'
    threads:
        all_cpus
    log:
        'output/logs/gstacks.log'
    singularity:
        stacks_container
    shell:
        'gstacks '
        '-I {params.stacks_i} '
        '-O {params.stacks_o} '
        '-M {input.popmap} '
        '-t {threads} '
        '--details '
        '&> {log}'

# 6 map the reads to the draft genome
rule markdup:
    input:
        'output/tmp/{individual}.sorted.bam'
    output:
        'output/050_stacks/bamfiles/{individual}.bam'
    log:
        'output/logs/markdup.{individual}.log'
    singularity:
        samtools
    shell:
        'samtools markdup '
        '-@ {all_cpus} '
        '-s '
        '{input} '
        '{output} '
        '2> {log}'

rule sort_sam:
    input:
        'output/tmp/{individual}.sam'
    output:
        pipe('output/tmp/{individual}.sorted.bam')
    log:
        'output/logs/sort_sam.{individual}.log'
    singularity:
        samtools
    shell:
        'samtools sort '
        '-O BAM '
        '--threads {all_cpus} '
        '{input} '
        '> {output} '
        '2> {log}'

rule map_indiv_reads:
    input:
        fq = process_reads('output/030_combined/{individual}.fq.gz'),
        ref = 'output/005_ref/ref.fasta'
    output:
        pipe('output/tmp/{individual}.sam')
    log:
        'output/logs/map_indiv_reads.{individual}.log'
    threads:
        all_cpus - 2
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t {threads} '
        '{input.ref} '
        '{input.fq} '
        '>> {output} '
        '2> {log}'

rule bwa_index:
    input:
        ref = ref
    output:
        expand('output/005_ref/ref.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
        ref = 'output/005_ref/ref.fasta'
    threads:
        1
    log:
        'output/logs/bwa_index.log'
    singularity:
        bwa_container
    shell:
        'cp {input.ref} {output.ref} ; '
        'bwa index '
        '-p {output.ref} '
        '{output.ref} '
        '2> {log}'

# generic rule for indexing vcf
rule sort_vcf:
    input:
        'output/{folder}/{file}.vcf'
    output:
        pipe('output/{folder}/{file}_sorted.vcf')
    singularity:
        samtools
    shell:
        'bcftools sort '
        '--temp-dir ' + tempfile.mkdtemp() + ' '
        '{input} '
        '>> {output} '
        '2> {log}'

rule index_vcf:
    input:
        'output/{folder}/{file}_sorted.vcf'
    output:
        gz = 'output/{folder}/{file}.vcf.gz',
        tbi = 'output/{folder}/{file}.vcf.gz.tbi'
    log:
        'output/logs/{folder}/{file}_index-vcf.log'
    wildcard_constraints:
        file = 'populations.snps'
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output.gz} 2> {log} '
        '; '
        'tabix -p vcf {output.gz} 2>> {log}'
