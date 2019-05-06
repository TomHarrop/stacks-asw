import pandas
import pickle
import multiprocessing

#############
# FUNCTIONS #
#############

def lookup_indiv(pickle_file, individual):
    with open(pickle_file, 'rb') as f:
        individual_i = pickle.load(f)
        sample_i = individual_i[individual]
        return(sample_i)

###########
# GLOBALS #
###########

reads_dir = 'data/raw_reads'
filtered_popmap = 'output/stacks_config/filtered_population_map.txt'

# containers
stacks_container = ('shub://TomHarrop/singularity-containers:stacks_2.3e')
stacks2beta_container = ('shub://TomHarrop/singularity-containers:stacks_2.0beta9'
                         '@bb2f9183318871f6228b51104056a2d0')
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'

#########
# SETUP #
#########

# read popmap
popmap = pandas.read_csv(filtered_popmap,
                         delimiter='\t',
                         header=None,
                         names=['individual', 'population'])
all_indivs = sorted(set(popmap['individual']))

#########
# RULES #
#########

rule target:
    input:
        'output/stacks_populations/r0/populations.snps.vcf',
        'output/combined_stats/individual_covstats_combined.csv'

subworkflow process_reads:
    snakefile: 'process_reads.Snakefile'


# # 12. filter the final catalog by r
rule populations:
    input:
        'output/map_to_genome/catalog.fa.gz',
        'output/map_to_genome/catalog.calls',
        map = filtered_popmap
    output:
        'output/stacks_populations/r0/populations.snps.vcf'
    params:
        stacks_dir = 'output/map_to_genome',
        outdir = 'output/stacks_populations/r0'
    threads:
        75
    log:
        'output/logs/populations_r0.log'
    singularity:
        stacks_container
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-O {params.outdir} '
        '-t {threads} '
        '-r 0 '
        '--genepop --vcf '
        '&> {log}'

# 11b map the catalog to the draft genome
rule integrate_alignments:
    input:
        fa = 'output/stacks_denovo/catalog.fa.gz',
        bam = 'output/map_to_genome/aln.bam',
    output:
        catalog = 'output/map_to_genome/catalog.fa.gz',
        tsv = 'output/map_to_genome/locus_coordinates.tsv',
        calls = 'output/map_to_genome/catalog.calls'
    params:
        stacks_dir = 'output/stacks_denovo',
        out_dir = 'output/map_to_genome'
    threads:
        1
    log:
        'output/logs/integrate_alignments.log'
    singularity:
        stacks_container
    shell:
        'stacks-integrate-alignments '
        '-P {params.stacks_dir} '
        '-B {input.bam} '
        '-O {params.out_dir} '
        '&> {log}' 

rule sam_to_bam:
    input:
        aln = 'output/map_to_genome/aln.sam'
    output:
        bam = 'output/map_to_genome/aln.bam',
    threads:
        1
    log:
        'output/logs/sam_to_bam.log'
    singularity:
        stacks_container
    shell:
        'samtools view '
        '--threads {threads} '
        '-O BAM '
        '-bh '
        '{input.aln} '
        '> {output.bam} '
        '2> {log}'


rule map_stacks_catalog:
    input:
        fa = 'output/stacks_denovo/catalog.fa.gz',
        index = expand('output/map_to_genome/draft_genome.fasta.{suffix}',
                       suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])

    output:
        temp('output/map_to_genome/aln.sam')
    params:
        prefix = 'output/map_to_genome/draft_genome.fasta',
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/map_stacks_catalog.log'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t {threads} '
        '{params.prefix} '
        '{input.fa} '
        '> {output} '
        '2> {log}'


rule bwa_index:
    input:
        'data/draft_genome.fasta'
    output:
        expand('output/map_to_genome/draft_genome.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix = 'output/map_to_genome/draft_genome.fasta'
    threads:
        1
    log:
        'output/logs/bwa_index.log'
    singularity:
        bwa_container
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input} '
        '2> {log}'


# # 11. generate final catalog
rule gstacks:
    input:
        expand('output/stacks_denovo/{individual}.matches.bam',
               individual=all_indivs),
        map = filtered_popmap
    output:
        'output/stacks_denovo/catalog.fa.gz',
        'output/stacks_denovo/catalog.calls'
    params:
        stacks_dir = 'output/stacks_denovo'
    threads:
        75
    log:
        'output/logs/gstacks.log'
    singularity:
        stacks_container
    shell:
        'gstacks '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-t {threads} '
        '&> {log}'

# # 10. convert matches to BAM
rule tsv2bam:
    input:
        expand('output/stacks_denovo/{individual}.matches.tsv.gz',
               individual=all_indivs),
        map = filtered_popmap
    output:
        expand('output/stacks_denovo/{individual}.matches.bam',
               individual=all_indivs)
    params:
        stacks_dir = 'output/stacks_denovo'
    threads:
        75
    log:
        'output/logs/tsv2bam.log'
    singularity:
            stacks_container
    shell:
        'tsv2bam '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-t {threads} '
        '&> {log}'

# # 9. match individuals to the catalog
rule sstacks:
    input:
        'output/stacks_denovo/catalog.tags.tsv.gz',
        map = filtered_popmap
    output:
        expand('output/stacks_denovo/{individual}.matches.tsv.gz',
               individual=all_indivs)
    params:
        stacks_dir = 'output/stacks_denovo'
    threads:
        75
    log:
        'output/logs/sstacks.log'
    singularity:
            stacks_container
    shell:
        'sstacks '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-p {threads} '
        '&> {log}'

# # 8. generate the catalog (takes ~ 1 week)
rule cstacks:
    input:
        expand('output/stacks_denovo/{individual}.alleles.tsv.gz',
               individual=all_indivs),
        expand('output/stacks_denovo/{individual}.snps.tsv.gz',
               individual=all_indivs),
        expand('output/stacks_denovo/{individual}.tags.tsv.gz',
               individual=all_indivs),
        map = filtered_popmap
    output:
        'output/stacks_denovo/catalog.tags.tsv.gz',
        'output/stacks_denovo/catalog.snps.tsv.gz',
        'output/stacks_denovo/catalog.alleles.tsv.gz'
    params:
        stacks_dir = 'output/stacks_denovo',
        n = '4'
    threads:
        75
    log:
        'output/logs/cstacks.log'
    singularity:
        stacks_container
    shell:
        'cstacks '
        '-p {threads} '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-n {params.n} '
        '&> {log}'

# 7d. combine individual coverage stats
rule individual_covstats_combined:
    input:
        expand('output/individual_covstats/{individual}.csv',
               individual=all_indivs)
    output:
        combined = 'output/combined_stats/individual_covstats_combined.csv'
    priority:
        1
    singularity:
        r_container
    script:
        'src/combine_csvs.R'

# # 7c. calculate coverage stats per individual
rule individual_covstats:
    input:
        tags_file = 'output/stacks_denovo/{individual}.tags.tsv.gz'
    output:
        covstats = 'output/individual_covstats/{individual}.csv'
    log:
        log = 'output/logs/individual_covstats/{individual}.log'
    threads:
        1
    singularity:
        r_container
    script:
        'src/calculate_mean_coverage.R'

# 7. assemble loci for individuals that passed the filter
rule ustacks:
    input:
        fq = process_reads('output/combined/{individual}.fq.gz'),
        pickle = process_reads('output/stacks_config/individual_i.p')
    params:
        wd = 'output/stacks_denovo',
        i = lambda wildcards, input:
            lookup_indiv(input.pickle, wildcards.individual),
        m = '3',
        M = '3'
    output:
        'output/stacks_denovo/{individual}.alleles.tsv.gz',
        'output/stacks_denovo/{individual}.snps.tsv.gz',
        'output/stacks_denovo/{individual}.tags.tsv.gz'
    threads:
        1
    log:
        'output/logs/ustacks_{individual}.log'
    singularity:
        stacks_container
    shell:
        'ustacks '
        '-p {threads} '
        '-t gzfastq '
        '-f {input.fq} '
        '-o {params.wd} '
        '-i {params.i} '
        '-m {params.m} '
        '-M {params.M} '
        '&> {log}'
