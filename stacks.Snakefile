import pandas
import pickle


#############
# FUNCTIONS #
#############

def read_sample_number(wildcards, input):
    with open(input.individual_i_pickle, 'rb') as f:
        individual_i = pickle.load(f)
        return(individual_i[wildcards.individual])


#########
# RULES #
#########

rule target:
    

subworkflow process_reads:
    snakefile: 'process_reads.Snakefile'


# 12. filter the final catalog by r
rule populations:
    input:
        'output/stacks_denovo/gstacks.fa.gz',
        'output/stacks_denovo/gstacks.vcf.gz',
        map = filtered_popmap
    output:
        'output/stacks_populations/r{r}/populations.sumstats_summary.tsv',
        'output/stacks_populations/r{r}/populations.markers.tsv',
        'output/stacks_populations/r{r}/populations.hapstats.tsv',
        'output/stacks_populations/r{r}/populations.sumstats.tsv',
        'output/stacks_populations/r{r}/populations.haplotypes.tsv',
        'output/stacks_populations/r{r}/populations.snps.genepop',
        'output/stacks_populations/r{r}/populations.snps.vcf'
    params:
        stacks_dir = 'output/stacks_denovo',
        outdir = 'output/stacks_populations/r{r}'
    threads:
        10
    log:
        'output/logs/populations_r{r}.log'
    singularity:
        stacks_container
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-O {params.outdir} '
        '-t {threads} '
        '-r {wildcards.r} '
        '--genepop --vcf '
        '&> {log}'

# 11. generate final catalog
rule gstacks:
    input:
        dynamic('output/stacks_denovo/{dyn_indiv3}.matches.bam'),
        map = filtered_popmap
    output:
        'output/stacks_denovo/gstacks.fa.gz',
        'output/stacks_denovo/gstacks.vcf.gz'
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

# 10. convert matches to BAM
rule tsv2bam:
    input:
        dynamic('output/stacks_denovo/{dyn_indiv2}.matches.tsv.gz'),
        map = filtered_popmap
    output:
        dynamic('output/stacks_denovo/{dyn_indiv3}.matches.bam')
    params:
        stacks_dir = 'output/stacks_denovo'
    threads:
        75
    log:
        'output/logs/tsv2bam{dyn_indiv3}.log'
    singularity:
            stacks_container
    shell:
        'tsv2bam '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-t {threads} '
        '&> {log}'

# 9. match individuals to the catalog
rule sstacks:
    input:
        'output/stacks_denovo/batch_1.catalog.tags.tsv.gz',
        'output/stacks_denovo/batch_1.catalog.snps.tsv.gz',
        'output/stacks_denovo/batch_1.catalog.alleles.tsv.gz',
        map = filtered_popmap
    output:
        dynamic('output/stacks_denovo/{dyn_indiv2}.matches.tsv.gz')
    params:
        stacks_dir = 'output/stacks_denovo'
    threads:
        75
    log:
        'output/logs/sstacks{dyn_indiv2}.log'
    singularity:
            stacks_container
    shell:
        'sstacks '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-p {threads} '
        '&> {log}'

# 8. generate the catalog (takes ~ 1 week)
rule cstacks:
    input:
        dynamic('output/stacks_denovo/{dyn_indiv}.alleles.tsv.gz'),
        dynamic('output/stacks_denovo/{dyn_indiv}.snps.tsv.gz'),
        dynamic('output/stacks_denovo/{dyn_indiv}.models.tsv.gz'),
        dynamic('output/stacks_denovo/{dyn_indiv}.tags.tsv.gz'),
        map = filtered_popmap
    output:
        'output/stacks_denovo/batch_1.catalog.tags.tsv.gz',
        'output/stacks_denovo/batch_1.catalog.snps.tsv.gz',
        'output/stacks_denovo/batch_1.catalog.alleles.tsv.gz'
    params:
        stacks_dir = 'output/stacks_denovo',
    threads:
        75
    log:
        'output/logs/cstacks.log'
    shell:
        'cstacks '
        '-p {threads} '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-n 3 '
        '&> {log}'

# 7d. combine individual coverage stats
rule individual_covstats_combined:
    input:
        dynamic('output/run_stats/individual_covstats/{dyn_indiv}.csv')
    output:
        combined = 'output/run_stats/individual_covstats_combined.csv'
    script:
        'src/combine_csvs.R'

# 7c. calculate coverage stats per individual
rule individual_covstats:
    input:
        tags_file = 'output/stacks_denovo/{dyn_indiv}.tags.tsv.gz'
    output:
        covstats = 'output/run_stats/individual_covstats/{dyn_indiv}.csv'
    log:
        log = 'output/logs/individual_covstats/{dyn_indiv}.log'
    threads:
        1
    script:
        'src/calculate_mean_coverage.R'

# 7b. combine individual assembly stats
rule individual_stats_combined:
    input:
        dynamic('output/run_stats/individual_stats/{dyn_indiv}.csv')
    output:
        combined = 'output/run_stats/individual_stats_combined.csv'
    script:
        'src/combine_csvs.R'

# 7a. calculate assembly stats per individual
rule individual_stats:
    input:
        alleles_file = 'output/stacks_denovo/{dyn_indiv}.alleles.tsv.gz',
        snps_file = 'output/stacks_denovo/{dyn_indiv}.snps.tsv.gz',
        tags_file = 'output/stacks_denovo/{dyn_indiv}.tags.tsv.gz'
    output:
        sample_stats = 'output/run_stats/individual_stats/{dyn_indiv}.csv'
    log:
        log = 'output/logs/individual_stats/{dyn_indiv}.log'
    threads:
        1
    script:
        'src/stacks_individual_stats.R'

# 7. assemble loci for individuals that passed the filter
rule ustacks:
    input:
        'output/run_stats/pass/{dyn_indiv}',
        individual_i_pickle = 'output/obj/individual_i.p'
    params:
        individual_i = read_sample_number,
        fastq = 'output/demux/{dyn_indiv}.fq.gz',
        wd = 'output/stacks_denovo'
    output:
        'output/stacks_denovo/{dyn_indiv}.alleles.tsv.gz',
        'output/stacks_denovo/{dyn_indiv}.snps.tsv.gz',
        'output/stacks_denovo/{dyn_indiv}.models.tsv.gz',
        'output/stacks_denovo/{dyn_indiv}.tags.tsv.gz'
    threads:
        15
    log:
        'output/logs/ustacks_{dyn_indiv}.log'
#    singularity:
#            stacks_container
    shell:
        # open the pickled dictionary and look up the sample_i
        with open(input.individual_i_pickle, 'rb') as f:
            individual_i = pickle.load(f)
        sample_i = individual_i[wildcards.individual]
        shell('ustacks '
              '-p {threads} '
              '-t gzfastq '
              '-f {params.fastq} '
              '-o {params.wd} '
              '-i {sample_i} '
              '-m 3 '
              '-M 3 '
              '&> {log}')


# 5. make a dictionary of sample:i for cstacks
rule enumerate_filtered_samples:
    input:
        map = process_reads('output/stacks_config/filtered_populations.txt')
    output:
        pickle = 'output/obj/individual_i.p'
    run:
        # read the filtered popmap
        my_popmap = pandas.read_csv(input.map,
                                    delimiter='\t',
                                    header=None)
        my_individuals = enumerate(sorted(set(my_popmap[0])))
        individual_i = {y: x for x, y in my_individuals}
        # pickle the individual_i dict for other rules to use
        with open(output.pickle, 'wb+') as f:
            pickle.dump(individual_i, f)
