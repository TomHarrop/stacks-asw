#!/usr/bin/env python3

import pandas
from pathlib import Path

###########
# GLOBALS #
###########

reads_dir = 'data/raw_reads'
filtered_popmap = 'output/stacks_config/filtered_population_map.txt'

#########
# SETUP #
#########

# singularity args
full_path = Path('.').resolve()
root_dir = ''.join(full_path.parts[0:2])
singularity_args = '-B {0}'.format(root_dir)

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
        'output/parameters/compare_defaults/optimised_samplestats_combined.csv'

subworkflow process_reads:
    snakefile: 'process_reads.Snakefile'

rule compare_defaults:
    input:
        'output/parameters/filtering/replicate_1_popmap.txt',
        'output/parameters/stats_Mm/samplestats_combined.csv',
        'output/parameters/stats_n/samplestats_combined.csv',
        process_reads(expand('output/combined/{individual}.fq.gz',
                             individual=all_indivs)),
        popmap = process_reads('output/stacks_config/filtered_population_map.txt')
    output:
        'output/parameters/compare_defaults/optimised_samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/parameters',
        indir = 'output/combined'
    log:
        'output/logs/parameters/compare_defaults.log'
    shell:
        'stacks_parameters '
        '--mode compare_defaults '
        '-M 3 -m 3 -n 4 '
        '-o {params.outdir} '
        '--individuals 20 '
        '--replicates 3 '
        '--threads {threads} '
        '--singularity_args \"{singularity_args}\" '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '


rule optim_n:
    input:
        'output/parameters/filtering/replicate_1_popmap.txt',
        'output/parameters/stats_Mm/samplestats_combined.csv',
        process_reads(expand('output/combined/{individual}.fq.gz',
                             individual=all_indivs)),
        popmap = process_reads('output/stacks_config/filtered_population_map.txt')
    output:
        'output/parameters/stats_n/samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/parameters',
        indir = 'output/combined'
    log:
        'output/logs/parameters/optim_n.log'
    shell:
        'stacks_parameters '
        '--mode optim_n '
        '-M 3 -m 3 '
        '-o {params.outdir} '
        '--individuals 20 '
        '--replicates 3 '
        '--threads {threads} '
        '--singularity_args \"{singularity_args}\" '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '


rule optim_mM:
    input:
        'output/parameters/filtering/replicate_1_popmap.txt',
        process_reads(expand('output/combined/{individual}.fq.gz',
                             individual=all_indivs)),
        popmap = process_reads('output/stacks_config/filtered_population_map.txt')
    output:
        'output/parameters/stats_Mm/samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/parameters',
        indir = 'output/combined'
    log:
        'output/logs/parameters/optim_mM.log'
    shell:
        'stacks_parameters '
        '--mode optim_Mm '
        '-o {params.outdir} '
        '--individuals 20 '
        '--replicates 3 '
        '--threads {threads} '
        '--singularity_args \"{singularity_args}\" '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '

rule optim_setup:
    input:
        process_reads(expand('output/combined/{individual}.fq.gz',
                             individual=all_indivs)),
        popmap = process_reads('output/stacks_config/filtered_population_map.txt')
    output:
        'output/parameters/filtering/replicate_1_popmap.txt'
    threads:
        50
    params:
        outdir = 'output/parameters',
        indir = 'output/combined'
    log:
        'output/logs/parameters/optim_setup.log'
    shell:
        'stacks_parameters '
        '--mode setup '
        '-o {params.outdir} '
        '--individuals 20 '
        '--replicates 3 '
        '--threads {threads} '
        '--singularity_args \"{singularity_args}\" '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '
