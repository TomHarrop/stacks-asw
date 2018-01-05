#!/usr/bin/env python3

import csv
import numpy
import os
import pandas
import pathlib
import pickle
import re
import shutil


#############
# FUNCTIONS #
#############

def get_full_path(binary):
    which = shutil.which(binary)
    # check if the binary exists
    if not which:
        raise EnvironmentError(
            'Dependency {0} not found in $PATH'.format(binary))
    # get the full path to binary
    binary_path = pathlib.Path(which).resolve()
    return str(binary_path)


def parse_key_and_write_config_files(key_file, outdir):
    '''Group key_file rows by 'Flowcell' and 'Lane', and write tsv of 'Barcode'
    and 'Sample' separately for each lane and flowcell to 'outdir' '''
    # generate dicts
    key_data = pandas.read_csv(key_file, delimiter='\t')
    grouped_key_data = key_data.groupby(['Flowcell', 'Lane'])
    # make output directory
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # write the key files by group
    for name, group in grouped_key_data:
        prefix = '_'.join([str(x) for x in name])
        config_file = os.path.join(outdir, '%s.config' % prefix)
        subset = group[['Barcode', 'Sample']]
        if len(subset) > 0:
            subset.to_csv(config_file,
                          sep='\t',
                          header=False,
                          index=False)
    # generate population map
    sample_to_population = {}
    for sample in key_data['Sample']:
        sample_to_population[sample] = re.sub('\d', '', sample).lower()
    # write the population map
    population_map = os.path.join(outdir, 'population_map.txt')
    with open(population_map, 'w') as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows([x, sample_to_population[x]]
                         for x in sample_to_population
                         if x in all_samples)


def touch(fname, mode=0o666, dir_fd=None, **kwargs):
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(f.fileno() if os.utime in os.supports_fd else fname,
                 dir_fd=None if os.supports_fd else dir_fd, **kwargs)


###########
# GLOBALS #
###########

# files and folders
key_file = 'data/SQ0003.txt'
reads_dir = 'data/raw_reads'
outdir = 'output'
filtered_popmap = 'output/stacks_config/filtered_populations.txt'

# filtering paramaters
r_values = list(str(x) for x in numpy.arange(0, 1.01, 0.1))

#########
# SETUP #
#########

# read key file
key_data = pandas.read_csv(key_file, delimiter='\t')
grouped_key_data = key_data.groupby(['Flowcell', 'Lane'])

# get a list of fastq files
all_read_files = []
read_dir_files = list((dirpath, filenames)
                      for (dirpath, dirnames, filenames)
                      in os.walk(reads_dir))

for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            all_read_files.append(os.path.join(dirpath, filename))

# get dicts of flowcell_lane:sample and sample:flowcell_lane
fc_lane_to_sample = {}
sample_to_fc_lane = {}
for name, group in grouped_key_data:
    fc_lane = '_'.join([str(x) for x in name])
    sample_list = list(group['Sample'])
    fc_lane_to_sample[fc_lane] = sample_list
    for sample in sample_list:
        sample_to_fc_lane[sample] = fc_lane

# get a list of samples
all_samples = sorted(set(x for x in sample_to_fc_lane.keys()
                         if any(list(sample_to_fc_lane[x] in y
                                     for y in all_read_files))))
all_fc_lanes = [x for x in fc_lane_to_sample
                if any([x in y for y in all_read_files])]


#########
# RULES #
#########

rule target:
    input:
        expand('output/stacks_populations/r{r}/populations.sumstats_summary.tsv',
               r=r_values),
        dynamic('output/run_stats/individual_stats/{dyn_indiv}.csv')
        # 'output/run_stats/population_stats_combined.csv'
        # 'output/run_stats/individual_stats_combined.csv'
        # 'output/run_stats/individual_covstats_combined.csv'

# extract per-flowcell/lane sample:barcode information
rule extract_barcode_config:
    input:
        key_file
    output:
        expand('output/stacks_config/{fc_lane}.config',
               fc_lane=fc_lane_to_sample.keys()),
        population_map = 'output/stacks_config/population_map.txt'
    run:
        parse_key_and_write_config_files(
            key_file,
            'output/stacks_config')

# for loop per fc_lane
for fc_lane in all_fc_lanes:
    rule:
        input:
            read_file = [x for x in all_read_files if fc_lane in x][0],
            config_file = 'output/stacks_config/{}.config'.format(fc_lane)
        output:
            expand('output/demux/{individual}.fq.gz',
                   individual=fc_lane_to_sample[fc_lane]),
        log:
            'output/logs/demux_{0}.log'.format(fc_lane)
        threads:
            1
        shell:
            'process_radtags '
            '-f {input.read_file} '
            '-i gzfastq -y gzfastq '
            '-b {input.config_file} '
            '-o output/demux '
            '-c -q '
            # '-r --barcode_dist_1 1 '    # rescue barcodes
            '-t 91 '                    # truncate output to 91 b
            '-w 0.1 '                   # window: approx. 9 bases
            '-s 15 '                    # minimum avg PHRED in window
            '--inline_null '
            '--renz_1 apeKI --renz_2 mspI '
            '&> {log}'

rule count_reads:
    input:
        tag_files = expand(
            'output/demux/{individual}.fq.gz',
            individual=all_samples)
    output:
        'output/run_stats/read_stats.txt'
    log:
        'output/logs/statswrapper.log'
    run:
        in_line = ','.join(input.tag_files)
        my_statswrapper = get_full_path('statswrapper.sh')
        shell('{my_statswrapper} '
              'in={in_line} '
              'out={output} '
              '2> {log}')

rule filter_samples:
    input:
        stats = 'output/run_stats/read_stats.txt',
        popmap = 'output/stacks_config/population_map.txt'
    params:
        sample_dir = 'output/demux',
    output:
        map = filtered_popmap,
        plot = 'output/run_stats/read_count_histogram.pdf',
        pop_counts = 'output/run_stats/individuals_per_population.csv'
    script:
        'src/filter_population_map.R'

rule enumerate_filtered_samples:
    input:
        map = filtered_popmap
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

rule select_filtered_samples:
    input:
        map = filtered_popmap
    output:
        dynamic('output/run_stats/pass/{dyn_indiv}')
    params:
        outdir = 'output/run_stats/pass'
    run:
        # read the filtered popmap 
        my_popmap = pandas.read_csv(input.map,
           delimiter='\t',
           header=None)
        # touch flag files
        for indiv in sorted(set(my_popmap[0])):
            my_path = os.path.join(params.outdir, indiv)
            touch(my_path)

rule ustacks:
    input:
        'output/run_stats/pass/{dyn_indiv}',
        individual_i_pickle = 'output/obj/individual_i.p'
    params:
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
    run:
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

# rule combine_individual_stats:
#     input:
#         dynamic('output/run_stats/individual_stats/{dyn_indiv}.csv')
#     output:
#         combined = 'output/run_stats/individual_stats_combined.csv'
#     script:
#         'src/combine_csvs.R'

# rule individual_covstats:
#     input:
#         tags_file = 'output/stacks_denovo/{dyn_indiv}.tags.tsv.gz'
#     output:
#         covstats = 'output/run_stats/individual_covstats/{dyn_indiv}.csv'
#     log:
#         log = 'output/logs/individual_covstats/{dyn_indiv}.log'
#     threads:
#         1
#     script:
#         'src/calculate_mean_coverage.R'

# rule combine_individual_covstats:
#     input:
#         dynamic('output/run_stats/individual_covstats/{dyn_indiv}.csv'),
#         map = filtered_popmap
#     output:
#         combined = 'output/run_stats/individual_covstats_combined.csv'
#     script:
#         'src/combine_csvs.R'

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
        'output/logs/sstacks.log'
    shell:
        'sstacks '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-p {threads} '
        '&> {log}'

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
        'output/logs/tsv2bam.log'
    shell:
        'tsv2bam '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-t {threads} '
        '&> {log}'

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
    shell:
        'gstacks '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-t {threads} '
        '&> {log}'

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
        'output/stacks_populations/r{r}/populations.haplotypes.tsv'
    params:
        stacks_dir = 'output/stacks_denovo',
        outdir = 'output/stacks_populations/r{r}'
    threads:
        15
    log:
        'output/logs/populations_r{r}.log'
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-O {params.outdir} '
        '-t {threads} '
        '-r {wildcards.r} '
        '&> {log}'

rule population_stats:
    input:
        sumstats = 'output/stacks_populations/r{r}/populations.sumstats.tsv',
        hapstats = 'output/stacks_populations/r{r}/populations.hapstats.tsv'
    output:
        pop_stats = 'output/run_stats/population_stats/{r}.csv'
    log:
        log = 'output/logs/population_stats/{r}.log'
    threads:
        1
    script:
        'src/stacks_population_stats.R'

rule combine_population_stats:
    input:
        expand('output/run_stats/population_stats/{r}.csv',
               r=r_values)
    output:
        combined = 'output/run_stats/population_stats_combined.csv'
    script:
        'src/combine_csvs.R'

