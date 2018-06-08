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

# filtering parameters
r_values = list(str(x) for x in numpy.arange(0, 1.01, 0.1))

# containers
stacks_container = 'shub://TomHarrop/singularity-containers:stacks_2.0b'


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
        filtered_popmap
  
# 4. filter the population map
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
    singularity:
        'shub://TomHarrop/singularity-containers:r_3.5.0'
    script:
        'src/filter_population_map.R'

# 3. check number of reads per sample to filter out low coverage indivs
rule count_reads:
    input:
        tag_files = expand(
            'output/filtering/kept/{individual}.fq.gz',
            individual=all_samples)
    params:
        inline = lambda wildcards, input: ','.join(input.tag_files)
    output:
        'output/run_stats/read_stats.txt'
    log:
        'output/logs/statswrapper.log'
    singularity:
        'shub://TomHarrop/singularity-containers:bbmap_38.00'
    shell:
        'statswrapper.sh '
        'in={params.inline} '
        'out={output} '
        '2> {log}'

# 2b. filter and truncate demuxed reads
rule trim_adaptors:
    input:
        'output/demux/{individual}.fq.gz'
    output:
        kept = 'output/filtering/kept/{individual}.fq.gz',
        discarded = 'output/filtering/discarded/{individual}.fq.gz',
        adaptor_stats = 'output/filtering/adaptor_stats/{individual}.txt',
        truncate_stats = 'output/filtering/truncate_stats/{individual}.txt',
        gc_hist = 'output/filtering/gc_hist/{individual}.txt'
    params:
        adaptors = 'data/adaptors/bbduk_adaptors_plus_AgR_common.fa'
    singularity:
        'shub://TomHarrop/singularity-containers:bbmap_38.00'
    log:
        adaptors = 'output/logs/filtering/{individual}_adaptors.txt',
        truncate = 'output/logs/filtering/{individual}_truncate.txt'
    threads:
        10
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input} '
        'interleaved=f '
        'out=stdout.fq '
        'stats={output.adaptor_stats} '
        'ref={params.adaptors} '
        'ktrim=r k=23 mink=11 hdist=1 '
        'findbestmatch=t '
        '2> {log.adaptors}'
        ' | '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fq '
        'interleaved=f '
        'outnonmatch={output.kept} '
        'outmatch={output.discarded} '
        'stats={output.truncate_stats} '
        'gchist={output.gc_hist} '
        'gcbins=auto '
        'overwrite=t '
        'forcetrimright=79 '
        'minlength=80 '
        'ziplevel=9 '
        '2> {log.truncate}'

# 2. for loop per fc_lane
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
        singularity:
            stacks_container
        shell:
            'process_radtags '
            '-f {input.read_file} '
            '-i gzfastq -y gzfastq '
            '-b {input.config_file} '
            '-o output/demux '
            '-c -q '
            # '-r --barcode_dist_1 1 '  # rescue barcodes
            '--barcode_dist_1 0 '       # don't allow bc mismatches
            # '-t 91 '                    # truncate output to 91 b
            '-w 0.1 '                   # window: approx. 9 bases
            '-s 15 '                    # minimum avg PHRED in window
            '--inline_null '
            '--renz_1 apeKI --renz_2 mspI '
            '&> {log}'

# 1. extract per-flowcell/lane sample:barcode information
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

