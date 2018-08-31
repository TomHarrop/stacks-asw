#!/usr/bin/env python3

import csv
import numpy
import os
import pandas
import pathlib2
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
    binary_path = pathlib2.Path(which).resolve()
    return str(binary_path)


def touch(fname, mode=0o666, dir_fd=None, **kwargs):
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(f.fileno() if os.utime in os.supports_fd else fname,
                 dir_fd=None if os.supports_fd else dir_fd, **kwargs)


###########
# GLOBALS #
###########

# files and folders
key_file = 'data/reads/SQ0727.txt'
filtered_popmap = 'output/stacks_config/filtered_populations.txt'

# filtering parameters
r_values = list(str(x) for x in numpy.arange(0, 1.01, 0.1))

# containers
stacks_container = 'shub://TomHarrop/singularity-containers:stacks_2.0b'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'

#########
# SETUP #
#########

# read key file
key_data = pandas.read_csv(key_file, delimiter='\t')

# remove spaces from mararoa-downs
key_data['sample'] = key_data['sample'].str.replace('\s', '-', regex=True)
key_data['sample'] = key_data['sample'].str.replace('.', '-', regex=False)

# add details for expected output
key_data['fc_lane'] = key_data[[
    'flowcell', 'lane']].astype(str).apply('_'.join, axis=1)
key_data['sample_fullname'] = key_data[[
    'sample',
    'flowcell',
    'lane',
    'row',
    'column']].astype(str).apply('_'.join, axis=1)

all_fc_lanes = sorted(set(key_data['fc_lane']))
all_fullnames = sorted(set(key_data['sample_fullname']))

# get a dict of fc to sample name
col_to_fcl = key_data.to_dict()['fc_lane']
col_to_sn = key_data.to_dict()['sample_fullname']

fc_name_to_sample_fullname = dict((k, []) for k in all_fc_lanes)

for key in col_to_sn:
    fc_name_to_sample_fullname[col_to_fcl[key]].append(col_to_sn[key])

# # get a list of fastq files
# all_read_files = []
# read_dir_files = list((dirpath, filenames)
#                       for (dirpath, dirnames, filenames)
#                       in os.walk(reads_dir))

# for dirpath, filenames in read_dir_files:
#     for filename in filenames:
#         if 'fastq.gz' in filename:
#             all_read_files.append(os.path.join(dirpath, filename))

# # get dicts of flowcell_lane:sample and sample:flowcell_lane
# fc_lane_to_sample = {}
# sample_to_fc_lane = {}
# for name, group in grouped_key_data:
#     fc_lane = '_'.join([str(x) for x in name])
#     sample_list = list(group['sample'])
#     fc_lane_to_sample[fc_lane] = sample_list
#     for sample in sample_list:
#         sample_to_fc_lane[sample] = fc_lane

# # get a list of samples
# all_samples = sorted(set(x for x in sample_to_fc_lane.keys()
#                          if any(list(sample_to_fc_lane[x] in y
#                                      for y in all_read_files))))
# all_fc_lanes = [x for x in fc_lane_to_sample
#                 if any([x in y for y in all_read_files])]


#########
# RULES #
#########

rule target:
    input:
        expand(
            'output/filtering/kept/{individual}.fq.gz',
            individual=all_fullnames)
  
# 4. filter the population map
# rule filter_samples:
#     input:
#         stats = 'output/run_stats/read_stats.txt',
#         popmap = 'output/stacks_config/population_map.txt'
#     params:
#         sample_dir = 'output/demux',
#     output:
#         map = filtered_popmap,
#         plot = 'output/run_stats/read_count_histogram.pdf',
#         pop_counts = 'output/run_stats/individuals_per_population.csv'
#     singularity:
#         'shub://TomHarrop/singularity-containers:r_3.5.0'
#     script:
#         'src/filter_population_map.R'

# # 3. check number of reads per sample to filter out low coverage indivs
# rule count_reads:
#     input:
#         tag_files = expand(
#             'output/filtering/kept/{individual}.fq.gz',
#             individual=all_samples)
#     params:
#         inline = lambda wildcards, input: ','.join(input.tag_files)
#     output:
#         'output/run_stats/read_stats.txt'
#     log:
#         'output/logs/statswrapper.log'
#     singularity:
#         'shub://TomHarrop/singularity-containers:bbmap_38.00'
#     shell:
#         'statswrapper.sh '
#         'in={params.inline} '
#         'out={output} '
#         '2> {log}'

# # 2b. filter and truncate demuxed reads
rule trim_adaptors:
    input:
        'output/demux/{individual}.fq.gz'
    output:
        kept = 'output/filtering/kept/{individual}.fq.gz',
        discarded = 'output/filtering/discarded/{individual}.fq.gz',
        adaptor_stats = 'output/filtering/adaptor_stats/{individual}.txt',
        truncate_stats = 'output/filtering/truncate_stats/{individual}.txt',
        gc_hist = 'output/filtering/gc_hist/{individual}.txt',
        lhist = 'output/filtering/lhist/{individual}.txt'
    params:
        adaptors = 'data/adaptors/bbduk_adaptors_plus_AgR_common.fa'
    singularity:
        bbduk_container
    log:
        adaptors = 'output/logs/filtering/{individual}_adaptors.txt',
        truncate = 'output/logs/filtering/{individual}_truncate.txt',
        lhist = 'output/logs/filtering/{individual}_lhist.txt'
    threads:
        1
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
        'reformat.sh '
        'in=stdin.fq '
        'out=stdout.fq '
        'lhist={output.lhist} '
        '2> {log.lhist}'
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
            read_file = 'data/reads/{0}_fastq.gz'.format(fc_lane),
            config_file = 'output/stacks_config/{0}.config'.format(fc_lane)
        output:
            expand(
                'output/demux/{individual}.fq.gz',
                individual=fc_name_to_sample_fullname[fc_lane])
        log:
            'output/logs/demux/{0}.log'.format(fc_lane)
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
rule write_config_files:
    input:
        key_file = key_file
    output:
        expand('output/stacks_config/{fc_lane}.config',
               fc_lane=all_fc_lanes),
        population_map = 'output/stacks_config/population_map.txt'
    params:
        outdir = 'output/stacks_config'
    log:
        'output/logs/write_config_files.log'
    singularity:
        r_container
    script:
        'src/write_config_files.R'

