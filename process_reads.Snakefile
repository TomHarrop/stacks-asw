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
all_individuals = sorted(set(key_data['sample']))

# get a dict of fc to sample name
col_to_fcl = key_data.to_dict()['fc_lane']
col_to_sn = key_data.to_dict()['sample_fullname']

fc_name_to_sample_fullname = dict((k, []) for k in all_fc_lanes)

for key in col_to_sn:
    fc_name_to_sample_fullname[col_to_fcl[key]].append(col_to_sn[key])

# get a dict of individual to sample_fullname
individual_to_sample_fullname = dict((k, []) for k in all_individuals)
for key in individual_to_sample_fullname:
    individual_to_sample_fullname[key] = sorted(
        set([x for x in all_fullnames
             if re.sub('_.*', '', x) == key]))

#########
# RULES #
#########

rule target:
    input:
        'output/combined_stats/reads.csv'
  
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

# 4. run reformat.sh to count reads and get a gc histogram
rule combine_individual_stats:
    input:
        read_stats = expand('output/individual_stats/reads/{individual}.txt',
                            individual=all_individuals),
        gc_stats = expand('output/individual_stats/gc_hist/{individual}.txt',
                          individual=all_individuals)
    output:
        read_stats = 'output/combined_stats/reads.csv',
        gc_stats = 'output/combined_stats/gc_stats.csv',
        gc_hist = 'output/combined_stats/gc_hist.csv'
    log:
        'output/logs/combine_individual_stats.log'
    threads:
        1
    singularity:
        r_container
    script:
        'src/combine_individual_stats.R'

rule count_reads:
    input:
        'output/combined/{individual}.fq.gz'
    output:
        reads = 'output/individual_stats/reads/{individual}.txt',
        gc = 'output/individual_stats/gc_hist/{individual}.txt'
    threads:
        1
    priority:
        1
    singularity:
        bbduk_container
    shell:
        'reformat.sh '
        'in={input} '
        'int=f '
        'out=/dev/null '
        'gchist={output.gc} '
        'gcbins=auto '
        '2> {output.reads}'

# 3 combine reads per-individual
rule combine_reads:
    input:
        lambda wildcards: expand(
            'output/filtering/kept/{sample_fullname}.fq.gz',
            sample_fullname=individual_to_sample_fullname[wildcards.individual])
    output:
        'output/combined/{individual}.fq.gz'
    shell:
        'cat {input} > {output}'

# # 2b. filter and truncate demuxed reads
rule trim_adaptors:
    input:
        'output/demux/{sample_fullname}.fq.gz'
    output:
        kept = 'output/filtering/kept/{sample_fullname}.fq.gz',
        discarded = 'output/filtering/discarded/{sample_fullname}.fq.gz',
        adaptor_stats = 'output/filtering/adaptor_stats/{sample_fullname}.txt',
        truncate_stats = ('output/filtering/truncate_stats/'
                          '{sample_fullname}.txt'),
        gc_hist = 'output/filtering/gc_hist/{sample_fullname}.txt',
        lhist = 'output/filtering/lhist/{sample_fullname}.txt'
    params:
        adaptors = 'data/adaptors/bbduk_adaptors_plus_AgR_common.fa'
    singularity:
        bbduk_container
    log:
        adaptors = 'output/logs/filtering/{sample_fullname}_adaptors.txt',
        truncate = 'output/logs/filtering/{sample_fullname}_truncate.txt',
        lhist = 'output/logs/filtering/{sample_fullname}_lhist.txt'
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
        'int=f '
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
                'output/demux/{sample_fullname}.fq.gz',
                sample_fullname=fc_name_to_sample_fullname[fc_lane])
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

