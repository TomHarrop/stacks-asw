#!/usr/bin/env python3

import pandas
import re

#############
# FUNCTIONS #
#############


###########
# GLOBALS #
###########

# files and folders
key_file = 'data/reads/SQ0727.txt'

# containers
stacks_container = ('shub://TomHarrop/variant-utils:stacks_2.53')
bbduk_container = 'shub://TomHarrop/seq-utils:bbmap_38.76'
r_container = 'shub://TomHarrop/r-containers:r_3.6.1'
pandas_container = 'shub://TomHarrop/py-containers:pandas_0.25.3'

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
        'output/000_config/filtered_population_map.txt',
        'output/000_config/individual_i.p'


# 5. make a dictionary of sample:i for cstacks
rule enumerate_filtered_samples:
    input:
        key_file = key_file
    output:
        pickle = 'output/000_config/individual_i.p'
    singularity:
        pandas_container
    script:
        'src/enumerate_filtered_samples.py'

# 4. filter the population map
rule filter_samples:
    input:
        popmap = 'output/000_config/population_map.txt',
        read_stats = 'output/040_stats/reads.csv',
        gc_stats = 'output/040_stats/gc_stats.csv'
    output:
        map = 'output/000_config/filtered_population_map.txt',
        plot = 'output/040_stats/read_count_histogram.pdf'
    log:
        'output/logs/filter_samples.log'
    singularity:
        r_container
    script:
        'src/filter_population_map.R'

# 4. run reformat.sh to count reads and get a gc histogram
rule combine_individual_stats:
    input:
        read_stats = expand('output/040_stats/reads/{individual}.txt',
                            individual=all_individuals),
        gc_stats = expand('output/040_stats/gc_hist/{individual}.txt',
                          individual=all_individuals)
    output:
        read_stats = 'output/040_stats/reads.csv',
        gc_stats = 'output/040_stats/gc_stats.csv',
        gc_hist = 'output/040_stats/gc_hist.csv'
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
        'output/030_combined/{individual}.fq.gz'
    output:
        reads = 'output/040_stats/reads/{individual}.txt',
        gc = 'output/040_stats/gc_hist/{individual}.txt'
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
            'output/020_filtered/kept/{sample_fullname}.fq.gz',
            sample_fullname=individual_to_sample_fullname[wildcards.individual])
    output:
        'output/030_combined/{individual}.fq.gz'
    shell:
        'cat {input} > {output}'

# # 2b. filter and truncate demuxed reads
rule trim_adaptors:
    input:
        'output/010_demux/{sample_fullname}.fq.gz'
    output:
        kept = 'output/020_filtered/kept/{sample_fullname}.fq.gz',
        discarded = 'output/020_filtered/discarded/{sample_fullname}.fq.gz',
        adaptor_stats = 'output/020_filtered/adaptor_stats/{sample_fullname}.txt',
        truncate_stats = ('output/020_filtered/truncate_stats/'
                          '{sample_fullname}.txt'),
        gc_hist = 'output/020_filtered/gc_hist/{sample_fullname}.txt',
        lhist = 'output/020_filtered/lhist/{sample_fullname}.txt'
    params:
        adaptors = 'data/adaptors/bbduk_adaptors_plus_AgR_common.fa'
    singularity:
        bbduk_container
    log:
        adaptors = 'output/logs/trim_adaptors.{sample_fullname}.adaptors.txt',
        truncate = 'output/logs/trim_adaptors.{sample_fullname}.truncate.txt',
        lhist = 'output/logs/trim_adaptors.{sample_fullname}.lhist.txt'
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
            config_file = 'output/000_config/{0}.config'.format(fc_lane)
        output:
            expand(
                'output/010_demux/{sample_fullname}.fq.gz',
                sample_fullname=fc_name_to_sample_fullname[fc_lane])
        log:
            'output/logs/demux.{0}.log'.format(fc_lane)
        threads:
            1
        singularity:
            stacks_container
        shell:
            'process_radtags '
            '-f {input.read_file} '
            '-i gzfastq -y gzfastq '
            '-b {input.config_file} '
            '-o output/010_demux '
            '-c -q '
            # '-r --barcode_dist_1 1 '  # rescue barcodes
            '--barcode_dist_1 0 '       # don't allow bc mismatches
            # '-t 91 '                  # truncate output to 91 b
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
        expand('output/000_config/{fc_lane}.config',
               fc_lane=all_fc_lanes),
        population_map = 'output/000_config/population_map.txt'
    params:
        outdir = 'output/000_config'
    log:
        'output/logs/write_config_files.log'
    singularity:
        r_container
    script:
        'src/write_config_files.R'

