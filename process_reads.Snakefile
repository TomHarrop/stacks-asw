#!/usr/bin/env python3

import pandas
import re
from pathlib import Path


#############
# FUNCTIONS #
#############

def aggregate_fullnames(wildcards):
    # we're not going to use wildcards, we're going to use the list of all
    # files. But we need to run `checkpoints` for snakemake's benefit
    for fc in all_fc_lanes:
        # check which fq files resulted from demuxing the fc_lane
        co = checkpoints.process_radtags.get(fc_lane=fc).output['fq']
    read_stats = snakemake.io.expand(
        'output/040_stats/reads/{individual}.txt',
        individual=all_individuals)
    gc_stats = snakemake.io.expand(
        'output/040_stats/gc_hist/{individual}.txt',
        individual=all_individuals)
    # return as a dict
    my_dict = {}
    my_dict['read_stats'] = read_stats
    my_dict['gc_stats'] = gc_stats
    return(my_dict)


def get_fq_files_for_indiv(wildcards):
    my_indiv = wildcards.individual
    if my_indiv in geo_individuals:
        my_key_file = geo_key_data
    elif my_indiv in para_individuals:
        my_key_file = para_key_data
    else:
        raise ValueError(f'wtf {my_indiv}')
    my_mask = my_key_file['sample'] == my_indiv
    my_df = my_key_file[my_mask]
    my_fullnames = sorted(set(my_df['sample_fullname']))
    my_fqs = snakemake.io.expand(
            'output/020_filtered/kept/{sample_fullname}.fq.gz',
            sample_fullname=my_fullnames)
    return(my_fqs)


def resolve_demuxed_file(wildcards):
    my_fullname = wildcards.sample_fullname
    if my_fullname in geo_fullnames:
        my_key_file = geo_key_data
    elif my_fullname in para_fullnames:
        my_key_file = para_key_data
    else:
        raise ValueError(f'wtf {my_fullname}')
    my_mask = my_key_file['sample_fullname'] == my_fullname
    my_df = my_key_file[my_mask]
    my_fc = sorted(set(my_df['fc_lane'].values))[0]
    my_fq = f'output/010_demux/{my_fc}/{my_fullname}.fq.gz'
    return(my_fq)


def resolve_read_file(wildcards):
    if wildcards.fc_lane in geo_fc_lanes:
        read_file = f'data/georeads/{wildcards.fc_lane}_fastq.gz'
        return {'read_file': read_file}
    elif wildcards.fc_lane in para_fc_lanes:
        fc_split = wildcards.fc_lane.split('_')
        read_file = next(Path('data/parareads/').glob(
            f'{fc_split[0]}*L*{fc_split[1]}*.fastq.gz')).as_posix()
        return {'read_file': read_file}
    else:
        raise ValueError(f'wtf {wildcards.fc_lane}')

###########
# GLOBALS #
###########

# files and folders
geo_key_file = 'data/georeads/SQ0727.txt'
para_key_file = 'data/para_key_data.csv'

# containers
bbduk_container = 'shub://TomHarrop/seq-utils:bbmap_38.76'
pandas_container = 'shub://TomHarrop/py-containers:pandas_0.25.3'
r_container = 'shub://TomHarrop/r-containers:r_3.6.1'
stacks_container = 'shub://TomHarrop/variant-utils:stacks_2.53'

#########
# SETUP #
#########

# geo samples

# read key file
geo_key_data = pandas.read_csv(geo_key_file, delimiter='\t')

# remove spaces from mararoa-downs and periods from Ophir
geo_key_data['sample'] = geo_key_data[
    'sample'].str.replace('\s', '-', regex=True)
geo_key_data['sample'] = geo_key_data[
    'sample'].str.replace('.', '-', regex=False)
geo_key_data['sample_name'] = geo_key_data['sample']

geo_key_data['sample'] = [
    f'geo_{x}' for x in geo_key_data['sample'].values]

# add details for expected output
geo_key_data['fc_lane'] = geo_key_data[[
    'flowcell', 'lane']].astype(str).apply('_'.join, axis=1)
geo_key_data['sample_fullname'] = geo_key_data[[
    'sample_name',
    'flowcell',
    'lane',
    'row',
    'column']].astype(str).apply('_'.join, axis=1)

geo_fc_lanes = sorted(set(geo_key_data['fc_lane']))
geo_fullnames = sorted(set(geo_key_data['sample_fullname']))
geo_individuals = sorted(set(geo_key_data['sample']))

# para samples

# read key file
para_key_data = pandas.read_csv(para_key_file, delimiter=',')

# add details for expected output
para_key_data['sample'] = para_key_data['sample_name']
para_key_data['fc_lane'] = para_key_data[[
    'key', 'lane']].astype(str).apply('_'.join, axis=1)
para_key_data['sample_fullname'] = para_key_data[[
    'sample',
    'key',
    'lane']].astype(str).apply('_'.join, axis=1)
para_key_data['sample'] = [
    f'para_{x}' for x in para_key_data['sample'].values]

para_fc_lanes = sorted(set(para_key_data['fc_lane']))
para_fullnames = sorted(set(para_key_data['sample_fullname']))
para_individuals = sorted(set(para_key_data['sample']))

# whole dataset
all_individuals = sorted(set(geo_individuals + para_individuals))
all_fc_lanes = sorted(set(geo_fc_lanes + para_fc_lanes))


#########
# RULES #
#########

wildcard_constraints:
    fc_lane = '|'.join(all_fc_lanes),
    individual = '|'.join(all_individuals)

rule target:
    input:
        expand('output/010_demux/{fc_lane}',
               fc_lane=all_fc_lanes),
        'output/000_config/filtered_population_map.txt'

# filter the population map
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

# count reads and get a gc histogram
rule combine_individual_stats:
    input:
        unpack(aggregate_fullnames)
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

# combine reads per-individual
rule combine_reads:
    input:
        get_fq_files_for_indiv
    output:
        'output/030_combined/{individual}.fq.gz'
    shell:
        'cat {input} > {output}'

# filter and truncate demuxed reads
rule trim_adaptors:
    input:
        resolve_demuxed_file
    output:
        kept = 'output/020_filtered/kept/{sample_fullname}.fq.gz',
        discarded = 'output/020_filtered/discarded/{sample_fullname}.fq.gz',
        adaptor_stats = 'output/020_filtered/adaptor_stats/{sample_fullname}.txt',
        truncate_stats = ('output/020_filtered/truncate_stats/'
                          '{sample_fullname}.txt'),
        gc_hist = 'output/020_filtered/gc_hist/{sample_fullname}.txt',
        lhist = 'output/020_filtered/lhist/{sample_fullname}.txt'
    params:
        adaptors = 'data/bbduk_adaptors_plus_AgR_common.fa'
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


# demux each fc_lane
checkpoint process_radtags:
    input:
        unpack(resolve_read_file),
        config_file = 'output/000_config/{fc_lane}.config'
    output:
        fq = directory('output/010_demux/{fc_lane}')
    log:
        'output/logs/demux.{fc_lane}.log'
    threads:
        1
    priority:
        100
    singularity:
        stacks_container
    shell:
        'mkdir -p {output.fq} ; '
        'process_radtags '
        '-f {input.read_file} '
        '-i gzfastq -y gzfastq '
        '-b {input.config_file} '
        '-o {output.fq} '
        '-c -q '
        # '-r --barcode_dist_1 1 '  # rescue barcodes
        '--barcode_dist_1 0 '       # don't allow bc mismatches
        # '-t 91 '                  # truncate output to 91 b
        '-w 0.1 '                   # window: approx. 9 bases
        '-s 15 '                    # minimum avg PHRED in window
        '--inline_null '
        '--renz_1 apeKI --renz_2 mspI '
        '&> {log} '
        # '|| true'                   # DOESN'T EXIT CLEANLY, WHYYYYY?

# extract per-flowcell/lane sample:barcode information
rule write_config_files:
    input:
        geo_key_file = geo_key_file,
        para_key_file = para_key_file
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
