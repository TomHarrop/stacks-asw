#!/usr/bin/env python3

import pandas
import re
from pathlib import Path


#############
# FUNCTIONS #
#############

def resolve_read_file(wildcards):
    if wildcards.fc_lane in geo_fc_lanes:
        return {'read_file':
                f'data/georeads/{wildcards.fc_lane}_fastq.gz'}
    elif wildcards.fc_lane in para_fc_lanes:
        fc_split = wildcards.fc_lane.split('_')
        return {'read_file':
                next(Path('data/parareads/').glob(
                    f'{fc_split[0]}*L*{fc_split[1]}*.fastq.gz'))}
    else:
        raise ValueError(f'wtf {wildcards.fc_lane}')


# def resolve_read_file(my_fc_lane):
#     print(f'function caught my_fc_lane: {my_fc_lane}')
#     if str(my_fc_lane) == '':
#         return {my_read_file: ''}
#     if my_fc_lane in geo_fc_lanes:
#         my_read_file = f'data/georeads/{my_fc_lane}_fastq.gz'
#         print(my_read_file)
#         return {'read_file': my_read_file}
#     elif my_fc_lane in para_fc_lanes:
#         fc_split = my_fc_lane.split('_')
#         my_read_file = next(Path('data/parareads/').glob(
#             f'{fc_split[0]}*L*{fc_split[1]}*.fastq.gz'))
#         print(my_read_file)
#         return {'read_file': my_read_file}
#     else:
#         raise ValueError(f'wtf {my_fc_lane}')


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

# remove spaces from mararoa-downs
geo_key_data['sample'] = geo_key_data[
    'sample'].str.replace('\s', '-', regex=True)
geo_key_data['sample'] = geo_key_data[
    'sample'].str.replace('.', '-', regex=False)

# add details for expected output
geo_key_data['fc_lane'] = geo_key_data[[
    'flowcell', 'lane']].astype(str).apply('_'.join, axis=1)
geo_key_data['sample_fullname'] = geo_key_data[[
    'sample',
    'flowcell',
    'lane',
    'row',
    'column']].astype(str).apply('_'.join, axis=1)

geo_fc_lanes = sorted(set(geo_key_data['fc_lane']))
geo_fullnames = sorted(set(geo_key_data['sample_fullname']))
geo_individuals = sorted(set(geo_key_data['sample']))

# get a dict of fc to sample name
col_to_fcl = geo_key_data.to_dict()['fc_lane']
col_to_sn = geo_key_data.to_dict()['sample_fullname']

geo_fc_name_to_sample_fullname = dict((k, []) for k in geo_fc_lanes)
for key in col_to_sn:
    geo_fc_name_to_sample_fullname[col_to_fcl[key]].append(col_to_sn[key])

# get a dict of individual to sample_fullname
geo_individual_to_sample_fullname = dict((k, []) for k in geo_individuals)
for key in geo_individual_to_sample_fullname:
    geo_individual_to_sample_fullname[key] = sorted(
        set([x for x in geo_fullnames
             if re.sub('_.*', '', x) == key]))

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

para_fc_lanes = sorted(set(para_key_data['fc_lane']))
para_fullnames = sorted(set(para_key_data['sample_fullname']))
para_individuals = sorted(set(para_key_data['sample']))

# get a dict of fc to sample name
pcol_to_fcl = para_key_data.to_dict()['fc_lane']
pcol_to_sn = para_key_data.to_dict()['sample_fullname']

para_fc_name_to_sample_fullname = dict((k, []) for k in para_fc_lanes)
for key in pcol_to_sn:
    para_fc_name_to_sample_fullname[pcol_to_fcl[key]].append(pcol_to_sn[key])

# get a dict of individual to sample_fullname
para_individual_to_sample_fullname = dict((k, []) for k in para_individuals)
for key in para_individual_to_sample_fullname:
    para_individual_to_sample_fullname[key] = sorted(
        set([x for x in para_fullnames
             if key in x]))

# all indivs
all_fc_lanes = sorted(set(para_fc_lanes + geo_fc_lanes))
fc_name_to_sample_fullname = {}
for key in geo_fc_name_to_sample_fullname:
    fc_name_to_sample_fullname[key] = geo_fc_name_to_sample_fullname[key]

for key in para_fc_name_to_sample_fullname:
    fc_name_to_sample_fullname[key] = para_fc_name_to_sample_fullname[key]

individual_to_sample_fullname = {}
for key in geo_individual_to_sample_fullname:
    new_key = f'geo_{key}'
    individual_to_sample_fullname[new_key] = geo_individual_to_sample_fullname[key]

for key in para_individual_to_sample_fullname:
    new_key = f'geo_{key}'
    individual_to_sample_fullname[new_key] = para_individual_to_sample_fullname[key]

all_individuals = sorted(set(individual_to_sample_fullname.keys()))


#########
# RULES #
#########

rule target:
    input:
        'output/000_config/filtered_population_map.txt'

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
def aggregate_fullnames(wildcards):
    read_stats = []
    gc_stats = []
    for fc in all_fc_lanes:
        co = checkpoints.process_radtags.get(fc_lane=fc).output['fq']
        fq_path = Path(co, '{sample_fullname}.fq.gz').as_posix()
        fq_wc = glob_wildcards(fq_path).sample_fullname
        fc_fq = expand(fq_path, sample_fullname=fq_wc)
        for individual in fq_wc:
            read_stats.append(f'output/040_stats/reads/{individual}.txt')
            gc_stats.append(f'output/040_stats/gc_hist/{individual}.txt')
    return {'read_stats': read_stats,
            'gc_stats': gc_stats}


rule combine_individual_stats:
    input:
        # read_stats = expand('output/040_stats/reads/{individual}.txt',
        #                     individual=all_individuals),
        # gc_stats = expand('output/040_stats/gc_hist/{individual}.txt',
        #                   individual=all_individuals)
        aggregate_fullnames
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
        'output/010_demux/{fc_lane}/{sample_fullname}.fq.gz'
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


# 2. for loop per fc_lane
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
        '&> {log}'

# 2. for loop per fc_lane
# for fc_lane in all_fc_lanes:
#     print(fc_lane)
#     rule:
#         input:
#             unpack(resolve_read_file),
#             config_file = 'output/000_config/{0}.config'.format(fc_lane)
#         output:
#             expand(
#                 'output/010_demux/{sample_fullname}.fq.gz',
#                 sample_fullname=fc_name_to_sample_fullname[fc_lane])
#         log:
#             'output/logs/demux.{0}.log'.format(fc_lane)
#         threads:
#             1
#         singularity:
#             stacks_container
#         shell:
#             'process_radtags '
#             '-f {input.read_file} '
#             '-i gzfastq -y gzfastq '
#             '-b {input.config_file} '
#             '-o output/010_demux '
#             '-c -q '
#             # '-r --barcode_dist_1 1 '  # rescue barcodes
#             '--barcode_dist_1 0 '       # don't allow bc mismatches
#             # '-t 91 '                  # truncate output to 91 b
#             '-w 0.1 '                   # window: approx. 9 bases
#             '-s 15 '                    # minimum avg PHRED in window
#             '--inline_null '
#             '--renz_1 apeKI --renz_2 mspI '
#             '&> {log}'

# 1. extract per-flowcell/lane sample:barcode information
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

