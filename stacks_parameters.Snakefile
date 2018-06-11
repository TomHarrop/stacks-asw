import os
import pandas
import pathlib


###########
# GLOBALS #
###########

reads_dir = 'data/raw_reads'
key_file = 'data/SQ0003.txt'

#########
# SETUP #
#########

singularity_options = "-B /{}".format(pathlib.Path().resolve().parts[1])

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
        process_reads(expand('output/filtering/kept/{individual}.fq.gz',
                             individual=all_samples)),
        popmap = process_reads('output/stacks_config/filtered_populations.txt')
    output:
        'output/parameters/compare_defaults/optimised_samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/parameters',
        indir = 'output/filtering/kept'
    log:
        'output/logs/parameters/compare_defaults.log'
    shell:
        'stacks_parameters '
        '--mode compare_defaults '
        '-M 3 -m 3 -n 4 '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 3 '
        '--threads {threads} '
        '--singularity_args \"{singularity_options}\" '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '


rule optim_n:
    input:
        'output/parameters/filtering/replicate_1_popmap.txt',
        'output/parameters/stats_Mm/samplestats_combined.csv',
        process_reads(expand('output/filtering/kept/{individual}.fq.gz',
                             individual=all_samples)),
        popmap = process_reads('output/stacks_config/filtered_populations.txt')
    output:
        'output/parameters/stats_n/samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/parameters',
        indir = 'output/filtering/kept'
    log:
        'output/logs/parameters/optim_n.log'
    shell:
        'stacks_parameters '
        '--mode optim_n '
        '-M 3 -m 3 '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 3 '
        '--threads {threads} '
        '--singularity_args \"{singularity_options}\" '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '


rule optim_mM:
    input:
        'output/parameters/filtering/replicate_1_popmap.txt',
        process_reads(expand('output/filtering/kept/{individual}.fq.gz',
                             individual=all_samples)),
        popmap = process_reads('output/stacks_config/filtered_populations.txt')
    output:
        'output/parameters/stats_Mm/samplestats_combined.csv'
    threads:
        50
    params:
        outdir = 'output/parameters',
        indir = 'output/filtering/kept'
    log:
        'output/logs/parameters/optim_mM.log'
    shell:
        'stacks_parameters '
        '--mode optim_Mm '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 3 '
        '--threads {threads} '
        '--singularity_args \"{singularity_options}\" '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '

rule optim_setup:
    input:
        process_reads(expand('output/filtering/kept/{individual}.fq.gz',
                             individual=all_samples)),
        popmap = process_reads('output/stacks_config/filtered_populations.txt')
    output:
        'output/parameters/filtering/replicate_1_popmap.txt'
    threads:
        50
    params:
        outdir = 'output/parameters',
        indir = 'output/filtering/kept'
    log:
        'output/logs/parameters/optim_setup.log'
    shell:
        'stacks_parameters '
        '--mode setup '
        '-o {params.outdir} '
        '--individuals 8 '
        '--replicates 3 '
        '--threads {threads} '
        '--singularity_args \"{singularity_options}\" '
        '{input.popmap} '
        '{params.indir} '
        '&> {log} '
