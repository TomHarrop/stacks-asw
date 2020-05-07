#!/usr/bin/env python3

from pathlib import Path


def resolve_path(x):
    return Path(x).resolve().as_posix()


funannotate = 'shub://TomHarrop/funannotate-singularity:funannotate_1.7.4'
bbduk_container = 'shub://TomHarrop/seq-utils:bbmap_38.76'


annot_contigs = [
    'contig_11164',     # 3 loci , 4 snps, 23 total
    'contig_13287',     # 1 locus, 3 snps, 10 total
    'contig_18336',     # 3 loci , 3 snps, 20 total
    'contig_2677',      # 2 loci , 3 snps, 22 total
    'contig_40523',     # 2 loci , 5 snps, 26 total
    'contig_8456',      # longest contig with sig SNPs
    'contig_3920'       # the test contig that worked interactively
    ]


rule target:
    input:
        expand('output/110_annotate/{contig}/training',
               contig=annot_contigs)


# try to predict
# rule funannotate_predict:
#     input:
#         'output/020_funannotate/training/funannotate_train.transcripts.gff3',
#         fasta = ('output/010_prepare/repeatmasker/'
#                  'asw-cleaned_sorted.fasta.masked'),
#         db = 'data/fundb_20200227',
#         trinity = 'data/Trinity.fasta'
#     output:
#         'output/020_funannotate/predict_results/ASW.gff3',
#         'output/020_funannotate/predict_results/ASW.mrna-transcripts.fa'
#     params:
#         fasta = lambda wildcards, input: resolve_path(input.fasta),
#         db = lambda wildcards, input: resolve_path(input.db),
#         wd = resolve_path('output/020_funannotate')
#     log:
#         'output/logs/funannotate_predict.log'
#     threads:
#         multiprocessing.cpu_count()
#     singularity:
#         funannotate
#     shell:
#         'cp /genemark/gm_key_64 ${{HOME}}/.gm_key ; '
#         'funannotate predict '
#         '-i {params.fasta} '
#         '-s ASW '
#         '--transcript_evidence {input.trinity} '
#         '-o {params.wd} '
#         '-d {params.db} '
#         '--cpus {threads} '
#         '--augustus_species lbonariensis '
#         '--optimize_augustus '
#         '--busco_seed_species tribolium2012 '
#         '--busco_db endopterygota '
#         '--organism other '
#         '--repeats2evm '
#         '--max_intronlen 10000 '
#         '&> {log}'

# run training algorithm
rule funannotate_train:
    input:
        fasta = 'output/110_annotate/{contig}/{contig}.fa.masked',
    output:
        directory('output/110_annotate/{contig}/training')
    params:
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        wd = resolve_path('output/020_funannotate'),
    log:
        'output/logs/funannotate_train.{contig}.log'
    threads:
        workflow.cores
    singularity:
        funannotate
    shell:
        'cp /genemark/gm_key_64 ${{HOME}}/.gm_key ; '
        'funannotate train '
        '--input {params.fasta} '
        '--out {params.wd} '
        '--max_intronlen 10000 '
        '--species ASW '
        '--cpus {threads} '
        '&> {log}'

# manually mask the assembly
rule rm_mask:
    input:
        cons = 'output/110_annotate/{contig}/consensi.fa',
        fasta = 'output/110_annotate/{contig}/{contig}.fa'
    output:
        'output/110_annotate/{contig}/{contig}.fa.masked'
    params:
        wd = resolve_path('output/110_annotate/{contig}'),
        lib = lambda wildcards, input: resolve_path(input.cons),
        fasta = lambda wildcards, input: resolve_path(input.fasta)
    log:
        resolve_path('output/logs/rm_mask.{contig}.log')
    threads:
        workflow.cores
    singularity:
        funannotate
    shell:
        'cd {params.wd} || exit 1 ; '
        'RepeatMasker '
        '-engine ncbi '
        '-pa {threads} '
        '-lib {params.lib} '
        '-dir {params.wd} '
        '-gccalc -xsmall -gff -html '
        '{params.fasta} '
        '&> {log}'


rule rm_model:
    input:
        'output/110_annotate/{contig}/{contig}.translation'
    output:
        'output/110_annotate/{contig}/families.stk',
        'output/110_annotate/{contig}/consensi.fa'
    params:
        wd = resolve_path('output/110_annotate/{contig}'),
    log:
        resolve_path('output/logs/rm_model.{contig}.log')
    threads:
        workflow.cores
    singularity:
        funannotate
    shell:
        'cd {params.wd} || exit 1 ; '
        'RepeatModeler '
        '-database {wildcards.contig} '
        '-engine ncbi '
        '-pa {threads} '
        '-dir {params.wd} '
        # '-recoverDir {params.wd} '
        '&> {log}'

rule rm_build:
    input:
        fasta = 'output/110_annotate/{contig}/{contig}.fa'
    output:
        'output/110_annotate/{contig}/{contig}.translation'
    params:
        wd = resolve_path('output/110_annotate/{contig}'),
        fasta = '{contig}.fa'
    log:
        resolve_path('output/logs/rm_build.{contig}.log')
    threads:
        workflow.cores
    singularity:
        funannotate
    shell:
        'cd {params.wd} || exit 1 ; '
        'BuildDatabase '
        '-name {wildcards.contig} '
        '-engine ncbi '
        '-dir {params.wd} '
        '&> {log} '


rule get_contig:
    input:
        'data/draft_genome.fasta'
    output:
        'output/110_annotate/{contig}/{contig}.fa'
    log:
        'output/logs/get_contig.{contig}.log'
    singularity:
        bbduk_container
    shell:
        'filterbyname.sh '
        'in={input} '
        'names={wildcards.contig} '
        'include=t '
        'out={output} '
        '2> {log}'
