#!/usr/bin/env python3

from pathlib import Path


#############
# FUNCTIONS #
#############

def resolve_path(path):
    return(Path(path).resolve().as_posix())


def resolve_filter_input(wildcards):
    if wildcards.mapped == 'mapped':
        return {
            'vcf': 'output/popgen/mapped/populations.vcf.gz',
            'tbi': 'output/popgen/mapped/populations.vcf.gz.tbi'
        }
    elif wildcards.mapped == 'denovo':
        return {
            'vcf': stacks(
                'output/stacks_populations/denovo/r0/populations.snps.vcf')
        }
    else:
        raise ValueError("WTF {wildcards.mapped}")


def stacks_mapping_resovler(wildcards):
    if wildcards.mapped == 'mapped':
        return {
            'catalog': stacks('output/map_to_genome/catalog.fa.gz'),
            'calls': stacks('output/map_to_genome/catalog.calls')
        }
    elif wildcards.mapped == 'denovo':
        return {
            'catalog': stacks('output/stacks_denovo/catalog.fa.gz'),
            'calls': stacks('output/stacks_denovo/catalog.calls')
        }
    else:
        raise ValueError("WTF {wildcards.mapped}")


###########
# GLOBALS #
###########

bioc_container = ('shub://TomHarrop/r-containers:bioconductor_3.10'
                  '@22b77812ec8211c7bbe29c9bbfc6dfba6a833982')
r_container = ('shub://TomHarrop/singularity-containers:r_3.6.0')
samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'
stacks252 = 'shub://TomHarrop/variant-utils:stacks_2.52'
# stacks2beta_container = ('shub://TomHarrop/'
#                          'singularity-containers:stacks_2.0beta9'
#                          '@bb2f9183318871f6228b51104056a2d0')
vcftools_container = 'shub://TomHarrop/variant-utils:vcftools_0.1.16'

# dict of extensions and arguments for vcftools
ext_to_arg = {
    'frq': 'freq2 --max-alleles 2',
    'idepth': 'depth',
    'ldepth.mean': 'site-mean-depth',
    'lqual': 'site-quality',
    'imiss': 'missing-indv',
    'lmiss': 'missing-site'}


#########
# RULES #
#########

subworkflow stacks:
    snakefile:
        'stacks.Snakefile'

subworkflow process_reads:
    snakefile: 'process_reads.Snakefile'


rule target:
    input:
        expand('output/popgen/{mapped}/dapc.pdf',
               mapped=['denovo', 'mapped']),
        expand('output/popgen/{mapped}/stacks_populations/populations.snps.vcf',
               mapped=['denovo', 'mapped']),
        expand('output/popgen/{mapped}/stacks_populations/fst_plot.pdf',
               mapped=['denovo', 'mapped'])

rule dapc:
    input:
        'output/popgen/{mapped}/stacks_populations/populations.snps.vcf'
    output:
        dapc_plot = 'output/popgen/{mapped}/dapc.pdf',
        pca_plot = 'output/popgen/{mapped}/pca.pdf',
        dapc_xv = 'output/popgen/{mapped}/dapc_xv.Rds'
    threads:
        1
    log:
        'output/logs/dapc.{mapped}.log'
    singularity:
        bioc_container
    script:
        'src/dapc.R'

rule plot_fst:
    input:
        fst = 'output/popgen/{mapped}/stacks_populations/populations.fst_summary.tsv'
    output:
        plot = 'output/popgen/{mapped}/stacks_populations/fst_plot.pdf'
    log:
        'output/logs/plot_fst.{mapped}.log'
    singularity:
        r_container
    script:
        'src/plot_fst.R'

rule populations:
    input:
        unpack(stacks_mapping_resovler),
        map = process_reads(
            'output/stacks_config/filtered_population_map.txt'),
        popmap = 'output/popgen/{mapped}/popmap.txt',
        whitelist = 'output/popgen/{mapped}/whitelist.txt'
    output:
        'output/popgen/{mapped}/stacks_populations/populations.snps.vcf',
        'output/popgen/{mapped}/stacks_populations/populations.fst_summary.tsv'
    params:
        stacks_dir = lambda wildcards, input:
            Path(input.catalog).parent,
        outdir = 'output/popgen/{mapped}/stacks_populations',
        smoothe = lambda wildcards, input:
            ('--fst-correction '
             '--smooth '
             '--bootstrap '
             '--bootstrap-wl ' + input.whitelist + ' '
             '--bootstrap-reps 1000 ') if wildcards.mapped == 'mapped' else ' '
    log:
        'output/logs/popgen/stacks_populations.{mapped}.log'
    singularity:
        stacks252
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-O {params.outdir} '
        '-W {input.whitelist} '
        '-r 0 '
        # '--vcf '
        # '--hwe '
        '--fstats '
        '{params.smoothe} '
        '&> {log}'

rule generate_whitelist:
    input:
        vcf = 'output/popgen/{mapped}/locusfilter.vcf',
        imiss = 'output/popgen/{mapped}/stats_locusfilter.imiss'
    params:
        imiss_rate = 0.2
    output:
        whitelist = 'output/popgen/{mapped}/whitelist.txt',
        popmap = 'output/popgen/{mapped}/popmap.txt'
    log:
        'output/logs/generate_whitelist.{mapped}.log'
    singularity:
        bioc_container
    script:
        'src/generate_whitelist.R'

# stats for filtering
rule stats_postfilter:
    input:
        vcf = 'output/popgen/{mapped}/locusfilter.vcf'
    output:
        'output/popgen/{mapped}/stats_locusfilter.{ext}'
    log:
        'output/logs/stats_postfilter_{ext}.{mapped}.log'
    params:
        wd = 'output/popgen/{mapped}',
        arg = lambda wildcards: ext_to_arg[wildcards.ext],
        vcf = lambda wildcards, input: resolve_path(input.vcf)
    singularity:
        vcftools_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'vcftools '
        '--vcf {params.vcf} '
        '--{params.arg} '
        '--out stats_locusfilter '
        '2> ' + resolve_path('{log}')

rule locusfilter:
    input:
        unpack(resolve_filter_input)
    output:
        'output/popgen/{mapped}/locusfilter.vcf'
    params:
        max_alleles = 2,
        min_maf = 0.05,
        f_missing = 0.2
    log:
        'output/logs/locusfilter.{mapped}.log'
    singularity:
        samtools
    shell:
        'bcftools view '
        '--max-alleles {params.max_alleles} '
        '--min-af {params.min_maf}:nonmajor '
        '--exclude "F_MISSING>{params.f_missing}" '
        '{input.vcf} '
        '> {output} '
        '2> {log}'

rule index_vcf:
    input:
        'output/popgen/{mapped}/populations.vcf.gz'
    output:
        'output/popgen/{mapped}/populations.vcf.gz.tbi'
    singularity:
        samtools
    shell:
        'tabix -p vcf {input}'

rule bgzip_vcf:
    input:
        'output/popgen/{mapped}/populations_sorted.vcf'
    output:
        'output/popgen/{mapped}/populations.vcf.gz'
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output}'

rule sort_vcf:      # segfaults on some computers
    input:
        'output/popgen/{mapped}/populations_header.vcf'
    output:
        temp('output/popgen/{mapped}/populations_sorted.vcf')
    singularity:
        samtools
    shell:
        'bcftools sort {input} > {output}'

rule add_vcf_header:
    input:
        vcf = stacks('output/stacks_populations/{mapped}/r0/populations.snps.vcf'),
        fai = 'output/map_to_genome/draft_genome.fasta.fai'
    output:
        temp('output/popgen/{mapped}/populations_header.vcf')
    singularity:
        samtools
    shell:
        'sed -e \'/#CHROM/,$d\' {input.vcf} > {output} ; '
        'awk \'{{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}}\' '
        '{input.fai} >> {output}  ; '
        'sed -n -e \'/#CHROM/,$p\' {input.vcf} >> {output}'

rule index_genome:
    input:
        'data/draft_genome.fasta'
    output:
        fa = 'output/map_to_genome/draft_genome.fasta',
        fai = 'output/map_to_genome/draft_genome.fasta.fai'
    singularity:
        samtools
    shell:
        'cp {input} {output.fa} ; samtools faidx {output.fa}'
