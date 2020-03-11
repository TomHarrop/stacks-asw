#!/usr/bin/env python3

from pathlib import Path


#############
# FUNCTIONS #
#############

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

bioc_container = ('shub://TomHarrop/singularity-containers:'
                  'bioconductor_3.9')
plink_container = 'shub://TomHarrop/singularity-containers:plink_1.09beta5'
stacks2beta_container = ('shub://TomHarrop/singularity-containers:stacks_2.0beta9'
                         '@bb2f9183318871f6228b51104056a2d0')
r_container = ('shub://TomHarrop/'
               'singularity-containers:r_3.6.0')

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
        'output/popgen/{mapped}/plink.raw'
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
        outdir = 'output/popgen/{mapped}/stacks_populations'
    log:
        'output/logs/popgen/stacks_populations.{mapped}.log'
    singularity:
        stacks2beta_container
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-O {params.outdir} '
        '-W {input.whitelist} '
        '-r 0 '
        '--genepop '
        '--plink '
        '--vcf '
        '--hwe '
        '--fstats '
        '--fasta_loci '
        '&> {log}'

rule generate_whitelist:
    input:
        plink = 'output/popgen/{mapped}/plink.raw'
    output:
        whitelist = 'output/popgen/{mapped}/whitelist.txt',
        popmap = 'output/popgen/{mapped}/popmap.txt'
    log:
        'output/logs/generate_whitelist.{mapped}.log'
    singularity:
        bioc_container
    script:
        'src/generate_whitelist.R'

rule convert_to_plink:
    input:
        'output/popgen/{mapped}/snps.ped',
        'output/popgen/{mapped}/snps.map'
    output:
        'output/popgen/{mapped}/plink.raw'
    params:
        workdir = 'output/{mapped}/popgen'
    threads:
        1
    singularity:
        plink_container
    shell:
        'cd {params.workdir} || exit 1 ; '
        'plink '
        '--ped snps.ped '
        '--map snps.map '
        '--recode A '
        '--allow-no-sex --allow-extra-chr --1 '
        '&> plink_log.txt'

rule filter_snps:
    input:
        'output/popgen/{mapped}/snps.gds'
    output:
        'output/popgen/{mapped}/snps.ped',
        'output/popgen/{mapped}/snps.map'
    params:
        maf = 0.05,
        missing_rate = 0.2,
        sample_missing_quantile = 0.8,
        ped_file = 'output/popgen/{mapped}/snps'
    threads:
        1
    log:
        'output/logs/filter_snps.{mapped}.log'
    singularity:
        bioc_container
    script:
        'src/filter_snps.R'

rule convert_to_gds:
    input:
        stacks('output/stacks_populations/{mapped}/r0/populations.snps.vcf')
    output:
        'output/popgen/{mapped}/snps.gds'
    threads:
        1
    log:
        'output/logs/convert_to_gds.{mapped}.log'
    singularity:
        bioc_container
    script:
        'src/convert_to_gds.R'
