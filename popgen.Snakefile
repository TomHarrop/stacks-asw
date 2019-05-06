#!/usr/bin/env python3

###########
# GLOBALS #
###########

bioc_container = ('shub://TomHarrop/singularity-containers:bioconductor_3.7')
plink_container = 'shub://TomHarrop/singularity-containers:plink_1.09beta5'

#########
# RULES #
#########

subworkflow stacks:
    snakefile:
        'stacks.Snakefile'

rule target:
    input:
        'output/popgen/dapc.pdf'

rule dapc:
    input:
        'output/popgen/plink.raw'
    output:
        dapc_plot = 'output/popgen/dapc.pdf',
        pca_plot = 'output/popgen/pca.pdf',
        dapc_xv = 'output/popgen/dapc_xv.Rds'
    threads:
        1
    log:
        'output/logs/dapc.log'
    singularity:
        bioc_container
    script:
        'src/dapc.R'

rule convert_to_plink:
    input:
        'output/popgen/snps.ped',
        'output/popgen/snps.map'
    output:
        'output/popgen/plink.raw'
    params:
        workdir = 'output/popgen'
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
        '--aec '
        '&> plink_log.txt'

rule filter_snps:
    input:
        'output/popgen/snps.gds'
    output:
        'output/popgen/snps.ped',
        'output/popgen/snps.map'
    params:
        maf = 0.05,
        missing_rate = 0.2,
        sample_missing_quantile = 0.8,
        ped_file = 'output/popgen/snps'
    threads:
        1
    log:
        'output/logs/filter_snps.log'
    singularity:
        bioc_container
    script:
        'src/filter_snps.R'

rule convert_to_gds:
    input:
        stacks('output/stacks_populations/r0/populations.snps.vcf')
    output:
        'output/popgen/snps.gds'
    threads:
        1
    log:
        'output/logs/convert_to_gds.log'
    singularity:
        bioc_container
    script:
        'src/convert_to_gds.R'
