#!/usr/bin/env python3

###########
# GLOBALS #
###########

bioc_container = ('shub://TomHarrop/singularity-containers:'
                  'bioconductor_3.9')
plink_container = 'shub://TomHarrop/singularity-containers:plink_1.09beta5'
stacks2beta_container = ('shub://TomHarrop/singularity-containers:stacks_2.0beta9'
                         '@bb2f9183318871f6228b51104056a2d0')

#########
# RULES #
#########

subworkflow stacks:
    snakefile:
        'stacks.Snakefile'

rule target:
    input:
        #'output/popgen/dapc.pdf'
        'output/popgen/stacks_populations/populations.snps.vcf'

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

rule populations:
    input:
        'output/stacks_denovo/catalog.fa.gz',
        'output/stacks_denovo/catalog.calls',
        popmap = 'output/popgen/popmap.txt',
        whitelist = 'output/popgen/whitelist.txt'
    output:
        'output/popgen/stacks_populations/populations.snps.vcf'
    params:
        stacks_dir = 'output/stacks_denovo',
        outdir = 'output/popgen/stacks_populations'
    threads:
        75
    log:
        'output/logs/popgen/stacks_populations.log'
    singularity:
        stacks2beta_container
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-O {params.outdir} '
        '-W {input.whitelist} '
        '-t {threads} '
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
        plink = 'output/popgen/plink.raw'
    output:
        whitelist = 'output/popgen/whitelist.txt',
        popmap = 'output/popgen/popmap.txt'
    log:
        'output/logs/generate_whitelist.log'
    singularity:
        bioc_container
    script:
        'src/generate_whitelist.R'

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
        '--allow-no-sex --allow-extra-chr --1 '
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
