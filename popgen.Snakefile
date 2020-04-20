#!/usr/bin/env python3

from pathlib import Path
import tempfile


#############
# FUNCTIONS #
#############

def resolve_path(path):
    return(Path(path).resolve().as_posix())


###########
# GLOBALS #
###########

bioc_container = ('shub://TomHarrop/r-containers:bioconductor_3.10'
                  '@22b77812ec8211c7bbe29c9bbfc6dfba6a833982')
r_container = 'shub://TomHarrop/singularity-containers:r_3.6.0'
samtools = 'shub://TomHarrop/align-utils:samtools_1.9'
stacks_container = 'shub://TomHarrop/variant-utils:stacks_2.53'
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
    snakefile: 'stacks.Snakefile'
    configfile: 'config.yaml'

subworkflow process_reads:
    snakefile: 'process_reads.Snakefile'
    configfile: 'config.yaml'

rule target:
    input:
        expand('output/070_populations/{popset}/populations.snps.vcf',
               popset=['geo', 'ns', 'para']),
        expand('output/070_populations/{popset}/dapc.pdf',
               popset=['geo', 'para']),
        expand('output/070_populations/{popset}/fst_plot.pdf',
               popset=['geo', 'para']),

rule dapc:
    input:
        'output/070_populations/{popset}/populations.snps.vcf'
    output:
        dapc_plot = 'output/070_populations/{popset}/dapc.pdf',
        pca_plot = 'output/070_populations/{popset}/pca.pdf',
        dapc_xv = 'output/070_populations/{popset}/dapc_xv.Rds'
    threads:
        1
    log:
        'output/logs/dapc.{popset}.log'
    singularity:
        bioc_container
    script:
        'src/dapc.R'

rule plot_fst:
    input:
        fst = ('output/070_populations/{popset}/'
               'populations.fst_summary.tsv')
    output:
        plot = 'output/070_populations/{popset}/fst_plot.pdf'
    log:
        'output/logs/plot_fst.{popset}.log'
    singularity:
        r_container
    script:
        'src/plot_fst.R'

rule populations:
    input:
        catalog = stacks('output/050_stacks/catalog.fa.gz'),
        calls = stacks('output/050_stacks/catalog.calls'),
        map = process_reads(
            'output/000_config/filtered_population_map.txt'),
        popmap = 'output/070_populations/{popset}/popmap.txt',
        whitelist = 'output/060_popgen/whitelist.txt'
    output:
        'output/070_populations/{popset}/populations.snps.vcf',
        'output/070_populations/{popset}/populations.fst_summary.tsv'
    params:
        stacks_dir = lambda wildcards, input:
            Path(input.catalog).parent,
        outdir = 'output/060_populations/{popset}',
        smoothe = (
            '--fst-correction '
            '--smooth '
            '--sigma 150000 '
            '--bootstrap '
            # '--bootstrap-wl ' + input.whitelist + ' '
            '--bootstrap-reps 1000 ')
    log:
        'output/logs/popgen_stacks_populations.{popset}.log'
    singularity:
        stacks_container
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.popmap} '
        '-O {params.outdir} '
        '-W {input.whitelist} '
        '-r 0 '
        '--vcf '
        # '--hwe '
        '--fstats '
        # '{params.smoothe} '
        '&> {log}'

rule generate_whitelist:
    input:
        vcf = 'output/060_popgen/populations_filtered.vcf',
        imiss = 'output/060_popgen/stats_locusfilter.imiss',
        fai = 'output/005_ref/ref.fasta.fai',
        para_data = 'data/para_sample_info.tsv'
    params:
        imiss_rate = 0.2
    output:
        whitelist = 'output/060_popgen/whitelist.txt',
        geo = 'output/070_populations/geo/popmap.txt',
        ns = 'output/070_populations/ns/popmap.txt',
        para = 'output/070_populations/para/popmap.txt'
    log:
        'output/logs/generate_whitelist.log'
    singularity:
        bioc_container
    script:
        'src/generate_whitelist.R'

# stats for filtering
rule stats_postfilter:
    input:
        vcf = 'output/060_popgen/populations_filtered.vcf'
    output:
        'output/060_popgen/stats_locusfilter.{ext}'
    log:
        resolve_path('output/logs/stats_postfilter_{ext}.log')
    params:
        wd = 'output/060_popgen',
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
        '2> {log}'

# filter the populations VCF
rule locusfilter:
    input:
        vcf = 'output/060_popgen/populations.vcf.gz',
        tbi = 'output/060_popgen/populations.vcf.gz.tbi'
    output:
        'output/060_popgen/populations_filtered.vcf'
    params:
        max_alleles = 2,
        min_maf = 0.05,
        f_missing = 0.2
    log:
        'output/logs/locusfilter.log'
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
        'output/060_popgen/populations.vcf.gz'
    output:
        'output/060_popgen/populations.vcf.gz.tbi'
    singularity:
        samtools
    shell:
        'tabix -p vcf {input}'

rule bgzip_vcf:
    input:
        'output/tmp/populations_sorted.vcf'
    output:
        'output/060_popgen/populations.vcf.gz'
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output}'

rule sort_vcf:      # segfaults on some computers
    input:
        'output/tmp/populations_header.vcf'
    output:
        temp('output/tmp/populations_sorted.vcf')
    log:
        'output/logs/sort_vcf.log'
    singularity:
        samtools
    shell:
        'bcftools sort '
        '--temp-dir ' + tempfile.mkdtemp() + ' '
        '{input} '
        '> {output} '
        '2> {log}'

rule add_vcf_header:
    input:
        vcf = stacks('output/050_stacks/populations/populations.snps.vcf'),
        fai = 'output/005_ref/ref.fasta.fai'
    output:
        temp('output/tmp/populations_header.vcf')
    singularity:
        samtools
    shell:
        'sed -e \'/#CHROM/,$d\' {input.vcf} > {output} ; '
        'awk \'{{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}}\' '
        '{input.fai} >> {output}  ; '
        'sed -n -e \'/#CHROM/,$p\' {input.vcf} >> {output}'

rule index_genome:
    input:
        stacks('output/005_ref/ref.fasta')
    output:
        fai = 'output/005_ref/ref.fasta.fai'
    singularity:
        samtools
    shell:
        'samtools faidx {input}'
