#!/usr/bin/env python3

import pandas
from pathlib import Path
import tempfile


#############
# FUNCTIONS #
#############

def aggregate_pops(wildcards):
    # get 'output/100_ehh/{popset}_pops'
    co = checkpoints.get_pop_indivs.get(
        popset=wildcards.popset).output[0]
    pop_path = Path(co, '{pop}.txt').as_posix()
    pops = glob_wildcards(pop_path).pop
    vcf_dict = {}
    for pop in pops:
        my_vcf_path = (
            'output/100_ehh/'
            f'{wildcards.popset}.{wildcards.pruned}/'
            f'{pop}.{{contigs}}.phased.vcf.gz')
        vcf_dict[pop] = snakemake.io.expand(
            my_vcf_path,
            contigs=longish_chrom)
    return vcf_dict


def resolve_path(path):
    return(Path(path).resolve().as_posix())


###########
# GLOBALS #
###########

bayescan = 'shub://TomHarrop/variant-utils:bayescan_2.1'
bioc_container = 'shub://TomHarrop/r-containers:bioconductor_3.11'
biopython = 'shub://TomHarrop/singularity-containers:biopython_1.73'
easysfs = 'shub://TomHarrop/variant-utils:easysfs_c2b26c5'
pgdspider = 'shub://MarissaLL/singularity-containers:pgdspider_2.1.1.5'
plink = 'shub://MarissaLL/singularity-containers:plink_1.9'
r_container = 'shub://TomHarrop/r-containers:r_3.6.3'
samtools = 'shub://TomHarrop/align-utils:samtools_1.9'
shapeit = 'shub://TomHarrop/variant-utils:shapeit_v2.r904'
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

# not complete, just testing
bayescan_sig_contigs = [
    'contig_11164',     # 3 loci , 4 snps, 23 total
    # 'contig_12006',   # 1 locus, 3 snps, 10 total (doesn't phase)
    'contig_13287',     # 1 locus, 3 snps, 10 total
    'contig_18336',     # 3 loci , 3 snps, 20 total
    'contig_2677',      # 2 loci , 3 snps, 22 total
    # 'contig_39072',   # 1 locus, 3 snps, 10 total (doesn't phase)
    'contig_40523',     # 2 loci , 5 snps, 26 total
    # 'scaffold_43207', # 1 locus, 4 snps, 11 total (doesn't phase)
    'contig_8456',      # longest contig with sig SNPs
    'contig_3920'       # the test contig that worked interactively
    ]

# chromosomes with > 20 SNPs for "phasing"
loci_per_chrom = 'data/loci_per_chrom.csv'
loci_chrom = pandas.read_csv(loci_per_chrom)
longish_chrom = sorted(set(loci_chrom.loc[loci_chrom['N'] > 20]['chrom']))

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
        # expand('output/070_populations/{popset}/populations.snps.vcf',
        #        popset=['geo', 'ns', 'para']),
        expand('output/060_popgen/dapc.{popset}.{pruned}.pdf',
               popset=['geo', 'para', 'rlpara'],
               pruned=['all', 'pruned']),
        # expand('output/070_populations/{popset}/fst_plot.pdf',
        #        popset=['geo', 'para']),
        expand('output/080_bayescan/{popset}.{pruned}/bayescan_qvals.pdf',
               popset=['geo', 'ns', 'para', 'rlpara', 'lpara', 'rpara'],
               pruned=['all', 'pruned']),
        expand(('output/090_demographics/{popset}.{pruned}/sfs/'
                'fastsimcoal2/populations_MSFS.obs'),
               popset=['ns'],
               pruned=['all', 'pruned']),
        'output/100_ehh/ns.all/xpehh.csv'    # not enough SNPs to phase pruned

# run extended haplotype homozygosity 
rule run_rehh:
    input:
        unpack(aggregate_pops),
        fai = 'output/005_ref/ref.fasta.fai'
    output:
        xpehh = 'output/100_ehh/{popset}.{pruned}/xpehh.csv',
        pdf = 'output/100_ehh/{popset}.{pruned}/xpehh.pdf'
    log:
        'output/logs/run_rehh.{popset}.{pruned}.log'
    singularity:
        bioc_container
    script:
        'src/run_rehh.R'

rule shapeit_haps:
    input:
        vcf = 'output/100_ehh/{popset}.{pruned}/{pop}.{contig}.vcf'
    output:
        haps = 'output/100_ehh/{popset}.{pruned}/{pop}.{contig}.haps',
        vcf = 'output/100_ehh/{popset}.{pruned}/{pop}.{contig}.phased.vcf'
    params:
        vcf = lambda wildcards, input: resolve_path(input.vcf),
        wd = 'output/100_ehh/{popset}.{pruned}',
    log:
        resolve_path(
            'output/logs/shapeit_haps.{popset}.{pruned}.{pop}.{contig}.log')
    singularity:
        shapeit
    shell:
        'cd {params.wd} || exit 1 ; '
        'shapeit '
        '--input-vcf {params.vcf} '
        '-O {wildcards.pop}.{wildcards.contig} '
        '-T {threads} '
        '--force '
        '&> {log} ; '
        'shapeit '
        '-convert '
        '--input-haps {wildcards.pop}.{wildcards.contig} '
        '--output-vcf {wildcards.pop}.{wildcards.contig}.phased.vcf '
        '&>> {log}'


rule pop_vcf:
    input:
        vcf = 'output/060_popgen/populations.{popset}.{pruned}.vcf.gz',
        popmap = 'output/100_ehh/{popset}_pops/{pop}.txt'
    output:
        'output/100_ehh/{popset}.{pruned}/{pop}.{contig}.vcf'
    log:
        'output/logs/pop_vcf.{popset}.{pruned}.{pop}.{contig}.log'
    singularity:
        samtools
    shell:
        'bcftools concat '
        '--rm-dups snps '
        '-a '
        '{input.vcf} '
        ' | '
        'bcftools view '
        '--regions {wildcards.contig} '
        '-S <( cut -f1 {input.popmap} ) '
        '- '
        '> {output} '
        '2> {log}'

checkpoint get_pop_indivs:
    input:
        popmap = 'output/070_populations/{popset}/popmap.txt',
    output:
        outdir = directory('output/100_ehh/{popset}_pops')
    log:
        'output/logs/get_pop_indivs.{popset}.log'
    singularity:
        r_container
    script:
        'src/get_pop_indivs.R'

# run fastsimcoal2 on the north-south populations
rule generate_sfs:
    input:
        vcf = 'output/060_popgen/populations.{popset}.{pruned}.vcf',
        popmap = 'output/070_populations/{popset}/popmap.txt',
        proj = 'output/090_demographics/{popset}.{pruned}/proj.txt'
    output:
        ('output/090_demographics/{popset}.{pruned}/sfs/'
         'fastsimcoal2/populations_MSFS.obs')
    params:
        wd = 'output/090_demographics/{popset}.{pruned}/sfs'
    log:
        'output/logs/generate_sfs.{popset}.{pruned}.log'
    singularity:
        easysfs
    shell:
        'easySFS.py '
        '-i {input.vcf} '
        '-p {input.popmap} '
        '-o {params.wd} '
        '-f '
        '--proj "$(cat {input.proj})" '
        '&> {log}'

rule get_best_proj:
    input:
        preview = 'output/090_demographics/{popset}.{pruned}/preview.txt'
    output:
        proj = 'output/090_demographics/{popset}.{pruned}/proj.txt'
    log:
        'output/logs/get_best_proj.{popset}.{pruned}.log'
    singularity:
        bioc_container
    script:
        'src/get_best_proj.R'

rule preview_projection:
    input:
        vcf = 'output/060_popgen/populations.{popset}.{pruned}.vcf',
        popmap = 'output/070_populations/{popset}/popmap.txt'
    output:
        'output/090_demographics/{popset}.{pruned}/preview.txt'
    log:
        'output/logs/preview_projection.{popset}.{pruned}.log'
    singularity:
        easysfs
    shell:
        'easySFS.py '
        '-i {input.vcf} '
        '-p {input.popmap} '
        '--preview '
        '> {output} '
        '2> {log}'

# bayescan
rule plot_bayescan:
    input:
        fst = 'output/080_bayescan/{popset}.{pruned}/bs/populations_fst.txt',
        vcf = 'output/060_popgen/populations.{popset}.{pruned}.vcf',
        fai = 'output/005_ref/ref.fasta.fai'
    output:
        pdf = 'output/080_bayescan/{popset}.{pruned}/bayescan_qvals.pdf'
    log:
        'output/logs/plot_bayescan.{popset}.{pruned}.log'
    singularity:
        r_container
    script:
        'src/plot_bayescan.R'

rule bayescan:
    input:
        geste = 'output/080_bayescan/{popset}.{pruned}/populations.geste'
    output:
        'output/080_bayescan/{popset}.{pruned}/bs/populations_fst.txt'
    params:
        outdir = 'output/080_bayescan/{popset}.{pruned}/bs',
        o = 'populations'
    log:
        'output/logs/bayescan.{popset}.{pruned}.log'
    threads:
        workflow.cores
    singularity:
        bayescan
    shell:
        'bayescan_2.1 '
        '{input.geste} '
        '-threads {threads} '
        '-od {params.outdir} '
        '-o {params.o} '
        '-pilot 5000 '
        '-nbp 20 '
        '-burn 15000 '
        '-n 30000 '
        '-thin 10 '
        '-pr_odds 500 '
        '-out_pilot '
        '-out_freq '
        '&> {log} '
        '|| true'       # nasty. but I need to see why bayescan is failing


rule convert_to_geste:
    input:
        spid = 'output/080_bayescan/{popset}.spid',
        vcf = 'output/060_popgen/populations.{popset}.{pruned}.vcf'
    output:
        geste = 'output/080_bayescan/{popset}.{pruned}/populations.geste'
    log:
        'output/logs/convert_to_geste.{popset}.{pruned}.log'
    threads:
        21
    singularity:
        pgdspider
    shell:
        'java -jar '
        '-XX:ParallelGCThreads={threads} '
        '/opt/pgdspider/PGDSpider2-cli.jar '
        '-inputfile {input.vcf} '
        '-inputformat VCF '
        '-outputfile {output.geste} '
        '-outputformat GESTE_BAYE_SCAN '
        '-spid {input.spid} '
        '&> {log}'

rule write_spid:
    input:
        spid = 'data/blank_spid.spid',
        popmap = 'output/070_populations/{popset}/popmap.txt'
    output:
        spid = 'output/080_bayescan/{popset}.spid'
    params:
        popmap = lambda wildcards, input: resolve_path(input.popmap)
    singularity:
        biopython
    script:
        'src/write_spid.py'

# plots
rule dapc:
    input:
        'output/060_popgen/populations.{popset}.{pruned}.vcf'
    output:
        dapc_plot = 'output/060_popgen/dapc.{popset}.{pruned}.pdf',
        pca_plot = 'output/060_popgen/pca.{popset}.{pruned}.pdf',
        dapc_xv = 'output/060_popgen/dapc_xv.{popset}.{pruned}.Rds'
    threads:
        1
    log:
        'output/logs/dapc.{popset}.{pruned}.log'
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

# get population stats from stacks
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
        outdir = 'output/070_populations/{popset}',
        smoothe = (
            '--fst-correction '
            '--smooth '
            '--sigma 150000 '
            '--bootstrap '
            # '--bootstrap-wl ' + input.whitelist + ' '
            '--bootstrap-reps 10000 ')
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
        '--hwe '
        '--fstats '
        '{params.smoothe} '
        '&> {log}'

# prune the vcfs
rule prune_vcf:
    input:
        vcf = 'output/060_popgen/populations.{popset}.all.vcf',
        prune = 'output/060_popgen/{popset}.prune.in'
    output:
        vcf = 'output/060_popgen/populations.{popset}.pruned.vcf'
    log:
        'output/logs/prune_vcf.{popset}.log'
    singularity:
        samtools
    shell:
        'bcftools view '
        '-i \'ID=@{input.prune}\' '
        '{input.vcf} '
        '> {output.vcf} '
        '2> {log}'


# get a set of LD-free SNPs
rule list_pruned_snps:
    input:
        vcf = 'output/060_popgen/populations.{popset}.all.vcf'
    output:
        'output/060_popgen/{popset}.prune.in'
    params:
        vcf = lambda wildcards, input: resolve_path(input.vcf),
        wd = 'output/060_popgen',
        indep = '50 10 0.1'     # 50 kb window, 10 SNPs, r2 < 0.1
    log:
        resolve_path('output/logs/list_pruned_snps.{popset}.log')
    singularity:
        plink
    shell:
        'cd {params.wd} || exit 1 ; '
        'plink '
        '--vcf {params.vcf} '
        '--double-id '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--indep-pairwise {params.indep} '
        '--out {wildcards.popset} '
        '&> {log}'

# populations is filtering out some of the SNPs that I want for bayescan /
# plots, use bcftools filter to make a VCF for that
rule filter_populations_vcf:
    input:
        vcf = 'output/060_popgen/populations.vcf.gz',
        tbi = 'output/060_popgen/populations.vcf.gz.tbi',
        popmap = 'output/070_populations/{popset}/popmap.txt'
    output:
        'output/060_popgen/populations.{popset}.all.vcf'
    params:
        min_maf = 0.05,
        f_missing = 0.2
    log:
        'output/logs/filter_populations_vcf.{popset}.log'
    singularity:
        samtools
    shell:
        'bcftools view '
        '--min-af {params.min_maf}:nonmajor '
        '--exclude "F_MISSING>{params.f_missing}" '
        '-S <( cut -f1 {input.popmap} ) '
        '{input.vcf} '
        '> {output} '
        '2> {log}'

# make separate popmaps for lincoln-ruakura para-nonpara comparisons
rule generate_rl_popmap:
    input:
        'output/070_populations/rlpara/popmap.txt'
    output:
        l = 'output/070_populations/lpara/popmap.txt',
        r = 'output/070_populations/rpara/popmap.txt'
    singularity:
        samtools
    shell:
        'grep "para_R" {input} > {output.r} ; '
        'grep "para_L" {input} > {output.l}'

rule generate_rlpara_popmap:
    input:
        'output/070_populations/para/popmap.txt'
    output:
        'output/070_populations/rlpara/popmap.txt'
    singularity:
        samtools
    shell:
        'grep -v "_poa" {input} '
        '| grep -v "invermay" '
        '| cut -f1,3 '
        '>{output}'

# use the locusfilter VCF to calculate per-indiv missingness
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

# stats for per-individual missingness filtering
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


# generic vcf index
rule generic_index_vcf:
    input:
        Path('{folder}', '{file}.vcf')
    wildcard_constraints:
        file = '(?!populations_sorted).*'
    output:
        gz = Path('{folder}', '{file}.vcf.gz'),
        tbi = Path('{folder}', '{file}.vcf.gz.tbi')
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output.gz} '
        '; '
        'tabix -p vcf {output.gz}'
