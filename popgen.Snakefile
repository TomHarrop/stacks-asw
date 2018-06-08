

# 13. re-run populations for PCA
rule populations_pca:
    input:
        'output/stacks_denovo/gstacks.fa.gz',
        'output/stacks_denovo/gstacks.vcf.gz',
        map = filtered_popmap
    output:
        'output/stacks_populations/for_pca/populations.sumstats_summary.tsv',
        'output/stacks_populations/for_pca/populations.markers.tsv',
        'output/stacks_populations/for_pca/populations.hapstats.tsv',
        'output/stacks_populations/for_pca/populations.sumstats.tsv',
        'output/stacks_populations/for_pca/populations.haplotypes.tsv',
        'output/stacks_populations/for_pca/populations.snps.genepop',
        'output/stacks_populations/for_pca/populations.snps.vcf'

    params:
        stacks_dir = 'output/stacks_denovo',
        outdir = 'output/stacks_populations/for_pca'
    threads:
        50
    log:
        'output/logs/populations_for_pca.log'
    singularity:
        stacks_container
    shell:
        'populations '
        '-P {params.stacks_dir} '
        '-M {input.map} '
        '-O {params.outdir} '
        '-t {threads} '
        '-r 0.8 -p 12 --min_maf 0.05 --max_obs_het 0.70 '
        '--genepop --vcf '
        '&> {log}'

# 12c. combine loci/SNP stats
rule combine_population_stats:
    input:
        expand('output/run_stats/population_stats/{r}.csv',
               r=r_values)
    output:
        combined = 'output/run_stats/population_stats_combined.csv'
    script:
        'src/combine_csvs.R'

# 12b. per-filter-run loci/SNP stats
rule population_stats:
    input:
        sumstats = 'output/stacks_populations/r{r}/populations.sumstats.tsv',
        hapstats = 'output/stacks_populations/r{r}/populations.hapstats.tsv'
    output:
        pop_stats = 'output/run_stats/population_stats/{r}.csv'
    log:
        log = 'output/logs/population_stats/{r}.log'
    threads:
        1
    script:
        'src/stacks_population_stats.R'

# 12a. rename the genepop file so adegenet will read it, French software FTW
rule genepop:
    input:
        expand('output/stacks_populations/r{r}/populations.snps.gen',
               r=r_values)

rule rename_genepop:
    input:
        'output/stacks_populations/r{r}/populations.snps.genepop'
    output:
        'output/stacks_populations/r{r}/populations.snps.gen'
    shell:
        'cp {input} {output}'
