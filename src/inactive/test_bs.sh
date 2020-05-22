
<( cut -f1 output/070_populations/ns/popmap.txt )

bcftools view \
	-S <( cut -f1 output/070_populations/ns/popmap.txt ) \
	-o test_bs/ns.vcf \
	output/060_popgen/populations.vcf.gz	

singularity exec \
pgdspider_2.1.1.5.sif \
java -jar /opt/pgdspider/PGDSpider2-cli.jar \
	-inputfile populations.snps.vcf \
	-inputformat VCF \
	-outputfile populations2.geste \
	-outputformat GESTE_BAYE_SCAN \
	-spid spid.spid

BayeScan2.1/source/bayescan_2.1 \
	populations.geste \
	-threads {threads} \
	-od bs1 \
	-o populations \
	-pilot 5000 \
	-nbp 20 \
	-burn 15000 \
	-n 30000 \
	-thin 10 \
	-pr_odds 500 \
	-out_pilot \
	-out_freq

bs1/populations_fst.txt