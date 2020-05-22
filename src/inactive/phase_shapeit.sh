
singularity exec \
	--writable-tmpfs \
	shub://TomHarrop/variant-utils:easysfs_c2b26c5 \
	easySFS.py \
	-i output/060_popgen/populations.ns.all.vcf \
	-p output/070_populations/ns/popmap.txt \
	-o test_sfs \
	--proj 22,22




###

sed -e '/#CHROM/,$d' \
	../output/070_populations/ns/populations.snps.vcf \
	> ns_header.vcf

awk '{{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}}' \
	../output/005_ref/ref.fasta.fai \
	>> ns_header.vcf

sed -n -e '/#CHROM/,$p' \
	../output/070_populations/ns/populations.snps.vcf \
	>> ns_header.vcf

bgzip ns_header.vcf
bcftools sort ns_header.vcf.gz > ns_sorted.vcf
bgzip ns_sorted.vcf
tabix ns_sorted.vcf.gz

bcftools view ns_sorted.vcf.gz --regions contig_3920 > ns_3920.vcf





./shapeit --input-vcf ns_3920.vcf -O ns_3920_phased -T 8 --force
./shapeit -convert --input-haps ns_3920_phased --output-vcf ns_3920_phased.vcf
bgzip ns_3920_phased.vcf
tabix ns_3920_phased.vcf.gz

bcftools query -l ns_3920_phased.vcf.gz > samples.txt
grep  "geo_Lincoln" samples.txt > lincoln.txt
bcftools view -S lincoln.txt -O z -o lincoln.vcf.gz ns_3920_phased.vcf.gz

# try with bcftools filtered VCF (more markers)
cp ../output/060_popgen/populations.ns.all.vcf ./
bgzip populations.ns.all.vcf
tabix populations.ns.all.vcf.gz

grep "North" ../output/070_populations/ns/popmap.txt \
	| cut -f 1 \
	> North.txt

grep "South" ../output/070_populations/ns/popmap.txt \
	| cut -f 1 \
	> South.txt

bcftools view \
	--regions contig_3920 \
	-S North.txt \
	populations.ns.all.vcf.gz > ns_3920_north.vcf

bcftools view \
	--regions contig_3920 \
	-S South.txt \
	populations.ns.all.vcf.gz > ns_3920_south.vcf


./shapeit --input-vcf ns_3920_north.vcf -O ns_3920_north_phased -T 8 --force
./shapeit -convert --input-haps ns_3920_north_phased --output-vcf ns_3920_north_phased.vcf
bgzip ns_3920_north_phased.vcf
tabix ns_3920_north_phased.vcf.gz

./shapeit --input-vcf ns_3920_south.vcf -O ns_3920_south_phased -T 8 --force
./shapeit -convert --input-haps ns_3920_south_phased --output-vcf ns_3920_south_phased.vcf
bgzip ns_3920_south_phased.vcf
tabix ns_3920_south_phased.vcf.gz


###

pgdspider






