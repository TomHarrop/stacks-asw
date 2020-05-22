#!/usr/bin/env Rscript

library(data.table)

vcf_file <- "output/060_popgen/populations.ns.all.vcf"
fai_file <- "output/005_ref/ref.fasta.fai"
bs_file <- "output/080_bayescan/ns.all/bs/populations_fst.txt"
xpehh_file <- "output/100_ehh/ns.all/xpehh.csv"

vcf_file <- snakemake@input[["vcf"]]
fai_file <- snakemake@input[["fai"]]
bs_file <- snakemake@input[["fst"]]

# get a list of loci to match with bs results
vcf_names <- c("chrom", "pos", "id", "ref", "alt", "qual")
loci <- fread(cmd = paste('grep -v "^#"', vcf_file),
              header = FALSE,
              col.names = vcf_names,
              select = 1:6)
loci[, snp_no := seq(1, .N)]

# read the FAI
fai_names <- c("Chr", "chr_length")
fai <- fread(fai_file, select = 1:2, col.names = fai_names)

# bayescan results
bs_names <- c("snp_no", "prob", "log10_PO", "qval", "alpha", "fst")
bs <- fread(bs_file,
            col.names = bs_names,
            skip = 1,
            header = FALSE)

# add locus info to bayescan results
bs_loci <- merge(bs, loci, by = "snp_no")
bs_loci[, nlog10_q := -log10(qval)]

bs_loci[, n_all := .N, by = chrom]
bs_loci[qval < 1e-2, bs_sig := .N, by = chrom]

setorder(bs_loci, -n_all, pos, na.last = TRUE)
sig_bs <- bs_loci[qval < 1e-2]

# xpehh results
xpehh <- fread(xpehh_file)
xpehh[, n_all := .N, by = CHR]
xpehh[LOGPVALUE > 4, xp_sig := .N, by = CHR]

setorder(xpehh, -n_all, POSITION, na.last = TRUE)
sig_xpehh <- xpehh[LOGPVALUE > 4]

# make tables
pn <- function(x){prettyNum(x, big.mark = ",")}
sp <- function(x){sprintf("%.2f", x)}

bs_tab <- unique(sig_bs[, .(contig = chrom, 
                            total_snps = n_all,
                            bs_snps = bs_sig,
                            bs_alpha = ifelse(bs_sig > 1,
                                              paste(sp(min(alpha)), sp(max(alpha)), sep = " -- "),
                                              sp(alpha)),
                            bs_region = ifelse(bs_sig > 1,
                                               paste(pn(min(pos)), pn(max(pos)), sep = " -- "),
                                               pn(pos))),
                        by = chrom],
                 by = "contig")
bs_tab[, chrom := NULL]

xp_tab <- unique(sig_xpehh[, .(contig = CHR, 
                               total_snps = n_all,
                               xp_snps = xp_sig,
                               xpehh = ifelse(xp_sig > 1,
                                              paste(sp(min(XPEHH_North_South)), sp(max(XPEHH_North_South)), sep = " -- "),
                                              sp(XPEHH_North_South)),
                               xp_region = ifelse(xp_sig > 1,
                                                  paste(pn(min(POSITION)), pn(max(POSITION)), sep = " -- "),
                                                  pn(POSITION))),
                           by = CHR],
                 by = "contig")
xp_tab[, CHR := NULL]

both_sig <- merge(bs_tab, xp_tab, all = TRUE)
both_sig[is.na(bs_snps), bs_snps := 0]
both_sig[is.na(xp_snps), xp_snps := 0]
both_sig[, total_sig := bs_snps + xp_snps]
both_sig[, both_sig := bs_snps > 0 & xp_snps > 0]

setorder(both_sig, -both_sig, -bs_snps, -xp_snps, na.last = TRUE)

# contig, number of sig snps, number of total snps, region covered



snp_table <- both_sig[, .(Contig = contig,
             `Total SNPs` = total_snps,
             `BayeScan SNPs` = bs_snps,
             `BayeScan region` = bs_region,
             `α` = bs_alpha,
             `XPEHH SNPs` = xp_snps,
             `XPEHH region` = xp_region,
             `XPEHH` = xpehh)]

# knitr::kable(snp_table,
#              align = paste(c("r", rep("l", 7)), collapse = ""),
#              format = "pandoc")

pander::pandoc.table(snp_table,
                     missing = "",
                     justify = paste(c("r", rep("l", 7)), collapse = ""),
                     style = "grid",
                     split.tables = Inf,
                     split.cells = c(Inf, 5, 9, Inf, Inf, 5, Inf, Inf))
