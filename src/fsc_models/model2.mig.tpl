//Parameters for the coalescence simulation program : simcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NPOP1
NPOP2
//Samples sizes and samples age
49
67
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration 0. Active until the first pop goes into bottleneck (looking backwards)
0.000 MIGRATE 
MIGRATE 0.000
//Migration 1
0.000 0.000
0.000 0.000
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
3 historical event
NENDBOT 0 0 0 NSCALE 0 1
SENDBOT 1 1 0 SSCALE 0 1
TDIV 0 1 1 RESIZE 0 1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 2.5e-8
