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
//Migration 0. On by default until the pops coalesce
0.000 MIGRATE
MIGRATE 0.000
//Migration 1
0.000 0.000 
0.000 0.000
//historical event: time, source, sink, migrants, new size rescale, growth rate, migr mat index
3 historical event
TDIV 1 0 1 DIVSCALE 0 1
ENDBOT 0 0 0 BOTSCALE 0 1 
STARTBOT 0 0 0 ANCSCALE 0 1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 2.5e-8
