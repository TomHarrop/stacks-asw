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
0
//historical event: time, source, sink, migrants, new size rescale, growth rate, migr mat index
4 historical event
NENDBOT 0 0 0 NBOTSCALE 0 0 
SENDBOT 1 1 0 SBOTSCALE 0 0 
TDIV 1 0 1 2 0 0
STARTBOT 0 0 0 ANCSCALE 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 2.5e-8

