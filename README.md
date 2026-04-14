# Usage of netsim.cpp
SIMULATION OF BIPARTITE GENE SHARING NETWORKS

Compilation:
g++ O3 netsim.cpp -o netsim

Input:
./netsim <alpha> <beta> <del> <ngene> <ngenome> <NRep>

Parameters:
- alpha: probability of functional innovation (adding a new gene to the network)
- beta: probability of emergence of a new species (adding a new genome to the network)
- del: deletion probability (epsilon)
- ngene: minimum number of genes to be reached before the simulation stops
- ngenome: minnimum number of genomes to be reached before the simulation stops
- NRep: number of replicates (independent simulations)

Each simulation stops when the number of genes and genomes are equal or greater than the provided values (an upper bound to the number of generations, MAXITER, is implemented just in case). The number of replicates is "NRep".

Results are recorded at two times:
 1. When the first node class (gene or genome) reaches its required value
 2. When both node classes have reached their required value

Besides the degree distributions, the program returns in stdout the following network-level statistics (one for each network):
- Time at which genes and genomes reached the required value
- Num genomes
- Num genes
- Num links
- Average genome degree
- Average gene degree
- Nestedness
In order to store that information, redirect stdout to a file.

NOTE: Modify MAXGENOME and MAXGEN before compilation according to the available memory and the required number of genes and genomes

WARNING:
- The nestedness calculation has been disabled because it is very slow. To calculate it, uncomment the relevant lines [l. 494, 512, and 532] in ProcessSingleRunOutput() before compilation.
- The code allows writing the final adjacency matrices (calls to function WriteAdjacencyMatrix() in main). This option is disabled but it can be uncommented if desired before compilation [lines 165 and 181]. Note that one file (one matrix) will be written for each simulation. For large NRep and network sizes, this means a lot of space.
- To obtain the nestedness, it is much more efficient to write the adjacency matrices and use matlab or other matrix-oriented software.
