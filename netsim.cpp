/* 07/04/26
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

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define MAXGEN 25000
#define MAXGENOME 20000
#define MAXITER 1000000

#define PI 3.141592654
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


//long seed = -1757948493;
long seed = - (long)time(NULL);	// Para los numeros aleatorios
long *idum = &seed;
double ran1(long *idum);



// Main functions
int SelectNodePrefAttach(int ddg[],int nl);			// Select a node with prob. proportional to its degree
int SelectNodeUnif(int ng);						// Select a node randomly
int CopyGenome(int adj[][MAXGENOME],int ddgen[],int ddgenome[],int idold,int idnew);	// Copy genome content to a new genome, return number of links made
double ProbSpeciation(int ki,double ps,double k0);	// Calculate speciation probability for a given genome
void SelectRandomLink(int *ig,int *iG,int adj[][MAXGENOME],int ddg[],int nl);	// Select indices of a random link, and assign to pointers ig and IG
void RemoveGen(int ig,int *ng_pt,int nG,int adj[][MAXGENOME],int ddg[]);		// Remove empty gene, replacing its content with the last non-empty gene and making Ngen--
void RemoveGenome(int iG,int *nG_pt,int ng,int adj[][MAXGENOME],int ddg[]);	// Idem, for genomes
void SingleRun(int *nG_pt,int *ng_pt,int *tt_pt,int ddgenome_pt[],int ddgen_pt[],double alpha,double del,double beta,int NumGene,int NumGenome,int Tmax,int *RunFlag_pt);
void ProcessSingleRunOutput(int *nG_,int *ng_,int ddgenome_[],int ddgene_[],int nr_,int degree_counts_genome_[],int degree_counts_gene_[],int flag);
void RecordResults(int degree_counts_genome[],int degree_counts_gene[],char outfileinfo[]);
void RecordTimes(int times[],char outfileinfo[]);
double CalculateNestedness(int *nG_,int *ng_,int *nl_,int ddgenome_[],int ddgene_[]);  // Calculate nestedness
void reset_adj(int adj[][MAXGENOME]);
void WriteAdjacencyMatrix(int *nG_,int *ng_,char outfileinfo[]);


// Global variables
int adj[MAXGEN][MAXGENOME] = {0};	// Adjacency matrix


int main (int argc,char *argv[])
{
  if (argc<7)
    {
    printf("Error: provide parameters alpha, beta, epsilon, NumGene, NumGenome, NRep\n");
    return -1;
    }

   // Define variables
   double alpha = 0.035;
   double beta = 0.025;
   double del = 0.;
   int tmax = MAXITER;
   int NRep = 10;
   int nr,ngenome,ngen,nlink,NumGene,NumGenome;
   int DDgen[MAXGEN] = {0};			// Degree distribution for genes
   int DDgenome[MAXGENOME] = {0};	// Degree distribution for genomes
   int degree_counts_genome1[MAXGENOME] = {0};
   int degree_counts_gene1[MAXGEN] = {0};
   int degree_counts_genome2[MAXGENOME] = {0};
   int degree_counts_gene2[MAXGEN] = {0};
   int totlink = 0;
   char outfileinfo[60] = "";
   int tt_single = 0;
   int times[2] = {0};
   int RunFlag = 0; // 0: first run, 1: reached number of genomes, 2: reached number of genes, 3: reached number of genes and genomes, -1: error (reached MAXGEN or MAXGENOME)

   // Parameter values
   alpha = atof(argv[1]);
   beta = atof(argv[2]);
   del = atof(argv[3]);
   NumGene = atoi(argv[4]);
   NumGenome = atoi(argv[5]);
   NRep = atoi(argv[6]);

   // Check that the gene and genome goals are below the maximum size for storage
   if (NumGene >= MAXGEN)
      {
      printf("Error: the number of required genes is too high (NumGene = %d; MAXGEN = %d)\n",NumGene,MAXGEN);
      return -1;
      }
   if (NumGenome >= MAXGENOME)
      {
      printf("Error: the number of required genomes is too high (NumGenome = %d; MAXGENOME = %d)\n",NumGenome,MAXGENOME);
      return -1;
      }

   // Print parameters and result headers
   printf("Parameter values: alpha = %f, beta = %f, del = %f, NumGene = %d, NumGenome = %d, seed = %ld\n",alpha,beta,del,NumGene,NumGenome,seed);
   //printf("Iter.\tNo. genomes\tNo. genes\tNo. links\t<k_G>\t<k_g>\tNestedness\n");

   // Run repeats
   for (nr=0; nr<NRep; nr++)
      {
      RunFlag = 0;
      //times[1] = 0;
      //times[2] = 0;
      printf("Simulating network %d\n",nr+1);
      sprintf(outfileinfo,"%s_%s_%s_%s_%s_%d",argv[1],argv[2],argv[3],argv[4],argv[5],nr);
      SingleRun(&ngenome,&ngen,&tt_single,DDgenome,DDgen,alpha,del,beta,NumGene,NumGenome,tmax,&RunFlag);
      // Process results of single run
      if (RunFlag == 1)
         {
         //times[1] = tt_single;
         //printf("Collecting intermediate data\n");
         ProcessSingleRunOutput(&ngenome,&ngen,DDgenome,DDgen,nr,degree_counts_genome1,degree_counts_gene1,RunFlag);
         SingleRun(&ngenome,&ngen,&tt_single,DDgenome,DDgen,alpha,del,beta,NumGene,NumGenome,tmax,&RunFlag);
         if (RunFlag == 3)
            {
            //times[2] = tt_single;
            //printf("Collecting final data\n");
            ProcessSingleRunOutput(&ngenome,&ngen,DDgenome,DDgen,nr,degree_counts_genome2,degree_counts_gene2,RunFlag);
     		//WriteAdjacencyMatrix(&ngenome,&ngen,outfileinfo);
            }
         else if (RunFlag == -1)
            printf("skipped\n");
         }
      else if (RunFlag == 2)
         {
         //times[2] = tt_single;
         //printf("Collecting intermediate data\n");
         ProcessSingleRunOutput(&ngenome,&ngen,DDgenome,DDgen,nr,degree_counts_genome2,degree_counts_gene2,RunFlag);
         SingleRun(&ngenome,&ngen,&tt_single,DDgenome,DDgen,alpha,del,beta,NumGene,NumGenome,tmax,&RunFlag);
         if (RunFlag == 3)
            {
            //times[1] = tt_single;
            //printf("Collecting final data\n");
            ProcessSingleRunOutput(&ngenome,&ngen,DDgenome,DDgen,nr,degree_counts_genome1,degree_counts_gene1,RunFlag);
     		//WriteAdjacencyMatrix(&ngenome,&ngen,outfileinfo);
            }
         else if (RunFlag == -1)
            printf("skipped\n");
         }
      else if (RunFlag == -1)
         printf("skipped\n");
      }

   // Print degree distributions
   sprintf(outfileinfo,"%s_%s_%s_%s_%s_%s_fixG",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
   RecordResults(degree_counts_genome1,degree_counts_gene1,outfileinfo);
   sprintf(outfileinfo,"%s_%s_%s_%s_%s_%s_fixg",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
   RecordResults(degree_counts_genome2,degree_counts_gene2,outfileinfo);
   // Print time to reach NumGenome and NumGene
   //sprintf(outfileinfo,"%s_%s_%s_%s_%s_%s",argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
   //RecordTimes(times,outfileinfo);

}



// Single run

void SingleRun(int *nG_pt,int *ng_pt,int *tt_pt,int ddgenome_pt[],int ddgen_pt[],double alpha,double del,double beta,int NumGene,int NumGenome,int Tmax,int *RunFlag_pt)
{
   // Main variables
   //int adj[MAXGEN][MAXGENOME] = {0};	// Adjacency matrix
   int ddgen[MAXGEN] = {0};			// Degree distribution for genes
   int ddgenome[MAXGENOME] = {0};	// Degree distribution for genomes
   int Ngen = 1;					// Number of genes
   int Ngenome = 1;					// Number of genomes
   int Nlink = 0;
   int tt = 0;
   int idgen, idgenome, idgenome0, dgi;
   double rnum, pspec;


   // Set initial state if RunFlag = 0
   if (*RunFlag_pt == 0)
      {
      reset_adj(adj);
      adj[0][0] = 1;		// Set initial adjacency matrix
      ddgen[0] = 1;
      ddgenome[0] = 1;
      Nlink = 1;
      }
   // Recover previous state if RunFlag = 1 or 2
   else if ((*RunFlag_pt == 1) || (*RunFlag_pt == 2))
      {
      Ngenome = *nG_pt;
      Ngen = *ng_pt;
      tt = *tt_pt;
      for (idgenome=0;idgenome<Ngenome;idgenome++)
         ddgenome[idgenome] = ddgenome_pt[idgenome];
      for (idgen=0;idgen<Ngen;idgen++)
         ddgen[idgen] = ddgen_pt[idgen];
      dgi = 0;
      for (idgenome=0;idgenome<Ngenome;idgenome++)
         {
         dgi = ddgenome[idgenome];
         Nlink += dgi;
         }
      }
   // Throw an error if RunFlag has another value
   else
      {
      printf("Error: invalid value of RunFlag.\n");
      return;
      }


   while (tt < Tmax)
      {
	  printf("%d\t%d\n",tt,Nlink);
      tt++;
      if ((Ngenome >= NumGenome) && (*RunFlag_pt != 1)) // If the number of genomes is reached for the first time, stop and adjust the RunFlag
         {
         //printf("Warning: requested number of genomes reached.\n");
         if (*RunFlag_pt == 0)
            *RunFlag_pt = 1; 		// Reached number of genomes
         else if (*RunFlag_pt == 2)
            *RunFlag_pt = 3;		// Reached number of genes and genomes
         break;
         }
      if ((Ngen >= NumGene) && (*RunFlag_pt != 2)) // If the number of genes is reached for the first time, stop and adjust the RunFlag
         {
         //printf("Warning: requested number of genes reached.\n");
         if (*RunFlag_pt == 0)
            *RunFlag_pt = 2; 		// Reached number of genomes
         else if (*RunFlag_pt == 1)
            *RunFlag_pt = 3;		// Reached number of genes and genomes
         break;
         }
      if ((Ngenome >= MAXGENOME) || (Ngen >= MAXGEN))
         {
         printf("Warning: maximum number of genes or genomes reached. Simulation terminated.\n");
         *RunFlag_pt = -1;
         break;
         }

      rnum = ran1(idum);
      // 1. HGT in every time step
      // Select gene for HGT
      idgen = SelectNodePrefAttach(ddgen,Nlink);
      // 1.1. Speciate with probability beta
      if (rnum < beta)
         {
         adj[idgen][Ngenome] = 1;
         ddgen[idgen] = ddgen[idgen] + 1;
         ddgenome[Ngenome] = 1;
         Nlink++;
         Ngenome++;
         rnum = rnum/beta;
         }
      // 1.2. HGT to existing genome
      else
         {
         // Select genome for HGT
         idgenome = SelectNodeUnif(Ngenome);
         // Create a new link in the bipartite network if the link did not already exist
         if (adj[idgen][idgenome] == 0)
            {
            adj[idgen][idgenome] = 1;
            ddgen[idgen] = ddgen[idgen] + 1;
            ddgenome[idgenome] = ddgenome[idgenome] + 1;
            Nlink++;
            rnum = (rnum-beta)/(1.-beta);
            }
         }
      // 2. Add new gene with probability alpha
      if (rnum < alpha)
         {
         rnum = rnum/alpha;
         // 2.1. Speciate with probability beta
         if (rnum < beta)
            {
            adj[Ngen][Ngenome] = 1;
            ddgen[Ngen] = 1;
            ddgenome[Ngenome] = 1;
            Nlink++;
            Ngen++;
            Ngenome++;
            }
         // 2.2. Add new gene to existing genome with probability 1 - beta
         else
            {
            idgenome = SelectNodeUnif(Ngenome);
            adj[Ngen][idgenome] = 1;
            ddgen[Ngen] = 1;
            ddgenome[idgenome] = ddgenome[idgenome] + 1;
            Nlink++;
            Ngen++;
            }
         }
      // 3. Delete a random link with probability del
      rnum = ran1(idum);
      if (rnum < del)
         {
         SelectRandomLink(&idgen,&idgenome,adj,ddgen,Nlink);
         adj[idgen][idgenome] = 0;
         ddgen[idgen] = ddgen[idgen] - 1;
         ddgenome[idgenome] = ddgenome[idgenome] - 1;
         Nlink--;
         // Take care of extinct genomes and/or genes (very inefficient as it is)
         if (ddgen[idgen] == 0)
            RemoveGen(idgen,&Ngen,Ngenome,adj,ddgen);
         if (ddgenome[idgenome] == 0)
            RemoveGenome(idgenome,&Ngenome,Ngen,adj,ddgenome);
	// Reset if all genomes go extinct
         if (Ngenome ==0 && Ngen == 0)
            {
            Ngenome = 1;
            Ngen = 1;
            Nlink = 1;
            adj[0][0] = 1;
            ddgenome[0] = 1;
            ddgen[0] = 1;
            }
         }
      }

   // Return results
   *ng_pt = Ngen;
   *nG_pt = Ngenome;
   *tt_pt = tt;
   for (idgenome=0;idgenome<Ngenome;idgenome++)
      ddgenome_pt[idgenome] = ddgenome[idgenome];
   for (idgen=0;idgen<Ngen;idgen++)
      ddgen_pt[idgen] = ddgen[idgen];

   /* Print genome degrees
   for (idgenome=0;idgenome<Ngenome;idgenome++)
      printf("Genome %d\t%d\n",idgenome,ddgenome[idgenome]);
*/

   // Print parameters and overall results
   //printf("Number of genomes: %d\nNumber of genes: %d\nNumber of links: %d\n",Ngenome,Ngen,Nlink);

}



// Auxiliary functions
void WriteAdjacencyMatrix(int *nG_,int *ng_,char outfileinfo[])
{
   FILE *xdata;
   int ii,jj;
   int ngenome = *nG_;
   int ngen = *ng_;
   char name[60] ="";

   printf("Saving results...\n");
   sprintf(name,"adj_%s.dat",outfileinfo);
   xdata=fopen(name,"w");

   for (ii=0;ii<ngen;ii++)
      {
      fprintf(xdata,"%d",adj[ii][0]);	// Print first element of the row
	  for (jj=1;jj<ngenome;jj++)
	     fprintf(xdata,"\t%d",adj[ii][jj]);
      fprintf(xdata,"\n");
      }
   fclose(xdata);
}

double CalculateNestedness(int *nG_,int *ng_,int *nl_,int ddgenome_[],int ddgene_[])
{
	int ngenome = *nG_;
	int ngen = *ng_;
	int nlink = *nl_;
	long int k2genome = 0;
	long int k2gen = 0;
	double kmgenome = double(nlink)/ngenome;
	double kmgen = double(nlink)/ngen;
	int ntot = ngen + ngenome;
	double nest0,nest = 0;
	int ii,jj,hh;
	int adjprod = 0;

	// Nestedness null model
	for (ii=0;ii<ngenome;ii++)
    	k2genome += ddgenome_[ii] * ddgenome_[ii];
    for (ii=0;ii<ngen;ii++)
    	k2gen += ddgene_[ii] * ddgene_[ii];
	nest0 = (double(k2genome) * ngen / ngenome) + (double(k2gen) * ngenome / ngen); // Numerator of null nestedness
	nest0 = nest0 / (kmgen * kmgenome * ntot * ntot); 		// Divide by denominator
	if (nest0 < 0)
		printf("!!! Negative nest0 = %f\t k2genome = %ld\t k2gen = %ld\n",nest0,k2genome,k2gen);

	// Empirical nestedness
	nest = 0;
    for (ii=0;ii<ngen;ii++) 	// Block gene x gene
    {
		for (jj=0;jj<ngen;jj++)
		{
			adjprod = 0;
			for (hh=0;hh<ngenome;hh++)
				adjprod += adj[ii][hh] * adj[jj][hh];
			nest += double(adjprod) / (ddgene_[ii] * ddgene_[jj]);
			if (nest < 0)
				printf("!!! Negative nest = %f\t ntot = %d, adjproc = %d\t ii = %d\t jj = %d\t In gene x gene block\n",nest,ntot,adjprod,ii,jj);
		}
	}
    for (ii=0;ii<ngenome;ii++) 	// Block genome x genome
    {
		for (jj=0;jj<ngenome;jj++)
		{
			adjprod = 0;
			for (hh=0;hh<ngen;hh++)
				adjprod += adj[hh][ii] * adj[hh][jj];
			nest += double(adjprod) / (ddgenome_[ii] * ddgenome_[jj]);
			if (nest < 0)
				printf("!!! Negative nest = %f\t ntot = %d, adjproc = %d\t ii = %d\t jj = %d\t In genome x gene block\n",nest,ntot,adjprod,ii,jj);
		}
	}
	nest = nest / (ntot * ntot);  // Normalize
	if (nest < 0)
		printf("!!! Negative nest = %f\t ntot = %d\n",nest,ntot);

	// Return relative nestedness
	return nest/nest0;
}

void ProcessSingleRunOutput(int *nG_,int *ng_,int ddgenome_[],int ddgene_[],int nr_,int degree_counts_genome_[],int degree_counts_gene_[],int flag)
{
   int nlink = 0;
   int dgi,ii;
   int ngenome = *nG_;
   int ngen = *ng_;
   double nest = 0;

   if (flag == -1)	// Initialize degree distributions
      {
      for (ii=0;ii<ngenome;ii++)
         ddgenome_[ii] = 0;
      for (ii=0;ii<ngen;ii++)
         ddgene_[ii] = 0;
      }

   if (flag == 1)	// Save data
      {
      for (ii=0;ii<ngenome;ii++)
         {
         dgi = ddgenome_[ii];
         nlink += dgi;
         degree_counts_genome_[dgi] = degree_counts_genome_[dgi] + 1;
         }
      for (ii=0;ii<ngen;ii++)
         {
         dgi = ddgene_[ii];
         degree_counts_gene_[dgi] = degree_counts_gene_[dgi] + 1;
         }
      //nest = CalculateNestedness(&ngenome,&ngen,&nlink,ddgenome_,ddgene_);
      // Print results of single run
      printf("%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n",nr_,ngenome,ngen,nlink,double(nlink)/ngenome,double(nlink)/ngen,nest);
      }

   if (flag == 2)	// Save data
      {
      for (ii=0;ii<ngenome;ii++)
         {
         dgi = ddgenome_[ii];
         nlink += dgi;
         degree_counts_genome_[dgi] = degree_counts_genome_[dgi] + 1;
         }
      for (ii=0;ii<ngen;ii++)
         {
         dgi = ddgene_[ii];
         degree_counts_gene_[dgi] = degree_counts_gene_[dgi] + 1;
         }
      //nest = CalculateNestedness(&ngenome,&ngen,&nlink,ddgenome_,ddgene_);
      // Print results of single run
      printf("%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n",nr_,ngenome,ngen,nlink,double(nlink)/ngenome,double(nlink)/ngen,nest);
      }

   if (flag == 3)	// Save data and initialize degree distributions
      {
      for (ii=0;ii<ngenome;ii++)
         {
         dgi = ddgenome_[ii];
         //ddgenome_[ii] = 0;
         nlink += dgi;
         degree_counts_genome_[dgi] = degree_counts_genome_[dgi] + 1;
         }
      for (ii=0;ii<ngen;ii++)
         {
         dgi = ddgene_[ii];
         //ddgene_[ii] = 0;
         degree_counts_gene_[dgi] = degree_counts_gene_[dgi] + 1;
         }
      //nest = CalculateNestedness(&ngenome,&ngen,&nlink,ddgenome_,ddgene_);
      for (ii=0;ii<ngenome;ii++)
      	ddgenome_[ii] = 0;
      for (ii=0;ii<ngen;ii++)
         ddgene_[ii] = 0;
      // Print results of single run
      printf("%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n",nr_,ngenome,ngen,nlink,double(nlink)/ngenome,double(nlink)/ngen,nest);
      }
}

void RecordResults(int degree_counts_genome[],int degree_counts_gene[],char outfileinfo[])
{
   FILE *xdata;
   int dgi,ii;
   int gtot = 0;
   int Gtot = 0;
   char name[60] ="";

   printf("Saving results...\n");

   for (ii=0;ii<MAXGEN;ii++)
      gtot += degree_counts_gene[ii];
   for (ii=0;ii<MAXGENOME;ii++)
      Gtot += degree_counts_genome[ii];

   sprintf(name,"ddgenome_%s.dat",outfileinfo);
   printf("Saving genome size distribution as %s\n",name);
   xdata=fopen(name,"w");
   fprintf(xdata,"#Genome size distribution\n");
   for (ii=0;ii<MAXGENOME;ii++)
      {
      dgi = degree_counts_genome[ii];
      if (dgi > 0)
         fprintf(xdata,"%d\t%f\n",ii,(double)dgi/Gtot);
      }
   fclose(xdata);

   sprintf(name,"ddgene_%s.dat",outfileinfo);
   printf("Saving gene abundance distribution as %s\n",name);
   xdata=fopen(name,"w");
   fprintf(xdata,"#Gene abundance distribution\n");
   for (ii=0;ii<MAXGEN;ii++)
      {
      dgi = degree_counts_gene[ii];
      if (dgi > 0)
         fprintf(xdata,"%d\t%f\n",ii,(double)dgi/gtot);
      }
   fclose(xdata);
}

void RecordTimes(int times[],char outfileinfo[])
{
   FILE *xdata;
   char name[60] ="";
   sprintf(name,"times_%s.txt",outfileinfo);
   xdata=fopen(name,"w");
   fprintf(xdata,"t_G\t%d\n",times[1]);
   fprintf(xdata,"t_g\t%d\n",times[2]);
   fclose(xdata);
}

void reset_adj(int adj[][MAXGENOME])
{
   int ii,jj;
   for (ii=0;ii<MAXGEN;ii++)
      {
      for (jj=0;jj<MAXGENOME;jj++)
         adj[ii][jj]=0;
      }
}

void RemoveGen(int ig,int *ng_pt,int nG,int adj[][MAXGENOME],int ddg[])
{
int ng = *ng_pt -1;
int iG;

ddg[ig] = ddg[ng];
ddg[ng] = 0;
for (iG=0;iG<nG;iG++)
   {
   adj[ig][iG] = adj[ng][iG];
   adj[ng][iG] = 0;
   }
*ng_pt = ng;
}

void RemoveGenome(int iG,int *nG_pt,int ng,int adj[][MAXGENOME],int ddg[])
{
int nG = *nG_pt -1;
int ig;

ddg[iG] = ddg[nG];
ddg[nG] = 0;
for (ig=0;ig<ng;ig++)
   {
   adj[ig][iG] = adj[ig][nG];
   adj[ig][nG] = 0;
   }
*nG_pt = nG;
}

void SelectRandomLink(int *ig,int *iG,int adj[][MAXGENOME],int ddg[],int nl)
{
int idgen = SelectNodePrefAttach(ddg,nl);
int nlg = ddg[idgen];
int ksum = 0;
int idgenome = 0;
double rnum;

do rnum  = ran1(idum);
while (rnum >= 1);
rnum = rnum * nlg;

do
	{
	ksum += adj[idgen][idgenome];
	idgenome++;
	}
while (ksum < rnum);
idgenome--;

if (adj[idgen][idgenome] != 1)
   {
    printf("Error: trying to delete an inexistent link\n");
    exit (EXIT_FAILURE);
   }

*iG = idgenome;
*ig = idgen;
return;
}

double ProbSpeciation(int ki,double ps,double k0)
{
double pspec = ps * exp(-double(ki)/k0);
return pspec;
}

int CopyGenome(int adj[][MAXGENOME],int ddgen[],int ddgenome[],int idold,int idnew)
{
int nl = 0;
int ig = 0;

while (nl < ddgenome[idold])
   {
   if (adj[ig][idold] == 1)
      {
      adj[ig][idnew] = 1;
      ddgen[ig] = ddgen[ig] + 1;
      nl++;
      }
   ig++;
   }

if (nl != ddgenome[idold])
   {
    printf("Error: inconsistency in the number of copied links\n");
    exit (EXIT_FAILURE);
    }
ddgenome[idnew] = nl;
return nl;
}


int SelectNodeUnif(int ng)
{
double rnum;
int gid;
do rnum  = ran1(idum);
while (rnum >= 1);
rnum = rnum*ng;
gid = (int) rnum;
return gid;
}

int SelectNodePrefAttach(int ddg[],int nl)
{
double rnum;
double ksum = 0;
int ix = 0;

do rnum  = ran1(idum);
while (rnum >= 1);
rnum = rnum*nl;

do
	{
	ksum += ddg[ix];
	ix++;
	}
while (ksum < rnum);
return (ix - 1);
}

double ran1(long *idum)
{
   int j;
   long k;
   static long iy=0;
   static long iv[NTAB];
   double temp;

   if (*idum <= 0 || !iy) {
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      for (j=NTAB+7; j>=0; j--) {
         k=(*idum)/IQ;
        *idum=IA*(*idum-k*IQ)-IR+k;
         if (*idum < 0) *idum += IM;
         if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if (*idum < 0) *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j] = *idum;
   if ((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}

