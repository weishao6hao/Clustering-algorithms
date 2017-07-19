#ifndef ADNESCLAVE3_H
#define ADNESCLAVE3_H

int SearchNextSampling_1(double **D_cluster,int N_instance,double Seuil,int Instance_start,int *busy,
			   int *groupe,int actif,int &Nbre);
int FindConnex_Sampling_1(int N_instance,struct Internal *I);
int ManageProbabiliteCluster(int **Liste,int N_instance,struct Internal *I,ParamFCNNrule *P);
/*
int SearchPartitionIndex(int pattern,int **liste, int n_partition, int *weight);
double RandIndex(int n_pattern,int **liste_1,int **liste_2,int n_partition_1,int n_partition_2,int *weight_1,int *weight_2);
double RandIndexFuzzy(int n_pattern,int **liste_1,int **liste_2,int n_partition_1,int n_partition_2,int *weight_1,
					  int *weight_2,int *ambigu_1,int *ambigu_2);
double AdjRandIndex( int n,int **liste_1,int **liste_2,int n_partition_1,int n_partition_2,int *weight_1,int *weight_2);
double  MutualInformation(int n,int **liste_1,int **liste_2,int n_partition_1,int n_partition_2,int *weight_1,int *weight_2);
double Cnp(int p,int n);
*/
#endif 
