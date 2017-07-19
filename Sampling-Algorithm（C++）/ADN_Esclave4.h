#ifndef ADNESCLAVE4_H
#define ADNESCLAVE4_H

double get(double **dist,int pattern, int *liste,int n,int *weight);

double ComparisonPartitionNuees(char *output,double **data,int n,int p,double **prototype,int *weightnnd,
								int n_proto,int *liste_instance_init,
								int n_clusterMax,double (*fdist)(double*,double*,int),int n_REP);

//double Silhouette_index(double **data,int n,int p,int n_cluster,int **liste,int *weight,int *weight_proto);

void GetWeight(int n_cluster, int *weight,int **liste,int *liste_proto,int *weight_proto,int n_proto,int **liste_weight);

//double Ball_index(double **data,int n,int p,int n_cluster,int **liste,int *weight);

int C_NueesDynamiquesWeighted(int N_iter,double Epsilon,double **data,int Nb_data,int *poids,int p,
							 double **cluster,int Nb_clust,int MinOccurence, 
							 double (*f_dist)(double *,double  *,int n), int *Cluster_weight,int Alea,int **);

int GenerePartitionFromCluster(double **pattern,int n,int p,double **cluster,int N_instance,
							   double (*fdist)(double  *,double  *,int n),
							   int **liste,int *weight,int *ambigu,double T);

/*
double ComparisonRANDSCAN(double **data,int n,int p,double **prototype,int *weightnnd,int n_proto,
								int n_cluster,double (*fdist)(double*,double*,int),int n_REP,
								double *Feat_data,double *Feat_proto,double dmax);

int GenerePartitionFromCluster_VlisteAvecFile(double **data,double **pattern,int n,int p,int **liste_in,int *weight_in,int N_cluster,int N_instance,
									  double (*fdist)(double  *,double  *,int n),
									  int **liste,int *weight,int *ambigu,double T);

void GenereListeFromCategory(int n,int *cat,int n_cluster,int **liste,int *weight);
int GenerePartitionFromCluster_Hierarchical(double **data,double **pattern,int n,int p,int **liste_in,int *weight_in,int N_cluster,int N_instance,
											double (*fdist)(double  *,double  *,int n),
											int **liste,int *weight,int *ambigu,double T);

											*/
#endif
