#ifndef DIDES_H
#define DIDES_H
int SamplingApproachP_DIDES(double **data,int n,int p,double Granularity,double **prototype,int &N_instance,int *weight, char *chemin,int *option);

// UtilDIDES
int FindThresholdNoise(int N_instance,int MinCard, struct Internal *I);
double FindThresholdDistance(int N_instance,double d_max,struct Internal *I);
int TraiteClusterOver(double **data,int n,double **data_int,int *Liste_active, int *weight_int,
					  int N_instance,int &pool,int Groupe_actif,double DensiteMin,struct Internal *I,ParamFCNNrule *Pg);
int TraiteClusterLess(int **Liste,int N_instance,int Groupe_actif,
					  int Nbp_Min,struct Internal *I,ParamFCNNrule *P);
int ModifieContenance_Sampling_1(int **Liste,int N_instance,struct Internal *I,ParamFCNNrule *P);
int SearchNextSampling_1(double **D_cluster,int N_instance,double Seuil,int Instance_start,int *busy,
			   int *groupe,int actif,int &Nbre);
int FindConnex_Sampling_1(int N_instance,struct Internal *I);
int ManageProbabiliteCluster(int **Liste,int N_instance,struct Internal *I,ParamFCNNrule *P);

#endif
