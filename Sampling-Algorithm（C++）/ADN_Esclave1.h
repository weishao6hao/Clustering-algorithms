#ifndef ADNESCLAVE1_H
#define ADNESCLAVE1_H

#include <stdio.h>

#ifndef MAX_CLUSTER
#define MAX_CLUSTER 100
#endif

#ifndef MAX_VALUE
#define MAX_VALUE 10000000000
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////
struct Tri
{
int *G_actif;
int *weight;
int *occ;
int *noise;
int *ptisole;
};


struct Internal
{
	int N_instance_Extreme;
	double D_distance_Extreme;
	int *disable;
	int *Nearest_proto;
	int *Nearest_between;
	int *Pre_individu;
	int *Liste_selectionne;
	int *free;
	int *rep;
	double Max_distance;
	double Treshold_groupe;
	double *Inertie_proto;
	double **D_cluster;
	double **data_proto;
	double *dmaxPrototype;
	double *Dist_limite;
	double *Nearest;
	int *Nearest_first;
	double *Pre_distance;
	int *Order_selection;
	int *Weight;
	int WorseWinner;
	double Ratio_massique;
	double *d1nn_proto;
	double Ave_d1nn_proto;
	double *C_proto;
	double CostP[20];
	double CostS[20];
	double *C_Oldproto;
	double *Inertie_1knn;
	double *Distribution_dMax;
	int C_stable;
	int *Distribution_Weight;
	int **Liste_pattern;
	double alpha_M,beta_M;
	double alpha_m,beta_m;
	struct Tri *Tr;
	int Seuil_instance_Ratio;
	int Seuil_instance_Inertie;
	int Nb_isole, Nb_bruit;
	int Ref_OK;
	int Occmin;
	int N_groupe;
	double Ref_distance;
	int Mass_Weight_Criterion;
	int Mass_OptionFinal;
};

struct Sousgroupe
{
int n_elem;
int *liste;
};

struct HPtC
{
int actif;
int selected;
double poids;
int npts;//nbre de points.
int p;
int *liste_poids;
int *liste_pattern;//liste de points
int closest;//cluster le plus proche.
double d_closest;//la distance la plus proche
double *G;
double Masse;
};

struct ParamFCNNrule
{
	double *Cost_n;
	double *d_1knn;//liste des max.
	double *d_1knnmin;//liste des min
	double *Ave_d1nn_proto;
	double **pattern;
	int *weight_1knn;//liste des poids des max.
	double *dav_1knn;//liste des moyennes des distances
	double *dec_1knn;//liste des moyennes des écart types des distances.
	double *Wt_1knn;//liste des estimés des moyennes des distances
	double *Wtnext_1knn;//liste des estimés des écart types des distances
	double *WI1_1knn;//liste des moyennes des poids.
	double *WI2_1knn;//liste des écart types des poids.
	double *Inertie_proto;
	double *II_g;
	double *II_1knn;
	double *R_massique;
	double Seuil_massique;
	double DensiteMin;
	double Thres_d_1knn;
	double Thres_Wt_1knn;
	double Thres_EcartW;
	double Thres_SequenceW;
	double Thres_SequenceD;
	double Ratio_distance;
	double Thres_FreqW;
	double Thres_RatioKmin;
	double Thres_RatioKmax;
	double Thres_EcartN;
	int NbPointsMin;
	int NbPointsNoise;
	double OccMin;
	double V_seuilH, V_seuilB;
	double C_seuilH, C_seuilB;
	double Coef;
	int SansStep2;
	int Auto;
	int Option_Alea_step2;
	int Option_Alea_Bruit;
	int Option_Alea_Cleaning_1;
	int Option_Alea_Cleaning_2;
	int Option_Assymptotique;
	int Init_Aleatoire;
	int n,p;
	int MaxCluster;
	int L_pattern;
	int Save;
	double Epsilon_Assymptote;
	double Epsilon_R;
	int Plage_R;
	double Dist_Noise;
	double Dist_Connexity;
	int Barycentre;
	int MaxClusterMemory;
	double Nb_distance[1000];
	double Nb_TriangularTest[1000];
	double Granularity;
	int Mass_WeightCriterion;
	int Mass_OptionFinal;
	int ActifMasse;
	double TradeoffDensite;
	double FactorDistance;
	int S2_Option_all;
	double S2_Beta;
	int Ratio_achieved;
	int Tha_achieved;
	double (*fdist)(double  *,double  *,int n);
};

ParamFCNNrule *InitParamFCNNrule(int n);
int UpdateWeight_DIDES(struct Internal *I,ParamFCNNrule *P,int N_prototype,int n_all,double **data);
int TestFinDIDES(struct Internal *I,ParamFCNNrule *P,int &N_instance);
double FindFromDistributionDIDES(double *dmax,int *Weight,int N_instance,int N_total,double SeuilExtreme);
int CalculInertieProto(double **data,int N_instance,struct Internal *I,ParamFCNNrule *P);
void SortPrototypes_DIDES(struct Internal *I,ParamFCNNrule *P,int &N_instance,int **Liste_pattern,int Min);
int ModeliseFunction(struct Internal *I,ParamFCNNrule *P,int N_prototype);
int RatioMassiqueAchieved(int N_instance,struct Internal *I,ParamFCNNrule *P);
void ProcessCentreGravitePrototype(double **data,int N_instance,struct Internal *I,ParamFCNNrule *P);
void Update1knnPrototypeMasse(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance);
int TestFin_DENDIS(struct Internal *I,ParamFCNNrule *P,int &N_instance);
int FindNewInstanceMassique(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance);
int SamplingApproach_DENDIS(double **data,int *Liste_instance,int &N_instance,int *weight,
					   int **Liste,double &dmax,ParamFCNNrule *P);
int SamplingApproachP_DENDIS(double **data,int n,int p,double Granularity,double **prototype,int &N_instance,int *weight,char*chemin,int *option);
int SamplingApproachP_PROBA(double **data,int n,int p,double Granularity,double **prototype,int &N_instance,
						int *Liste_instance,int *weight, char*chemin,int *option);	
int SamplingApproach_DIDES(double **data,int *Liste_instance,int &N_instance,int *weight,int **Liste,double &dmax,ParamFCNNrule *P);
int SamplingApproachP_DIDES(double **data,int n,int p,double Granularity,double **prototype,int &N_instance,int *weight,char*chemin,int *option=NULL);
int InitProgramme(double **data,struct Internal *I,ParamFCNNrule *P);
void InitBoucle(struct Internal *I,int N_instanceOld,int New_instance);
int UpdateWithGravity(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance);
int FindNewInstance_DIDES(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance);
int FindNewInstance_PROBA(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance);
void ProcessCriteria(int N_instance,struct Internal *I,ParamFCNNrule *P);
void UpdateAlgorithm(struct Internal *I,int N_instance,int Data_Newinstance);
int TestFin(struct Internal *I,ParamFCNNrule *P,int &N_instance);
int CalculInertieProto(double **data,int N_instance,struct Internal *I,ParamFCNNrule *P);
int InertieAchieved(int N_instance,struct Internal *I,ParamFCNNrule *P);
void AttributionPattern(int N_instance,int **Liste_pattern,struct Internal *I,ParamFCNNrule *P);
void Update1knnPrototype(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance);
double FindFromDistribution(double *dmax,int *Weight,int N_instance,int SeuilBruit,double SeuilExtreme);
void LibereMemoireAlgorithme(struct Internal *I,ParamFCNNrule *P);
int CalculCost(struct Internal *I,ParamFCNNrule *P,int N_instance);
#endif
