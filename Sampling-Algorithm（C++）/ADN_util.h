#ifndef ADNUTIL_H
#define ADNUTIL_H

int EcritFichier1D(char *fichier,int n,double *data);
int LoadFichier1D(char *fichier,int n,double *data);
void GenereMatriceDistance(double **data,int n,int p,int *liste,double **dist,double (*f_dist)(double *,double  *,int n));
double CalculCostEnd(double **data,int n,int p,double **prototype,int n_proto,double (*f_dist)(double *,double  *,int n));
double FindFromDistributionDIDES(double *dmax,int *Weight,int N_instance,int N_total,double SeuilExtreme);
int EcritFile(char *fichier,int Open,double *vecteur,int n);
int SelectinstanceWithListe(double **pattern,int p,int *liste_det,int n_inst,int Min,int *weightG,double **data);
int RegressionExponentielle(double *Pattern,int n_pattern,double &b,double &a);
int UpdateWithGravity_DENDIS(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance);
void SortPrototypes_DENDIS(struct Internal *I,ParamFCNNrule *P,int &N_instance,int **Liste_pattern,int Min);
void CalculusStatVectorBool(double *pattern,int N_pattern,int *booleen,double &Moy,
				  double &Ecart,double &Min,double &Max);
void CalculusStatVector(double *pattern,int N_pattern,double &Moy,double &Ecart,double &Min,double &Max);
double InertieGenAleaListe(double **data,int n,int p,int n_output);
void GenereMatriceDistance(double **data,int n,int p,int *liste,double **dist,double (*f_dist)(double *,double  *,int n));
int DuplicateData(double **data,int n,int p,double **data_out,int k_facteur);
int N_clusterDBSCAN(int *cat,double **prototype,int &n_proto,int p,int &bruit,double**,int *,int *weight_in,int *weight_out);
void NormaliseMinMax(double **data,int n,int p);
int LitFichier2DS(char *fichier,int n,int p,double **data);
int ReturnColonne(char *buffer);
int ExtraitDimensionFile(char *fichier,int *ligne,int *colonne);
int GenAleaListe(int n,int *liste,int n_output);
void Delete1D(int *x);
void Delete1D(double *x);
/*inline*/ void Delete2D(double **x,int n);
/*inline*/ double **Reservmemdouble(int n,int p);
/*inline*/ int **Reservmemint(int n,int p);
/*inline*/ void Delete2D(int **x,int n);
int *Reservememint(int n);
void GenAlea(double **data,int n,int p,double **output,int n_output);
int max(int a,int b);
int min(int a,int b);
double fmax(double a,double b);
double fmin(double a,double b);
int SelectinstanceWithListe(double **pattern,int p,int *liste_det,int n_inst,int Min,int *weightG,double **data);
void SortWeight(int Ngeneration,int **weight,int Min,int *ReelWeight,double **Radius,double *Indice);
void InitTri(struct Tri *Tr,int n);
void LibereTri(struct Tri *Tr);
int InegalTriangular(double **D_cluster,double D_nearest,int Cluster_actuel, int Cluster_candidat);
int SamplingApproach_1(double **data,int *Liste_instance,int &N_instance,int *weight,int **Liste,ParamFCNNrule *P);
void DisablePrototype(double *d_max,double *dist_proto,int N_proto,int *disable);
double **Reservmemdouble(int n,int p);
int tri_decroissant_int(const void *a,const void *b);
int tri_decroissant_double(const void *a,const void *b);
int tri_decroissant(const void *a,const void *b);
int tri_croissant(const void *a,const void *b);
int tri_croissant_int(const void *a,const void *b);
double Dist_eucl(double *p1,double *p2,int taille);
double d_Eucl(double *i1, double *i2, int n);
int EcritFichier2D(char *fichier,int n,int p,double **data);
void CalculusStatVector(double *pattern,int N_pattern,double &Moy,double &Ecart);
void NormaliseEcartType(double **data,int n,int p);



#endif 
