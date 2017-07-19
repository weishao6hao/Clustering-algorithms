
#include "ADN_Esclave1.h"
#include "ADN_Esclave2.h"
#include "ADN_util.h"


void InitSingleParamFCNNrule(ParamFCNNrule *P)
{
	P->SansStep2 = 0;
	P->TradeoffDensite = 2;
	P->FactorDistance = 1;
	P->Thres_RatioKmin = 0.9;
	P->Thres_RatioKmax = 0.9;
	P->Thres_d_1knn =  15.;
	P->Thres_Wt_1knn = 10.;//Percentange
	P->Thres_SequenceD = 15;
	P->Thres_SequenceW = 5;
	P->Thres_EcartN = 90;
	P->Thres_EcartW= 5.;//Delta des inerties du step t à t+1.
	P->Thres_FreqW = 0.9;
	P->Ratio_distance = 0.5;
	P->Coef = 0.8;
	P->V_seuilH = 50;
	P->C_seuilH = 2;
	P->V_seuilB = 20;
	P->C_seuilB = 1;
	P->Option_Alea_Cleaning_1 = 0;
	P->Option_Alea_Cleaning_2 = 0;
	P->Option_Alea_Bruit = 0;
	P->Option_Alea_step2 = 1;
	P->Init_Aleatoire = 1;
	P->Seuil_massique = 0.2;
	P->Option_Assymptotique = 0;
	P->DensiteMin=10;
	P->NbPointsMin=1;
	P->Granularity = 0.01;
	P->NbPointsNoise=0;
	P->OccMin = 10;
	P->Auto = 0;
	P->L_pattern = 1;
	P->Save = 0;
	P->p = 1;
	P->Epsilon_Assymptote = 0.001;
	P->Epsilon_R = 0.01;
	P->Plage_R = 5;
	P->Dist_Noise = 1.25;
	P->Dist_Connexity = 1.1;
	P->Barycentre = 0;
	P->Save = 0;
	P->Option_Alea_Bruit = 1;
	P->Mass_WeightCriterion = 1;
	P->Mass_OptionFinal = 1;
	P->ActifMasse = 0;
	P->S2_Option_all = 0;
	P->S2_Beta = 1.5;
	P->Ratio_achieved = 0;
	P->Tha_achieved = 0;
}

void InitSingleNull(ParamFCNNrule *P)
{
	P->Ave_d1nn_proto  = NULL;
	P->d_1knn  = NULL;
	P->d_1knnmin  = NULL;
	P->weight_1knn  = NULL;
	P->dav_1knn  = NULL;
	P->dec_1knn  = NULL;
	P->Wt_1knn  = NULL;
	P->Wtnext_1knn  = NULL;
	P->WI1_1knn  = NULL;
	P->WI2_1knn  = NULL;
	P->II_g  = NULL;
	P->II_1knn  = NULL;
	P->R_massique  = NULL;
	P->Inertie_proto  = NULL;
	P->pattern = NULL;
	P->Cost_n = NULL;
}	


ParamFCNNrule *InitParamFCNNrule(int n)
{
struct ParamFCNNrule *P;
int pmax = 20;

	P = new ParamFCNNrule;

	P->n = n;
	P->fdist = d_Eucl;
	P->Ave_d1nn_proto = new double[n];
	P->d_1knn = new double[n];
	P->d_1knnmin = new double[n];//liste des min
	P->weight_1knn = new int[n];//liste des poids des max.
	P->dav_1knn = new double[n];//liste des moyennes des distances
	P->dec_1knn = new double[n];//liste des moyennes des écart types des distances.
	P->Wt_1knn = new double[n];//liste des estimés des moyennes des distances
	P->Wtnext_1knn = new double[n];//liste des écart types des poids.
	P->WI1_1knn = new double[n];//liste des estimés des écart types des distances
	P->WI2_1knn = new double[n];//liste des moyennes des poids.
	P->II_g = new double[n];//liste des estimés des moyennes des poids.
	P->II_1knn = new double[n];//liste des estimés des écart types des poids.
	P->R_massique = new double[n];
	P->Inertie_proto = new double[n];
	P->Cost_n = new double[n];
	P->pattern = Reservmemdouble(n,pmax);
	P->MaxClusterMemory = n;

	for (int i =0;i<n;i++)
	{
		P->Cost_n[i] = MAX_VALUE;
		P->d_1knn[i] =0;
		P->d_1knnmin[i] = 0;
		P->weight_1knn[i] = 0;
		P->dav_1knn[i] = 0;
		P->dec_1knn[i] = 0;
		P->Wt_1knn[i] = 0;
		P->Wtnext_1knn[i] = 0;
		P->WI1_1knn[i] = 0;
		P->WI2_1knn[i] = 0;
		P->II_g[i] = 0;
		P->II_1knn[i] = 0;
		P->Ave_d1nn_proto[i] = 0;
		P->Inertie_proto[i] = 0;
		P->R_massique[i] = 0;
	}
	P->MaxCluster = n;
	
	for (int i=0;i<1000;i++)
	{
		P->Nb_distance[i] = 0;
		P->Nb_TriangularTest[i] = 0;
	}

	return P;
}

void LibereParamFCNNrule(struct ParamFCNNrule *P)
{
	if (P->d_1knn) delete [] P->d_1knn;
	if (P->d_1knnmin)   delete [] P->d_1knnmin;
	if (P->weight_1knn) delete [] P->weight_1knn;
	if (P->dav_1knn) delete [] P->dav_1knn;
	if (P->dec_1knn) delete [] P->dec_1knn;
	if (P->Wt_1knn) delete [] P->Wt_1knn;
	if (P->Wtnext_1knn) delete [] P->Wtnext_1knn;
	if (P->WI1_1knn) delete [] P->WI1_1knn;
	if (P->WI2_1knn) delete [] P->WI2_1knn;
	if (P->II_g) delete [] P->II_g;
	if (P->II_1knn) delete [] P->II_1knn;
	if (P->R_massique) delete [] P->R_massique;
	if (P->Ave_d1nn_proto) delete [] P->Ave_d1nn_proto;
	if (P->Inertie_proto) delete [] P->Inertie_proto;
	if (P->Cost_n) delete [] P->Cost_n;
	Delete2D(P->pattern,P->MaxClusterMemory);
	
//	delete P;

}
