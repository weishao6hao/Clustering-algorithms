#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "ADN_Esclave1.h"
#include "ADN_util.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////
int InitMemory(struct Internal *I,ParamFCNNrule *P)
{
	I->dmaxPrototype = new double[P->MaxCluster+1];
	I->disable = new int[P->MaxCluster+1];
	I->Dist_limite = new double[max(P->n,P->MaxCluster+1)];
	I->rep = new int[max(P->n,P->MaxCluster+1)];
	I->free = new int[max(P->n,P->MaxCluster+1)];
	I->Nearest = new double[max(P->n,P->MaxCluster+1)];
	I->Pre_individu = new int[max(P->n,P->MaxCluster+1)];
	I->Pre_distance = new double[max(P->n,P->MaxCluster+1)];
	I->Liste_selectionne = new int[max(P->n,P->MaxCluster+1)];
	I->Nearest_proto =  new int[max(P->n,P->MaxCluster+1)];
	I->Weight = new int[max(P->n,P->MaxCluster+1)];
	I->Inertie_proto = new double[P->MaxCluster+1];
	I->d1nn_proto = new double[P->MaxCluster+1];
	I->Inertie_1knn = new double[P->MaxCluster+1];
	I->Distribution_dMax = new double[P->MaxCluster+1];
	I->Distribution_Weight = new int[P->MaxCluster+1];
	I->Order_selection = new int[max(P->n,P->MaxCluster+1)];
	I->D_cluster = Reservmemdouble(P->MaxCluster+1,P->MaxCluster+1);
	I->data_proto = Reservmemdouble(P->n,MAX_CLUSTER);
	I->Nearest_between = new int[P->MaxCluster+1];
	I->Nearest_first = new int[P->MaxCluster+1];
	I->C_proto = new double[P->p];
	I->C_Oldproto = new double[P->p];
	I->Mass_OptionFinal = P->Mass_OptionFinal;
	I->Mass_Weight_Criterion = P->Mass_WeightCriterion;
	I->Tr = new struct Tri;
	InitTri(I->Tr,P->MaxCluster+1);

	for (int i=0;i<P->MaxCluster+1;i++)
		I->Weight[i] = 0;
	
	return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void LibereMemoireAlgorithme(struct Internal *I,ParamFCNNrule *P)
{
	delete [] I->dmaxPrototype;
	delete []I->disable;
	delete []I->Dist_limite;
	delete []I->rep;
	delete []I->free;
	delete []I->Nearest;
	delete []I->Pre_individu; 
	delete []I->Pre_distance; 
	delete []I->Liste_selectionne;
	delete []I->Nearest_proto;
	delete []I->Weight;
	delete []I->Inertie_proto;
	delete []I->Order_selection;
	delete []I->C_proto;
	delete []I->C_Oldproto;
	delete []I->d1nn_proto;
	delete []I->Inertie_1knn;
	delete []I->Distribution_dMax;
	delete []I->Distribution_Weight;
	delete []I->Nearest_between;
	delete []I->Nearest_first; 
	Delete2D(I->D_cluster,P->MaxCluster+1);
	Delete2D(I->data_proto,P->n);
	LibereTri(I->Tr);
	delete I->Tr;	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int InitParametres(double **data,struct Internal *I,ParamFCNNrule *P)
{
int i,j,winner;
double vector[1000];
double d_min = 100000;
double dist;
	
	for (i=0;i<P->MaxCluster;i++)
		for (j=0;j<P->MaxCluster;j++)
			I->D_cluster[i][j] = 0;
	
	for (i=0;i<P->n;i++)
		{I->free[i]=1; I->Nearest[i]=MAX_VALUE; I->Pre_distance[i]=-MAX_VALUE; I->Nearest_proto[i]=0;}//on initialise avec le cluster 0. 
	
	for (i=0;i<P->MaxCluster;i++)
		{I->dmaxPrototype[i] = MAX_VALUE; I->disable[i] = 0;	}

	winner = 0;
	if (P->Init_Aleatoire == 0)
	{	
		for (j=0;j<P->p;j++)
			vector[j] = MAX_VALUE;

		//on cherche le min...
		for (i=0;i<P->n;i++)
		{
			for (j=0;j<P->p;j++)
				if (vector[j] > data[i][j])
					vector[j] = data[i][j];		
		}
		for (i=0;i<P->n;i++)
		{	
			dist = P->fdist(vector,data[i],P->p);
			if (dist < d_min)
			{
				d_min = dist;
				winner = i;
			}
		}
		I->Nearest[winner] = 0;
	}
	else
	{
		winner = rand()%(P->n-1);
		I->Nearest[winner] = 0;
	}
	I->Weight[0] = P->n;
	
	I->Nearest_first[0] = 0;
	I->Pre_individu[0] = winner;
	I->free[winner] = 0;
	I->Liste_selectionne[0] = winner;
	I->Order_selection[winner] = 0;

	for (i=0;i<P->n;i++)//pour tous les éléments de T.
		I->Nearest_proto[i] = 0;	

	I->C_stable = 0;
	I->Ref_distance = -1;
	I->Ref_OK = 0;
	I->Nb_isole = 0;
	I->Nb_bruit = 0;
	I->N_instance_Extreme = 10000;

return winner;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
int TestFin_DENDIS(struct Internal *I,ParamFCNNrule *P,int &N_instance)
{

	if (N_instance == P->n-1)
		return 1;

	if (N_instance >= P->MaxCluster)
		return 1;

	if (P->Auto == 1)//on utilise la distance...
	{
		if (I->Mass_Weight_Criterion == 0 && P->SansStep2 == 1)
			return 1;

		if (P->d_1knn[N_instance - 1] < I->Ref_distance)
			return 1;
	}
		
	return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
void Update1knnPrototypeMasse(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance)
{
int i;
double Moy,Ecart,Min,Max;
int *booleen;

	Moy = 0;
	//Update des distances entre clusters.
	P->dav_1knn[N_instance - 1] = 0;

	booleen = new int[N_instance];
	for (i=0;i<N_instance;i++)
	{
		if (I->Weight[i] > P->Granularity*0.1*P->n)
			booleen[i] = 1;
		else
			booleen[i] = 0;
		CalculusStatVectorBool(I->Pre_distance,N_instance,booleen,Moy,Ecart,Min,Max);
	}

//	I->Ref_distance = Moy -Ecart;
	delete []booleen;

	for (i=0;i<N_instance;i++)
	{
		I->D_cluster[i][N_instance] = P->fdist(data[I->Liste_selectionne[i]],data[I->Pre_individu[I->WorseWinner]],P->p); 
		I->D_cluster[N_instance][i] = I->D_cluster[i][N_instance];
//		P->dav_1knn[N_instance - 1] += I->Pre_distance[i]; 
	}
//	P->dav_1knn[N_instance - 1] /= (double)(N_instance);
	P->dav_1knn[N_instance - 1] = Moy;//la moyenne dénuée du bruit.


	//Calcul des 1 plus proches voisins entre prototype.
	if (P->Auto == 1 || P->Save==1)
	{
		if (N_instance == 1)//donc 2...avec le nouveau!
			I->Ave_d1nn_proto = I->Pre_distance[I->WorseWinner];
		else
			I->Ave_d1nn_proto = (I->Ave_d1nn_proto*(double)(N_instance) - I->d1nn_proto[I->WorseWinner] + 2*I->Pre_distance[I->WorseWinner]) / (double)(N_instance + 1);

		I->d1nn_proto[I->WorseWinner] = I->Pre_distance[I->WorseWinner];//le nouveau prototype est forcément plus proche de son prototype d'origine: les 2 distances min sont modifiées.
		I->d1nn_proto[N_instance] = I->Pre_distance[I->WorseWinner];
		P->Ave_d1nn_proto[N_instance] = I->Ave_d1nn_proto;
	}

	I->Inertie_1knn[N_instance - 1]  = I->Ave_d1nn_proto;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
void UpdateInstance(int N_instance,int Groupe_actif,int *Liste_selectionne,int &pool,struct Internal *I)
{
int pattern,j;

	for (j=0;j<N_instance;j++)
	{
		if (I->Tr->G_actif[j] == Groupe_actif)
		{
			pattern = I->Liste_selectionne[j];
			Liste_selectionne[pool] = pattern;
			pool++;
		}
	}

}
//////////////////////////////////////////////////////////////////////////////////////////////////////
int InitProgramme(double **data,struct Internal *I,ParamFCNNrule *P)
{
int First_instance=0;

	InitMemory(I,P);
//	I->Max_distance = InertieGenAleaListe(data,P->n,P->p,max(P->n/100,500));
	I->Max_distance = 1;
	First_instance = InitParametres(data,I,P);
	I->Ref_distance = -1;
	return First_instance;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void InitBoucle(struct Internal *I,int N_instanceOld,int New_instance)
{
int i;

	I->Liste_selectionne[N_instanceOld] = New_instance;

	if (N_instanceOld > 0)//On disable certains prototyppes au sens de l'inégalité triangulaire...
		DisablePrototype(I->Pre_distance,I->D_cluster[N_instanceOld],N_instanceOld+1,I->disable);
	
	for (i=0;i<N_instanceOld+1;i++)//on initialise la distance pour l'ensemble DS.
	{	I->Pre_distance[i] = 0; I->Pre_individu[i] = -1; I->Dist_limite[i] = 0; }

}
//////////////////////////////////////////////////////////////////////////////////////////////////////
int UpdateWithGravity(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance)
{
int i,j,winner;

	for (i=0;i<N_instance;i++)//pour tous les éléments de T.
	{
		for (j=0;j<P->p;j++)
			P->pattern[i][j] = 0;
	}

	for (i=0;i<P->n;i++)//pour tous les éléments de T.
	{
		winner = I->Nearest_proto[i];//l'instance la plus proche est mémorisé...on garde en mémoire l'ancien.
		for (j=0;j<P->p;j++)
			P->pattern[winner][j] = P->pattern[winner][j] + data[i][j]; 
	}


	for (i=0;i<N_instance;i++)
		for (j=0;j<P->p;j++)
			P->pattern[i][j] /= (double)(I->Weight[i]);

return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int FindNewInstance_DIDES(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance)
{
int i;
int winner;
double dist, DistWinner = 0;
int NewInstance;
int Data_NewInstance;

	NewInstance = N_instance - 1;
	Data_NewInstance = I->Liste_selectionne[NewInstance];
	I->Nearest_proto[Data_NewInstance] = NewInstance;

	if (N_instance > 1)
	{
		I->Weight[N_instance-1] = 1;//il y a donc 1 individu dans la nouvelle instance
		I->Weight[I->WorseWinner] = I->Weight[I->WorseWinner] - 1; 	// 1 de moins dans l'ancienne.
	}

	I->WorseWinner = -1;
	for (i=0;i<P->n;i++)//pour tous les éléments de T.
	{
		if (I->free[i] == 1)//s'il est libre.
		{
			winner = I->Nearest_proto[i];//l'instance la plus proche est mémorisé...on garde en mémoire l'ancien.
			if (I->disable[I->Nearest_proto[i]] == 0)//On ne peut pas sur la base du prototype se décider.
			{
//				P->Nb_TriangularTest[N_instance - 1] += 1;
				//si la distance au proto actuel est supérieur à la distance limite au nouveau
					//si d(I->Nearest_proto[i],NewInstance) > 2*I->Nearest[i]) retour 0, le pattern ne peut pas appartenir au nouveau prototype. Cas contraire, on ne sait pas.
					if (InegalTriangular(I->D_cluster,I->Nearest[i],I->Nearest_proto[i],NewInstance))
					{
						dist = P->fdist(data[i],data[Data_NewInstance],P->p);
//						P->Nb_distance[N_instance - 1] += 1;
						if (dist < I->Nearest[i])//plus proche d'un prototype...
						{				
							I->Nearest[i] = dist;//on update le Farest.					
							I->Weight[winner] = I->Weight[winner] - 1;//on diminue de 1 l'instance impactée.
							I->Weight[NewInstance] = I->Weight[NewInstance] + 1; //on augmente de 1 l'instance impactée.
							I->Nearest_proto[i] = NewInstance;//on update le prototype...
							winner = NewInstance;//le prototype gagnant.	
						}
//					}//fin du I->Nearest[i]<=I->Dist_limite[winner]
				}//fin du for if.
			}//on regarde son plus proche.
			
			if (I->Nearest[i] >= I->Pre_distance[winner])//le plus éloigné du plus proche
			{
				I->Pre_distance[winner] = I->Nearest[i];//la distance entre le prototype et son plus éloigné.
				I->Pre_individu[winner] = i;// i devient le plus eloigné des plus proches!		
				if (I->Pre_distance[winner] > DistWinner)//le prototype ayant le plus éloigné des individus.
				{
					I->WorseWinner = winner;//le nouveau prototype.
					DistWinner = I->Pre_distance[winner];
				}
			}//fin du if nearest[i].
		}//fin du if free.
	}//fin du for i
	
	if (I->WorseWinner == -1)
		return -1;
	else
		return I->Pre_individu[I->WorseWinner];
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
int FindNewInstance_PROBA(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance)
{ 
int i;
int winner;
double dist;
int NewInstance;
int Data_NewInstance,maxW;
double p1, p2;
double MaxScore,maxD; 
	
	NewInstance = N_instance - 1;
	Data_NewInstance = I->Liste_selectionne[NewInstance];
	I->Nearest_proto[Data_NewInstance] = NewInstance;

	if (N_instance > 1)
	{
		I->Weight[N_instance-1] = 1;//il y a donc 1 individu dans la nouvelle instance
		I->Weight[I->WorseWinner] = I->Weight[I->WorseWinner] - 1; 	// 1 de moins dans l'ancienne.
	}

	I->WorseWinner = -1;
	for (i=0;i<P->n;i++)//pour tous les éléments de T.
	{
		if (I->free[i] == 1)//s'il est libre.
		{
			winner = I->Nearest_proto[i];//l'instance la plus proche est mémorisé...on garde en mémoire l'ancien.
			if (I->disable[I->Nearest_proto[i]] == 0)//On ne peut pas sur la base du prototype se décider.
			{
				//si la distance au proto actuel est supérieur à la distance limite au nouveau
					//si d(I->Nearest_proto[i],NewInstance) > 2*I->Nearest[i]) retour 0, le pattern ne peut pas appartenir au nouveau prototype. Cas contraire, on ne sait pas.
					if (InegalTriangular(I->D_cluster,I->Nearest[i],I->Nearest_proto[i],NewInstance))
					{
						dist = P->fdist(data[i],data[Data_NewInstance],P->p);
						if (dist < I->Nearest[i])//plus proche d'un prototype...
						{				
							I->Nearest[i] = dist;//on update le Farest.					
							I->Weight[winner] = I->Weight[winner] - 1;//on diminue de 1 l'instance impactée.
							I->Weight[NewInstance] = I->Weight[NewInstance] + 1; //on augmente de 1 l'instance impactée.
							I->Nearest_proto[i] = NewInstance;//on update le prototype...
							winner = NewInstance;//le prototype gagnant.	
						}
//					}//fin du I->Nearest[i]<=I->Dist_limite[winner]
				}//fin du for if.
			}//on regarde son plus proche.

			if (I->Nearest[i] >= I->Pre_distance[winner])//le plus éloigné du plus proche
			{
				I->Pre_distance[winner] = I->Nearest[i];//la distance entre le prototype et son plus éloigné.
				I->Pre_individu[winner] = i;// i devient le plus eloigné des plus proches!		

			}//fin du if nearest[i].
	
		}//fin du if free.
	}//fin du for i
	
	maxW=0;
	maxD = 0;
	for (i=0;i<N_instance;i++)
	{
		if (I->Weight[i] >= maxW)
			maxW = I->Weight[i];
		if (I->Pre_distance[i] >=maxD)
			maxD = I->Pre_distance[i];
	}

	MaxScore = 0;
	for (i=0;i<N_instance;i++)
	{
		p1 = (double)I->Weight[i]/(double)maxW;
		p2 = (double)I->Pre_distance[i]/(double)maxD;
		if (p1*p2 > MaxScore)
		{
			MaxScore = p1*p2;
			I->WorseWinner = i;
		}
	}
//	printf("I=%d MaxScore=%lf %d\n",i,MaxScore,I->Pre_individu[I->WorseWinner]);
	if (I->WorseWinner == -1)
		return -1;
	else
		return I->Pre_individu[I->WorseWinner];
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
int CalculCost(struct Internal *I,ParamFCNNrule *P,int N_instance)
{
int i;
double  CostC;
double Dist;

	CostC = 0;
	Dist = 0;
	for (i=0;i<N_instance;i++)
	{
		CostC = CostC + I->Weight[i]*(I->Pre_distance[i]);
		Dist += I->Pre_distance[i];
	}


	CostC /= (double)(P->n*I->Max_distance);
	P->Cost_n[N_instance - 1] = CostC;

	if (CostC<P->Granularity)
		return 1;
	else 
		return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void Update1knnPrototype(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance)
{
int i,j;
double min_d=1000;


	if (N_instance == 1)//donc 2...avec le nouveau!
	{
		I->Ave_d1nn_proto = I->Pre_distance[I->WorseWinner];
//		P->dav_1knn[N_instance-1] = I->Pre_distance[I->WorseWinner];
	}
	else
	{
		I->Ave_d1nn_proto = (I->Ave_d1nn_proto*(double)(N_instance) - I->d1nn_proto[I->WorseWinner] + 
 			                 2*I->Pre_distance[I->WorseWinner]) / (double)(N_instance + 1);

//		P->dav_1knn[N_instance-1] = (P->dav_1knn[N_instance-2]*(double)(N_instance-1) - I->d1nn_proto[I->WorseWinner] + 
// 									2*I->Pre_distance[I->WorseWinner]) / (double)(N_instance); 
	}
	
	
	//Update des distances entre clusters.		
	for (i=0;i<N_instance;i++)
	{
		I->D_cluster[i][N_instance] = P->fdist(data[I->Liste_selectionne[i]],data[I->Pre_individu[I->WorseWinner]],P->p); 
		I->D_cluster[N_instance][i] = I->D_cluster[i][N_instance];
	}

	//Calcul des 1 plus proches voisins entre prototype.
	if (P->Auto == 1 || P->Save==1)
	{
		
		I->Nearest_between[I->WorseWinner] = N_instance-1;//le plus proche prototype devient l'instance sélectionné.
		I->Nearest_between[N_instance-1]  = I->WorseWinner;//l'instance selectionné devient le plus proche prototype.
		
		if (N_instance < MAX_CLUSTER)//chaque instance est son plus proche voisin...
			I->Nearest_first[N_instance] = N_instance;
		else
		{
			min_d = 1000;
			for (j=0;j<MAX_CLUSTER;j++)//On raisonne avec un nombre de CLUSTER  = 3.
			{
				if (I->D_cluster[N_instance][j] < min_d)
				{
					min_d = I->D_cluster[N_instance][j];
					I->Nearest_first[N_instance] = j;
				}
			}
		}

		if (N_instance <= MAX_CLUSTER)
		{
			for (i=0;i<P->n;i++) //pour chaque data on update sa distance la plus proche avec le nombre d'instance.
				I->data_proto[i][N_instance-1] = I->Nearest[i];
		}

	//	printf("%lf\n",I->data_proto[0][N_instance-1]);
		
		I->d1nn_proto[I->WorseWinner] = I->Pre_distance[I->WorseWinner];//le nouveau prototype est forcément plus proche de son prototype d'origine: les 2 distances min sont modifiées.
		I->d1nn_proto[N_instance] = I->Pre_distance[I->WorseWinner];
		
	}

	I->Inertie_1knn[N_instance - 1]  = I->Ave_d1nn_proto;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
void AttributionPattern(int N_instance,int **Liste_pattern,struct Internal *I,ParamFCNNrule *P)
{
int i,j;
int *w_int;
int winner;

	if (P->L_pattern == 0)
		return;

	w_int = new int[P->n];
	
	for (i=0;i<P->n;i++)
		w_int[i] = 0;
	
	for (j=0;j<P->n;j++)
	{
		winner = I->Nearest_proto[j];
		w_int[winner]++;
	}

	for (i=0;i<N_instance;i++)//On a un pattern non attribué à ce niveau.
	{
		if (Liste_pattern[i]!=NULL)
		{
			delete []Liste_pattern[i];
			Liste_pattern[i] = NULL;
		}
	
		if (w_int[i] > 0)
			Liste_pattern[i] = new int[w_int[i]];
		w_int[i] = 0;
	}
	
	for (j=0;j<P->n;j++)
	{
		winner = I->Nearest_proto[j];
		Liste_pattern[winner][w_int[winner]] = j;
		w_int[winner]++;
	}
	delete []w_int;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
void ProcessCriteria(int N_instance,struct Internal *I,ParamFCNNrule *P)
{
int i;
double dmin, dmax;
double M;

//	if (P->Auto == 0 && P->Save==0)
//		return;
	if (P->Auto == 0)
		return;

	dmin = 100000;
	dmax = 0;
	M=0;
	for (i=0;i<N_instance;i++)
	{
		if (I->Pre_distance[i] < dmin)
			if (I->Weight[i] > P->NbPointsNoise)
				dmin = I->Pre_distance[i];
		
		if (I->Pre_distance[i] > dmax)
			dmax = I->Pre_distance[i];

		M+=I->Pre_distance[i];
	}

	P->dav_1knn[N_instance - 1] = M/(double)(N_instance);
	P->d_1knn[N_instance - 1] = dmax;
	P->d_1knnmin[N_instance - 1] = dmin;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////
int TestFin(struct Internal *I,ParamFCNNrule *P,int &N_instance)
{
	if (N_instance >= (P->MaxCluster) || N_instance == P->n-1)
		return 1;
	
	if (P->Auto == 1)
	{
		if (P->d_1knn[N_instance - 1] < I->Ref_distance)
			return 1;
	}

	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void ProcessCentreGravitePrototype(double **data,int N_instance,struct Internal *I,ParamFCNNrule *P)
{
int i;
static double Error;
double Epsilon;

	if ((P->Auto==0) && (P->Save==0))
		return;

	if (I->C_stable == 1)
		return;

	if (N_instance == 1)
	{
		Error = 0;
		for (i=0;i<P->p;i++)
		{
			I->C_Oldproto[i] = I->C_proto[i] = data[I->Liste_selectionne[0]][i];
			Error += fabs(I->C_Oldproto[i]);
		}
		Error = Error / (100*(double)(P->p));
	
	}
	else
	{
		Epsilon = 0;
		for (i=0;i<P->p;i++)
		{
			I->C_proto[i] = (I->C_proto[i]*(double)(N_instance - 1) + data[I->Liste_selectionne[N_instance - 1]][i])/(double)(N_instance);		
			Epsilon += fabs(I->C_proto[i] - I->C_Oldproto[i]);
		}

		for (i=0;i<P->p;i++)
			I->C_Oldproto[i] = I->C_proto[i];

		Epsilon = Epsilon / (double) (P->p);
		if (Epsilon < Error)
			I->C_stable = 1;
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////
void UpdateAlgorithm(struct Internal *I,int N_instance,int Data_Newinstance)
{
	if (Data_Newinstance != -1)
	{
		I->Nearest[Data_Newinstance] = 0;
		I->free[Data_Newinstance] = 0;//on l'enlève de la liste.
		I->Liste_selectionne[N_instance] = Data_Newinstance;//on met dans la liste.
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
int FindNewInstanceMassique(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance)
{ 
int i,winner;
double dist, DistWinner = 0;
int NewInstance;
int Data_NewInstance;
int maxWeight,winWeight;
//double Min,Max,Moy,Ecart;
int OK_fin=0;
double maxmax = 0;
double RatioW;

	NewInstance = N_instance - 1;
	Data_NewInstance = I->Liste_selectionne[NewInstance];
	I->Nearest_proto[Data_NewInstance] = NewInstance;

	if (N_instance > 1 && I->WorseWinner!=-1)
	{
		I->Weight[N_instance-1] = 1;//il y a donc 1 individu dans la nouvelle instance
		I->Weight[I->WorseWinner] = I->Weight[I->WorseWinner] - 1; 	// 1 de moins dans l'ancienne.
	}

	I->WorseWinner = -1;

	for (i=0;i<P->n;i++)//pour tous les éléments de T.
	{
		if (I->free[i] == 1)//s'il est libre.
		{
			winner = I->Nearest_proto[i];//l'instance la plus proche est mémorisé...on garde en mémoire l'ancien.
			if (I->disable[I->Nearest_proto[i]] == 0)//On ne peut pas sur la base du prototype se décider.
			{
				if (InegalTriangular(I->D_cluster,I->Nearest[i],I->Nearest_proto[i],NewInstance))
				{
					dist = P->fdist(data[i],data[Data_NewInstance],P->p);
					if (dist < I->Nearest[i])//plus proche d'un prototype...
					{				
						I->Nearest[i] = dist;//on update le Farest.					
						I->Weight[winner] = I->Weight[winner] - 1;//on diminue de 1 l'instance impactée.
						I->Weight[NewInstance] = I->Weight[NewInstance] + 1; //on augmente de 1 l'instance impactée.
						I->Nearest_proto[i] = NewInstance;//on update le prototype...
						winner = NewInstance;//le prototype gagnant.	
					}
				}//fin du for if.
			}//on regarde son plus proche.
			
			if (I->Nearest[i] >= I->Pre_distance[winner])//le plus éloigné du plus proche
			{
				I->Pre_distance[winner] = I->Nearest[i];//la distance entre le prototype et son plus éloigné.
				I->Pre_individu[winner] = i;// i devient le plus eloigné des plus proches!		
			}//fin du if nearest[i].
		}//fin du if free.
	}//fin du for i

	//on calcule le max max.
	for (i=0;i<N_instance;i++)
		if (I->Pre_distance[i] > maxmax)
		{
			maxmax = I->Pre_distance[i];
			winner = i;
		}

	if (I->Mass_Weight_Criterion==1)
	{		
		winWeight = maxWeight = 0;
		OK_fin = 0;
		for (i=0;i<N_instance;i++)
		{
			//On conditionne avec les distances.
			if (N_instance >= 2)
			{
				if (I->Weight[i] > maxWeight)//recherche du MAX
				{
					RatioW =  1 + P->TradeoffDensite*(double)(I->Weight[i])/(double)(P->NbPointsMin);
					if  (I->Pre_distance[i]*RatioW>=P->dav_1knn[N_instance - 2])//Moyenne à l'n-2. 
					{
						maxWeight = I->Weight[i];//le poids maximum.
						winWeight = i;//le numéro du prototype gagnant
						OK_fin = 1;
					}
				}
			}
			else
			{
				if (I->Weight[i] > maxWeight)//aucune condition sur les distances.
				{
					maxWeight = I->Weight[i];//le poids maximum.
					winWeight = i;//le numéro du prototype gagnant
					OK_fin = 1;
				}
			}
		}	

		if (maxWeight <= P->NbPointsMin || OK_fin==0)
		{
			I->Ref_distance = P->FactorDistance * P->dav_1knn[N_instance - 2];
			
			I->Mass_Weight_Criterion = 0;			
//			printf("Stop=%dinst Ref=%lf\n",N_instance,I->Ref_distance);
		}

		if (OK_fin == 0)//il n'y en a plus...
			return -2;
		
		NewInstance = I->Liste_selectionne[winWeight];//on va chercher l'individu le plus éloigné.
		I->WorseWinner = winWeight;//numero de l'instance.
	}
	else
	{
		DistWinner = 0;
		for (i=0;i<N_instance;i++)
		{
			if (I->Pre_distance[i] > DistWinner)//le prototype ayant le plus éloigné des individus.
			{
				I->WorseWinner = i;//le nouveau prototype.
				DistWinner = I->Pre_distance[I->WorseWinner];
			}
		}
		winWeight = I->WorseWinner;
//		printf("Worse=%d\n",I->WorseWinner);
	}

	if (I->WorseWinner == -1)
		return -1;
	else
		return I->Pre_individu[winWeight];
}
///////////////////////////////////////////////////////////////////////////////////
int UpdateWeight_DIDES(struct Internal *I,ParamFCNNrule *P,int N_prototype,int n_all,double **data)
{
int i,j,winner;
double dist, dmax;
double W_total;

	if (N_prototype >= n_all)//on ne fait rien! 
		return 1;

	W_total = 0;
	for(i=N_prototype;i<n_all;i++)
		W_total += I->Weight[i];

	winner = 0;
	for (i=0;i<P->n;i++)
	{
		if (I->Nearest_proto[i] >= N_prototype)//il faut placer ce pattern.
		{
			dmax = 100000000;
			for (j=0;j<N_prototype;j++)
			{
				dist = P->fdist(data[i],data[I->Liste_selectionne[j]],P->p);
				if (dist < dmax)
				{
					dmax = dist;
					winner = j;
				}
			}
			I->Nearest_proto[i] = winner;
			I->Weight[winner] += 1;
			if (I->Pre_distance[winner] < dmax)
				I->Pre_distance[winner] = dmax;
		}
	}

	return 1;	
}
///////////////////////////////////////////////////////////////////////////////////
int ModeliseFunction(struct Internal *I,ParamFCNNrule *P,int N_prototype)
{
int N_Instance_Extreme;
double value_Extreme,value_min;
int i;

	RegressionExponentielle(P->d_1knn,N_prototype,I->alpha_M,I->beta_M);//Modélisation de la courbe des max

	N_Instance_Extreme = (int)(-1. / (1 - exp((log (1-P->Epsilon_Assymptote)/I->beta_M))));//nombre d'instance pour ne plus avoir de variations...
	value_Extreme = I->alpha_M * pow(N_Instance_Extreme,I->beta_M);//valeur extreme du dmin.

	RegressionExponentielle(P->d_1knnmin,N_prototype,I->alpha_m,I->beta_m);//Modélisation de la courbe des min...
	value_min = I->alpha_m * pow(N_prototype,I->beta_m);
	I->N_instance_Extreme = N_Instance_Extreme;
	I->D_distance_Extreme = value_Extreme;
	I->Ref_OK = 1;
//	seuilMassique = (P->Granularity * (double)(P->n));
	if (P->Option_Assymptotique == 1)//A rajouter le lissage.
	{
		I->Ref_distance = FindFromDistributionDIDES(I->Distribution_dMax,I->Weight,I->Seuil_instance_Ratio,P->n,value_Extreme);
		if (P->d_1knn[N_prototype-1]<I->Ref_distance)//la dmax est inférieure au seuil, il faut arrêter avant.
		{
			for (i=0;i<N_prototype;i++)
			{
				if (P->d_1knn[i]<I->Ref_distance)
					break;
			}
			return i+1;
		}
	}
	else
	{
		I->Ref_distance = value_min;
	}

return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FindFromDistribution(double *dmax,int *Weight,int N_instance,int SeuilBruit,double SeuilExtreme)
{
int i;
double dmin;
double Moy, Ecart;
int N_Pattern_reel;

	Moy = 0;
	Ecart = 0;
	N_Pattern_reel = 0;
	dmin = 10000;

	for (i=0;i<N_instance;i++)
	{
		if (Weight[i] > SeuilBruit && dmax[i]>SeuilExtreme)
		{
			Moy += dmax[i];
			N_Pattern_reel++;
			if (dmax[i] < dmin)
				dmin = dmax[i];
		}
	}
	
	if (N_Pattern_reel > 0)
	{
		Moy = Moy / (double)(N_Pattern_reel);
	
		for (i=0;i<N_instance;i++)
			if (Weight[i] > SeuilBruit && dmax[i]>SeuilExtreme)
				Ecart = Ecart + (Moy - dmax[i]) * (Moy - dmax[i]);

		Ecart = sqrt(Ecart / (double)(N_Pattern_reel));
	}
	
	dmin = fmax(dmin,Moy - 2*Ecart);

return dmin;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FindFromDistributionDIDES(double *dmax,int *Weight,int N_instance,int N_total,double SeuilExtreme)
{
int i;
double pp, s_pp;

	double M_local=0;
	int nr=0;
	s_pp = 0;
	for (i=0;i<N_instance;i++)
	{
		if (dmax[i]>=SeuilExtreme)
		{
			pp = 1. - Weight[i]/(double)(N_total);
			M_local += pp * dmax[i];
			s_pp += pp;
			nr++;
		}
	}

	if (nr>0)
		M_local /= (double)(s_pp);
	else
		M_local = SeuilExtreme;


	return M_local;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int Test_endMassique(int N_instance,struct Internal *I,ParamFCNNrule *P)
{
int i;
int NbpointNoise;

	if (P->Save == 0 && P->Auto == 0)
		return 0;

	NbpointNoise = 0;
	I->Ratio_massique = 0;
	for (i=0;i<N_instance;i++)
	{
		if ( I->Weight[i] > P->NbPointsNoise)
		{
			if (I->Weight[i] <= P->NbPointsMin)				
				I->Ratio_massique = I->Ratio_massique + I->Weight[i];
		}
		else
			NbpointNoise++;
	}
		
	I->Ratio_massique = I->Ratio_massique / (double)(P->n-NbpointNoise);
	P->R_massique[N_instance - 1] = I->Ratio_massique;
	
	//C'est la condition.
	if (I->Ratio_massique >= P->Seuil_massique || N_instance >= I->N_instance_Extreme)
	{
		I->Seuil_instance_Ratio = N_instance;
		for (i=0;i<N_instance;i++)
		{
			I->Distribution_dMax[i] = I->Pre_distance[i];
			I->Distribution_Weight[i] = I->Weight[i];
		}
		return 1;
	}

return 0;
}

int RatioMassiqueAchieved(int N_instance,struct Internal *I,ParamFCNNrule *P)
{
int i;
int NbpointNoise;
double Ratio_Threshold;

	if (P->Save == 0 && P->Auto == 0)
		return 0;
	
	if (I->Ref_OK == 1)
		return 0;

	if (P->Ratio_achieved == 0)
	{
		NbpointNoise = 0; 	I->Ratio_massique = 0;
		for (i=0;i<N_instance;i++)
		{
			if ( I->Weight[i] > P->NbPointsNoise)
			{
				if (I->Weight[i] <= P->OccMin)
					I->Ratio_massique = I->Ratio_massique + I->Weight[i];
			}
		else
			NbpointNoise++;
		}
		I->Ratio_massique = I->Ratio_massique / (double)(P->n-NbpointNoise);
		P->R_massique[N_instance - 1] = I->Ratio_massique;
		//C'est la condition.
		if (I->Ratio_massique >= P->Seuil_massique)
		{
			P->Ratio_achieved = 1;
			I->Seuil_instance_Ratio = N_instance;
			for (i=0;i<N_instance;i++)
			{
				I->Distribution_dMax[i] = I->Pre_distance[i];
				I->Distribution_Weight[i] = I->Weight[i];
			}
			if (P->Tha_achieved == 1)
				return 1;
		}
	}

	if (P->Tha_achieved == 0)
	{
		NbpointNoise = 0; Ratio_Threshold = 0;
		for (i=0;i<N_instance;i++)
		{
			if ( I->Weight[i] > P->NbPointsNoise)
			{
				if (I->Weight[i] <= P->NbPointsMin)
					Ratio_Threshold = Ratio_Threshold + I->Weight[i];
			}
		else
			NbpointNoise++;
		}
		Ratio_Threshold = Ratio_Threshold / (double)(P->n-NbpointNoise);
		if (Ratio_Threshold >= P->Seuil_massique)
		{
			P->Tha_achieved = 1;
			if (P->Ratio_achieved == 1)
				return 1;
		}
	}

	
	

	
	return 0;


}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SortPrototypes_DIDES(struct Internal *I,ParamFCNNrule *P,int &N_instance,int **Liste_pattern,int Min)
{
int i,j;
int *Weight_tempon;
int **Liste_tempon;
int *SelectionIndividu;
int compteur = 0;

	if (P->Auto == 0)
		return;

	SelectionIndividu = new int[N_instance];
	Liste_tempon = new int*[N_instance];
	Weight_tempon = new int[N_instance];
	
	
	compteur = 0;
	for (i=0;i<N_instance;i++)
	{
		Liste_tempon[i] = NULL;
		if (I->Weight[i]>Min)
		{
			if (I->Tr->noise[i] != 2)
			{
				Weight_tempon[compteur] = I->Weight[i];
				Liste_tempon[compteur] = new int[I->Weight[i]];
				for (j=0;j<I->Weight[i];j++)
					Liste_tempon[compteur][j] = Liste_pattern[i][j];//MEMOIRE!!
				SelectionIndividu[compteur] = I->Liste_selectionne[i];							
				compteur++;
			}
		}
		if (Liste_pattern[i])
		{
			delete []Liste_pattern[i];
			Liste_pattern[i] = NULL;
		}
	}
	
	for (i=0;i<compteur;i++)
	{	
		Liste_pattern[i] = new int[Weight_tempon[i]];

		for (j=0;j<Weight_tempon[i];j++)
			Liste_pattern[i][j] = Liste_tempon[i][j];
		I->Liste_selectionne[i] = SelectionIndividu[i];
		I->Weight[i]=Weight_tempon[i];
	}	

	delete []Weight_tempon;
	delete []SelectionIndividu;
	for (i=0;i<N_instance;i++)
	{
		if (Liste_tempon[i])
			delete []Liste_tempon[i];
	}
	delete []Liste_tempon;
	
	N_instance = compteur;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int CalculInertieProto(double **data,int N_instance,struct Internal *I,ParamFCNNrule *P)
{
int i;
double Inertie;
double dist;

	if (P->Auto == 0 && P->Save==0)
		return 0;

	if (I->Ref_OK == 1)
		return 0;

	if (I->C_stable ==0)
	{
		Inertie = 0;
		for (i=0;i<N_instance;i++)
		{
			dist = P->fdist(I->C_proto,data[I->Liste_selectionne[i]],P->p);
			Inertie = Inertie + dist*dist; 
		}
		Inertie = Inertie / (double)(N_instance);
	}
	else
	{
		dist = P->fdist(I->C_proto,data[I->Liste_selectionne[N_instance - 1]],P->p);
		Inertie = (I->Inertie_proto[N_instance - 2] * (double)(N_instance-1) + dist*dist) / (double)(N_instance);
	}

	if (P->Auto == 1 || P->Save==1)
	{
		I->Inertie_proto[N_instance - 1] = Inertie;
		P->Inertie_proto[N_instance - 1] = Inertie;
	}

	return 1;
}

int InertieAchieved(int N_instance,struct Internal *I,ParamFCNNrule *P)
{
int i;
double R_1,R_2;
static double I_1knn[5];
static double I_proto[5];
	
	if (P->Save == 0 && P->Auto == 0)
		return 0;

	R_1 = R_2 = 0;

	if (N_instance >= (2*P->Plage_R))
	{
		for (i=1;i<P->Plage_R;i++)
		{
			I_1knn[i] = I_1knn[i-1];
			I_proto[i] = I_proto[i-1];
		}
		I_1knn[0]  = fmin(I->Inertie_proto[N_instance - 1],I->Inertie_proto[N_instance - 2])/fmax(I->Inertie_proto[N_instance - 1],I->Inertie_proto[N_instance - 2]);
		I_proto[0] = fmin(I->Inertie_1knn[N_instance - 1],I->Inertie_1knn[N_instance - 2])/fmax(I->Inertie_1knn[N_instance - 1],I->Inertie_1knn[N_instance - 2]);
		
		R_1 = R_2 = 0;
		for (i=0;i<P->Plage_R;i++)
		{
			R_1 = R_1 + I_1knn[i];
			R_2 = R_2 + I_proto[i];
		}
		R_1 = R_1 / (P->Plage_R);
		R_2 = R_2 / (P->Plage_R);
	}
	else
	{
		for (i=1;i<P->Plage_R;i++)
		{
			I_1knn[i]  = rand()%(100);
			I_proto[i] = rand()%(100);
		}
	}

	if (R_1 > (1-P->Epsilon_R) && R_2 > (1-P->Epsilon_R))
	{
		I->Seuil_instance_Inertie = N_instance;
		return 1;
	}
	else
		return 0;

}


