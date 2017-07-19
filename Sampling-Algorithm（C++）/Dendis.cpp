

#include "ADN_Esclave1.h"
#include "ADN_Esclave2.h"
#include "ADN_util.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int SamplingApproachP_DENDIS(double **data,int n,int p,double Granularity,double **prototype,int &N_instance,int *weight, char *chemin,int *option)	
{
int *Liste_instance;
int **Liste;
ParamFCNNrule *P;
int MaxCluster,i,N_prototype;
double dmax;
char Output[500];

	MaxCluster = N_instance;
	Liste_instance = new int[MaxCluster];
	Liste = new int *[MaxCluster];
	for (i=0;i<MaxCluster;i++)
		Liste[i] = NULL;

	P = InitParamFCNNrule(n);
	InitSingleParamFCNNrule(P);	
	P->p = p;
	P->n = n;
	
	P->Option_Alea_Cleaning_2 = 0;
	P->Option_Alea_Cleaning_1 = 0;
	P->SansStep2 = 1;
	P->Init_Aleatoire = option[1];
	P->Option_Assymptotique = 1;
	P->Auto = 1;
	P->L_pattern = 1;
	P->Save = 0;
	P->Dist_Connexity = 1;//lié à la moyenne
	P->Barycentre = 1;
	P->Granularity = Granularity;
	P->TradeoffDensite = 1;
	P->FactorDistance = 1;
	P->NbPointsMin = max(5,(int)(double)(n)*P->Granularity); //Agit sur le nombre de prototypes/
	P->OccMin = (int)((double)(n) * P->Granularity);
	P->NbPointsNoise = max(1,(int)(double)(n)/1000);//bruit
	P->DensiteMin = (int)((double)(n)*P->Granularity);//le nombre de patterns par occurence.
	P->n = n; P->p = p; 
	P->MaxCluster = MaxCluster;

//	type = 1; // on raisonne sur les distances sur type = 0, sur les points si type 1.
//	OptionPreliminaire = 0;
	P->S2_Beta = 1;
	N_instance=0;

	N_prototype = SamplingApproach_DENDIS(data,Liste_instance,N_instance,weight,Liste,dmax,P);
	N_prototype = SelectinstanceWithListe(data,p,Liste_instance,N_prototype,0,NULL,prototype);

	if (option[2] == 1)
	{
		if (chemin!=NULL)
			sprintf(Output,"%sdmaxDENDIS.txt",chemin);
		else
			sprintf(Output,"dmaxDENDIS.txt");
		EcritFile(Output,1,P->d_1knn,N_prototype);
		if (chemin!=NULL)
			sprintf(Output,"%sdavDENDIS.txt",chemin);
		else
			sprintf(Output,"davDENDIS.txt");
		EcritFile(Output,1,P->dav_1knn,N_prototype);
		if (chemin!=NULL)
			sprintf(Output,"%sCostDENDIS.txt",chemin);
		else
			sprintf(Output,"CostDENDIS.txt");
		EcritFile(Output,1,P->Cost_n,N_prototype);

	}

	Delete2D(Liste,MaxCluster);
	LibereParamFCNNrule(P);
	delete P;

	delete []Liste_instance;

	return N_prototype;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int SamplingApproach_DENDIS(double **data,int *Liste_instance,int &N_instance,int *weight,
					   int **Liste,double &dmax,ParamFCNNrule *P)
{
struct Internal *I;
int terminate;
int New_instance;

	if (P==NULL) 
		return 0;
	
	I = new Internal;
	//on initialise les data et on extrait le premier pattern soit aléatoire soit par min min.
	New_instance = InitProgramme(data,I,P);//New instance est le 1ere prototype, un des n data...

	N_instance = 0;

	terminate = 0;
	N_instance = 0;
	while (!terminate)
	{
		if (New_instance != -2)//vient de l'étape inférieure où on a stoppé...
		{

			N_instance++;//on augmente ici car on commence avec 1 instance.
			InitBoucle(I,N_instance-1,New_instance);//on initialise certains paramètres et surtout on disable les prototypes (cf. partie optimisation) dont les patterns ne peuvent pas être le futur prototype.
		}
		New_instance = FindNewInstanceMassique(data,I,P,N_instance);//détermine le nouveau prototype et on réajuste les poids sur N_instance - 1. Le coeur du programme
		if (New_instance >= 0)//caution
		{
			CalculCost(I,P,N_instance);//optionnel ici.
			ProcessCriteria(N_instance,I,P);//on déduit  le max(max) et min(max) à partir des infos de la ligne précédente. (on a la distribution des max sur s)
			UpdateAlgorithm(I,N_instance,New_instance);//l'instance sélectionnée ne peut plus être sélectionné..elle est marquée. On stocke le prototype sélectionnée dans le vecteur des prototypes.
			terminate = TestFin_DENDIS(I,P,N_instance);//test de la fin suivant le critère d'occurence et si mode Auto sur la distance de référence.
			if (terminate == 0)
				Update1knnPrototypeMasse(data,I,P,N_instance);//Update des 1knn entre prototype et la moyenne. En mode auto, on a tjs besoin cf inégalité triangulaire
		}
		
		if (terminate == 1)// Terminé suivant les critères paramétrés, on va dégraisser ou engraisser le nombre de prototypes.
		{
			AttributionPattern(N_instance,Liste,I,P);//on remplit les patterns attachés à chaque instance. Aucun calcul de distance. Chaque pattern connaît son prototype.
			SortPrototypes_DENDIS(I,P,N_instance,Liste,P->NbPointsNoise);//on enlève les protoypes avec le LABEL NOISE.
		}
	}

	if (P->Barycentre == 1)//dans ce cas on prend le barycentre..C'est une option inutilisée. Au lieu de prendre le proto on prend le barycentre des patterns attachés.
		UpdateWithGravity_DENDIS(data,I,P,N_instance);
	
	
	for (int i=0;i<N_instance;i++)//aspect informatique uniquement
	{	
		Liste_instance[i] = I->Liste_selectionne[i];
		weight[i] = I->Weight[i];
	}

	LibereMemoireAlgorithme(I,P);
	delete I;
	
	return N_instance;
}
