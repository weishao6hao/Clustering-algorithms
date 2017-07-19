#include "math.h"
#include "time.h"
#include <stdlib.h>

#include "ADN_Esclave1.h"
#include "ADN_Esclave2.h"
#include "ADN_Esclave3.h"
#include "ADN_util.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int SamplingApproachP_DIDES(double **data,int n,int p,double Granularity,double **prototype,int &N_instance,int *weight, char *chemin,int *option)	
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
	P->MaxCluster = MaxCluster;
	P->Granularity = Granularity;
	P->NbPointsNoise = (int)(0.1 * P->Granularity* (double)(n));//toujours plus petit que la granularité
	P->NbPointsNoise = max(1,P->NbPointsNoise);//bruit
	P->Barycentre = 0;
	P->Init_Aleatoire = 0;
	P->Mass_OptionFinal = 0;
	P->Option_Alea_Cleaning_2 = option[0];
	P->Init_Aleatoire = option[1];
	P->Option_Assymptotique = 1;
	P->Auto = 1;
	P->L_pattern = 1;
	P->Dist_Noise = 1.2;
	P->OccMin = ((double)(n) * P->Granularity);
	P->NbPointsMin = (int)(0.01 * (double)(n));
	P->Seuil_massique = 0.2;
	P->Mass_WeightCriterion = 0;//On ne veut pas prendre en compte la masse si 0.
	P->DensiteMin = (int)((double)(n)*P->Granularity);//le nombre de patterns par occurence.

	N_prototype = SamplingApproach_DIDES(data,Liste_instance,N_instance,weight,Liste,dmax,P);
	N_prototype = SelectinstanceWithListe(data,p,Liste_instance,N_prototype,0,NULL,prototype);

	if (option[2] == 1)
	{
		if (chemin!=NULL)
			sprintf(Output,"%sdmaxDIDES.txt",chemin);
		else
			sprintf(Output,"dmaxDIDES.txt");
		EcritFile(Output,1,P->d_1knn,N_prototype);
		if (chemin!=NULL)
			sprintf(Output,"%sdavDIDES.txt",chemin);
		else
			sprintf(Output,"davDIDES.txt");
		EcritFile(Output,1,P->dav_1knn,N_prototype);
		if (chemin!=NULL)
			sprintf(Output,"%sCostDIDES.txt",chemin);
		else
			sprintf(Output,"CostDIDES.txt");
		EcritFile(Output,1,P->Cost_n,N_prototype);
	}

	Delete2D(Liste,MaxCluster);
	LibereParamFCNNrule(P);
	delete P;
	delete []Liste_instance;

	return N_prototype;
}

/*
data: input des données
Liste_instance: la liste des prototypes sélectionnés parmi les data
N_instance: le nombre d'instance
weight: le poids associé à chaque instance
dmax: max(max) final
P: c'est une structure qui permet de guider l'algorithme (basic version, sauvegarde,...).
*/
int SamplingApproach_DIDES(double **data,int *Liste_instance,int &N_instance,int *weight,int **Liste,double &dmax,ParamFCNNrule *P)
{
struct Internal *I;
int terminate;
int Stop_RatioMassique, Stop_Inertie;
int New_instance, Out, N_prototype;

	if (P==NULL) 
		return 0;
	
	I = new Internal;
	//on initialise les data et on extrait le premier pattern soit aléatoire soit par min min.
	New_instance = InitProgramme(data,I,P);//New instance est le 1ere prototype, un des n data...

	Stop_RatioMassique = Stop_Inertie = 0; N_instance = 0; N_prototype = 0;

	terminate = 0;
	N_instance = 0;
	while (!terminate)
	{
		N_instance++;//on augmente ici car on commence avec 1 instance.
		InitBoucle(I,N_instance-1,New_instance);//on initialise certains paramètres et surtout on disable les prototypes (cf. partie optimisation) dont les patterns ne peuvent pas être le futur prototype.		
		New_instance = FindNewInstance_DIDES(data,I,P,N_instance);//détermine le nouveau prototype et on réajuste les poids sur N_instance - 1. Le coeur du programme

		if (New_instance != -1)//caution
		{
			CalculCost(I,P,N_instance);//optionnel ici.
			ProcessCriteria(N_instance,I,P);//on déduit  le max(max) et min(max) à partir des infos de la ligne précédente. (on a la distribution des max sur s)
			UpdateAlgorithm(I,N_instance,New_instance);//l'instance sélectionnée ne peut plus être sélectionné..elle est marquée. On stocke le prototype sélectionnée dans le vecteur des prototypes.			
			terminate = TestFinDIDES(I,P,N_instance);//test de la fin suivant le critère d'occurence et si mode Auto sur la distance de référence.
			if (terminate == 0)//ce n'est pas terminé.
			{
				Update1knnPrototype(data,I,P,N_instance);//Update des 1knn entre prototype et la moyenne. En mode auto, on a tjs besoin cf inégalité triangulaire
				if	(Stop_RatioMassique == 0)//c'est calculé que si non atteint.
					Stop_RatioMassique = RatioMassiqueAchieved(N_instance,I,P);	//calcul du mass ratio et test de la limite sur ce critère							
				
				if (Stop_RatioMassique == 1)//les 2 sont atteints. Le Ratio est bon et l'inertie est stable.
				{
					if (I->Ref_OK == 0)
					{
						Out = ModeliseFunction(I,P,N_instance);//modélisation par la fonction puissance et calcul de la distance de référence qui va servir de threshold.			
						if (Out!=0)
						{
							N_prototype = N_instance;//on est allé plus loin que nécessaire...
							N_instance = Out;
							terminate = TestFinDIDES(I,P,N_instance);//c'est termine.
						}
					}
				}
			}
			else
				terminate = 1;
		}
		else
			terminate = 1;
		
		if (terminate == 1)// Terminé suivant les critères paramétrés.
		{				
			UpdateWeight_DIDES(I,P,N_instance,N_prototype,data);

			if (P->d_1knn)
				dmax = P->d_1knn[N_instance - 1];
			if (Liste)
				AttributionPattern(N_instance,Liste,I,P);//on remplit les patterns attachés à chaque instance. Aucun calcul de distance. Chaque pattern connaît son prototype.
			
			if (P->Option_Alea_Cleaning_2 == 1 || P->Option_Alea_Cleaning_1==1)//si traitement spécifique: Bruit + Bias probabilité	(option high densite OU elagage)
				N_instance = ManageProbabiliteCluster(Liste,N_instance,I,P);//on va gérer les high densité et la partie enlevé (c'est un paramètre) d'élagage	
			
			if (Liste!=NULL)
				SortPrototypes_DIDES(I,P,N_instance,Liste,0);//on enlève les protoypes avec le LABEL NOISE.				
		}

	}

	for (int i=0;i<N_instance;i++)//aspect informatique uniquement
	{	
		Liste_instance[i] = I->Liste_selectionne[i];
		weight[i] = I->Weight[i];
	}

	if (P->Barycentre == 1)//dans ce cas on prend le barycentre..C'est une option inutilisée. Au lieu de prendre le proto on prend le barycentre des patterns attachés.
		UpdateWithGravity(data,I,P,N_instance);
	
	LibereMemoireAlgorithme(I,P);
	delete I;
	
	return N_instance;
}


//  *****************   UtilDIDES  *****************************

int TraiteClusterOver(double **data,int n,double **data_int,int *Liste_active, int *weight_int,
					  int N_instance,int Groupe_actif,double DensiteMin,struct Internal *I,ParamFCNNrule *Pg)
{
int compteur,j,k;
int nbpattern,Output;
int pattern;
//double dmax;
int **Liste_voisins;
int n_instance;
int *fct_inverse;
int *Liste_candidat;
int *TListe, *TPoids;

	fct_inverse = new int[n];
	Liste_candidat = new int[n];
	compteur = 0;		
	TListe = new int[n];
	TPoids = new int[n];

	//on va récupérer les data correspondantes.
	for (j=0;j<N_instance;j++)
	{
		if (I->Tr->G_actif[j] == Groupe_actif)
		{
			Liste_candidat[compteur] = I->Liste_selectionne[j];//le vrai i.
			fct_inverse[I->Liste_selectionne[j]] = j;
			for (k=0;k<Pg->p;k++)
				data_int[compteur][k] = data[Liste_candidat[compteur]][k];
			compteur++;
		}
	}
	
	nbpattern = (int)( (double)(I->Tr->weight[Groupe_actif])/DensiteMin);
	nbpattern = max(nbpattern, I->Occmin);			
	Pg->n = compteur; 
	Pg->MaxCluster=nbpattern;
	Liste_voisins = Reservmemint(nbpattern,compteur);
	
	//Il faut remettre!
//	Output = SamplingApproach_1(data_int,Liste_active,Output,weight_int,Liste_voisins,dmax,Pg);
Output=0;	
	for (j=0;j<Output;j++)//les références sélectionnées.
	{
		pattern = Liste_candidat[Liste_active[j]];
		TListe[j] = pattern;
		TPoids[j] = 0;
		for (int k=0;k<weight_int[j];k++)//on calcule son poids...
		{
			n_instance = Liste_voisins[j][k];
			n_instance = Liste_candidat[n_instance];//c'est le vrai pattern.
			TPoids[j] = TPoids[j] +  I->Weight[fct_inverse[n_instance]];
		}
	}

	for (j=0;j<compteur;j++)//le nombre de 
	{
		n_instance = fct_inverse[Liste_candidat[j]];
		I->Weight[n_instance]=0;
	}

	Delete2D(Liste_voisins,nbpattern);
	delete []Liste_candidat;
	delete []fct_inverse;
	delete []TListe;
	delete []TPoids;

	return 1;
}

int AleatoireWeight(int *liste,int n,struct Internal *I)
{
int i, j, Max, reel;
int *occ;
int random;
int value;

	Max = 0;
	for (i=0;i<n;i++)
	{
		if (I->Weight[liste[i]] > Max)
			Max = I->Weight[liste[i]];
	}
	occ = new int[n*Max];

	reel = 0;
	for (i=0;i<n;i++)
		for (j=0;j<I->Weight[liste[i]];j++)
		{
			occ[reel] = i;
			reel++;
		}
	random = rand()%(reel-1);
	value = occ[random];
	delete []occ;

return value;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
int TraiteClusterLess(int **Liste,int N_instance,int &pool,int Groupe_actif,
					  int Nbp_Min,struct Internal *I,ParamFCNNrule *P)
{
int i,j;
int *Liste_pattern,*Liste_instance;
int Ref, Pattern;
int index,NewPoids;
int poids_total,PoidsPartiel;
int Nbre_associe,Nbre_instance_ajoute;
double Ratio;

	srand (time (NULL));
/*
	poids_total = I->Tr->weight[Groupe_actif];
	Nbre_associe = poids_total / Nbp_Min;//nombre de patterns par prototypes.
	Nbre_instance_ajoute = Nbre_associe-I->Tr->occ[Groupe_actif];//le nombre de prototypes à rajouter.
*/	
	poids_total = I->Tr->weight[Groupe_actif];
	Nbre_associe = I->Tr->weight[Groupe_actif] / Nbp_Min;//nombre de patterns par prototypes.
	Nbre_instance_ajoute = Nbp_Min - I->Tr->occ[Groupe_actif];
	if (Nbre_associe == 0 || Nbre_instance_ajoute<=0)
		return 1;


	Liste_pattern = new int[P->n];
	Liste_instance = new int[P->n];

	//ici on regarde les instances attachées au même groupe.
	index = 0;
	for (i=0;i<N_instance;i++)
		if (I->Tr->G_actif[i] == Groupe_actif)
		{
			Liste_instance[index] = i;
			index++;
		}

//	return 1;
	for (i=0;i<Nbre_instance_ajoute;i++) //On veut Nbp_Min en tout instance donc on en selectionne Nbp_min - les existantes dans le groupe.
	{
		if (I->Tr->occ[Groupe_actif] > 1)// A AMELIORER!!!!
//			Ref = rand()%(int)(I->Tr->occ[Groupe_actif]-1);//on extrait aléatoirement l'instance dans le groupe, c'est le numéro/
			Ref = AleatoireWeight(Liste_instance,index,I);
		else
			Ref = 0;//c'est le numéro du pattern attaché...

		Ref = Liste_instance[Ref];//le numéro d'instance comme découvert par l'algorithme.
		
		if (I->Weight[Ref] > 1)
			Pattern = rand()%(I->Weight[Ref]-1);//on extrait un de ces patterns attachés dans la liste.
		else
			Pattern = 0;

		I->Liste_selectionne[pool] = Liste[Ref][Pattern];//on place la nouvelle instance dans la liste.
		
		if (i==Nbp_Min - 1)
			Nbre_associe = poids_total - ((Nbp_Min-1)*Nbre_associe);

		if (Nbre_associe < 1) continue;

		if (Liste[pool]!=NULL)
			delete []Liste[pool];
		Liste[pool] = new int[Nbre_associe];
		
		for (j=0;j<Nbre_associe;j++)
			Liste[pool][j] = I->Liste_selectionne[pool];//on met le même pattern dans la liste attaché.	

		I->Weight[pool] = Nbre_associe;//on augment le poids de cette instance
		pool++;//on augmente le nombre d'instance.
	}

	NewPoids = poids_total - Nbre_associe*Nbre_instance_ajoute;//on a rajouté des instances qui ont un certain poids, on les déduit du poids total.
	Ratio = (double)(NewPoids)/(double)(poids_total);//on diluer chacune des instances présentes.
	PoidsPartiel = Nbre_instance_ajoute*Nbre_associe;//le poids à partager avec les prototypes de références.

	for (i=0;i<I->Tr->occ[Groupe_actif];i++)//on parcours les instances du groupes.
	{
		Ref = Liste_instance[i];//dans le groupe.
		if (i<I->Tr->occ[Groupe_actif]-1)
		{
			I->Weight[Ref] = (int)((double)(I->Weight[Ref]) * Ratio);
			I->Weight[Ref] = max(1,I->Weight[Ref]);
			PoidsPartiel += I->Weight[Ref];
		}
		else
		{
			I->Weight[Ref] = poids_total - PoidsPartiel;
			I->Weight[Ref] = max(1,I->Weight[Ref]);
		}
	}

	delete []Liste_pattern;
	delete []Liste_instance;

return 1;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ModifieContenance_Sampling_1(int **Liste,int N_instance,struct Internal *I,ParamFCNNrule *P)
{
int i;
int *Liste_active, *weight_int;
int pool=0;
//ParamFCNNrule *Pg;
double s, Ecart, Min, Max;
double s_O, Ecart_O, Min_O, Max_O;
double densite[1000], occ[1000];
int *booleen;


	weight_int = new int[P->MaxCluster];
	booleen = new int[P->MaxCluster];
	Liste_active = new int[P->MaxCluster];
//	Pg = new ParamFCNNrule;
//	InitSingleParamFCNNrule(Pg);
	/*
	Pg->Save=0;
	Pg->n = P->n;
	Pg->p = P->p;
	Pg->Auto=0;
	Pg->fdist = P->fdist;
	*/	
	//	DensiteMin = P->DensiteMin; 
	I->Occmin = (int)P->OccMin;
	pool = N_instance;
	
	for (i=0;i<I->N_groupe;i++)//pour chaque groupe.
	{
		if (I->Tr->weight[i] > (int)(P->Granularity*(double)(P->n)))
			booleen[i] = 1;
		else
			booleen[i] = 0;

		densite[i] = ((double)I->Tr->weight[i]/(double)I->Tr->occ[i]);//nombre d'individus par instance.		
		occ[i] = I->Tr->occ[i];
	}

	CalculusStatVectorBool(densite,I->N_groupe,booleen,s,Ecart,Min,Max);//on regarde les densites.	
	CalculusStatVectorBool(occ,I->N_groupe,booleen,s_O, Ecart_O, Min_O, Max_O);//on regarde les occurences.
	
	//on le cas des Sousreprésentés en nombre d'instance: on augmente le nombre d'instance...
	if (P->Option_Alea_Cleaning_2 == 1)
	{
		for (i=0;i<I->N_groupe;i++)//pour chaque groupe.
		{
			if (booleen[i]==0)
				continue;
			double RR = (double)(I->Tr->weight[i])/((double)(P->Granularity)*(double)(P->n));			
			double OccMin = (double)(I->Tr->weight[i]) / s;
			OccMin = fmin(OccMin,RR);
//			OccMin = fmax(OccMin, RR);//On regarde ceux qui sont trés dense.
			if (densite[i] > s)//le poids divisé par le nombre d'occurence.
			{		
//					if (I->Tr->occ[i]<=MinPattern)//On teste s'il y en a pas assez de représentativité
//					{
						TraiteClusterLess(Liste,N_instance,pool,i,(int)OccMin,I,P);
//					}
			}
		}
	}
	

	//on va diminuer le nombre d'instance
//	if (P->Option_Alea_Cleaning_1 == 1)
//	{
		//on va d'abord traiter le cas des Surreprésentés en nombre d'instance
//		for (i=0;i<I->N_groupe;i++)//pour chaque groupe.
//		{
//			Densite = ((double)I->Tr->weight[i]/(double)I->Tr->occ[i]);//nombre d'individus par instance.		
//			if (Densite < DensiteMin && I->Tr->occ[i]>I->Occmin)//On teste s'il y en a trop.	
//				TraiteClusterOverNew(data,Liste,N_instance,i,DensiteMin,I,P);

//		}
//	}

	N_instance = pool;	

//	LibereParamFCNNrule(Pg);
//	delete Pg;
	delete []weight_int;
	delete []Liste_active;	
	delete []booleen;

	return N_instance;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
int SearchNextSampling_1(double **D_cluster,int N_instance,double Seuil,int Instance_start,int *busy,
			   int *groupe,int actif,int &Nbre)
{
int i;

	if (Nbre == 0 && busy[Instance_start]==0)//Si non intégré on initialise.
	{
		Nbre = 1;
		busy[Instance_start] = 1;
		groupe[Instance_start] = actif;
	}

	if (groupe[Instance_start]<actif)//On l'intégère que si elle n'est pas déjà intégrée....
		return 0;
	
	for (i=0;i<N_instance;i++)
	{
		if (busy[i]==0)//il ne faut pas qu'il soit dans un groupe inférieure et qu'il soit libre. busy garantie les 2.
			if (D_cluster[Instance_start][i] <= Seuil) //on le regarde s'il n'est pas déjà dans un autre groupe.
			{
				busy[i] = 1;
				Nbre = Nbre + 1;
				groupe[i] = actif;
				SearchNextSampling_1(D_cluster,N_instance,Seuil,i,busy,groupe,actif,Nbre);
			}
	}

return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int FindConnex_Sampling_1(int N_instance,struct Internal *I)
{
int i;
int *Busy;
int Groupe_actif, Groupe_actifWithoutNoise,Nbre;
int occ_instance;

	Busy = new int[N_instance];
	Groupe_actif = 0;
	occ_instance = 0;
	for (i=0;i<N_instance;i++)
	{
		I->Tr->G_actif[i] = -1;
		I->Tr->weight[i] =0;
		I->Tr->occ[i] = 0;
		Busy[i] = 0;
	}

	for (i=0;i<N_instance;i++)
	{
		if (I->Tr->noise[i] == 1)//détecter comme bruit donc c'est un groupe à part.
		{
			Busy[i] = 1;
			occ_instance++;
		}
	}

	Groupe_actifWithoutNoise = 0;
	for (i=0;i<N_instance;i++)
	{
		Nbre = 0;
		SearchNextSampling_1(I->D_cluster,N_instance,I->Treshold_groupe,i,Busy,I->Tr->G_actif,Groupe_actif,Nbre);
		if (Nbre > 0) 
		{
			Groupe_actif++;
			Groupe_actifWithoutNoise++;
		}
		occ_instance = occ_instance + Nbre;
		
		if (occ_instance == N_instance)
			break;
	}

	I->N_groupe = Groupe_actif;
	
	//on somme les poids respectifs de chaque groupe.
	for (i=0;i<N_instance;i++)
	{
		if (I->Tr->G_actif[i] != -1)
		{
			I->Tr->weight[I->Tr->G_actif[i]] += I->Weight[i];
			I->Tr->occ[I->Tr->G_actif[i]] += 1;//on augmente le nombre d'occurence.
		}
	}

//	printf("Nb Groupe = %d\n",Groupe_actifWithoutNoise);

	delete []Busy;	

	return Groupe_actif;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int ManageProbabiliteCluster(int **Liste,int N_instance,struct Internal *I,ParamFCNNrule *P)
{
int i;
//int Nb_groupe;
double Ave_withnoise=0;
int N_Consistant=0;
	

	if (P->Auto==0)
		return N_instance;
	
	I->Nb_bruit = 0;
	if (P->Option_Alea_Bruit == 1)//on repère les instances bruités.
	{
		Ave_withnoise = 0;
		N_Consistant = 0;
		for (i=0;i<N_instance;i++)
		{
			if (I->Weight[i]>P->NbPointsNoise)//point bruit
			{
				Ave_withnoise = Ave_withnoise + I->d1nn_proto[i];
				N_Consistant++;
			}
		}
		if (N_Consistant > 0)		
			Ave_withnoise = (double)(Ave_withnoise) / (double)(N_Consistant); 
		else
			Ave_withnoise = 0;

		for (i=0;i<N_instance;i++)//on balaye les instances.
		{
			if (I->d1nn_proto[i]>=2*I->Treshold_groupe)//point bruit.)//on compare la distance de son pattern le plus éloigné à la moyenne des 1knn entre prototypes.
			{
				if (I->Weight[i]<=P->NbPointsNoise)//point bruit.
				{
					I->Tr->noise[i] = 1;//affecté comme bruit.	
					I->Nb_bruit++;//le nombre de points bruit est augmenté.
				}//fin du for if
			}//fin du for if
		}//fin du for i

	}//fin du for P.
//	printf("%d \n",I->Nb_bruit);


	FindConnex_Sampling_1(N_instance,I);
	N_instance = ModifieContenance_Sampling_1(Liste,N_instance,I,P);


return N_instance;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int FindThresholdNoise(int N_instance,int MinCard, struct Internal *I)
{
int i;
//double Moy,Ecart,Min,Max;
//double *D_selectionne;
//int *weight;
int reel;
int s=0;

	s = 0;
	reel = 0;
	for (i=0;i<N_instance;i++)
	{
		if (I->Weight[i] < MinCard && I->Pre_distance[i]>I->D_distance_Extreme)
		{
			s+=I->Weight[i];
			reel++;
		}
	}

	if (reel > 0)
		s = (int)((double)(s) / (double)(reel));
	else
		s = 1;

return s;
}

double FindThresholdDistance(int N_instance,double d_max,struct Internal *I)
{
int i;
double Moy,Ecart,Min,Max;
double *D_selectionne, threshold;
int *weight;
int reel;

	D_selectionne = new double[N_instance];
	weight = new int[N_instance];
	reel = 0;
	for (i=0;i<N_instance;i++)
	{
		if (I->Pre_distance[i] < 2*d_max && I->Pre_distance[i]>I->D_distance_Extreme)
		{
			D_selectionne[reel] = I->d1nn_proto[i];
			weight[reel] = I->Weight[i];
			reel++;
		}
	}

	if (reel > 0)
	{
		CalculusStatVector(D_selectionne,reel,Moy,Ecart,Min,Max);
		threshold = Moy;
	}
	else
		threshold = I->D_distance_Extreme;

		
	delete []D_selectionne;
	delete []weight;
	return threshold;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int TestFinDIDES(struct Internal *I,ParamFCNNrule *P,int &N_instance)
{
int MinCard;

	if (N_instance >= (P->MaxCluster) || N_instance == P->n-1)
		return 1;
	
	if (P->Auto == 1 && I->Ref_OK==1)
	{
		if (P->d_1knn[N_instance - 1] < I->Ref_distance)
		{
			MinCard = (int)(P->Granularity * (double)(P->n));
			I->Treshold_groupe = FindThresholdDistance(N_instance,P->d_1knn[N_instance - 1],I);
			P->NbPointsNoise = FindThresholdNoise(N_instance,MinCard,I);
			return 1;
		}
	}

	return 0;
}

