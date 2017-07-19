#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ADN_Esclave1.h"
#include "ADN_Esclave2.h"
#include "ADN_Esclave3.h"
#include "ADN_Esclave4.h"
//#include "dbscan.h"
#include "ADN_util.h"

//Pour chaque point on calcule 
//1. Distance moyenne dans son cluster. =>  (a)
//2. Distance moyenne avec tous les autres clusters. d_exter1,2...n
//3. On garde la plus petite distance moyenne. b = min (d_exter1,...)
//4. (a-b)/fabs(a,b);

double get(double **dist,int pattern, int *liste,int n,int *weight)
{
int i;
double s, reel;
double *P;

	P = new double[n];

	reel = 0;
	s=0;

	if (weight == NULL)
	{
		for (i=0;i<n;i++)
			P[i] = 1;
	}
	else
	{
		for (i=0;i<n;i++)
			P[i] = weight[liste[i]];
	}

	for (i=0;i<n;i++)
	{
		if (liste[i] != pattern)
		{
			s += P[i]*dist[pattern][liste[i]];
			reel+=P[i];
		}
	}
	
	delete P;

	if (reel != 0)
		return s/=(double)(reel);
	else
		return 0;
}
/*
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Silhouette_index(double **data,int n,int p,int n_cluster,int **liste,int *weight,int *weight_proto)
{
int i,j,k;
double **dist;
double a, b, si,inter, V;
int pattern_i;
long Masse = 0;

	if (n <= n_cluster) 
		return 0;

	dist = Reservmemdouble(n,n);
	if (!dist) return 0;

	GenereMatriceDistance(data,n,p,NULL,dist,d_Eucl);
	si = 0;
	for (i=0;i<n_cluster;i++)
	{
		for (j=0;j<weight[i];j++)
		{
			pattern_i=liste[i][j];
			a = get(dist,pattern_i,liste[i],weight[i],weight_proto);
			b = 10000; 
			for (k=0;k<n_cluster;k++)
				if (k!=i)
				{
					inter = get(dist,pattern_i,liste[k],weight[k],weight_proto);
					if (inter < b)
						b = inter;
				}

			if (fmax(a,b)!=0)
				V = (b-a)/fmax(b,a);
			else
				V = 0;

			if (weight_proto!=NULL)
			{
				si += weight_proto[pattern_i]*V;
				Masse += weight_proto[pattern_i];
			}
			else
			{
				si += V;
				Masse += 1;
			}
		}
	}

	si /= (double)(Masse);
	Delete2D(dist,n);
	
	return si;
}


double Ball_index(double **data,int n,int p,int n_cluster,int **liste,int *weight,int *weight_cluster)
{
int i,j,k;
int pattern,*P;
double **Gravity;
double Cost, dist;
int Masse;

	P = new int[n];
	Gravity = Reservmemdouble(n_cluster,p);
	for (i=0;i<n_cluster;i++)
		for (k=0;k<p;k++)
			Gravity[i][k] = 0;
	
	if (weight_cluster == NULL)
	{
		for (i=0;i<n_cluster;i++)
			P[i] = 1;
	}
	else
	{
		for (i=0;i<n_cluster;i++)
			P[i] = weight_cluster[i];
	}

	for (i=0;i<n_cluster;i++)
	{
		Masse = 0;
		if (weight[i] == 0) continue;
		for (j=0;j<weight[i];j++)
		{
			pattern = liste[i][j];
			for (k=0;k<p;k++)
				Gravity[i][k] += P[i]*data[pattern][k];
			Masse += P[i];
		}
		
		if (Masse != 0)
		{
			for (k=0;k<p;k++)
				Gravity[i][k] /= (double)(Masse);
		}
		
	}

	Cost=0;
	for (i=0;i<n_cluster;i++)
	{
		if (weight[i] == 0) continue;
		dist = 0;
		Masse = 0;
		for (j=0;j<weight[i];j++)
		{
			pattern = liste[i][j];
			for (k=0;k<p;k++)
				dist += P[i]*(Gravity[i][k] - data[pattern][k])*(Gravity[i][k] - data[pattern][k]);
			Masse += P[i];
		}
		if (Masse == 0) continue;
		dist = dist / (double)(Masse);
		Cost += dist;
	}
	Cost /= (double)(n_cluster);

	Delete2D(Gravity,n_cluster);
	delete P;
return Cost;
}
*/
/*
ncluster:
weight: le poids de chaque cluster.
n_proto, liste_proto,weight_proto: la liste des prototypes, le nombre et le poids.
liste_weight: la liste des poids de chaque instance consideré dans la cluster.
*/
void GetWeight(int n_cluster, int *weight,int **liste,int *weight_proto,int n_proto,int **liste_weight)
{
int i,j,k;
int pattern;


	for (i=0;i<n_cluster;i++)
	{
		if (liste_weight[i] != NULL)
			delete liste_weight[i];
		
		liste_weight[i] = new int[weight[i]];
		
		for (j=0;j<weight[i];j++)
		{
			pattern = liste[i][j];
			for (k=0;k<n_proto;k++)
				if (pattern == k)
				{
					liste_weight[i][j] = weight_proto[k];
					break;
				}
		}
	}
}
/*
double ComparisonPartitionNuees(char *output,double **data,int n,int p,double **prototype,int *weightnnd,
								int n_proto,int *liste_instance_init,
								int n_clusterMax,double (*fdist)(double*,double*,int),int n_REP)
{
double **cluster_prototype, **cluster_data;
int *ambigu_proto, *ambigu_data;
int *weight_algoP, *weight_algoD;
int **liste_proto,**liste_data;
int alea;
double T = 0.99;
double RI, RI_moy,Ball_data,Ball_moy;
double Ball_proto, Ball_pmoy;
double Sil_index, Sil_indexp;
double Sil_index_pmoy,Sil_index_moy;
double Ball_proto_min, Ball_proto_max;
int n_cluster, n_cluster_d,n_cluster_p;
int **liste_ind_cluster,**liste_weight;
FILE *fich;

	Ball_proto= Ball_data = Sil_index = Sil_indexp = 0;
	cluster_prototype = Reservmemdouble(n_clusterMax,p);
	if (cluster_prototype == NULL) return -1;
	cluster_data = Reservmemdouble(n_clusterMax,p);
	if (cluster_data == NULL) return -1;

	ambigu_proto = new int[n];
	ambigu_data = new int[n];
	liste_ind_cluster = new int*[n_clusterMax];
	liste_weight = new int*[n_clusterMax];
	liste_proto = new int*[n_clusterMax];
	liste_data = new int*[n_clusterMax];
	weight_algoP = new int[n_clusterMax];
	weight_algoD = new int[n_clusterMax];


	for (int i=0;i<n_clusterMax;i++)
	{  liste_weight[i] =  liste_ind_cluster[i] = liste_proto[i] = liste_data[i] = NULL; }
	

	fich = fopen(output,"w");
	for (int C=2; C<n_clusterMax;C++)
	{
		n_cluster = C;
		RI_moy = 0;
		Ball_pmoy = Ball_moy = 0;
		Sil_index_pmoy = 0;
		Sil_index_moy = 0;
		Ball_proto_min = 100000;
		Ball_proto_max = 0;
		for (int REP=0;REP<n_REP;REP++)
		{
			RI = Ball_data = Ball_proto = Sil_index = Sil_indexp = 0;
			if (n_proto <= n_cluster) continue;
			//initialisation.
			srand(time(NULL));
			for (int i=0;i<n_cluster;i++)
			{
				alea = rand()%(n-1);
				for (int j=0;j<p;j++)
					cluster_prototype[i][j] = cluster_data[i][j] = data[alea][j];
			}

			n_cluster_p = C_NueesDynamiquesWeighted(100,0.01,prototype,n_proto,weightnnd,p,cluster_prototype,n_cluster,1,fdist,weight_algoP,0,liste_ind_cluster);//aléatoire			
			Ball_proto = Ball_index(prototype,n_proto,p,n_cluster_p,liste_ind_cluster,weight_algoP,weightnnd);
			Sil_indexp =	Silhouette_index(prototype,n_proto,p,n_cluster_p,liste_ind_cluster,weight_algoP,weightnnd);
			GenerePartitionFromCluster(data,n,p,cluster_prototype,n_cluster_p,fdist,liste_proto,weight_algoP,ambigu_proto,T);

			n_cluster_d = C_NueesDynamiquesWeighted(100,0.01,data,n,NULL,p,cluster_data,n_cluster,1,fdist,weight_algoD,0,liste_ind_cluster);//aléatoire			
			Ball_data = Ball_index(data,n,p,n_cluster_d,liste_ind_cluster,weight_algoD,NULL);
			Sil_index = Silhouette_index(data,n,p,n_cluster_d,liste_ind_cluster,weight_algoD,NULL);		
			GenerePartitionFromCluster(data,n,p,cluster_data,n_cluster_d,d_Eucl,liste_data,weight_algoD,ambigu_data,T);	
			//printf("(%d,%d)",n_cluster_p,n_cluster_d);
		
			if (Ball_proto<Ball_proto_min)
				Ball_proto_min = Ball_proto;
			if (Ball_proto>Ball_proto_max)
				Ball_proto_max = Ball_proto;

			//liste_proto (dans data), weight_algoP (cluster), weightnnd est FAUX
	//		RI = RandIndexFuzzy(n,liste_proto,liste_data,n_cluster_p,n_cluster_d,weight_algoP,weight_algoD,ambigu_proto,ambigu_data);
	//		RI = AdjRandIndex(n,liste_proto,liste_data,n_cluster_p,n_cluster_d,weight_algoP,weight_algoD);
	//		RI = MutualInformation(n,liste_proto,liste_data,n_cluster_p,n_cluster_d,weight_algoP,weight_algoD);
			RI_moy += RI;
			Ball_moy  += Ball_data;
			Ball_pmoy += Ball_proto;
			Sil_index_moy += Sil_index;
			Sil_index_pmoy += Sil_indexp;
		}
	
		Sil_index_pmoy /= (double)(n_REP);
		Sil_index_moy /= (double)(n_REP);
		Ball_moy /= (double)(n_REP);
		Ball_pmoy /= (double)(n_REP);
		RI_moy /= (double)(n_REP);
		printf("%d Ins:%d:I=%5.3lf\tBall=(%5.3lf,%5.3lf),=Ext(%5.3lf,%5.3lf),Sil=(%5.3lf,%5.3lf)\n",C,n_proto,RI_moy,Ball_moy,Ball_pmoy,Ball_proto_min,Ball_proto_max,Sil_index_moy,Sil_index_pmoy);
		fprintf(fich,"%d\t%d\t%d\t%5.3lf\t%5.3lf\t%5.3lf\t%5.3lf\t%5.3lf\n",n,C,n_proto,RI_moy,Ball_moy,Ball_pmoy,Sil_index_moy,Sil_index_pmoy);

	}

	fclose(fich);
	
	Delete2D(cluster_prototype,n_clusterMax);
	Delete2D(cluster_data,n_clusterMax);
	Delete2D(liste_data,n_clusterMax);
	Delete2D(liste_proto,n_clusterMax);
	Delete2D(liste_ind_cluster,n_clusterMax);
	Delete2D(liste_weight,n_clusterMax);
	delete ambigu_proto;
	delete ambigu_data;
	delete weight_algoP;
	delete weight_algoD;

	return (RI_moy)/(double)(n_REP);
}
*/
void GenereListeFromCategory(int n,int *cat,int n_cluster,int **liste,int *weight)
{
int i, classe;

	for (i=0;i<n_cluster;i++)
		weight[i] = 0;

	for (i=0;i<n;i++)
		weight[cat[i]]++;

	for (i=0;i<n_cluster;i++)
		liste[i] = new int [weight[i]];

	for (i=0;i<n_cluster;i++)
		weight[i] = 0;

	for (i=0;i<n;i++)
	{
		classe = cat[i];
		liste[classe][weight[classe]] = i;
		weight[classe]++;
	}
}
/*
double ComparisonRANDSCAN(double **data,int n,int p,double **prototype,int *weightnnd,int n_proto,
								int n_cluster,double (*fdist)(double*,double*,int),int n_REP,
								double *Feat_data,double *Feat_proto, double dmax)
{
double **cluster_prototype, **cluster_data;
double **prototrie, **datatrie;
int *ambigu_proto, *ambigu_data;
int *weight_algoP, *weight_algoD;
int *weight_algoP_A, *weight_algoD_A,*weight_trie;
int **liste_proto_A,**liste_data_A;
int **liste_proto,**liste_data;
int *cat_prototrie, *cat_datatrie;
int alea,i,OK_cluster;
double T = 1;
double RI, RI_moy,Ball_data,Ball_moy;
double Ball_proto, Ball_pmoy, r_proto;
int *cat_proto, *cat_data, bruit;
int Ncluster_d,Ncluster_p,n_trie;
double dist, r,n_ancien,n_ancien_proto, r_base;

	RI_moy = 0;
	Ball_pmoy = Ball_moy = 0;

	cat_proto = new int[n];
	cat_prototrie = new int[n];
	cat_data = new int[n];
	cat_datatrie = new int[n];
	weight_algoP = new int[n_proto];
	weight_algoD = new int[n_cluster];
	weight_algoP_A = new int[n_proto];
	weight_trie = new int[n_proto];
	ambigu_proto = new int[n];
	ambigu_data = new int[n];
	cluster_prototype = Reservmemdouble(n_cluster,p);
	cluster_data = Reservmemdouble(n_cluster,p);
	liste_proto = new int*[n_cluster];
	liste_proto_A = new int*[n_proto];
	liste_data = new int*[n_cluster];
	prototrie = Reservmemdouble(n_proto,p);
	datatrie = Reservmemdouble(n,p);


	for (int i=0;i<n_cluster;i++)
	{  liste_proto[i] = liste_data[i] = NULL; }
	for (int i=0;i<n_proto;i++)
	{  liste_proto_A[i] = NULL; }
	for (i=0;i<n;i++)
		ambigu_data[i] = 0;

	n_REP = 1;
	n_ancien = n;
	n_ancien_proto = n_proto;
	Ball_data=0;
	bruit = 1;

	for (int REP=0;REP<n_REP;REP++)
	{
		//initialisation. 
		srand(time(NULL));

		dist = 0.01;
		
		for (int z=1;z<50;z++)
		{
			OK_cluster = 0;
			n = n_ancien;
			n_proto = n_ancien_proto;
			r_proto = (1)/(1+dist);//dist = 1/1+r => dist + dist*r = 1; r = (1-dist)/dist.
			//r = 0.7;
			LaunchDBSCAN(prototype,n_proto,p,r_proto,0,cat_proto);
			Ncluster_p = N_clusterDBSCAN(cat_proto,prototype,n_proto,p,bruit,prototrie,cat_prototrie,weightnnd,weight_trie);
			Ball_proto = Ball_index(prototrie,n_trie,p,Ncluster_p,liste_proto,weight_algoP,NULL);//_A le poids des prototypes, P poids data

			printf("Ncluster = %d bruit = %d\n",Ncluster_p,bruit);
			if (Ncluster_p==1) break;

			if (Ncluster_p >=2 && Ncluster_p<20 && bruit<0.1*(double)(n))
			{
				OK_cluster = 1;
				r_base = r_proto;
				//
				GenereListeFromCategory(n_proto,cat_prototrie,Ncluster_p,liste_proto_A,weight_algoP_A);
				

//	n_cluster_p = C_NueesDynamiquesWeighted(100,0.01,prototype,n_proto,weightnnd,p,cluster_prototype,n_cluster,1,fdist,weight_algoP,0,liste_ind_cluster);//aléatoire			
//			Ball_proto = Ball_index(prototype,n_proto,p,n_cluster_p,liste_ind_cluster,weight_algoP,weightnnd);
		


//				printf("Ncluster = %d bruit = %d Ball=%5.2lf\n",Ncluster_p,bruit,Ball_proto);

				double pas=0.05;
				for (int j=0;j<=14;j++)
				{
					r = (1)/(1+(0.3+(double)(j)*pas)*dist);
					LaunchDBSCAN(data,n,p,r,0,cat_data);	//DBSCAN sur les data.
					n_trie = n;
					Ncluster_d = N_clusterDBSCAN(cat_data,data,n_trie,p,bruit,datatrie, cat_datatrie,NULL,NULL);//on sort des datatries et des categories tries.
//					printf("n=%d, ntrie=%d\n",n,n_trie);
					
					if (Ncluster_d >= 2 && abs(Ncluster_d-Ncluster_p)<=1)
					{
						//on part des prototypes et on recherche pour chaque cluster les patterns attractés.
						GenerePartitionFromCluster_Hierarchical(datatrie,prototrie,n_trie,p,liste_proto_A,weight_algoP_A,Ncluster_p,n_proto,d_Eucl,
											       liste_proto,weight_algoP,ambigu_proto,T);
						GenereListeFromCategory(n_trie,cat_datatrie,Ncluster_d,liste_data,weight_algoD);
								
						Ball_data = Ball_index(datatrie,n_trie,p,Ncluster_d,liste_data,weight_algoD,NULL);//_A le poids des prototypes, P poids data		
//						printf("DATA:d = %lf\tbruit=%d Ncluster=%d \n",dist/2,bruit,Ncluster_d);
					//	RI = RandIndexFuzzy(n_trie,liste_proto,liste_data,Ncluster_p,Ncluster_d,weight_algoP,weight_algoD,ambigu_proto,ambigu_data);
					//	RI = AdjRandIndex(n,liste_proto,liste_data,Ncluster_p,Ncluster_d,weight_algoP,weight_algoD);
						
						printf("C=(%d,%d) R(=%5.2lf,%5.2lf) RI:%5.2lf  Ball:(%5.2lf,%5.2lf)\n",Ncluster_p, Ncluster_d,r_proto, r,RI,Ball_proto,Ball_data);
					}
				}
			}
			if (OK_cluster==1)
			{
				dist+=0.005;
				break;
			}
			if (OK_cluster == 0)
				dist+=0.08;

			if (OK_cluster==0 && Ncluster_p<=100)
				dist+=0.005;
		}

	}
	
	n = n_ancien;
	n_proto = n_ancien_proto;
	Delete2D(cluster_prototype,n_cluster);
	Delete2D(cluster_data,n_cluster);
	Delete2D(liste_data,n_cluster);
	Delete2D(liste_proto,n_cluster);
	Delete2D(liste_proto_A,n_ancien_proto);
	Delete2D(prototrie,n_ancien_proto);
	Delete2D(datatrie,n);
	delete cat_prototrie;
	delete cat_datatrie;
	delete ambigu_proto;
	delete ambigu_data;
	delete weight_algoP;
	delete weight_algoD;
	delete weight_algoP_A;
	delete weight_trie;
	delete cat_proto;
	delete cat_data;

	return (RI_moy)/(double)(n_REP);
}
*/

int C_NueesDynamiquesWeighted(int N_iter,double Epsilon,double **data,int Nb_data,int *poids,int p,
							 double **cluster,int Nb_clust,int MinOccurence, 
							 double (*f_dist)(double *,double  *,int n), int *Cluster_weight,int Alea,int **Liste)
{
//int N_iter;
double Error;
int i_clust;
int i_dat;
int win;
double min_dist;
double dist;
double **Cluster_temporaire;
int i,j;
int boucle = 0;
int Nb_cluster_valide;
int RandomP;
int S_poids;
int PoidsCluster[100];
int *Belong;

	Error = 10000;
	S_poids = 0;
	Belong = new int[Nb_data];

	if (poids != NULL)
	{
		for (i=0;i<Nb_data;i++)
			S_poids += poids[i];
	}
	else
		S_poids = Nb_data; 
	
	if (Alea == 1)
	{
		for (i=0;i<Nb_clust;i++)
		{
			RandomP = rand()%(Nb_data-1);
			for (j=0;j<p;j++)
				cluster[i][j] = data[RandomP][j];
		}
	}

	if (N_iter == 0)
		return Nb_clust;

	Cluster_temporaire = Reservmemdouble(Nb_clust,p);
	
	if (Liste!=NULL)
	{
		for (i=0;i<Nb_clust;i++)
		{
			if (Liste[i]!=NULL) 
			{
				delete Liste[i];	
				Liste[i] = NULL;
			}
			PoidsCluster[i] = 0;
		}
	}

	while(boucle < N_iter || (Error > Epsilon) )
	{
		//initialisation du cluster temporaire à 0.
		for (i_clust = 0; i_clust < Nb_clust;i_clust++)
		{
			for (j=0;j<p;j++)
				Cluster_temporaire[i_clust][j] = 0;
		}

		//initialisation des poids à 0
		for (i_clust = 0; i_clust < Nb_clust;i_clust++)
			Cluster_weight[i_clust] = 0;

		//calcul des distances des data aux clusters.		
		for (i_dat = 0; i_dat < Nb_data;i_dat++)
		{
			min_dist = 10000;
			win = -1;
			for (i_clust = 0; i_clust < Nb_clust;i_clust++)
			{
				dist = f_dist(data[i_dat],cluster[i_clust],p);
				if (dist < min_dist)
				{
					min_dist = dist;
					win = i_clust;
				}
			}
	
			if (win != -1)
			{
				//on affecte la data i_dat au cluster win.
				if (poids)
				{
					for (j=0;j<p;j++)
						Cluster_temporaire[win][j] = Cluster_temporaire[win][j] + (double)(poids[i_dat])*data[i_dat][j];
				
					//on augmente le poids du cluster win.
					Cluster_weight[win] += poids[i_dat];
			
				}
				else
				{
					for (j=0;j<p;j++)
						Cluster_temporaire[win][j] = Cluster_temporaire[win][j] + data[i_dat][j];
					//on augmente le poids du cluster win.
					Cluster_weight[win] +=1; 
				}
			}
		}
		
		Error = 0;
		for (i_clust = 0; i_clust < Nb_clust;i_clust++)
		{
			if (Cluster_weight[i_clust] != 0)
			{
				for (j=0;j<p;j++)
					Cluster_temporaire[i_clust][j] = Cluster_temporaire[i_clust][j] / (double)(Cluster_weight[i_clust]);
			
				for (j=0;j<p;j++)
					Error = Error + (Cluster_temporaire[i_clust][j] - cluster[i_clust][j]) * (Cluster_temporaire[i_clust][j] - cluster[i_clust][j]);
				
				for (j=0;j<p;j++)
					cluster[i_clust][j] = Cluster_temporaire[i_clust][j];				
								
			}
/*
			else
			{
				RandomP = rand()%(Nb_data-1);
				for (j=0;j<p;j++)
					cluster[i_clust][j] = data[RandomP][j];
			}
*/
		}
		boucle++;
		
	}

	Nb_cluster_valide = Nb_clust;

	for (i_clust = 0; i_clust < Nb_cluster_valide;i_clust++)
	{
		if (Cluster_weight[i_clust] < max(1,MinOccurence))
		{
			for (int k=i_clust;k<Nb_cluster_valide-1;k++)
			{
				for (j=0;j<p;j++)
					cluster[k][j] = cluster[k+1][j];
				Cluster_weight[k] = Cluster_weight[k+1];
			}		
			Nb_cluster_valide--;
		}
	}
	
	Delete2D(Cluster_temporaire,Nb_clust);

	//ICI PB
	
	for (i=0;i<Nb_cluster_valide;i++)	
		PoidsCluster[i] = 0;
		
	for (i_dat = 0; i_dat < Nb_data;i_dat++)
	{
		min_dist = 10000;
		win = -1;
		for (i_clust = 0; i_clust < Nb_cluster_valide;i_clust++)
		{
			dist = f_dist(data[i_dat],cluster[i_clust],p);
			if (dist < min_dist)
			{
				min_dist = dist;
				win = i_clust;
			}
		}
	
		Belong[i_dat] = win;
		PoidsCluster[win] += 1;			
	}

	if (Liste!=NULL)
	{
		for (i_clust = 0; i_clust < Nb_cluster_valide;i_clust++)
		{
			if (Cluster_weight[i_clust])				
				Liste[i_clust] = new int[PoidsCluster[i_clust]];
		}
	
		for (i_clust = 0; i_clust < Nb_cluster_valide;i_clust++)
			Cluster_weight[i_clust] = PoidsCluster[i_clust]=0;

		for (i_dat = 0; i_dat < Nb_data;i_dat++)
		{
			win = Belong[i_dat];
			Liste[win][PoidsCluster[win]] = i_dat;				
			PoidsCluster[win]++;											
		}
	}

	for (i_clust = 0; i_clust < Nb_cluster_valide;i_clust++)
		Cluster_weight[i_clust] = PoidsCluster[i_clust];
	

	delete Belong;
	return Nb_cluster_valide;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//Pattern,n,p: la base initiale
//cluster: les clusters trouvés.
//N_instance: la taille du clustering
//
//liste: liste des individus pour chaque cluster.
//weight: le nombre d'individus par cluster

int GenerePartitionFromCluster(double **pattern,int n,int p,double **cluster,int N_instance,
							   double (*fdist)(double  *,double  *,int n),
							   int **liste,int *weight,int *ambigu,double T)
{
int i, j, winner;
double DistWinner;
double *dist;
int Nb_ambigu = 0;
int NbInt;
int W,Pattern;
FILE *fich;

	
	//FAIRE UNE SAUVEGARDE INTERMEDIAIRE...
	fich = fopen("FileInt","w");
	dist = new double[N_instance]; 
	
	for (i=0;i<N_instance;i++)//on regarde les S actifs.	
		weight[i] = 0;

	for (i=0;i<n;i++)//pour tous les éléments de T.
	{
		DistWinner = 10000000;
		ambigu[i] = 0;
		winner = -1;
		for (j=0;j<N_instance;j++)//on regarde les S actifs.	
		{
			dist[j] = fdist(pattern[i],cluster[j],p);
			if (dist[j] < DistWinner)//plus proche d'un prototype...
			{				
				DistWinner = dist[j];
				winner = j;//le prototype gagnant.	
			}
		}
			
		if (DistWinner != 0)
		{	
			for (j=0;j<N_instance;j++)
				if (j!=winner && (DistWinner/dist[j])>T)
					ambigu[i] = 1;
		}

		if (ambigu[i] == 1)
			Nb_ambigu++;

		fprintf(fich,"%d\t%d\t",winner,i);
		weight[winner] = weight[winner] + 1; 
	}
	fclose(fich);

		NbInt = 0;
		for (i=0;i<N_instance;i++)
		{						
			if (liste[i]!=NULL)
			{
				delete liste[i];
				liste[i] = NULL;
			}
			if (weight[i] > 0)
				liste[i] = new int[weight[i]];
			NbInt += weight[i];
			weight[i] = 0;
		}
	
		fich = fopen("FileInt","r");
		for (i=0;i<NbInt;i++)
		{
			fscanf(fich,"%d%d",&W,&Pattern);
			liste[W][weight[W]] = Pattern;
			weight[W]+=1;
		}
		fclose(fich);

delete dist;


return 	Nb_ambigu;
}

int GenerePartitionFromCluster_Hierarchical(double **data,double **pattern,int n,int p,int **liste_in,int *weight_in,int N_cluster,int N_instance,
											double (*fdist)(double  *,double  *,int n),
											int **liste,int *weight,int *ambigu,double T)
{
int i, j, k, winner;
double DistWinner;
double *dist;
int Nb_ambigu = 0;
double **cluster;
int *identite;
int pool;

	winner = 0;
	identite = new int[n];
	cluster = Reservmemdouble(N_instance,p);

	pool = 0;
	for (i=0;i<N_cluster;i++)//pour chaque cluster
	{
		for (j=0;j<weight_in[i];j++)//un nombre de prototype.
		{
			identite[pool] = i;
			for (k=0;k<p;k++)
				cluster[pool][k] = pattern[liste_in[i][j]][k];
			pool++;
		}
	}//CLUSTER est la matrice des instances.

	dist = new double[N_instance]; 
	
	for (i=0;i<N_cluster;i++)//on regarde les S actifs.	
		weight[i] = 0;

	for (i=0;i<n;i++)//pour tous les éléments de T.
	{
		DistWinner = 10000;
		ambigu[i] = 0;
		for (j=0;j<N_instance;j++)//on regarde les S actifs.	
		{
			dist[j] = fdist(data[i],cluster[j],p);
			
			if (dist[j] < DistWinner)//plus proche d'un prototype...
			{				
				DistWinner = dist[j];
				winner = identite[j];//le prototype gagnant.	
			}
		}

		if (DistWinner != 0)
		{	
			for (j=0;j<N_instance;j++)
			{
				if (identite[j]!=winner && (DistWinner/dist[j])>T)
					ambigu[i] = 1;
			}
		}

		if (ambigu[i] == 1)
			Nb_ambigu++;

		weight[winner] = weight[winner] + 1; 
	}//On CALCULE le poids de chaque cluster.
	
	//On reserve la mémoire exactement nécessaire.
	for (i=0;i<N_cluster;i++)
	{
		if (liste[i]!=NULL)
			liste[i] = NULL;
		if (weight[i]>0)
			liste[i] = new int[weight[i]];
	}

	for (i=0;i<N_cluster;i++)//on regarde les S actifs.	
		weight[i] = 0;

	for (i=0;i<n;i++)//pour tous les éléments de T.
	{
		DistWinner = 10000;
		ambigu[i] = 0;
		for (j=0;j<N_instance;j++)//on regarde les S actifs.	
		{
			dist[j] = fdist(data[i],cluster[j],p);
			
			if (dist[j] < DistWinner)//plus proche d'un prototype...
			{				
				DistWinner = dist[j];
				winner = identite[j];//le prototype gagnant.	
			}
		}

		if (DistWinner != 0)
		{	
			for (j=0;j<N_instance;j++)
			{
				if (identite[j]!=winner && (DistWinner/dist[j])>T)
					ambigu[i] = 1;
			}
		}

		if (ambigu[i] == 1)
			Nb_ambigu++;

		liste[winner][weight[winner]] = i;
		weight[winner] = weight[winner] + 1; 
	}

	delete identite;
	delete dist;
	Delete2D(cluster,N_instance);

return 	Nb_ambigu;
}
