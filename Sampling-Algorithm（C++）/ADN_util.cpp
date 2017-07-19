#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "math.h"
#include "time.h"
#include "ADN_Esclave1.h"
#include "ADN_util.h"

extern int *G_1Dx,**Gint_2dy;
extern double *G_1Ddx,**G_2Dy;


int EcritFichier1D(char *fichier,int n,double *data)
{
int i;
FILE *fich;

	fich = fopen(fichier,"w");
	if (!fich) return 0;

	for (i=0;i<n;i++)
		fprintf(fich,"%6.4lf\n",data[i]);
	
	fclose(fich);
	return 1;
}

int LoadFichier1D(char *fichier,int n,double *data)
{
int i;
FILE *fich;

	fich = fopen(fichier,"r");
	if (!fich) return 0;

	for (i=0;i<n;i++)
		fscanf(fich,"%lf",&(data[i]));
	
	fclose(fich);
	return 1;
}

void GenereMatriceDistance(double **data,int n,int p,int *liste,double **dist,double (*f_dist)(double *,double  *,int n))
{
int i,j;

	for (i=0;i<n;i++)
		dist[i][i] = 0;

	if (liste!=NULL)
	{
		for (i=0;i<n;i++)
			for (j=i+1;j<n;j++)
			{
				dist[i][j] = f_dist(data[liste[i]], data[liste[j]],p);
				dist[j][i] = dist[i][j];
			}
	}
	else
	{
		for (i=0;i<n;i++)
			for (j=i+1;j<n;j++)
			{
				dist[i][j] = f_dist(data[i], data[j],p);
				dist[j][i] = dist[i][j];
			}
	}
}
double CalculCostEnd(double **data,int n,int p,double **prototype,int n_proto,double (*f_dist)(double *,double  *,int n))
{
int i,j, winner;
double *Max_dist, Min_local;
int *weight;
double d_calc, cost;

	Max_dist = new double[n_proto];
	weight = new int[n_proto];
	for (i=0;i<n_proto;i++)
	{
		Max_dist[i] = 0;
		weight[i] = 0;
	}

	winner = 0;
	for (i=0;i<n;i++)
	{
		Min_local = MAX_VALUE;
		for (j=0;j<n_proto;j++)
		{
			d_calc = f_dist(data[i],prototype[j],p);
			if (d_calc < Min_local)
			{
				Min_local = d_calc;
				winner = j;
			}
		}
		weight[winner] += 1;
		if (Max_dist[winner] < Min_local)
			Max_dist[winner] = Min_local;
	}

	cost = 0;
	for (j=0;j<n_proto;j++)
		cost += (double)(weight[j]) * Max_dist[j];

	cost = cost / (double)(n);
	
	delete [] Max_dist;
	delete [] weight;

	return cost;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////:
int EcritFile(char *fichier,int Open,double *vecteur,int n)
{
int i;
FILE *fich;

	if (Open == 1)
		fich = fopen(fichier,"w");
	else
		fich = fopen(fichier,"a");

	if (!fich) 
		return 0;

	for (i=0;i<n;i++)
		fprintf(fich,"%lf\n",vecteur[i]);
	
	fprintf(fich,"\n");

	fclose(fich);
	return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////:
// Y = b*x¨a
int RegressionExponentielle(double *Pattern,int n_pattern,double &b,double &a)
{
int i;
double *Y;
double *X;
double Av_x, Av_y;
double d_x, d_y;
double d_xx, d_xy;

	Y = new double[n_pattern];
	X = new double[n_pattern];

	for (i=0;i<n_pattern;i++)
		if (Pattern[i] != 0)
			Y[i] = log(Pattern[i]);
		else
			Y[i] = 0;

	for (i=0;i<n_pattern;i++)
		X[i] = log((double)(i+1));

	Av_x = Av_y = 0;
	for (i=0;i<n_pattern;i++)
	{
		Av_x = Av_x + X[i];
		Av_y = Av_y + Y[i];
	}

	Av_x /= (double)(n_pattern);
	Av_y /= (double)(n_pattern);
	
	d_xx = d_xy = 0;
	for (i=0;i<n_pattern;i++)
	{
		d_x = X[i] - Av_x;
		d_y = Y[i] - Av_y;
		d_xx = d_xx + d_x * d_x;
		d_xy = d_xy + d_x * d_y;
	}

	a = d_xy/d_xx;
	b = Av_y - a*Av_x;
	b = exp(b);
	
	delete [] Y;
	delete [] X;
return 1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int UpdateWithGravity_DENDIS(double **data,struct Internal *I,ParamFCNNrule *P, int N_instance)
{
int i,j,winner;
int *I_mod;
int **liste;
int *masse, pattern_actif;
double dist_min, dist_calc;

	I_mod = new int[N_instance];
	liste = new int*[N_instance];
	masse = new int[N_instance];

	for (i=0;i<N_instance;i++)//pour tous les éléments de T.
	{
		liste[i] = NULL;
		masse[i] = 0;
		if (I->Pre_distance[i] > P->S2_Beta * P->dav_1knn[N_instance - 2])
			I_mod[i] = 1;
		else
			I_mod[i] = 0;

		for (j=0;j<P->p;j++)
			P->pattern[i][j] = 0;
	}

	if (P->S2_Option_all == 1)
	{
		for (i=0;i<P->n;i++)//pour tous les éléments de T.
		{
			winner = I->Nearest_proto[i];//l'instance la plus proche est mémorisé...on garde en mémoire l'ancien.
			for (j=0;j<P->p;j++)
				P->pattern[winner][j] = P->pattern[winner][j] + data[i][j]; 
		}

		for (i=0;i<N_instance;i++)
			for (j=0;j<P->p;j++)
				P->pattern[i][j] /= (double)(I->Weight[i]);
	}
	else // on ne regarde que ceux qui sont exentrés.
	{
		for (i=0;i<P->n;i++)//pour tous les éléments de T.
		{
			winner = I->Nearest_proto[i];//l'instance la plus proche est mémorisé...on garde en mémoire l'ancien.
			if (I_mod[winner] == 1)
			{
				if (liste[winner] == NULL)
					liste[winner] = new int[I->Weight[winner]];
				
				liste[winner][masse[winner]] = i;
				masse[winner] += 1;
				
				for (j=0;j<P->p;j++)
					P->pattern[winner][j] = P->pattern[winner][j] + data[i][j]; 				
			}
		}
		
		//on calcule le barycentre pour les instances concernés.
		for (i=0;i<N_instance;i++)
			if (I_mod[i] == 1)
			{
				//le nouveau centre.
				for (j=0;j<P->p;j++)
					P->pattern[i][j] /= (double)(I->Weight[i]);
//				printf("AVANT%5.2lf %5.2lf\n",data[I->Liste_selectionne[i]][0],data[I->Liste_selectionne[i]][1]); 
//				printf("BARYCENTRE%5.2lf %5.2lf\n",P->pattern[i][0],P->pattern[i][1]); 
				dist_min = 1000;
				for (j=0;j<masse[i];j++)//pour les patterns de ce prototype.
				{
					pattern_actif = liste[i][j];//le pattern.
					dist_calc = P->fdist(data[pattern_actif],P->pattern[i],P->p);
					if (dist_calc < dist_min)
					{
						dist_min = dist_calc;
						I->Liste_selectionne[i] = pattern_actif;
					}
				}				
//				printf("APRES%5.2lf %5.2lf\n",data[I->Liste_selectionne[i]][0],data[I->Liste_selectionne[i]][1]); 
			}		
	}

	//on calcule la moyenne des voisinages.

	
	//On élimine le BRUIT.
	for (i=0;i<N_instance;i++)
	{
		if (I->Weight[i] <= 0.1*P->Granularity*(double)P->n)
			I->Weight[i] = 0;

		if (I->Weight[i] <= (0.5*P->Granularity*(double)P->n) && I->Pre_distance[i] >= P->S2_Beta*P->dav_1knn[N_instance - 2])
			I->Weight[i] = 0;
	}

	Delete2D(liste,N_instance);
	delete [] masse;
	delete [] I_mod;

return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SortPrototypes_DENDIS(struct Internal *I,ParamFCNNrule *P,int &N_instance,int **Liste_pattern,int Min)
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
		if (I->Weight[i]>0)
		{
			Weight_tempon[compteur] = I->Weight[i];
			Liste_tempon[compteur] = new int[I->Weight[i]];
			for (j=0;j<I->Weight[i];j++)
				Liste_tempon[compteur][j] = Liste_pattern[i][j];//MEMOIRE!!
			SelectionIndividu[compteur] = I->Liste_selectionne[i];							
			compteur++;
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
	delete [] Weight_tempon;
	delete [] SelectionIndividu;

	for (i=0;i<N_instance;i++)
	{
		if (Liste_tempon[i])
			delete []Liste_tempon[i];
	}
	delete []Liste_tempon;
	
	N_instance = compteur;
}

void CalculusStatVectorBool(double *pattern,int N_pattern,int *booleen,double &Moy,
				  double &Ecart,double &Min,double &Max)
{
int i;
int Occurence = 0;

	Ecart = 0; 
	Moy = 0;
	Min = 10000;
	Max = -10000;
	

	for (i=0;i<N_pattern;i++)
		if (booleen[i] == 1)
		{
			Moy = Moy + pattern[i];
			
			if (pattern[i] < Min)
				Min = pattern[i];
			if (pattern[i] > Max)
				Max = pattern[i];

			Occurence++;
		}	
	
	if (Occurence != 0)
		Moy = Moy / (double)(Occurence);
	else
		Moy = 0;

	for (i=0;i<N_pattern;i++)
		if (booleen[i] == 1)
			Ecart = Ecart + (Moy - pattern[i]) * (Moy - pattern[i]);
	
	if (Occurence != 0)
		Ecart = sqrt(Ecart / (double)(Occurence));
	else 
		Ecart = 0;
	
}

void CalculusStatVector(double *pattern,int N_pattern,double &Moy,double &Ecart,double &Min,double &Max)
{
int i;
int Occurence = 0;

	Ecart = 0; 
	Moy = 0;
	Min = 10000;
	Max = -10000;
	

	for (i=0;i<N_pattern;i++)
		{
			Moy = Moy + pattern[i];
			
			if (pattern[i] < Min)
				Min = pattern[i];
			if (pattern[i] > Max)
				Max = pattern[i];

			Occurence++;
		}	
	
	if (Occurence != 0)
		Moy = Moy / (double)(Occurence);
	else
		Moy = 0;

	for (i=0;i<N_pattern;i++)
			Ecart = Ecart + (Moy - pattern[i]) * (Moy - pattern[i]);
	
	if (Occurence != 0)
		Ecart = sqrt(Ecart / (double)(Occurence));
	else 
		Ecart = 0;
	
}

double InertieGenAleaListe(double **data,int n,int p,int n_output)
{
int i,j;
int random;
double **data_test;
double moyenne[100];
double Inertie;

	srand(time(NULL)); // initialisation de rand

	data_test = Reservmemdouble(n_output,p);
	if (n < n_output) 
		return 0;

	for (i=0;i<p;i++)
		moyenne[i] = 0;

	for (i=0;i<n_output;i++)
	{
		random = rand()%(n-1);
		for (j=0;j<p;j++)
		{
			data_test[i][j] = data[random][j]; 
			moyenne[j] += data_test[i][j];
		}
	}

	for (j=0;j<p;j++) 
		moyenne[j] /= (double)(n_output);

	Inertie = 0;
	for (i=0;i<n_output;i++)
		Inertie += d_Eucl(moyenne,data_test[i],p); 
	Inertie = Inertie / (double)(n_output);

	Delete2D(data_test,n_output);
	
	return Inertie;
}

void NormaliseEcartType(double **data,int n,int p)
{
int i,j;
double *vector;
double mean, ecart;

	vector = new double[n];	
	for (j=0;j<p;j++)
	{
		for (i=0;i<n;i++)
			vector[i] = data[i][j];

		CalculusStatVector(vector,n,mean,ecart);			
		for (i=0;i<n;i++)
			data[i][j] = (data[i][j] - mean) / (ecart);
	}

delete []vector;
}

void CalculusStatVector(double *pattern,int N_pattern,double &Moy,double &Ecart)
{
int i;
	
	Ecart = 0; 
	Moy = 0;
	
	for (i=0;i<N_pattern;i++)
		Moy = Moy + pattern[i];
	
	
	Moy = Moy / (double)(N_pattern);
	
	for (i=0;i<N_pattern;i++)
		Ecart = Ecart + (Moy - pattern[i]) * (Moy - pattern[i]);

	Ecart = sqrt(Ecart / (double)(N_pattern));
	
}


void NormaliseMinMax(double **data,int n,int p)
{
int i,j;
double Vmin, Vmax;

	
	for (j=0;j<p;j++)
	{
		Vmin=100000;
		Vmax=-10000;
		for (i=0;i<n;i++)
		{
			if (data[i][j] < Vmin)
				Vmin = data[i][j];
			if (data[i][j] > Vmax)
				Vmax = data[i][j];
		}
		for (i=0;i<n;i++)
			data[i][j] = (data[i][j] - Vmin) / (Vmax - Vmin);
	}

}

int LitFichier2DS(char *fichier,int n,int p,double **data)
{
int i,j;
double val;
FILE *fich;

	fich = fopen(fichier,"r");
	if (!fich) return 0;

	for (i=0;i<n;i++)
	{
		for (j=0;j<p;j++)
		{
			fscanf(fich,"%lf ",&val);
			data[i][j] = val;
		}

//		fscanf(fich,"%d ",&cat[i]);
	}

	fclose(fich);
	return 1;
}


int ReturnColonne(char *buffer)
{
int length;
char *read;
char *var;
char cretour = '\0';
char espace = 32;
char tabulation = 9;
int fin_ligne = 0;
int deb, fin, i, compteur;
double variable;
int p;

	if (buffer ==NULL)
		return 0;

	read = new char[10000];
	var = new char[1000];
	p = 0;
	length = strlen(buffer);
	strcpy(read,buffer);
	read[length] = cretour;
	compteur = 0;

	while (read[compteur]==espace)
		compteur++;

	deb = compteur;
	while(fin_ligne==0)
	{
		while (read[compteur]!=espace && read[compteur]!=tabulation && read[compteur] != cretour)
			compteur++;
		
		fin = compteur;
		
		for (i=deb;i<fin;i++)
			var[i-deb] = read[i];
		var[fin-deb] = cretour; 

		sscanf(var,"%lf",&variable);
		p++;
		compteur = fin+1;

		while (read[compteur]==espace || read[compteur]==tabulation)
			compteur++;

		if (read[compteur-1] == cretour || compteur>=length-1)
		{
			delete [] read;
			delete [] var;
			return p;
		}
		deb = compteur;
	}

	delete [] read;
	delete [] var;

return p;
}

int ExtraitDimensionFile(char *fichier,int *ligne,int *colonne)
{
FILE *fich;
int n;
int p;
char *buffer=NULL;
char *read=NULL;
int oldp;

	fich = fopen(fichier,"r");
	if (!fich) 
		return -1;

	oldp=0;
	n = 0;
	p = 0;
	read = new char[3000];
	buffer = new char[10001];

	while(fgets(buffer,10000,fich) && strlen(buffer)>= 2)
	{
		p = ReturnColonne(buffer);
		if (n>1 && p!=oldp)
			break;
		n++;
		oldp = p;
	}
			
	*ligne = n;
	*colonne = p;

	delete []read;
	delete []buffer;

	fclose(fich);
	
	return n;
}

int EcritFichier2D(char *fichier,int n,int p,double **data)
{
int i,j;
FILE *fich;

	fich = fopen(fichier,"w");
	if (!fich) return 0;

	for (i=0;i<n;i++)
	{
		for (j=0;j<p;j++)
			fprintf(fich,"%6.4lf\t",data[i][j]);

		fprintf(fich,"\n");
	}

	fclose(fich);
	return 1;
}

/*
soit P1 un prototype existant et P2 le nouveau prototype.
On connait d(P1,P2) et on connaît aussi dmax(P1)
On veut savoir si tous les individus de P1 ne peuvent pas être voisins de P.
C'est vrai si d(P1M) < 1/2 * d(P1,P2)

On veut savoir pour tous les protoypes. disable[i] = 1 si c'est le cas.
*/
void DisablePrototype(double *d_max,double *dist_proto,int N_proto,int *disable)
{
int i;

	for (i=0;i<N_proto - 1;i++)
	{
		if (d_max[i] < 0.5 * dist_proto[i])
			disable[i] = 1;
		else
			disable[i] = 0;
	}
}

//La distance entre les différents clusters,
int InegalTriangular(double **D_cluster,double D_nearest,int Cluster_actuel, int Cluster_candidat)
{
	if (D_cluster[Cluster_actuel][Cluster_candidat] >= 2*D_nearest) 
		return 0;
	else
		return 1;
}

void InitTri(struct Tri *Tr,int n)
{
	Tr->G_actif = new int[n];
	Tr->occ = new int[n];
	Tr->weight = new int[n];
	Tr->noise = new int[n];
	Tr->ptisole = new int[n];
	
	for (int i=0;i<n;i++)
	{
		Tr->G_actif[i] = 0;
		Tr->occ[i] = 0;
		Tr->weight[i] = 0;
		Tr->noise[i] = 0;
		Tr->ptisole[i] = 0;
	
	}
}


void LibereTri(struct Tri *Tr)
{
	delete []Tr->G_actif;
	delete []Tr->occ;
	delete []Tr->weight;
	delete []Tr->noise;
	delete []Tr->ptisole;
}

double fmax(double a,double b)
{
	if (a > b)
		return a;
	else
		return b;
}

double fmin(double a,double b)
{
	if (a < b)
		return a;
	else
		return b;
}

int max(int a,int b)
{
	if (a > b)
		return a;
	else
		return b;
}

int min(int a,int b)
{
	if (a < b)
		return a;
	else
		return b;
}


void Delete1D(int *x)
{
	if (x == NULL)
		return;

	G_1Dx = x;
	
	delete []G_1Dx;

	G_1Dx = NULL;
}

void Delete1D(double *x)
{
	if (x == NULL)
		return;

	G_1Ddx = x;
	
	delete []G_1Dx;

	G_1Dx = NULL;

}

/*inline*/ void Delete2D(double **x,int n)
{
int i;
	
	if (x==NULL) 
		return;

	G_2Dy = x;

	for (i=0;i<n;i++)
		if (G_2Dy[i])
			delete [] G_2Dy[i];

	delete [] G_2Dy;

	G_2Dy = NULL;

}

/*inline*/ double **Reservmemdouble(int n,int p)
{
int i;
double **data;

	data = new double*[n];
	if (!data) return NULL;

	for (i=0;i<n;i++)
	{
		data[i] = new double [p];
		
		if (data[i] == NULL)
		{
			for(int j=0;j<i-1;j++)
				delete []data[j];
			delete []data;
			return NULL;
		}
	}

return data;
}


/*inline*/ int **Reservmemint(int n,int p)
{
int i;
int **data;

	data = new int*[n];
	if (!data) return NULL;

	for (i=0;i<n;i++)
		data[i] = new int [p];

return data;
}

/*inline*/ void Delete2D(int **x,int n)
{
int i;
	
	if (x==NULL) 
		return;

	Gint_2dy = x;

	for (i=0;i<n;i++)	
	{
		if (Gint_2dy[i]!=NULL) 
		{
			delete []Gint_2dy[i];
			Gint_2dy[i] = NULL;
		}
	}
	delete []Gint_2dy;

	Gint_2dy = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int *Reservememint(int n)
{

int *categ;
	
	categ = new int[n];
	if (!categ) 
		return NULL;

return categ;
}

int GenAleaListe(int n,int *liste,int n_output)
{
int i;
int n_actif;
int random;

	srand(time(NULL)); // initialisation de rand
	n_actif = 0;

	if (n < n_output) 
		return 0;
	
	for (i=0;i<n_output;i++)
		liste[n_actif++] = random = rand()%(n);

	return 1;
}

void GenAlea(double **data,int n,int p,double **output,int n_output)
{
int i,j;
int n_actif;
int random;
double **tampon;
int ninput;

	n_actif = 0;
	if (n < n_output) 
		return;

	tampon = Reservmemdouble(n,p);

	for (i=0;i<n;i++)
		for (j=0;j<p;j++)
			tampon[i][j] = data[i][j];

	ninput = n;

	while (n_actif < n_output)
	{
		random = rand()% (ninput);

		for (j=0;j<p;j++)
			output[n_actif][j] = tampon[random][j];

		for (i=random;i<n_actif-1;i++)
			for (j=0;j<p;j++)
				tampon[i][j] = tampon[i+1][j];

		n_actif++;
		ninput--;
	}

	Delete2D(tampon,n);
}

void SortWeight(int Ngeneration,int **weight,int Min,int *ReelWeight,double **Radius,double *Indice)
{
int i,j;
double s;

	for (i=0;i<Ngeneration;i++)
	{
		ReelWeight[i] = 0;
		Indice[i] = 0;
		s = 0;
		for (j=0;j<Ngeneration;j++)
			if (weight[i][j] >= Min)
			{
				ReelWeight[i] += 1;
				Indice[i] = Indice[i] + weight[i][j] * Radius[i][j];
				s = s + weight[i][j];
			}
		Indice[i] = Indice[i] / s;
	}
}

int tri_decroissant_int(const void *a,const void *b)
{ 
	if(*(int*)(a)<=*(int*)(b)) 
		return(1); 
	else 
		return(-1);
}

int tri_decroissant_double(const void *a,const void *b)
{ 
	if(*(double*)(a)<=*(double*)(b)) 
		return(1); 
	else 
		return(-1);
}

int tri_decroissant(const void *a,const void *b)
{ 
	if(*(double*)(a)<=*(double*)(b)) 
		return(1); 
	else 
		return(-1);
}

/////////////////////////////////////////////////////////////////
int tri_croissant(const void *a,const void *b)
{ 
	if(*(double*)(a)<=*(double*)(b)) 
		return(-1); 
	else 
		return(1);
}

/////////////////////////////////////////////////////////////////
int tri_croissant_int(const void *a,const void *b)
{ 
	if(*(int*)(a)<=*(int*)(b)) 
		return(-1); 
	else 
		return(1);
}

/////////////////////////////////////////////////////////////////
double Dist_eucl(double *p1,double *p2,int taille)
{
int i;
double distance = 0;	

	for (i=0;i<taille;i++)
		distance = distance + ((p1[i] - p2[i])*(p1[i] - p2[i]));

	if (distance == 0)
		return 0.;
	else
		return sqrt(distance);
}

///////////////////////////////////////////////////////////
double d_Eucl(double *i1, double *i2, int n)
{
int i;
double d = 0;

	for (i=0;i<n;i++)
		d = d + (i1[i] - i2[i])*(i1[i] - i2[i]);
	
	if (d!= 0.0)
	{
		d = sqrt(d);
		return d;
	}
	else
		return 0;
}

///////////////////////////////////////////////////////////
