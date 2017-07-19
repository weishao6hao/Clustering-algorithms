

#include "math.h"
#include "ADN_Esclave1.h"
#include "ADN_Esclave2.h"
#include "ADN_Esclave3.h"
#include "ADN_util.h"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int SelectinstanceWithListe(double **pattern,int p,int *liste_det,int n_inst,int Min,int *weightG,double **data)
{
int i,j;
int Nb=0;

	if (weightG != NULL)
	{
		for (i=0;i<n_inst;i++)
		{
			if (weightG[i] >= Min)
			{
				for (j=0;j<p;j++)
					data[Nb][j] =  pattern[liste_det[i]][j]; 
				weightG[Nb] = weightG[i];
				Nb++;
			}
		}
	}
	else
	{
		for (i=0;i<n_inst;i++)
		{
			for (j=0;j<p;j++)
				data[Nb][j] =  pattern[liste_det[i]][j]; 
			Nb++;
		}
	}

	return Nb;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////:
/*
input: les data

*/
int GenerePartitionFromCluster_VlisteAvecFile(double **data,double **pattern,int n,int p,
											  int **liste_in,int *weight_in,int N_cluster,int N_instance,
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
	
	identite = new int[n];
	cluster = Reservmemdouble(N_instance,p);
	winner = 0;
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
	}

	dist = new double[N_instance]; 
	
	for (i=0;i<N_instance;i++)//on regarde les S actifs.	
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

//		liste[winner][weight[winner]] = i;
		weight[winner] = weight[winner] + 1; 
	}
	
	for (i=0;i<N_instance;i++)
	{
		if (liste[i]!=NULL)
			liste[i] = NULL;
		liste[i] = new int[weight[i]];
	}

	for (i=0;i<N_instance;i++)//on regarde les S actifs.	
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



////////////////////////////////////////////////////////////////////////////////////
int SearchPartitionIndex(int pattern,int **liste, int n_partition, int *weight)
{
int i,j;

	for (i=0;i<n_partition;i++)
		for (j=0;j<weight[i];j++)
		{
			if (liste[i][j] == pattern)
				return i;
		}

return -1;
}
/*
////////////////////////////////////////////////////////////////////////////////////
double RandIndex(int n_pattern,int **liste_1,int **liste_2,int n_partition_1,int n_partition_2,int *weight_1,int *weight_2)
{
int i,j;
int *NumPart_1;
int *NumPart_2;
int x_1, x_2, y_1, y_2;
double a,b,c,d;
double Index;

	NumPart_1 = new int[n_pattern];
	NumPart_2 = new int[n_pattern];

	for (i=0;i<n_pattern;i++)
	{
		NumPart_1[i] = SearchPartitionIndex(i,liste_1,n_partition_1,weight_1);
		NumPart_2[i] = SearchPartitionIndex(i,liste_2,n_partition_2,weight_2);
		if (NumPart_1[i] == -1) 
		{
			delete NumPart_1;
			delete NumPart_2;
			return -1;
		}

		if (NumPart_2[i] == -1) 
		{
			delete NumPart_1;
			delete NumPart_2;
			return -1;
		}
	}

	a = b = c = d = 0;
	for (i=0;i<n_pattern;i++)
	{
		x_1 = NumPart_1[i]; 
		x_2 = NumPart_2[i];
		for (j=i+1;j<n_pattern;j++)
		{
			if (i!=j)
			{	
				y_1 = NumPart_1[j];
				y_2 = NumPart_2[j];
				if ( (x_1 == y_1) && (x_2 == y_2))
					a++;
		
				if ( (x_1 != y_1) && (x_2 == y_2))
					b++;
		
				if ( (x_1 == y_1) && (x_2 != y_2))
					c++;
		
				if ( (x_1 != y_1) && (x_2 != y_2))
					d++;
			}
		}

	
	}
	
	delete NumPart_1;
	delete NumPart_2;
	
	long double Denominateur = ( (long double)(n_pattern) * (long double)(n_pattern - 1)/2.0);
	Index = (double)( a/long double (Denominateur) + d/long double (Denominateur));

	return Index;
}


double  MutualInformation(int n,int **liste_1,int **liste_2,int n_partition_1,int n_partition_2,int *weight_1,int *weight_2)
{
int i,j,k,l;
double R_index, X_i,X_j;
int **Matrix;
double Adj_R_index, Exp_index,Ratio;
double *Rl, *Rc;
double H_i, H_j;

	Rl = new double[n_partition_1];
	Rc = new double[n_partition_2];

	Matrix = Reservmemint(n_partition_1,n_partition_2);

	for (i=0;i<n_partition_1;i++)
		for (j=i;j<n_partition_2;j++)
			Matrix[i][j] = 0;

	for (i=0;i<n_partition_1;i++)
		for (j=i;j<n_partition_2;j++)
		{
			for (k=0;k<weight_1[i];k++)
				for (l=0;l<weight_2[j];l++)
					if (liste_1[i][k] == liste_2[j][l])
						Matrix[i][j] += 1;
			Matrix[j][i] = Matrix[i][j];
		}

	H_i = 0;
	for (i=0;i<n_partition_1;i++)
	{
		X_i = 0;
		for (j=0;j<n_partition_2;j++)
			X_i += Matrix[i][j];
		Rl[i] = X_i/(double)(n);
		if (Rl[i]!=0)
			H_i += Rl[i]*log(Rl[i]);
	}

	H_j=0;
	for (j=0;j<n_partition_2;j++)
	{
		X_j = 0;
		for (i=0;i<n_partition_1;i++)
			X_j += Matrix[i][j];
		Rc[j] = X_j/(double)(n);
		if (Rc[j]!=0)
			H_j += Rc[j]*log(Rc[j]);
	}

	R_index = 0;
	double fij;
	for (i=0;i<n_partition_1;i++)
		for (j=0;j<n_partition_2;j++)
		{
			fij = (double)Matrix[i][j]/(double)(n);
			if (Rl[i]!=0 && Rc[j]!=0 && fij!=0)
				R_index += fij * log(fij/(Rl[i]*Rc[j]));
		}
	Delete2D(Matrix,n_partition_1);
	
	delete Rl;
	delete Rc;

	return R_index/(-0.5*(H_i+H_j));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Cnp(int p,int n)
{
int i;
double Factor = 100000;
double s1, s2;

	if (p>n) return 0;
	
	s1 = (double)(n)/Factor;
	for (i=n-1;i>(n-p);i--)
		s1 = s1 * ((double)(i)/Factor);  
	
	s2 = (double)(p)/Factor;
	for (i=p-1;i>1;i--)
		s2 = s2 * ((double)(i)/Factor);
	
return (s1/s2);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double AdjRandIndex( int n,int **liste_1,int **liste_2,int n_partition_1,int n_partition_2,int *weight_1,int *weight_2)
{
int i,j,k,l;
double R_index, X_i,X_j,R_i, R_j;
int **Matrix;
double Adj_R_index, Exp_index,Ratio;

	Matrix = Reservmemint(n_partition_1,n_partition_2);
	for (i=0;i<n_partition_1;i++)
		for (j=0;j<n_partition_2;j++)
			Matrix[i][j] = 0;

	for (i=0;i<n_partition_1;i++)
		for (j=0;j<n_partition_2;j++)
		{
			for (k=0;k<weight_1[i];k++)
				for (l=0;l<weight_2[j];l++)
				{
					if (liste_1[i][k] == liste_2[j][l])
						Matrix[i][j] += 1;
				}

//			if (j<n_partition_1 && i<n_partition_2)
//				Matrix[j][i] = Matrix[i][j];
		}

	

	R_index = 0;
	for (i=0;i<n_partition_1;i++)
		for (j=0;j<n_partition_2;j++)
			R_index += Cnp(2,Matrix[i][j]);

	R_i = 0;
	for (i=0;i<n_partition_1;i++)
	{
		X_i = 0;
		for (j=0;j<n_partition_2;j++)
			X_i += Matrix[i][j];
		R_i += Cnp(2,X_i);
	}

	R_j = 0;
	for (j=0;j<n_partition_2;j++)
	{
		X_j = 0;
		for (i=0;i<n_partition_1;i++)
			X_j += Matrix[i][j];
		R_j += Cnp(2,X_j);
	}

	Ratio = Cnp(2,n);
	Exp_index = (R_i*R_j)/Ratio;

	Adj_R_index = (R_index - Exp_index) / (0.5*(R_i+R_j) - Exp_index); 
	
	Delete2D(Matrix,n_partition_1);

	return Adj_R_index;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double RandIndexFuzzy(int n_pattern,int **liste_1,int **liste_2,int n_partition_1,int n_partition_2,int *weight_1,
					  int *weight_2,int *ambigu_1,int *ambigu_2)
{
int i,j;
int *NumPart_1;
int *NumPart_2;
int x_1, x_2, y_1, y_2;
double a,b,c,d;
double Index;
long double n_pattern_reel;

	NumPart_1 = new int[n_pattern];
	NumPart_2 = new int[n_pattern];

	for (i=0;i<n_pattern;i++)
	{
		NumPart_1[i] = SearchPartitionIndex(i,liste_1,n_partition_1,weight_1);
		NumPart_2[i] = SearchPartitionIndex(i,liste_2,n_partition_2,weight_2);
		if (NumPart_1[i] == -1) 
		{
			delete NumPart_1;
			delete NumPart_2;
			return -1;
		}

		if (NumPart_2[i] == -1) 
		{
			delete NumPart_1;
			delete NumPart_2;
			return -1;
		}
	}

	a = b = c = d = 0;
	n_pattern_reel = 0;
	for (i=0;i<n_pattern;i++)
		if (ambigu_1[i] == 0 && ambigu_2[i] == 0)
		{
		x_1 = NumPart_1[i]; 
		x_2 = NumPart_2[i];
		for (j=i+1;j<n_pattern;j++)
		{
			if (i!=j && ambigu_1[j] == 0 && ambigu_2[j] == 0)
			{
				n_pattern_reel++;
				y_1 = NumPart_1[j];
				y_2 = NumPart_2[j];
				if ( (x_1 == y_1) && (x_2 == y_2))
					a++;
		
				if ( (x_1 != y_1) && (x_2 == y_2))
					b++;
		
				if ( (x_1 == y_1) && (x_2 != y_2))
					c++;
		
				if ( (x_1 != y_1) && (x_2 != y_2))
					d++;
			}
		}

	
	}
	
	delete NumPart_1;
	delete NumPart_2;
	
//	long double Denominateur = ( (long double)(n_pattern_reel) * (long double)(n_pattern_reel - 1)/2.0);
	long double Denominateur = n_pattern_reel;
	if (Denominateur == 0) 
	{
		printf("probleme dans index fuzzy\n");
		return 0;
	}
	
	Index = (double)( a/long double (Denominateur) + d/long double (Denominateur));

	return Index;
}
*/
