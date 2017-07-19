// ProjetPro.cpp : définit le point d'entrée pour l'application console.
//



#define PROTRAS 0
#define DENDIS 1
#define DIDES 2


#include "stdio.h"
#include <string.h>
#include <stdlib.h>
#include "time.h"
#include "ADN_Esclave1.h"
#include "ADN_Esclave2.h"
#include "ADN_Esclave3.h"
#include "ADN_util.h"

/////////////////////////////
int *G_1Dx,**Gint_2dy;
double *G_1Ddx,**G_2Dy;
/////////////////////////////

int ExecuteSampling(char *folder,char *input,int Methode,double Epsilon,int *option);

int main(int argc, char ** argv)
//******************************
{
  if(argc < 4)
    {
      printf("\nThe program %s needs 3 arguments : ", argv[0]);
      printf("\n\t* the file name;");   
      printf("\n\t* the granularity.") ;
      printf("\n\t* the sampling method DENDIS (1) or DIDES (2)") ;
      printf("\nOption : ");
      printf("\n\t* -d directory for input and output data (default: current)");
      printf("\n\t* -s save the results in three files (0 or 1), default: 0 (unsaved)");
      printf("\n\t* -i random (1) or deterministic (0) init, default:  1");
      printf("\n\t* -p DIDES post-processing (0 or 1), default:1  (with)");
      printf("\n");
      return 1;
    }

char *folder = NULL;
int Methode;
double Epsilon;
int option[3];

Epsilon = strtod(argv[2], NULL);
Methode = atoi(argv[3]);

 if(Methode > 2 || Methode < 0)
   {
     printf("\nInvalid method argument (%d), set to DENDIS\n", Methode);
     Methode = 1;
   }
	
option[0] = 1;//DIDES post-processing
option[1] = 1;//random init
option[2] = 0;//Save results 3 files.

folder = NULL; // Empty string for current directory

    if(argc > 4) 
      {
	for(int i = 4; i < argc; i++)
	  {
	    if(! strncmp("-p", argv[i], 2) ) option[0] = atoi(argv[i + 1]);
	    else if(! strncmp("-i", argv[i], 2) ) option[1] = atoi(argv[i + 1]);
	    else if(! strncmp("-s", argv[i], 2) ) option[2] = atoi(argv[i + 1]);
	    else if(! strncmp("-d", argv[i], 2) ) folder = argv[i + 1];
	  }
      }

    printf("Data:  %s, gr:  %f, method; %d \n", argv[1], Epsilon, Methode);
    printf("Options, 0: %d, 1: %d, 2: %d, f: %s\n", option[0], option[1],option[2],folder);

	if (ExecuteSampling(folder,argv[1],Methode,Epsilon,option)==1)
		printf("done");
	else
		printf("error");
	return 0;
}

int ExecuteSampling(char *folder,char *input,int Methode,double Epsilon,int *option)
{
double **data;
int n, p;
double **prototype;
int N_instance, *weight;
int npoint, pmax, pused, MaxCluster = 1000;
char fsort[500],cheminoutput[200];
int *Liste_instance;
double cost;
void *pointeur=NULL;

	MaxCluster = 2000;
	npoint = 100000;
	pmax = 100;
	pused = 10;

	
	data = Reservmemdouble(npoint,pmax);
	if (!data) return 0;

	prototype = Reservmemdouble(MaxCluster,pmax);
	if (!prototype){Delete2D(data,npoint); return 0;}

	weight = new int[MaxCluster];
	if (!weight) { Delete2D(data,npoint); Delete2D(prototype,MaxCluster); return 0;}

	Liste_instance = new int[MaxCluster];
	if (!Liste_instance) { Delete2D(data,npoint); Delete2D(prototype,MaxCluster); delete []weight; return 0;}

	if (folder)
		sprintf(fsort,"%s%s",folder,input);
	else
		sprintf(fsort,"%s",input);

	if (option[2] == 1)
	{
		if (folder)
		{
			sprintf(cheminoutput,"%s",folder);
			pointeur = cheminoutput;
		}
		else
			pointeur = NULL;
	}

	if (ExtraitDimensionFile(fsort,&n,&p)!=-1)//Read the file and extract dimensions.
	{
		if (n==0 || p==0) 
			return 0; 			
		n = min(n,npoint);//possibilité la taille des patterns
		LitFichier2DS(fsort,n,p,data);//lecture fichier
		p = min(p,pused);//possibilité de limiter la taille de l'espace		
		NormaliseEcartType(data,n,p);//normalisation des data
		if (folder)
			sprintf(fsort,"%snorm%s",folder,input);
		else
			sprintf(fsort,"norm%s",input);
		EcritFichier2D(fsort,n,p,data);
	}
	else
		return 0;
	
	N_instance = MaxCluster - 1;//En paramètre un nombre d'instance MAX
	switch(Methode)
	{
		#ifdef DENDIS
		case DENDIS:
		SamplingApproachP_DENDIS(data,n,p,Epsilon,prototype,N_instance,weight,(char *)pointeur,option);
		cost = CalculCostEnd(data,n,p,prototype,N_instance,d_Eucl);
		printf("Sampling done DENDIS(%lf):(%d,%d);%d,%lf\n",Epsilon,n,p,N_instance,cost);
		if (folder)
			sprintf(fsort,"%soutDENDIS.txt",folder);
		else
			sprintf(fsort,"outDENDIS.txt");
		EcritFichier2D(fsort,N_instance,p,prototype);
		break;
		#endif
	
		#ifdef DIDES
		case DIDES:
		SamplingApproachP_DIDES(data,n,p,Epsilon,prototype,N_instance,weight,(char *)pointeur,option);
		cost = CalculCostEnd(data,n,p,prototype,N_instance,d_Eucl);
		printf("Sampling done DIDES(%lf):(%d,%d);%d,%lf\n",Epsilon,n,p,N_instance,cost);
		if (folder)
			sprintf(fsort,"%soutDIDES.txt",folder);
		else
			sprintf(fsort,"outDIDES.txt");
		EcritFichier2D(fsort,N_instance,p,prototype);
		break;
		#endif

		default: break;
	}
	

	Delete2D(data,npoint);
	Delete2D(prototype,MaxCluster);
	delete [] Liste_instance;
	delete [] weight;
	return 1;
}
