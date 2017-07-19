input file: example data (only this format).

USING a fonction.

int example(void)
{
char folder[500];
char input[500];
int Methode=1;
double Epsilon=0.001;
int option[3];
size_t memsize;

	option[0] = 1;
	option[1] = 0;
	option[2] = 0;

	sprintf(folder,"c:\\mydata\\essai\\");
	sprintf(input,"0synthetic.txt");

	
	while(!kbhit())
	{
		printf("GO %d\t", rand()%100);
//		printf("GO %d(%lu)\t",i,(unsigned long)getMemorySize( ));
		ExecuteSampling(folder,input,Methode,Epsilon,option);
//		printf("(%lu)\n",(unsigned long)getMemorySize( ));
	}
}

USING the Cmd.

Scenario 1:
samplingapplication 0synthetic.txt 0.02 1 => Execute in the same folder without option (exe + input file + granularity + method)

Scenario 2:
samplingapplication 0synthetic.txt 0.02 1 -s 1=> Execute in the same folder with the save option (3 files)

Scenario 3:
samplingapplication 0synthetic.txt 0.02 1 -d c:\mydata\essai\\ =>  (0synthetic.txt inside the folder "c:\mydata\essai\\") 


