//	Galois Field 2^8

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <malloc.h>

#define		W	8
#define		N	10
#define		M	6

long int fsize, fsize2;
int Device[N+M];
int AA[N][N];		//	Matrix to recalculate data at recovery time
unsigned char buf[1024]={'\0'};

FILE *InFile, *OutFiles[N+M], *RecoFile;

void OpenInAndOutFiles(void)
{
	char inputfile[80], outputfile[80];

	printf("\n\nInput Filename :: "); gets(inputfile);
	printf("\n\nOutput Filename :: "); gets(outputfile);

	if((InFile = fopen(inputfile, "rb")) == NULL)
	{
		printf("\n\nFile can't open."); exit(0);
	}

	fseek (InFile, 0, SEEK_END);   // non-portable
    fsize=ftell (InFile);

    printf("\n\nSize of file = %ld", fsize);

    fseek ( InFile , 0 , SEEK_SET );

    //rewind (InFile);

	char filename[][20]={"Dff1.txt", "Dff2.txt", "Dff3.txt",
						"Dff4.txt", "Dff5.txt", "Dff6.txt",
						"Dff7.txt", "Dff8.txt", "Dff9.txt", "Dff10.txt",
						"Dff11.txt", "Dff12.txt", "Dff13.txt",
						"Dff14.txt", "Dff15.txt", "Dff16.txt",
						"Dff17.txt", "Dff18.txt", "Dff19.txt", "Dff20.txt"};

	int ch=32;
	for(int i=0; i<N+M; i++)
	{
		if((OutFiles[i] = fopen(filename[i], "wb")) == NULL)
		{
			printf("\n\nFile can't create."); exit(0);
		}

		fprintf(OutFiles[i], "%ld", fsize);
		fputc(ch, OutFiles[i]);
	}

	if((RecoFile = fopen(outputfile,"wb")) == NULL)
	{
		printf("\n\nCan't create recovered file.");
		exit(0);
	}
}


unsigned int prim_poly_4 = 023;
unsigned int prim_poly_8 = 0435;
unsigned int prim_poly_16 = 0210013;
unsigned short *gflog, *gfilog;

void setup_tables(void)
{
	unsigned int b, log, x_to_w, prim_poly, b1;

	//FILE *logfile;

	switch(W) {
		case 4: prim_poly = prim_poly_4; break;
		case 8: prim_poly = prim_poly_8; break;
		case 16: prim_poly = prim_poly_16; break;
		//default: return -1;
	}

	x_to_w = 1 << W;

	gflog = (unsigned short *) malloc (sizeof(unsigned short) * x_to_w);
	gfilog = (unsigned short *) malloc (sizeof(unsigned short) * x_to_w);

	b = 1;

	for (log = 0; log < x_to_w-1; log++) {

		gflog[b] = (unsigned short) log;  	//printf("%u,",log);
		gfilog[log] = (unsigned short) b;	//printf("%u\n",b);

		b = b << 1;

		if (b & x_to_w) b = (b ^ prim_poly)&(x_to_w-1);

		
}

unsigned short AddSub(unsigned short a, unsigned short b)		//	if a & b are same then result will be 0
{
	return a^b;
}

#define NW (1 << W) /* In other words, NW equals 2 to the w-th power */

unsigned short gMult(unsigned short a, unsigned short b)
{
	int sum_log;

	if (a == 0 || b == 0) return 0;

	sum_log = gflog[a] + gflog[b];

	if (sum_log >= NW-1) sum_log -= NW-1;

	return gfilog[sum_log];
}

int gDiv(int a, int b)
{
	int diff_log;

	if (a == 0) return 0;
	if (b == 0) return -1; 			/* Can’t divide by 0 */

	diff_log = gflog[a] - gflog[b];

	if (diff_log < 0) diff_log += NW-1;

	return gfilog[diff_log];
}

int gPower(int i, int j)
{
	int result=i;

	if(j==0) return 1;

	for(int k = 2; k <= j; k++) result = gMult(result, i);

	return result;
}

int matrixF[M][N], matrixI[N][N]={{0}};

void setup_matrixFI(void)
{
	int i, j;

	for(i=1; i<=M; i++)
	{
		for(j=1; j<=N; j++) {
			matrixF[i-1][j-1] = gPower(j, i-1);
			//printf("%5d", matrixF[i-1][j-1]);
		}
		//printf("\n\n");
	}

	for(i=0; i<N; i++) matrixI[i][i] = 1;


}




void gInverse(int aa[][N], int bb[][N])		//	find inverse of a matrix over GF
{
	int i, j, k;
	int x;

	for(i=0; i<N; i++)
		for(j=0; j<N; j++) {
                if(i==j) bb[i][j]=1.0;
                else bb[i][j]=0.0;
		}

	for(k=0; k<N; k++){


	x=aa[k][k];		//	this element needs to be 1;

	if(x) {
		for(j=0; j<N; j++){
			aa[k][j] = gDiv(aa[k][j], x); 		//	aa[k][j] /= x;
			bb[k][j] = gDiv(bb[k][j], x);		//	bb[k][j] /= x;
		}
	}

	else {		//	Find a non zero element at aa[i][k] position
		i=k+1;
		while(!aa[i][k]) i++;

		if(i>=N) {
			printf("Inversion not possible, i=%d, k=%d\n\n", i, k);

			for(int i=0; i<N; i++) {
				for(int j=0; j<N; j++) printf("%10d", aa[i][j]);

				printf("\n\n");
			}

			exit(0);
		}

		//	Now interchange row i & row k

		int t;

		for(j=0;j<N;j++) {

			t=aa[i][j]; aa[i][j] = aa[k][j]; aa[k][j] = t;

			t=bb[i][j]; bb[i][j] = bb[k][j]; bb[k][j] = t;
		}

		x=aa[k][k];		//	this element needs to be 1;

		for(j=0; j<N; j++){
			aa[k][j] = gDiv(aa[k][j], x); 		//	aa[k][j] /= x;
			bb[k][j] = gDiv(bb[k][j], x);		//	bb[k][j] /= x;
		}
	}

	for(i=0; i<N; i++){

		if(i != k) {

			x = aa[i][k];

			for(j=0; j<N; j++){

				int tt;

				tt = gMult(aa[k][j], x);
				aa[i][j] = AddSub(aa[i][j],tt);		//	aa[i][j] -= aa[k][j]*x;

				tt = gMult(bb[k][j], x);
				bb[i][j] = AddSub(bb[i][j],tt);		//	bb[i][j] -= bb[k][j]*x;

			}
		}
	}

	}

}

void find_Checksum(void)
{
	int i, k, j;

	for(i=0; i<M; i++)
	{
		k=0;

		for(j=0; j<N; j++)
		{
			k = AddSub(k, gMult(Device[j], matrixF[i][j]));
		}
		

		Device[N+i] = k;
		//printf("%d", Device[N+i]);		//	Device[N+ indicates checksum bytes
	}
}


void FragmentFile(void)
{
	int ch;
	int i, tx, buf_loc;
	//long int file_ch=0;


	printf("\n\nWriting data to files...\n\n");


	i=0;

	while( (tx = fread(buf, 1, 1024, InFile)) > 0)
	{
		buf_loc=0;

		while(buf_loc < tx) {


		Device[i%N] = (int)buf[buf_loc];	i++; buf_loc++;	//file_ch++;

			if( (i%N == 0) && i) 	//	Device full
			{
				find_Checksum();	//	calculate checksum for device.

				for(int j=0; j<N+M; j++) {
					//printf("%c", Device[j]);
					fputc(Device[j], OutFiles[j]);
				}
				i=0;
			}
		}

		
	}

	//printf("\n\n%d bytes are Written...\n\n", file_ch);

	if(i) {			//	padding needed

		while(i<N) { Device[i++] = 0; printf("#"); }

		find_Checksum();	//	calculate checksum for device.

		for(int j=0; j<N+M; j++) {
			//printf("%c", Device[j]);
			fputc(Device[j], OutFiles[j]);
		}
	}

	fclose(InFile);
	for(int j=0; j<N+M; j++) fclose(OutFiles[j]);
}


void RecoverSystem(int failList[], int nn)
{


	int dc[N+M]={0}, i;		//	initial data and checksum all good

	for(i=0; i<nn; i++) dc[failList[i]-1]=1;	// 1 means failed block

	printf("\n\nGood blocks are...\n");

	for(i=0; i<N+M; i++) if(!dc[i]) printf("%d, ", i+1);	//	show good blocks

	int E[N];

	int k=0;
	i=0;

	while(k<N) {

		while(dc[i]) i++;

		printf("\n\n%d", i+1);

		if(i<N) {				//	Finding A' and E'
			for(int j=0; j<N; j++) AA[k][j] = matrixI[i][j];
			//E[k] = Device[i];
		}
		else {
			for(int j=0; j<N; j++) AA[k][j] = matrixF[i-N][j];
			//E[k] = Device[i];	//	E[k] = Checksum[i-N];
		}

		k++; i++;
	}

	char filename[][20]={"Dff1.txt", "Dff2.txt", "Dff3.txt",
						"Dff4.txt", "Dff5.txt", "Dff6.txt",
						"Dff7.txt", "Dff8.txt", "Dff9.txt", "Dff10.txt",
						"Dff11.txt", "Dff12.txt", "Dff13.txt",
						"Dff14.txt", "Dff15.txt", "Dff16.txt",
						"Dff17.txt", "Dff18.txt", "Dff19.txt", "Dff20.txt"};

	printf("\n\nFile is recovering from following files...\n");

	int goodFiles[N];
	i=0;
	for(int j=0; j<N;j++) {
		while(dc[i])i++; goodFiles[j] = i;
		i++;
		printf("%d, ", goodFiles[j]+1);
	}

	printf("\n\nNow opening those good files...");

	for(i=0; i<N; i++)
	{
		if((OutFiles[i] = fopen(filename[goodFiles[i]], "rb")) == NULL)
		{
			printf("\n\nCan't open file to recover...");
			exit(0);
		}
	}

	
	int BB[N][N];

	gInverse(AA, BB);

	int loop=0, j;

	long int file2_ch=0;

	//	getting the first element from each file as fsize

	int ch;
	for(i=0; i<N; i++) {
		fscanf(OutFiles[i],"%d", &fsize2);
		ch = fgetc(OutFiles[i]);
	}

	printf("\n\nTarget file should have %d bytes.", fsize2);

	//getch();

	for(i=0; i<N; i++) E[i] = fgetc(OutFiles[i]);

	//while(E[0] != EOF)

	while(file2_ch<fsize2)

	{

	//printf("\n\nRecovered Device...\n\n");

	for(i=0; i<N && file2_ch<fsize2; i++)
	{
		Device[i]=0;
		for(j=0; j<N; j++) Device[i] = AddSub(Device[i], gMult(BB[i][j], E[j]));

		fputc(Device[i], RecoFile); file2_ch++;
	}

	for(i=0; i<N; i++) E[i] = fgetc(OutFiles[i]);
	}

	for(i=0; i<N; i++) fclose(OutFiles[i]);

	fclose(RecoFile);

}




int main(void)
{
	int x, y;

	setup_tables();

	
	setup_matrixFI();

	OpenInAndOutFiles();

	FragmentFile();


	int nn;

	printf("\n\nThere are total %d files of Data and Checksum.", N+M);
	printf("\n\nHow many files corrupted? : ");
	scanf("%d", &nn);

	if(nn>M) {
		printf("\n\nRecover is not possible.");
		exit(0);
	}

	int failList[M];

	printf("\nEnter list(from 1 to %d) of failed files :: ", N+M);

	for(int i=0; i<nn; i++) scanf("%d", &failList[i]);

	RecoverSystem(failList, nn);



	free(gflog);
	free(gfilog);

	return 0;

}


