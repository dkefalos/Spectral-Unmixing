/*  SAD: Calculate Spectral Angle Distance of a set Signatures with an Image
    Copyright (C) 2014  Dimitris Kefalos

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

#define FP_TYPE double
FP_TYPE mySQRT(FP_TYPE);
#define myABS(x) (((x)<0)?-(x):(x))
#define mySQR(x) ((x)*(x))
#define EPSILON 1.e-19
#define ROTATE(a,i,j,k,l) {g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);}

int readtxt(FILE*, FILE*, int, int, float**);
int SAD_txt(FILE*, int, int, int, float**, float**);

int compare (const void * x, const void * y){

double da = *(const double *)x;
double db = *(const double *)y;

if (da == db)
return 0;

return (da > db) ? -1 : 1;
}

int main(int argc, char *argv[]) {

FILE     *ftxt1, *report, *ftxt2;
int      bands, nsign1, nsign2, i, nthreads;
float    **a, **b;

if(argc!=7){
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * SAD:      Calculate Spectral Angle Distance of a set Signatures       *");
    puts(" *           with another set of Signatures                              *");
    puts(" *                                                                       *");
    puts(" *         -   Input parameters (at command line) :                      *");
    puts(" *                 -        First Text File Path                         *");
    puts(" *                 -        Second Text File Path (Endmembers)           *");
    puts(" *                 -        Number of Signatures in First set            *");
    puts(" *                 -        Number of Signatures in Second set           *");
    puts(" *                 -        Number of Bands                              *");
    puts(" *                 -        Output File Path                             *");
    puts(" *                                                                       *");
    puts(" *         -   Works on the file types supported by GDAL                 *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * SAD         by Dimitris Kefalos  (dkefalos@gmail.com)                 *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    system("PAUSE");
    exit(1);
    }

/* Argument Handling */

nsign1=atoi(argv[3]);
nsign2=atoi(argv[4]);
bands=atoi(argv[5]);

if (nsign1<0){
    puts("Please give a valid number of signatures");
    exit(-1);
}

if (nsign2<0){
    puts("Please give a valid number of signatures");
    exit(-1);
}

if (bands<0){
    puts("Please give a valid number of bands");
    exit(-1);
}

/* File Opening */

if ( (report=fopen(argv[6] , "wt")) == NULL){
    printf("\nCannot open output file.\n");
    exit (-1);
    }

if ((ftxt1=fopen(argv[1], "rt")) == NULL){
    printf("\nCannot open text file 1.\n");
    exit (-1);
    }

if ((ftxt2=fopen(argv[2], "rt")) == NULL){
    printf("\nCannot open text file 2.\n");
    exit (-1);
    }

if ( (a = malloc(nsign1*sizeof(float*))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

for (i=0; i< nsign1; i++){
    if ( (a[i] = malloc(bands*sizeof(float))) == NULL){
        puts("\nBad Memory Allocation.\n");
        return(-1);
        }
    }

if ( (b = malloc(nsign2*sizeof(float*))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

for (i=0; i< nsign2; i++){
    if ( (b[i] = malloc(bands*sizeof(float))) == NULL){
        puts("\nBad Memory Allocation.\n");
        return(-1);
        }
    }

/* Read Signature Text File */
fprintf(report, "Signatures :\n");
readtxt(ftxt1, report, bands, nsign1, a);
fprintf(report, "Endmembers :\n");
readtxt(ftxt2, report, bands, nsign2, b);

/* Calculate SADs */

SAD_txt(report, nsign1, nsign2, bands, a, b);
puts("SADtxt Executed");

/* Closing */

for (i=0; i<nsign2; i++){
    free(b[i]);
}
free(b);
for (i=0; i<nsign1; i++){
    free(a[i]);
}
free(a);
fflush(NULL);
if (fclose(ftxt2)!=0){
   puts("Error in txt2 closing");
   }

if (fclose(ftxt1)!=0){
   puts("Error in txt1 closing");
   }

if (fclose(report)!=0){
   puts("Error in report closing");
   }

return(0);
}

int readtxt(FILE* txtfile, FILE* report, int bands, int nendmembers, float** a){

int      i, j;
double   tmp;
char     ignore[256];

fgets(ignore, sizeof(ignore), txtfile);
for (i=0; i<bands; i++){
    for(j=0; j<nendmembers; j++){
        if (fscanf(txtfile, "%lf ", &tmp)!=1){
            puts("not found");
            }
        a[j][i]=tmp;
        }
    fscanf(txtfile, "\n");
    }

for (i=0; i<bands; i++){
    for(j=0; j<nendmembers; j++){
        fprintf(report, "%f ", a[j][i]);
        }
    fprintf(report, "\n");
    }
fprintf(report, "\n");

return (0);
}

int SAD_txt(FILE* report, int nsign1, int nsign2, int bands, float** a, float** b){

int     i, j, k, l;
double  **pro, *sum, **sqr, **sad;

/* Memory Allocation */

if ((sum=calloc(nsign1,sizeof(double)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((pro=calloc(nsign1,sizeof(double*)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((sqr=calloc(nsign1,sizeof(double*)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((sad=malloc(nsign1*sizeof(double*)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

for (i=0; i<nsign1; i++){
    if ((pro[i]=calloc(nsign2,sizeof(double)))==NULL){
        puts("Error in Memory Allocation");
        exit(1);
        }

    if ((sqr[i]=calloc(nsign2,sizeof(double)))==NULL){
        puts("Error in Memory Allocation");
        exit(1);
        }

    if ((sad[i]=malloc(nsign2*sizeof(double)))==NULL){
        puts("Error in Memory Allocation");
        exit(1);
        }
    }

/* SAD Calculation */

for (i=0; i<nsign1; i++){
    for (j=0; j<nsign2; j++){
        for (k=0; k<bands; k++){
            pro[i][j]+=a[i][k]*b[j][k];
            sqr[i][j]+=b[j][k]*b[j][k];
            }
        }
    for(l=0; l<bands; l++){
        sum[i]+=a[i][l]*a[i][l];
        }
    }

fprintf(report, "\npro :\n");
for (i=0; i<nsign1; i++){
    for (j=0; j<nsign2; j++){
        fprintf(report, " %lf", pro[i][j]);
        }
    fprintf(report, "\n");
    }

fprintf(report, "\nsqr :\n");
for (i=0; i<nsign1; i++){
    for (j=0; j<nsign2; j++){
        fprintf(report, " %lf", sqr[i][j]);
        }
    fprintf(report, "\n");
    }

fprintf(report, "\nsum :\n");
for (i=0; i<nsign1; i++){
    fprintf(report, "%lf\n", sum[i]);
    }

for (i=0; i<nsign1; i++){
    for (j=0; j<nsign2; j++){
            sad[i][j]=acos((pro[i][j])/((sqrt(sum[i]))*(sqrt(sqr[i][j]))));
        }
    }

/* File Writing */

fprintf(report, "\nSADs :\n");
for (i=0; i<nsign1; i++){
    for (j=0; j<nsign2; j++){
        fprintf(report, " %lf", sad[i][j]);
        }
    fprintf(report, "\n");
    }

/* Closing SAD_txt */

for (i=0; i<nsign1; i++){
    free(sad[i]);
    free(sqr[i]);
    free(pro[i]);
    }
free(sad);
free(sqr);
free(pro);
free(sum);

return (0);
}
