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
int SAD_image(FILE*, FILE*, FILE*, int, int, int, int, float**);

int compare (const void * x, const void * y){

double da = *(const double *)x;
double db = *(const double *)y;

if (da == db)
return 0;

return (da > db) ? -1 : 1;
}

int main(int argc, char *argv[]) {

FILE     *fin, *report, *ftxt, *fsad;
int      rows, cols, bands, nendmembers, i, nthreads;
char     ersfile[100], binfile[100], tiffile[100], sysstring[200], txtfile[100], sadfile[100];
float    **a;

if(argc!=7){
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * SAD:      Calculate Spectral Angle Distance of a set Signatures       *");
    puts(" *           with an Image                                               *");
    puts(" *                                                                       *");
    puts(" *         -   Input parameters (at command line) :                      *");
    puts(" *                 -        Image filename                               *");
    puts(" *                 -        Text file Path                               *");
    puts(" *                 -        Number of Rows                               *");
    puts(" *                 -        Number of Columns                            *");
    puts(" *                 -        Number of Bands                              *");
    puts(" *                 -        Number of Signatures                         *");
    puts(" *                                                                       *");
    puts(" *         -   Works on the file types supported by GDAL                 *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * Output : - Report :    Imagename_SAD_report.txt                       *");
    puts(" *          - SAD Image : Imagename_SAD                                  *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * SAD         by Dimitris Kefalos  (dkefalos@gmail.com)                 *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    system("PAUSE");
    exit(1);
    }

/* File Manipulation */

strcpy(tiffile, argv[1]);
i=(strlen(tiffile)-4);
memcpy(binfile, tiffile, strlen(tiffile)-4);
binfile[i]='\0';
i=(strlen(binfile));
memcpy(ersfile, binfile, strlen(binfile));
ersfile[i]='\0';
sprintf(ersfile, "%s.ers", binfile);
memcpy(txtfile, binfile, strlen(binfile));
txtfile[i]='\0';
sprintf(txtfile, "%s_SAD_report.txt", binfile);
i=(strlen(binfile));
memcpy(sadfile, binfile, strlen(binfile));
sadfile[i]='\0';
sprintf(sadfile, "%s_SAD", binfile);

/* Argument Handling */

rows=atoi(argv[3]);
cols=atoi(argv[4]);
bands=atoi(argv[5]);
nendmembers=atoi(argv[6]);

if (rows<0){
    puts("Please give a valid number of rows");
    exit(-1);
}

if (cols<0){
    puts("Please give a valid number of columns");
    exit(-1);
}

if (bands<0){
    puts("Please give a valid number of bands");
    exit(-1);
}

if (nendmembers<0){
    puts("Please give a valid number of endmembers");
    exit(-1);
}

/* GDAL Call */

sprintf(sysstring, "gdal_translate -ot Float32 -of ERS %s %s", tiffile, ersfile);
printf("\n");
printf(sysstring);
printf("\n");
if ((i=system(sysstring))!=0 && (i=system(sysstring))!=1){
    puts("\nGDAL did not run correctly");
    exit(-1);
    }
puts("GDAL has Executed");

/* File Opening */

if ( (fin=fopen(binfile , "rb")) == NULL){
    printf("\nCannot open %s file.\n", binfile);
    exit (-1);
    }

if ( (fsad=fopen(sadfile , "wb+")) == NULL){
    printf("\nCannot open %s file.\n", sadfile);
    exit (-1);
    }

if ( (report=fopen(txtfile , "wt")) == NULL){
    printf("\nCannot open %s file.\n", txtfile);
    exit (-1);
    }

if ((ftxt=fopen(argv[2], "rt")) == NULL){
    printf("\nCannot open text file.\n");
    exit (-1);
    }

if ( (a = malloc(nendmembers*sizeof(float*))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

for (i=0; i< nendmembers; i++){
    if ( (a[i] = malloc(bands*sizeof(float))) == NULL){
        puts("\nBad Memory Allocation.\n");
        return(-1);
        }
    }

/* Read Signature Text File */

readtxt(ftxt, report, bands, nendmembers, a);

/* Create SAD Image */

SAD_image(fin, fsad, report, rows, cols, bands, nendmembers, a);

/* Closing */

for (i=0; i<nendmembers; i++){
    free(a[i]);
}
free(a);
fflush(NULL);
if (fclose(ftxt)!=0){
   puts("Error in txt closing");
   }

if (fclose(report)!=0){
   puts("Error in report closing");
   }

if (fclose(fsad)!=0){
   puts("Error in fin closing");
   }

if (fclose(fin)!=0){
   puts("Error in fin closing");
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

fprintf(report, "Signatures :\n");
for (i=0; i<bands; i++){
    for(j=0; j<nendmembers; j++){
        fprintf(report, "%f ", a[j][i]);
        }
    fprintf(report, "\n");
    }

return (0);
}

int SAD_image(FILE* fin, FILE* fsad, FILE* report, int rows, int cols, int bands, int nendmembers, float** a){

int     i, j, k, l;
float   *buffer, *buffer2;
double  **pro, *sum, **sqr;

/* Memory Allocation */

if ((buffer=malloc(cols*bands*sizeof(float)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((buffer2=malloc(cols*nendmembers*sizeof(float)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((sum=malloc(nendmembers*sizeof(double)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((pro=malloc(cols*sizeof(double*)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((sqr=malloc(cols*sizeof(double*)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

for (i=0; i<cols; i++){
    if ((pro[i]=malloc(nendmembers*sizeof(double)))==NULL){
        puts("Error in Memory Allocation");
        exit(1);
        }

    if ((sqr[i]=malloc(nendmembers*sizeof(double)))==NULL){
        puts("Error in Memory Allocation");
        exit(1);
        }
    }

/* SAD Calculation and File Writing */

for (i=0; i<rows; i++){
    fread(buffer, sizeof(float), cols*bands, fin);
    printf("Calculating SADs at row %d of %d\r", i+1, rows);
    for (j=0; j<cols; j++){
        for (k=0; k<nendmembers; k++){
            pro[j][k]=0.0;
            sum[k]=0.0;
            sqr[j][k]=0.0;
            for (l=0; l<bands; l++){
                pro[j][k]+=a[k][l]*buffer[l*cols+j];
                sum[k]+=a[k][l]*a[k][l];
                sqr[j][k]+=buffer[l*cols+j]*buffer[l*cols+j];
                }
            }
        }
    for (j=0; j<nendmembers; j++){
        for (k=0; k<cols; k++){
            buffer2[j*cols+k]=acos((pro[k][j])/((sqrt(sum[j]))*(sqrt(sqr[k][j]))));
            }
        }
    fwrite(buffer2, sizeof(float), nendmembers*cols, fsad);
    }

/* Closing SAD_image */

for (i=0; i<cols; i++){
    free(sqr[i]);
    free(pro[i]);
    }
free(sqr);
free(pro);
free(sum);
free(buffer2);
free(buffer);

return (0);
}
