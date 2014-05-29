/*  MNF - Performs MNF transformation on hyperspectral Image
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
# include <omp.h>
# include <time.h>

#define FP_TYPE double
FP_TYPE mySQRT(FP_TYPE);
#define myABS(x) (((x)<0)?-(x):(x))
#define mySQR(x) ((x)*(x))
#define EPSILON 1.e-19
#define ROTATE(a,i,j,k,l) {g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);}

void jacobi(double **, int, double *, double **, int *);
void eigsort(double *, double **, int);
int MNF(FILE*, FILE*, FILE*, int, int, int, int, int, int);
int Matrix_inversion(FP_TYPE*, FP_TYPE*, int);
int Cholesky_decomposition(FP_TYPE *, FP_TYPE *, int );
int NND(FILE*, FILE*, int, int, int, int);
int MRTBM(FILE*, FILE*, FILE*, int, int, int, int);

int compare (const void * x, const void * y){

double da = *(const double *)x;
double db = *(const double *)y;

if (da == db)
return 0;

return (da > db) ? -1 : 1;
}

int main(int argc, char *argv[]) {

FILE     *fin, *report, *fsignal;
int      rows, cols, bands, i, j, k, m, nthreads, opt, nvalue, npixels;
char     ersfile[100], binfile[100], tiffile[100], sysstring[200], txtfile[100], mnffile[100];
float    *buffer;

if(argc!=7){
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * MNF:   Perform MNF Transformation on Hyperspectral Image              *");
    puts(" *                                                                       *");
    puts(" *         -   Input parameters (at command line) :                      *");
    puts(" *                 -        Image filename                               *");
    puts(" *                 -        Rows                                         *");
    puts(" *                 -        Columns                                      *");
    puts(" *                 -        Bands                                        *");
    puts(" *                 -        Noise Estimation Method                      *");
    puts(" *                          1 for NND                                    *");
    puts(" *                          2 for MRTBM                                  *");
    puts(" *                 -        Null Cell Value                              *");
    puts(" *                                                                       *");
    puts(" *         -   Works on the file types supported by GDAL                 *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * Output : - Report :    Imagename_MNF_report.txt                       *");
    puts(" *          - MNF Image : Imagename_MNF                                  *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * Unmixing    by Dimitris Kefalos  (dkefalos@gmail.com)                 *");
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
sprintf(txtfile, "%s_MNF_report.txt", binfile);
i=(strlen(binfile));
memcpy(mnffile, binfile, strlen(binfile));
mnffile[i]='\0';
sprintf(mnffile, "%s_MNF", binfile);

/* Argument Handling */

rows=atoi(argv[2]);
cols=atoi(argv[3]);
bands=atoi(argv[4]);
opt=atoi(argv[5]);
nvalue=atoi(argv[6]);

if (opt!=1 && opt!=2){
    puts("Please choose valid noise estimation method");
    exit(-1);
}

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

/* GDAL Call */

sprintf(sysstring, "gdal_translate -ot Float32 -of ERS %s %s", tiffile, ersfile);
printf("\n");
printf(sysstring);
printf("\n");
if ((i=system(sysstring))!=0 && (i=system(sysstring))!=1 ){
    printf("\nGDAL did not run correctly %d", i);
    exit(-1);
    }
puts("GDAL has Executed");

/* File Opening */

if ( (fin=fopen(binfile , "rb")) == NULL){
    printf("\nCannot open %s file.\n", binfile);
    exit (-1);
    }

if ( (report=fopen(txtfile , "wt")) == NULL){
    printf("\nCannot open %s file.\n", txtfile);
    exit (-1);
    }

if ((fsignal=fopen(mnffile , "wb+")) == NULL){
    printf("\nCannot open MNF file.\n");
    exit (-1);
    }

if ( (buffer = malloc(cols*bands*sizeof(float))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

nthreads=omp_get_max_threads();
printf("Number of Threads : %d\n\n", nthreads);

/* Masked Pixels Count */

npixels=0;
for (i=0; i<rows; i++){
    fread(buffer, cols*bands, sizeof(float), fin);
    for (j=0; j<cols; j++){
        m=0;
        for (k=0; k<bands; k++){
            if (buffer[k*cols+j]==nvalue){
                m++;
                }
            }
        if (m==bands){
            npixels++;
            }
        }
    }
rewind(fin);
printf("\nMasked pixels=%d\n", npixels);

/* Transformation */

MNF(fin, fsignal, report, rows, cols, bands, opt, nvalue, npixels);

/* Closing */

fflush(NULL);
free(buffer);
if (fclose(fsignal)!=0){
   puts("Error in fsignal closing");
   }

if (fclose(report)!=0){
   puts("Error in report closing");
   }

if (fclose(fin)!=0){
   puts("Error in fin closing");
   }

return(0);
}

int MNF(FILE *fin, FILE *fsignal, FILE* report, int rows, int cols, int bands, int opt, int nvalue, int npixels){

FILE    *fnoise, *fpca1;
int     i, j, k, l, nrot;
float   *buffer, *buffer2, *sum;
double  *Nmean, **Ncov, *Neigval, **Neigvec, **sumv, tstart, tend;

/* File Opening */

if ((fpca1=tmpfile()) == NULL){
    printf("\nCannot open signal file.\n");
    exit (-1);
    }

if ((fnoise=tmpfile()) == NULL){
    printf("\nCannot open noise file.\n");
    exit (-1);
    }

/* Memory Allocation */

if ((buffer2=malloc(cols*bands*sizeof(float)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((buffer=malloc(cols*bands*sizeof(float)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((sum=calloc(bands,sizeof(float)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ( (sumv=calloc(bands,sizeof(double*))) == NULL ){
    puts("Error in Memory Allocation.");
    exit(-1);
    }

if ((Nmean=malloc(bands*sizeof(double)))== NULL){
    puts("Error in Nmean Memory Allocation");
	exit(-1);
    }

if ( (Ncov=malloc(bands*sizeof(double*))) == NULL ){
    puts("Error in Ncov Memory Allocation.");
    exit(-1);
    }

if ( (Neigval=malloc(bands*sizeof(double))) == NULL ){
    puts("Error in Neigenval Memory Allocation.");
    exit(-1);
    }

if ( (Neigvec=malloc(bands*sizeof(double*))) == NULL ){
    puts("Error in Neigenvec Memory Allocation.");
    exit(-1);
    }

for (i=0; i<bands; i++){
    if ( (Ncov[i]=malloc(bands*sizeof(double))) == NULL ){
        puts("Error in (Ncov[i]) Memory Allocation.");
        exit(-1);
        }
    if ( (Neigvec[i]=malloc(bands*sizeof(double))) == NULL ){
        puts("Error in (Neigenvec[i]) Memory Allocation.");
        exit(-1);
        }
    if ((sumv[i]=calloc(bands,sizeof(double))) == NULL){
        puts("Error in Memory Allocation.");
        exit(-1);
        }
}

/* Noise Estimation Option */

if (opt==1){
    NND(fin, fnoise, rows, cols, bands, nvalue);
    puts("Noise calculated using NND Method");
}

if (opt==2){
    MRTBM(fin, fnoise, report, rows, cols, bands, nvalue);
    puts("Noise calculated using MRTBM Method");
}

/* Noise Sum Calculation for Every Band */

rewind(fnoise);
for(i=0; i<rows; i++){
    fread(buffer, sizeof(float), cols*bands, fnoise);
    for(j=0; j<bands; j++){
        for(k=0; k<cols; k++){
            if (buffer[j*cols+k]!=nvalue){
                sum[j]+=buffer[j*cols+k];
                }
            }
        }
    }

/* Noise Mean Calculation for Every Band */

for (i=0; i<bands; i++){
    Nmean[i]=sum[i]/(rows*cols-npixels);
    }

fprintf(report, "\nNoise Mean :\n");
for (i=0; i<bands; i++){
    fprintf(report, "Band %d : %10.5lf\n", i+1, Nmean[i]);
    }

/* sumv Calculation for Every Band */

rewind(fnoise);
tstart=omp_get_wtime();

# pragma omp parallel default(none) private(i, j, k, l) shared(rows, cols, bands, fnoise, buffer, sumv, Nmean, nvalue)
{
for (i=0; i<rows; i++){
    # pragma omp single
    {
    fread (buffer, sizeof(float), cols*bands, fnoise);
    printf("Calculating Image Statistics at row %d of %d\r", i+1, rows);
    }
    # pragma omp for schedule(dynamic)
    for (j=0; j<bands; j++){
        for (k=0; k<bands; k++){
            for (l=0; l<cols; l++){
                if (buffer[(j*cols)+l]!=nvalue){
                    sumv[j][k]+=(buffer[(j*cols)+l]-Nmean[j])*(buffer[(k*cols)+l]-Nmean[k]);
                    }
                }
            }
        }
    }
}

tend=omp_get_wtime();
printf("\nDuration of Statistics Parallel region : %2.1lf sec\n", tend-tstart);

fprintf(report, "\nsumv:\n");
for (i=0; i<bands; i++){
    fprintf (report, "Band %d", i+1);
    for (j=0; j<bands; j++){
        fprintf(report, " %10.5lf", sumv[i][j]);
        }
    fprintf(report, "\n");
    }

/* Noise Variance-Covariance Matrix Calculation */

for (i=0; i<bands; i++){
    for (j=0; j<bands; j++){
        Ncov[i][j]=sumv[i][j]/(rows*cols-npixels);
        }
    }

fprintf(report, "\nNoise Variance-Covariance Matrix:\n");
for (i=0; i<bands; i++){
    fprintf (report, "Band %d", i+1);
    for (j=0; j<bands; j++){
        fprintf(report, " %10.5lf", sqrt(fabs(Ncov[i][j])));
        }
    fprintf(report, "\n");
    }

/* Noise Eigenvectors Calculation and Sorting */

jacobi(Ncov, bands, Neigval, Neigvec, &nrot);
eigsort(Neigval, Neigvec, bands);
fflush(NULL);

fprintf(report, "\nNoise Eigenvectors:\n");
for (i=0; i<bands; i++){
    for (j=0; j<bands; j++){
        fprintf(report,"%10.5f ", Neigvec[i][j]);
        }
    fprintf(report, "\n");
    }

fprintf(report,"\nNoise Eigenvalues :\n");
for (i=0; i<bands; i++){
    fprintf(report," %10.5lf\n", Neigval[i]);
    }

for (i=0; i<bands; i++){
    sum[i]=0.0;
    for (j=0; j<bands; j++){
        sumv[i][j]=0.0;
        }
    }

/* Sum Calculation for Every Band */

rewind(fin);
for(i=0; i<rows; i++){
    fread(buffer, sizeof(float), cols*bands, fin);
    for(j=0; j<bands; j++){
        for(k=0; k<cols; k++){
            if (buffer[j*cols+k]!=nvalue){
                sum[j]+=buffer[j*cols+k];
                }
            }
        }
    }

/* Mean Calculation for Every Band */

for (i=0; i<bands; i++){
    Nmean[i]=sum[i]/(rows*cols-npixels);
    }

fprintf(report, "\nData Mean :\n");
for (i=0; i<bands; i++){
    fprintf(report, "Band %d : %10.5lf\n", i+1, Nmean[i]);
    }

/* White Noise Image Creation */

rewind(fin);
tstart=omp_get_wtime();
# pragma omp parallel default(none) private (i, j, k , l) shared(buffer, buffer2, rows, cols, bands, Neigvec, fin, fpca1, Nmean, Ncov, nvalue)
{
for (i=0; i<rows; i++){
    # pragma omp single
    {
    fread(buffer, sizeof(float), cols*bands, fin);
    printf("Calculating White Noise Image at row %d of %d\r", i+1, rows);
    }
    # pragma omp for schedule(dynamic)
    for (j=0; j<bands; j++){
        for (k=0; k<cols; k++){
            buffer2[j*cols+k]=0.0;
            for(l=0; l<bands; l++){
                if (buffer[l*cols+k]!=nvalue){
                    buffer2[j*cols+k]+=1.*(buffer[(l*cols)+k]-Nmean[l])*(Neigvec[l][j]);
                    }
                else {
                    buffer2[j*cols+k]=nvalue;
                    }
                }
            if (buffer2[j*cols+k]!=nvalue){
                buffer2[j*cols+k]=buffer2[j*cols+k]/sqrt(fabs(Ncov[j][j]));
                }
            }
        }
    # pragma omp single
    {
    fwrite(buffer2,sizeof(float),bands*cols,fpca1);
    }
    }
}
tend=omp_get_wtime();
printf("\nDuration of Noise Whitening Parallel region : %2.1lf sec\n", tend-tstart);

for (i=0; i<bands; i++){
    sum[i]=0.0;
    for (j=0; j<bands; j++){
        sumv[i][j]=0.0;
        }
    }

/* White Noise Image Sum Calculation for Every Band */

rewind(fpca1);
for(i=0; i<rows; i++){
    fread(buffer, sizeof(float), cols*bands, fpca1);
    for(j=0; j<bands; j++){
        for(k=0; k<cols; k++){
            if (buffer[j*cols+k]!=nvalue){
                sum[j]+=buffer[j*cols+k];
                }
            }
        }
    }

/* White Noise Image Mean Calculation for Every Band */

for (i=0; i<bands; i++){
    Nmean[i]=sum[i]/(rows*cols-npixels);
    }

fprintf(report, "\nWhite Noise Mean :\n");
for (i=0; i<bands; i++){
    fprintf(report, "Band %d : %10.5lf\n", i+1, Nmean[i]);
    }

/* White Noise Image sumv Calculation for Every Band */

rewind(fpca1);
tstart=omp_get_wtime();
# pragma omp parallel default(none) private(i, j, k, l) shared(buffer, fpca1, Nmean, sumv, rows, cols, bands, nvalue)
{
for (i=0; i<rows; i++){
    # pragma omp single
    {
    fread (buffer, sizeof(float), cols*bands, fpca1);
    printf("Calculating White Noise Image Statistics at row %d of %d\r", i+1, rows);
    }
    # pragma omp for schedule(dynamic)
    for (j=0; j<bands; j++){
        for (k=0; k<bands; k++){
            for (l=0; l<cols; l++){
                if (buffer[(j*cols)+l]!=nvalue){
                    sumv[j][k]+=(buffer[(j*cols)+l]-Nmean[j])*(buffer[(k*cols)+l]-Nmean[k]);
                    }
                }
            }
        }
    }
}

tend=omp_get_wtime();
printf("\nDuration of Statistics Parallel region : %2.1lf sec\n", tend-tstart);

fprintf(report, "\nWhite sumv:\n");
for (i=0; i<bands; i++){
    fprintf (report, "Band %d", i+1);
    for (j=0; j<bands; j++){
        fprintf(report, " %10.5lf", sumv[i][j]);
        }
    fprintf(report, "\n");
    }

/* White Noise Image Covariance Matrix Calculation */

for (i=0; i<bands; i++){
    for (j=0; j<bands; j++){
        Ncov[i][j]=sumv[i][j]/(rows*cols-npixels);
        }
    }

fprintf(report, "\nWhite Noise Variance-Covariance Matrix:\n");
for (i=0; i<bands; i++){
    fprintf (report, "Band %d", i+1);
    for (j=0; j<bands; j++){
        fprintf(report, " %15.5lf", Ncov[i][j]);
        }
    fprintf(report, "\n");
    }

/* White Noise Image Eigenvectors and Sorting */

jacobi(Ncov, bands, Neigval, Neigvec, &nrot);
eigsort(Neigval, Neigvec, bands);

fprintf(report, "\nWhite Noise Eigenvectors:\n");
for (i=0; i<bands; i++){
    for (j=0; j<bands; j++){
        fprintf(report,"%15.5f ", Neigvec[i][j]);
        }
    fprintf(report, "\n");
    }

fflush(NULL);

/* MNF Image Creation */

tstart=omp_get_wtime();
rewind(fpca1);
# pragma omp parallel default(none) private(i, j, k, l) shared(buffer, buffer2, rows, cols, bands, fpca1, fsignal, Neigvec, Nmean, nvalue)
{
for (i=0; i<rows; i++){
    # pragma omp single
    {
    fread(buffer, sizeof(float), cols*bands, fpca1);
    printf("Calculating Pr. Components of Wh. Noise Image at row %d of %d\r", i+1, rows);
    }
    # pragma omp for schedule(dynamic)
    for (j=0; j<bands; j++){
        for (k=0; k<cols; k++){
            buffer2[j*cols+k]=0.0;
            for(l=0; l<bands; l++){
                if (buffer[l*cols+k]!=nvalue){
                    buffer2[j*cols+k]+=1.*(buffer[(l*cols)+k]-Nmean[l])*(Neigvec[l][j]);
                    }
                else {
                    buffer2[j*cols+k]=nvalue;
                    }
                }
            }
        }
    # pragma omp single
    {
    fwrite(buffer2,sizeof(float),bands*cols,fsignal);
    }
    }
}
tend=omp_get_wtime();
printf("\nDuration of Wh. Noise PCA Parallel region : %2.1lf sec", tend-tstart);
fflush(NULL);

/* Closing MNF */

for(i=0; i<bands; i++){
    free(sumv[i]);
    free(Neigvec[i]);
    free(Ncov[i]);
    }
free(Neigvec);
free(Neigval);
free(Ncov);
free(Nmean);
free(sumv);
free(sum);
free(buffer);
free(buffer2);
if (fclose(fpca1)!=0){
   puts("Error in fpca1 closing");
   }
if (fclose(fnoise)!=0){
   puts("Error in fnoise closing");
   }

puts("\nMNF Executed");
return (0);
}

int NND(FILE* fin, FILE* fnoise, int rows, int cols, int bands, int nvalue){

int     i, j, k;
float   *buffer, *buffer2, *noise1, dif[2];

/* Memory Allocation */

if ((buffer2=malloc(cols*bands*sizeof(float)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((buffer=malloc(cols*bands*sizeof(float)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((noise1=malloc(cols*bands*sizeof(float)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

/* Noise Estimation of the first row */

for (i=0; i<bands*cols; i++){
    noise1[i]=0.0;
    }
fwrite(noise1 , sizeof(float), cols*bands, fnoise);

// Noise Estimation of the rest

rewind(fin);
fread(buffer2, sizeof(float), cols*bands,fin);

for (i=1; i<rows; i++){
    fread(buffer2, sizeof(float), cols*bands, fin);
    for (j=0; j<bands; j++){
        for (k=1; k<cols-1; k++){
            if (buffer[j*cols+k]== nvalue || buffer2[j*cols+k]==nvalue || buffer2[j*cols+k+1]==nvalue){
                noise1[j*cols+k]=nvalue;
                }
            else {
                dif[0]=buffer2[j*cols+k]-buffer2[j*cols+k+1];
                dif[1]=buffer2[j*cols+k]-buffer[j*cols+k];
                noise1[j*cols+k]=(dif[0]+dif[1])/2;
                }
            }
        }
    fwrite(noise1, sizeof(float), cols*bands, fnoise);
    memcpy(buffer, buffer2, cols*bands*sizeof(float));
    }

/* NND Closing */

free(noise1);
free(buffer);
free(buffer2);

return (0);
}

int MRTBM(FILE* fin, FILE* fnoise, FILE* report, int rows, int cols, int bands, int nvalue){

int     i, j, k, l;
float   *noise2, *buffer, **Z, *z;
double  **zTz, *zTzl, *zTzl_inv, **zTz_inv, **R, *b, *c, tstart, tend, tmp;

/* Memory Allocation */

if ((noise2=malloc(cols*rows*sizeof(float)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ((buffer=malloc(cols*bands*sizeof(float)))==NULL){
    puts("Error in Memory Allocation");
	exit(1);
    }

if ( (zTz=calloc(bands,sizeof(double*))) == NULL ){
    puts("Error in zTz Memory Allocation.");
    exit(-1);
    }

if ( (zTzl=malloc(bands*bands*sizeof(double))) == NULL ){
     puts("Error in zTzl Memory Allocation.");
     exit(-1);
     }

if ( (zTzl_inv=malloc(bands*bands*sizeof(double))) == NULL ){
     puts("Error in zTzl_inv Memory Allocation.");
     exit(-1);
     }

if ( (zTz_inv=malloc(bands*sizeof(double*))) == NULL ){
     puts("Error in zTz_inv Memory Allocation.");
     exit(-1);
     }

if ( (R=calloc((bands-1),sizeof(double*))) == NULL ){
     puts("Error in R Memory Allocation.");
     exit(-1);
     }

if ( (b=calloc((bands-1),sizeof(double))) == NULL ){
     puts("Error in b Memory Allocation.");
     exit(-1);
     }

if ( (c=calloc(rows*cols,sizeof(double))) == NULL ){
     puts("Error in c Memory Allocation.");
     exit(-1);
     }

if ( (z=calloc(rows*cols,sizeof(float))) == NULL ){
     puts("Error in z Memory Allocation.");
     exit(-1);
     }

if ( (Z=malloc(rows*cols*sizeof(float*))) == NULL ){
     puts("Error in c Memory Allocation.");
     exit(-1);
     }

for (i=0; i<bands; i++){
    if ( (zTz[i]=calloc(bands,sizeof(double))) == NULL ){
       puts("Error in zTz[i] Memory Allocation.");
       exit(-1);
       }
    if ( (zTz_inv[i]=malloc(bands*sizeof(double))) == NULL ){
       puts("Error in zTz_inv[i] Memory Allocation.");
       exit(-1);
       }
    }

for( i=0; i<rows*cols; i++){
    if ( (Z[i]=malloc((bands-1)*sizeof(float))) == NULL ){
        puts("Error in Z[i] Memory Allocation.");
        exit(-1);
        }
    }
for (i=0; i<bands-1; i++){
    if ( (R[i]=calloc(rows*cols,sizeof(double))) == NULL ){
       puts("Error in R[i] Memory Allocation.");
       exit(-1);
        }
}

tstart=omp_get_wtime();
for (i=0; i<bands; i++){
    printf("Calcul. Noise using Mul. Regression at band %d of %d", i+1, bands);

/* Z Matrix Formation */

    for (j=0; j<i; j++){
        fseek(fin, j*cols*sizeof(float), SEEK_SET);
        for (k=0; k<rows; k++){
            fread(buffer, sizeof(float), cols, fin);
            for (l=0; l<cols; l++){
                Z[k*cols+l][j]=buffer[l];
                }
            fseek(fin, (bands-1)*cols*sizeof(float), SEEK_CUR);
            }
        }

    for (j=i; j<bands-1; j++){
        fseek(fin, (j+1)*cols*sizeof(float), SEEK_SET);
        for (k=0; k<rows; k++){
            fread(buffer, sizeof(float), cols, fin);
            for (l=0; l<cols; l++){
                Z[k*cols+l][j]=buffer[l];
                }
            fseek(fin, (bands-1)*cols*sizeof(float), SEEK_CUR);
            }
        }

    fseek(fin, i*cols*sizeof(float), SEEK_SET);
    for (j=0; j<rows; j++){
        fread(buffer, sizeof(float), cols, fin);
        for (k=0; k<cols; k++){
            z[j*cols+k]=buffer[k];
            }
        fseek(fin, (bands-1)*cols*sizeof(float), SEEK_CUR);
        }

    for (j=0; j<bands-1; j++){
        for (k=0; k<bands-1; k++){
            zTz[j][k]=0.0;
        }
    }

/* zTz Matrix Calculation */

//# pragma omp parallel for default(none) private(j, k, l, tmp) shared(rows, cols, bands, zTz, Z)
    for (j=0; j<bands-1; j++){
        for (k=0; k<bands-1; k++){
            tmp=0.0;
            for (l=0; l<rows*cols; l++){
                tmp+=Z[l][k]*Z[l][j];
                }
            zTz[j][k]=tmp;
            }
        }

/* zTz Matrix Inversion */

    for (j=0; j<bands-1; j++){
        for (k=0; k<bands-1; k++){
            zTzl[j*(bands-1)+k]=zTz[j][k];
            }
        }

    Matrix_inversion(zTzl, zTzl_inv, bands-1);

    for (j=0; j<bands-1; j++){
        for (k=0; k<bands-1; k++){
            zTz_inv[j][k]=zTzl_inv[j*(bands-1)+k];
            }
        }

/* R Matrix Calculation */

    for (j=0; j<bands-1; j++){
        for (k=0; k<rows*cols; k++){
            tmp=0.0;
            for (l=0; l<bands-1; l++){
                tmp+=zTz_inv[j][l]*Z[k][l];
                }
            R[j][k]=tmp;
            }
        }

/* Regression Factors Calculation */

    for (j=0; j<bands-1; j++){
        tmp=0.0;
        for (k=0; k<rows*cols; k++){
            tmp+=R[j][k]*z[k];
            }
        b[j]=tmp;
        }

/* c Matrix Calculation */

    for (j=0; j<rows*cols; j++){
        tmp=0.0;
        for (k=0; k<bands-1; k++){
            tmp+=Z[j][k]*b[k];
            }
        c[j]=tmp;
        }

/* Noise Calculation */

    for (j=0; j<rows*cols; j++){
        if (z[j]!=nvalue){
            noise2[j]=(float)(z[j]-c[j]);
            }
        else {
            noise2[j]=nvalue;
            }
        }

fprintf(report, "\nnoise=\n");
for (j=0; j<400; j++){
    fprintf(report, "%f ", noise2[j]);
    }
fprintf(report, "\n");

/* Noise Image Writing */

    if (fseek(fnoise, i*cols*sizeof(float), SEEK_SET) != 0){
       puts ("Bad File Rewind");
       }

    for (j=0; j<rows; j++){
        fwrite(noise2+j*cols, sizeof(float), cols, fnoise);
        if (fseek ( fnoise, cols*(bands-1)*sizeof(float), SEEK_CUR) != 0){
           puts( "Bad File Control");
           }
        fflush(NULL);
        }

tend=omp_get_wtime();
printf(" Duration: %2.2lf min\r", (tend-tstart)/60.0);
    }

/* MRTBM Closing */

for(i=0; i<bands-1; i++){
    free(R[i]);
    }
for(i=0; i<rows*cols; i++){
    free(Z[i]);
    }
for(i=0; i<bands; i++){
    free(zTz_inv[i]);
    free(zTz[i]);
    }
free(Z);
free(z);
free(c);
free(b);
free(R);
free(zTz_inv);
free(zTzl_inv);
free(zTzl);
free(zTz);
free(buffer);
free(noise2);

return(0);
}

int Matrix_inversion(FP_TYPE *ata,FP_TYPE *ata_inv,int n) {
    FP_TYPE *tmpa, *b;
    int *l;
    int r_value;
    register int i, j, k;

    tmpa=malloc(n*n*sizeof(FP_TYPE));
    b=malloc(n*sizeof(FP_TYPE));
    l=malloc(n*sizeof(int));
    r_value=1;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            for (k=0; k<n; k++)
                tmpa[j*n+k]=ata[j*n+k];
            b[j]=0.0;
        }
        b[i]=1.0;
        r_value*=Cholesky_decomposition(tmpa, b, n);
        for (j=0; j<n; j++)
            ata_inv[j*n+i]=b[j];
    }
    free(l);
    free(b);
    free(tmpa);
    return r_value;
}

int Cholesky_decomposition(FP_TYPE *a, FP_TYPE *b, int n) {
    FP_TYPE s;
    register int i, j, k;

    for (i=0; i<n; i++) {
        s=0.0;
        if (i)
            for (j=0; j<i; j++)
                s+=mySQR(a[i*n+j]);
        a[i*n+i]=mySQRT(a[i*n+i]-s);
        if (myABS(a[i*n+i])<EPSILON)
            return 0;
        if (i<n-1)
            for (j=i+1; j<n; j++) {
                s=0.0;
                if (i)
                    for (k=0; k<i; k++)
                        s+=a[i*n+k]*a[j*n+k];
                a[j*n+i]=(a[j*n+i]-s)/a[i*n+i];
            }
    }
    b[0]/=a[0*n+0];
    for (k=1; k<n; k++) {
        s=0.0;
        for (j=0; j<k; j++)
            s+=a[k*n+j]*b[j];
        b[k]=(b[k]-s)/a[k*n+k];
    }
    b[n-1]/=a[(n-1)*n+n-1];
    for (i=0; i<n-1; i++) {
        s=0.0;
        for (j=n-i-1; j<n; j++)
            s+=a[j*n+n-i-2]*b[j];
        b[n-i-2]=(b[n-i-2]-s)/a[(n-i-2)*n+n-i-2];
    }
    return 1;
}

FP_TYPE mySQRT(FP_TYPE x) {
    FP_TYPE px,a;

    a=px=x;
    x/=2.0;
    while (myABS(x-px)>EPSILON) {
        px=x;
        x=(px+a/px)/2;
    }
    return x;
}

/* Jacobi function computes all eigenvalues and eigenvectors of a real symmetric matrix
a[1..n][1..n]. On output, elements of a above the diagonal are destroyed.
d[1..n] returns the eigenvalues of a.	v[1..n][1..n] is a matrix whose
columns contain, on output, the normalized eigenvectors of a.
nrot returns the number of Jacobi rotations that were required. */

void jacobi(double **a, int n, double *d, double **v, int *nrot){

int j,iq,ip,i;
double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

if ((b=(double *)malloc(n*sizeof(double)))==NULL){
    puts("Jacobi reports: Error in Memory Allocation");
	exit(1);
    }

if ((z=(double *)malloc(n*sizeof(double)))==NULL){
    puts("Jacobi reports: Error in Memory Allocation");
    exit(1);
    }
for(ip=0;ip<n;ip++){
	for(iq=0;iq<n;iq++)
		v[ip][iq]=0.0;
	v[ip][ip]=1.0;
	}

for(ip=0;ip<n;ip++){
	d[ip]=b[ip]=a[ip][ip];
	z[ip]=0.0;
	}

*nrot=0;

for(i=1;i<=50;i++){
	sm=0.0;
	for(ip=0;ip<n-1;ip++){
		for(iq=ip+1;iq<n;iq++)
			sm+=fabs(a[ip][iq]);
		}
	if(fabs(sm)<1.0e-14){
		free(z);
		free(b);
		return;
		}
	if(i<4)
        tresh=0.2*sm/(n*n);
	else
		tresh=0.0;
	for(ip=0;ip<n-1;ip++){
		for(iq=ip+1;iq<n;iq++){
			g=100.0*fabs(a[ip][iq]);
			if (i>4&&(fabs(d[ip])+g)==fabs(d[ip])&&(fabs(d[iq])+g)==fabs(d[iq]))
				a[ip][iq]=0.0;
			else if(fabs(a[ip][iq])>tresh){
				h=d[iq]-d[ip];
                if((fabs(h)+g)==fabs(h))
                    t=(a[ip][iq])/h;
                else{
                    theta=0.5*h/(a[ip][iq]);
                    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                    if (theta<0.0)
                        t=-t;
                    }
                c=1.0/sqrt(1+t*t);
                s=t*c;
                tau=s/(1.0+c);
                h=t*a[ip][iq];
                z[ip]-=h;
                z[iq]+=h;
                d[ip]-=h;
                d[iq]+=h;
                a[ip][iq]=0.0;
		    	for(j=0;j<=ip-1;j++){
					ROTATE(a,j,ip,j,iq)
					}
				for(j=ip+1;j<=iq-1;j++){
                    ROTATE(a,ip,j,j,iq)
			    	}
				for(j=iq+1;j<n;j++){
                    ROTATE(a,ip,j,iq,j)
                    }
			    for(j=0;j<n;j++){
					ROTATE(v,j,ip,j,iq)
					}
			    ++(*nrot);
				}
		    }
        }
    for(ip=0;ip<n;ip++){
		b[ip]+=z[ip];
		d[ip]=b[ip];
		z[ip]=0.0;
		}
    }
}

/* Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output
from jacobi routine, eigensort function sorts the eigenvalues into descending order,
and rearranges the columns of v correspondingly. The method is straight insertion. */

void eigsort(double d[],double **v,int n){

int k,j,i;
double p;

for(i=0;i<n-1;i++){
	p=d[k=i];
	for(j=i+1;j<n;j++)
		if(d[j]>=p)
			p=d[k=j];
	if(k!=i){
		d[k]=d[i];
		d[i]=p;
		for(j=0;j<n;j++){
			p=v[j][i];
			v[j][i]=v[j][k];
			v[j][k]=p;
			}
		}
    }
}
