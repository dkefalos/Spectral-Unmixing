/*  SEE-E - Performs Simple Endmember Extraction Enchanced on hyperspectral Images
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
int PCA(FILE*, FILE*, int, int, int, int, FILE*, int, int);
int image_add(FILE*, FILE*, FILE*, int, int, int, FILE*, int);
int SEE_E(FILE*, FILE*, FILE*, FILE*, FILE*, int, int , int, int, double**, int, int);
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

FILE     *fin, *report, *report2, *report3, *fsignal, *fplus;
int      rows, cols, bands, nendmembers, i, j, k, m, nthreads, opt, nvalue, npixels;
char     ersfile[100], binfile[100], tiffile[100], sysstring[200], txtfile[100], txtfile2[100], txtfile3[100];
float    *buffer;
double   **a;

if(argc!=8){
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * SEE-E:   Perform SEE-E on Hyperspectral Images                        *");
    puts(" *                                                                       *");
    puts(" *         -   Input parameters (at command line) :                      *");
    puts(" *                 -        Image filename                               *");
    puts(" *                 -        Number of Rows                               *");
    puts(" *                 -        Number of Columns                            *");
    puts(" *                 -        Number of Bands                              *");
    puts(" *                 -        Number of Endmembers                         *");
    puts(" *                 -        Noise Estimation Method                      *");
    puts(" *                          1 for MNF with NND                           *");
    puts(" *                          2 for MNF with MRTBM                         *");
    puts(" *                          3 for PCA instead of MNF                     *");
    puts(" *                 -        Null Cell Value                              *");
    puts(" *                                                                       *");
    puts(" *         -   Works on the file types supported by GDAL                 *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * Output : - Report :                   Imagename_SEE-E_report.txt      *");
    puts(" *          - Endmembers File :          Imagename_endmembers.txt        *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * SEE-E       by Dimitris Kefalos  (dkefalos@gmail.com)                 *");
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
sprintf(txtfile, "%s_SEE-E_report.txt", binfile);
i=(strlen(binfile));
memcpy(txtfile2, binfile, strlen(binfile));
txtfile2[i]='\0';
sprintf(txtfile2, "%s_endmembers.txt", binfile);
memcpy(txtfile3, binfile, strlen(binfile));
txtfile3[i]='\0';
sprintf(txtfile3, "%s_SEE-E_endmembers_spectral.txt", binfile);

/* Argument Handling */

rows=atoi(argv[2]);
cols=atoi(argv[3]);
bands=atoi(argv[4]);
nendmembers=atoi(argv[5]);
opt=atoi(argv[6]);
nvalue=atoi(argv[7]);

if (opt!=1 && opt!=2 && opt!=3){
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

if (nendmembers<0){
    puts("Please give a valid number of endmembers");
    exit(-1);
}

/* GDAL Call */

sprintf(sysstring, "gdal_translate -ot Float32 -of ERS %s %s", tiffile, ersfile);
printf("\n");
printf(sysstring);
printf("\n");
if ((i=system(sysstring))!=0 && (i=system(sysstring))!=1 && (i=system(sysstring))!=256 ){
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

if ( (report2=fopen(txtfile2 , "wt")) == NULL){
    printf("\nCannot open %s file.\n", txtfile2);
    exit (-1);
    }

if ( (report3=fopen(txtfile3 , "wt")) == NULL){
    printf("\nCannot open %s file.\n", txtfile3);
    exit (-1);
    }

if ((fsignal=tmpfile()) == NULL){
    printf("\nCannot open signal file.\n");
    exit (-1);
    }

if ( (fplus=tmpfile()) == NULL){
    printf("\nCannot open PCA1 file.\n");
    exit (-1);
    }

if ( (a = malloc(bands*sizeof(double*))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

if ( (buffer = malloc(cols*bands*sizeof(float))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

for (i=0; i< bands; i++){
    if ( (a[i] = malloc(nendmembers*sizeof(double))) == NULL){
        puts("\nBad Memory Allocation.\n");
        return(-1);
        }
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

/* Transformation Option */

if (opt==1 || opt ==2){
    MNF(fin, fsignal, report, rows, cols, bands, opt, nvalue, npixels);
    }
else if (opt==3){
    puts("Executing PCA:");
    PCA(fin, report, rows, cols, bands, bands, fsignal, nvalue, npixels);
    }

/* Endmember Spectal Signature Extraction */

image_add(fin, fsignal, report, rows, cols, bands, fplus, nvalue);
SEE_E(fplus, fin, report, report2, report3, rows, cols, bands, nendmembers, a, nvalue, npixels);
puts("\nSEE-E Executed");

/* Closing */

for (i=0; i<bands; i++){
    free(a[i]);
}
free(buffer);
free(a);
fflush(NULL);
if (fclose(fplus)!=0){
   puts("Error in fplus closing");
   }

if (fclose(fsignal)!=0){
   puts("Error in fsignal closing");
   }

if (fclose(report2)!=0){
   puts("Error in report closing");
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

int image_add(FILE* fin, FILE* fsignal, FILE* report, int rows, int cols, int bands, FILE* fplus, int nvalue){

FILE    *fpca2;
int     maxrow, maxcol, i, j;
float   *buffer, *buffer2, *buffer3, *mpv, smax;

/* Memory Allocacion */

if ((fpca2=tmpfile()) == NULL){
    printf("\nCannot open pca2 file.\n");
    exit (-1);
    }

if ( (buffer = malloc(cols*sizeof(float))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

if ( (buffer2 = malloc(cols*bands*sizeof(float))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

if ( (buffer3 = malloc(3*cols*bands*sizeof(float))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

if ( (mpv = malloc(bands*sizeof(float))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

/* Maximum Projected Value */

smax=-500000.0;
rewind(fsignal);
for (i=0; i<rows; i++){
    fread(buffer, sizeof(float), cols, fsignal);
    for (j=0; j<cols; j++){
        if (buffer[j]!=nvalue){
            if (smax<buffer[j]){
                smax=buffer[j];
                maxrow=i;
                maxcol=j;
                }
            }
        }
    fseek(fsignal, cols*(bands-1)*sizeof(float), SEEK_CUR);
    }

fprintf(report, "\nsmax= %5.5lf maxrow= %d maxcol= %d", smax, maxrow, maxcol);

/* Maximum Projected Value Signature */

fseek(fin, maxrow*cols*bands*sizeof(float), SEEK_SET);
fread(buffer2, sizeof(float), cols*bands, fin);
for(i=0; i<bands; i++){
    mpv[i]=buffer2[i*cols+maxcol];
    }

/* Image Add Preparation */

for (i=0; i<bands; i++){
	for (j=0; j<cols*3; j++){
		buffer3[i*cols*3+j]=mpv[i];
		}
	}

/* Image Add Writing */

rewind(fin);
for (i=0 ; i<rows; i++){
    for (j=0; j<bands; j++){
        fread(buffer2, sizeof(float), cols, fin);
        fwrite(buffer2, sizeof(float), cols, fplus);
        fwrite(buffer3+j*cols*3, sizeof(float), cols*3, fplus);
        }
    }
puts("Image add Finished\n");

/* Closing Image Add */

free(mpv);
free(buffer3);
free(buffer2);
free(buffer);
if (fclose(fpca2)!=0){
   puts("Error in fpca2 closing");
   }
return (0);
}

int SEE_E(FILE* fplus, FILE* fin, FILE* report, FILE* report2, FILE* report3, int rows, int cols, int bands, int nendmembers, double **a, int nvalue, int npixels){

FILE     *fpca4;
int      *coordsr, *coordsc, *endm, c, i, j, k, l, *tempint;
float    *hmin, *hmax, *buffer, *buffer2, **ssignature, threshold, **dif;
double   *sum, **pro, **sqr, **sad, *ssad, ssadmin, ssadmax;

/* File Opening */

if ((fpca4=tmpfile()) == NULL){
    printf("\nCannot open pca4 file.\n");
    exit (-1);
    }
fflush(NULL);

/* Memory Allocation */

if ( (hmin = malloc((nendmembers-1)*sizeof(float))) == NULL){
    puts("\nBad hmin Memory Allocation.\n");
    return(-1);
    }

if ( (hmax = malloc((nendmembers-1)*sizeof(float))) == NULL){
    puts("\nBad hmax Memory Allocation.\n");
    return(-1);
    }

if ( (coordsr = malloc(2*(nendmembers-1)*sizeof(int))) == NULL){
    puts("\nBad coordsr Memory Allocation.\n");
    return(-1);
    }

if ( (coordsc = malloc(2*(nendmembers-1)*sizeof(int))) == NULL){
    puts("\nBad coordsc Memory Allocation.\n");
    return(-1);
    }

if ( (ssignature = malloc(2*(nendmembers-1)*sizeof(float*))) == NULL){
    puts("\nBad ssingnature Memory Allocation.\n");
    return(-1);
    }

if ( (sum = calloc(2*(nendmembers-1),sizeof(double))) == NULL){
    puts("\nBad sum Memory Allocation.\n");
    return(-1);
    }

if ( (pro = calloc(2*(nendmembers-1),sizeof(double*))) == NULL){
    puts("\nBad pro Memory Allocation.\n");
    return(-1);
    }

if ( (sqr = calloc(2*(nendmembers-1),sizeof(double*))) == NULL){
    puts("\nBad sqr Memory Allocation.\n");
    return(-1);
    }

if ( (buffer = malloc(cols*4*(nendmembers-1)*sizeof(float))) == NULL){
    puts("\nBad buffer Memory Allocation.\n");
    return(-1);
    }

if ( (buffer2 = malloc(cols*bands*sizeof(float))) == NULL){
    puts("\nBad buffer2 Memory Allocation.\n");
    return(-1);
    }

if ( (ssad = calloc(2*(nendmembers-1),sizeof(double))) == NULL){
    puts("\nBad ssad Memory Allocation.\n");
    return(-1);
    }

if ( (sad = calloc(2*(nendmembers-1),sizeof(double*))) == NULL){
    puts("\nBad sad Memory Allocation.\n");
    return(-1);
    }

if ( (dif = calloc(2*(nendmembers-1),sizeof(float*))) == NULL){
    puts("\nBad dif Memory Allocation.\n");
    return(-1);
    }

if ( (endm = calloc(2*(nendmembers-1),sizeof(int))) == NULL){
    puts("\nBad end Memory Allocation.\n");
    return(-1);
    }

if ( (tempint = calloc(2*(nendmembers-1),sizeof(int))) == NULL){
    puts("\nBad tempint Memory Allocation.\n");
    return(-1);
    }

for (i=0; i<(2*(nendmembers-1)); i++){
    if ( (ssignature[i] = malloc(bands*sizeof(float))) == NULL){
        puts("\nBad ssignature[i] Memory Allocation.\n");
        return(-1);
        }

    if ( (pro[i] = calloc(2*(nendmembers-1),sizeof(double))) == NULL){
        puts("\nBad pro[i] Memory Allocation.\n");
        return(-1);
        }

    if ( (sqr[i] = calloc(2*(nendmembers-1),sizeof(double))) == NULL){
        puts("\nBad sqr[i] Memory Allocation.\n");
        return(-1);
        }
     if ( (sad[i] = calloc(2*(nendmembers-1),sizeof(double))) == NULL){
        puts("\nBad sad[i] Memory Allocation.\n");
        return(-1);
        }
    if ( (dif[i] = calloc(2*(nendmembers-1),sizeof(float))) == NULL){
        puts("\nBad dif[i] Memory Allocation.\n");
        return(-1);
        }
    }

/* PCA on Enchanced Image */

fflush(NULL);
puts("PCA on Enchanced Image");
PCA(fplus, report, rows, cols*4, bands, nendmembers-1, fpca4, nvalue, npixels);
fflush(NULL);

/* Min-Max Calculation for Every Band */

for(i=0; i<(nendmembers-1); i++){
    hmin[i]=500000.;
    hmax[i]=-500000.;
    }

//fseek(fpca4, cols*(nendmembers-1)*4*sizeof(float), SEEK_SET);
rewind(fpca4);
rewind(fin);
for (i=0; i<rows; i++){
	fread(buffer, sizeof(float), cols*(nendmembers-1)*4, fpca4);
	fread(buffer2, sizeof(float), cols*bands, fin);
	for (j=0; j<(nendmembers-1); j++){
		for (k=0; k<cols; k++){
		    if (buffer2[j*cols*4+k]!=nvalue){
                if (buffer[j*cols*4+k]<hmin[j]){
                    hmin[j]=buffer[j*cols*4+k];
                    coordsr[j]=i;
                    coordsc[j]=k;
                    }
                if (buffer[j*cols*4+k]>hmax[j]){
                    hmax[j]=buffer[j*cols*4+k];
                    coordsr[(nendmembers-1)+j]=i;
                    coordsc[(nendmembers-1)+j]=k;
                    }
                }
			}
		}
	}

for (i=0; i<nendmembers-1; i++){
    fprintf(report, " %10.5lf  %10.5lf\n", hmin[i], hmax[i]);
    }

/* Endmember Canditate Signature Reading */

rewind(fin);
for (i=0; i<2*(nendmembers-1); i++){
    fseek(fin, (long)coordsr[i]*cols*bands*sizeof(float), SEEK_SET);
    fread(buffer2, sizeof(float), cols*bands, fin);
    for (j=0; j<bands; j++){
        ssignature[i][j]=buffer2[j*cols+coordsc[i]];
        }
    }

fprintf(report, "\nEndmember Candidates:");
for (i=0; i<(nendmembers-1); i++){
    fprintf(report, "\n%d: row=%d, col=%d :", i+1, coordsr[i]+1, coordsc[i]+1);
    for(j=0; j<bands; j++){
        fprintf(report, " %f", ssignature[i][j]);
        }
    }
for (i=(nendmembers-1); i<2*(nendmembers-1); i++){
    fprintf(report, "\n%d: row=%d, col=%d :", i+1, coordsr[i]+1, coordsc[i]+1);
    for(j=0; j<bands; j++){
        fprintf(report, " %f", ssignature[i][j]);
        }
    }

/* SAD Calculation for Every Endmember */

for (i=0; i<2*(nendmembers-1); i++){
    for (j=0; j<2*(nendmembers-1); j++){
        for (k=0; k<bands; k++){
            if (i!=j){
                pro[i][j]+=ssignature[i][k]*ssignature[j][k];
                sqr[i][j]+=ssignature[j][k]*ssignature[j][k];
                }
            }
        }
    for (l=0; l<bands; l++){
        sum[i]+=ssignature[i][l]*ssignature[i][l];
        }
    }

fprintf(report, "\npro :\n");
for (i=0; i<2*(nendmembers-1); i++){
    for (j=0; j<2*(nendmembers-1); j++){
        fprintf(report, " %lf", pro[i][j]);
        }
    fprintf(report, "\n");
    }

fprintf(report, "\nsqr :\n");
for (i=0; i<2*(nendmembers-1); i++){
    for (j=0; j<2*(nendmembers-1); j++){
        fprintf(report, " %lf", sqr[i][j]);
        }
    fprintf(report, "\n");
    }

fprintf(report, "\nsum :\n");
for (i=0; i<2*(nendmembers-1); i++){
    fprintf(report, "%lf\n", sum[i]);
    }

for (i=0; i<2*(nendmembers-1); i++){
    for (j=0; j<2*(nendmembers-1); j++){
        if (i!=j){
            sad[i][j]=acos((pro[i][j])/((sqrt(sum[i]))*(sqrt(sqr[i][j]))));
            }
        }
    }

for (i=0; i<2*(nendmembers-1); i++){
    if (sum[i]==0){
        for (j=0; j<2*(nendmembers-1); j++){
            sad[i][j]=0;
            sad[j][i]=0;
            }
        }
    }

fprintf(report, "\nSADs:\n");
for (i=0; i<2*(nendmembers-1); i++){
    fprintf(report, "%d: ", i+1);
    for (j=0; j<2*(nendmembers-1); j++){
        fprintf(report, " %10.5lf", sad[i][j]);
        }
    fprintf(report, "\n");
    }

/* SAD Adding */

for (i=0; i<2*(nendmembers-1); i++){
    for (j=0; j<2*(nendmembers-1); j++){
        if (i!=j){
            ssad[i]+=sad[i][j];
            }
        }
    }

fprintf(report, "\nsSADs:\n");
for (i=0; i<2*(nendmembers-1); i++){
    fprintf(report, "ssad[%d] = %15.5lf\n", i, ssad[i]);
    }

/* SAD Normalization */

ssadmax=0.;
ssadmin=5000.;
for (i=0; i<2*(nendmembers-1); i++){
    if (ssadmax<ssad[i]){
        ssadmax=ssad[i];
        }
    if (ssadmin>ssad[i]){
        ssadmin=ssad[i];
        }
    }

for (i=0; i<2*(nendmembers-1); i++){
    ssad[i]= (ssad[i]-ssadmin)/(ssadmax-ssadmin);
    }

fprintf(report, "\nNormalized sSADs:\n");
for (i=0; i<2*(nendmembers-1); i++){
    fprintf(report, "ssad[%d] = %15.5lf\n", i, ssad[i]);
    }

/* SAD Differencing */

for (i=0; i<2*(nendmembers-1); i++){
    for (j=i+1; j<2*(nendmembers-1); j++){
        dif[i][j]=fabs(ssad[i]-ssad[j]);
        }
    }

fprintf(report, " \nDifT:\n");
for (i=0; i<2*(nendmembers-1); i++){
    for (j=0; j<2*(nendmembers-1); j++){
        fprintf(report, "dif[%d][%d]= %10.5lf     ", j, i, dif[j][i]);
        }
    fprintf(report, " \n");
    }

/* Endmember Selection */

fflush(NULL);
threshold=0.0;
c=1000;
while (threshold<0.1 && c>nendmembers){
    c=0;
    threshold= threshold+0.00001;
    for (i=0; i<2*(nendmembers-1); i++){
        endm[i]=0;
        for (j=i+1; j<2*(nendmembers-1); j++){
            if (dif[i][j]<threshold){
                tempint[i]++;
                }
            }
        if (tempint[i]==0){
            c++;
            endm[i]=1;
            }
        }
    }

printf("Extracted %d Endmembers\n", c);
k=0;
fprintf(report, "\nEndmembers :\n");
for (i=0; i<2*(nendmembers-1); i++){
    if (endm[i]!=0){
        fprintf(report, "\nEndmember %d: at row %d, col %d:\t", i+1, coordsr[i]+1, coordsc[i]+1);
        for (j=0; j<bands; j++){
            fprintf(report, "\t%5.5lf", ssignature[i][j]);
            a[j][k]=ssignature[i][j];
            }
        k++;
        }
    }

fprintf(report2, "Endmembers :\n");
for (i=0; i<nendmembers; i++){
    for (j=0; j<bands; j++){
        fprintf(report2, "%f ", a[j][i]);
        }
    fprintf(report2, "\n");
    }

for (i=0; i<bands; i++){
    for (j=0; j<nendmembers; j++){
        fprintf(report3, "\t%f" ,(a[i][j]));
        }
    fprintf(report3, "\n");
}

/* SEE Closing */

if (fclose(fpca4)!=0){
   puts("Error in fpca4 closing");
   }

for (i=0; i<(2*(nendmembers-1)); i++){
    free(dif[i]);
    free(sad[i]);
    free(sqr[i]);
    free(pro[i]);
    free(ssignature[i]);
    }
free(tempint);
free(endm);
free(dif);
free(sad);
free(ssad);
free(buffer2);
free(buffer);
free(sqr);
free(pro);
free(sum);
free(ssignature);
free(coordsc);
free(coordsr);
free(hmax);
free(hmin);

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

int PCA(FILE* fin, FILE* report, int rows, int cols, int bands, int noc, FILE * fpca, int nvalue, int npixels){

int              i, j, k, l, nrot;
double           *sum, **sumv, **cov, *mean, *eigval, **eigvec, tstart, tend;
float            *pc, *buffer;

if (noc==-1){
    noc=bands;
    }

/* Memory Allocation */

if ( (buffer=malloc(bands*cols*sizeof(float))) == NULL ){
    puts("Error in (buffer) Memory Allocation.");
    exit(-1);
    }

if ( (sum=calloc(bands,sizeof(double))) == NULL ){
    puts("Error in (sum) Memory Allocation.");
    exit(-1);
    }

if ( (sumv=calloc(bands,sizeof(double*))) == NULL ){
    puts("Error in (sumv) Memory Allocation.");
    exit(-1);
    }

if ( (mean=calloc(bands,sizeof(double))) == NULL ){
    puts("Error in (mean) Memory Allocation.");
    exit(-1);
    }

if ( (cov=malloc(bands*sizeof(double*))) == NULL ){
    puts("Error in (cov) Memory Allocation.");
    exit(-1);
    }

if ( (eigval=malloc(bands*sizeof(double))) == NULL ){
    puts("Error in (eigval) Memory Allocation.");
    exit(-1);
    }

if ( (eigvec=malloc(bands*sizeof(double*))) == NULL ){
    puts("Error in (eigvec) Memory Allocation.");
    exit(-1);
    }

if ( (pc=calloc(noc*cols,sizeof(float))) == NULL ){
    puts("Error in (pc) Memory Allocation.");
    exit(-1);
    }

for (i=0; i<bands; i++){
    if ((sumv[i]=calloc(bands,sizeof(double))) == NULL){
        puts("Error in (sumv[i]) Memory Allocation.");
        exit(-1);
        }
    if ((cov[i]=malloc(bands*sizeof(double))) == NULL){
        puts("Error in (cov[i]) Memory Allocation.");
        exit(-1);
        }
    if ( (eigvec[i]=malloc(bands*sizeof(double))) == NULL ){
        puts("Error in (eigvec[i]) Memory Allocation.");
        exit(-1);
        }
    }

rewind(fin);

/* Sum Calculation for Every Band*/

for (i=0; i<rows; i++){
    fread (buffer, sizeof(float), cols*bands, fin);
    for (j=0; j<bands; j++){
        for (k=0; k<cols; k++){
            if (buffer[j*cols+k]!=nvalue){
                sum[j]+=buffer[(j*cols)+k];
                }
            }
        }
    }

fprintf(report,"\nPCA :\n");
fprintf(report,"\nSum :\n");
for (i=0; i<bands; i++){
    fprintf(report,"Band %d : %10.5lf\n", i+1, sum[i]);
    }

/* Mean Calculation for Every Band */

for (i=0; i<bands; i++){
    mean[i]=sum[i]/(rows*cols-npixels);
    }

fprintf(report,"\nMean :\n");
for (i=0; i<bands; i++){
    fprintf(report,"Band %d : %10.5lf\n", i+1, mean[i]);
    }

rewind(fin);
/* sumv Calculation for Every Band*/

tstart=omp_get_wtime();
# pragma omp parallel default(none) private(i, j, k, l) shared(buffer, fin, mean, sumv, rows, cols, bands, nvalue)
{
for (i=0; i<rows; i++){
    # pragma omp single
    {
    fread (buffer, sizeof(float), cols*bands, fin);
    printf("Calculating Statistics at row %d of %d\r", i+1, rows);
    }
    # pragma omp for schedule(dynamic)
    for (j=0; j<bands; j++){
        for (k=0; k<bands; k++){
            for (l=0; l<cols; l++){
                if (buffer[j*cols+l]!=nvalue){
                    sumv[j][k]+=(buffer[(j*cols)+l]-mean[j])*(buffer[(k*cols)+l]-mean[k]);
                    }
                }
            }
        }
    }
}

tend=omp_get_wtime();
printf("\nDuration of Statistics Parallel region : %2.1lf sec\n", tend-tstart);

/* Variance-Covariance Matrix Calculation */

for (i=0; i<bands; i++){
    for (j=0; j<bands; j++){
        cov[i][j]=sumv[i][j]/(rows*cols-npixels);
        }
    }

fprintf(report,"\nVariance-Covariance Matrix:\n");
for (i=0; i<bands; i++){
    fprintf(report, "Band %d", i+1);
    for (j=0; j<bands; j++){
        fprintf(report," %15.5lf", cov[i][j]);
        }
    fprintf(report,"\n");
    }
fflush(NULL);
/* Eigenvalues and Eigenvectors Calculation and Sorting */

jacobi(cov, bands, eigval, eigvec, &nrot);
eigsort(eigval, eigvec, bands);

//fprintf(report,"\nNumber of Jacobi rotations : %d\n", nrot);

fprintf(report,"\nEigenvalues :\n");
for (i=0; i<bands; i++){
    fprintf(report," %10.5lf\n", eigval[i]);
    }

fprintf(report,"\nEigenvectors :\n");
for (i=0; i<bands; i++){
    for (j=0; j<bands; j++){
        fprintf(report," %10.5lf", eigvec[i][j]);
        }
    fprintf(report,"\n");
    }


/* PC Calculation */

rewind(fin);
tstart=omp_get_wtime();
# pragma omp parallel default (none) private(i, j, k, l) shared(fin, buffer, rows, cols, bands, noc, pc, eigvec, fpca, mean, nvalue)
{
for (i=0; i<rows; i++){
    # pragma omp single
    {
    fread(buffer, sizeof(float), cols*bands, fin);
    printf("Calculating Pr. Components at row %d of %d\r", i+1, rows);
    }
    # pragma omp for schedule (dynamic)
    for (j=0; j<noc; j++){
        for (k=0; k<cols; k++){
            pc[j*cols+k]=0.0;
            for(l=0; l<bands; l++){
                if (buffer[l*cols+k]!=nvalue){
                    pc[j*cols+k]+=1.*(buffer[(l*cols)+k]-mean[j])*(eigvec[l][j]);
                    }
                else {
                    pc[j*cols+k]=nvalue;
                    }
                }
            }
        }
    # pragma omp single
    {
    fwrite(pc, sizeof(float), noc*cols, fpca);
    }
    }
}

tend=omp_get_wtime();
printf("\nDuration of PCA Parallel region : %2.1lf sec\n", tend-tstart);

/* Closing */

fflush(NULL);
for (i=0; i<bands; i++){
    free(sumv[i]);
    free(cov[i]);
    free(eigvec[i]);
    }
free(pc);
free(eigvec);
free(eigval);
free(cov);
free(mean);
free(sumv);
free(sum);
free(buffer);

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
