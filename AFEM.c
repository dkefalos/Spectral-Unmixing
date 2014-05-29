/*  AFEM - Performs Abutance Fraction Estimation of given Endmembers of hyperspectral Images
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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Endmember Text File must follow the format:

    Endmembers :
    54 58 56 35......
    54 54 58 48......

    */

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

int readtxt(FILE*, FILE*, int, int, double**);
int Unmixing(FILE*, FILE*, int, int, int, int , char*, double**, int);
int Matrix_inversion(FP_TYPE*, FP_TYPE*, int);
int Cholesky_decomposition(FP_TYPE *, FP_TYPE *, int );

int compare (const void * x, const void * y){

double da = *(const double *)x;
double db = *(const double *)y;

if (da == db)
return 0;

return (da > db) ? -1 : 1;
}

int main(int argc, char *argv[]) {

FILE     *fin, *report, *ftxt;
int      rows, cols, bands, nendmembers, i, j, k, m, nthreads, nvalue, npixels;
char     ersfile[100], binfile[100], tiffile[100], sysstring[200], txtfile[100], enddir[100];
float    *buffer;
double   **a;

if(argc!=8){
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * AFEM:   Abuntance Fraction Estimation of the Endmembers of a Image    *");
    puts(" *                                                                       *");
    puts(" *         -   Input parameters (at command line) :                      *");
    puts(" *                 -        Image filename                               *");
    puts(" *                 -        Endmember Text file Path                     *");
    puts(" *                 -        Number of Rows                               *");
    puts(" *                 -        Number of Columns                            *");
    puts(" *                 -        Number of Bands                              *");
    puts(" *                 -        Number of Endmembers                         *");
    puts(" *                 -        Null Cell Value                              *");
    puts(" *                                                                       *");
    puts(" *         -   Works on the file types supported by GDAL                 *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * Output : - Report :                   Imagename_AEM_report.txt        *");
    puts(" *          - Abundance Fraction Image : Imagename_abuntance_fraction    *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * AFEM         by Dimitris Kefalos  (dkefalos@gmail.com)                *");
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
sprintf(txtfile, "%s_AEM_report.txt", binfile);
i=(strlen(binfile));
memcpy(enddir, binfile, strlen(binfile));
enddir[i]='\0';
sprintf(enddir, "%s_abundance_fraction", binfile);

/* Argument Handling */

rows=atoi(argv[3]);
cols=atoi(argv[4]);
bands=atoi(argv[5]);
nendmembers=atoi(argv[6]);
nvalue=atoi(argv[7]);

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

if ((ftxt=fopen(argv[2], "rt")) == NULL){
    printf("\nCannot open text file.\n");
    exit (-1);
    }

if ( (a = malloc(bands*sizeof(double*))) == NULL){
    puts("\nBad Memory Allocation.\n");
    return(-1);
    }

if ( (buffer = malloc(bands*cols*sizeof(float))) == NULL){
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

/* Read Endmember Text File */

readtxt(ftxt, report, bands, nendmembers, a);

/* Endmember Abundance Estimation */

Unmixing(fin, report, rows, cols, bands, nendmembers, enddir, a, nvalue);
puts("\nUnmixing Executed");

/* Closing */

for (i=0; i<bands; i++){
    free(a[i]);
}
free(buffer);
free(a);
fflush(NULL);
if (fclose(ftxt)!=0){
   puts("Error in txt closing");
   }

if (fclose(report)!=0){
   puts("Error in report closing");
   }

if (fclose(fin)!=0){
   puts("Error in fin closing");
   }

return(0);
}

int readtxt(FILE* txtfile, FILE* report, int bands, int nendmembers, double** a){

int      i, j;
double   tmp;
char     ignore[256];

fgets(ignore, sizeof(ignore), txtfile);
for (i=0; i<nendmembers; i++){
    for(j=0; j<bands; j++){
        if (fscanf(txtfile, "%lf ", &tmp)!=1){
            puts("not found");
            }
        a[j][i]=tmp;
        }
    fscanf(txtfile, "\n");
    }

for (i=0; i<nendmembers; i++){
    for(j=0; j<bands; j++){
        fprintf(report, "%f ", a[j][i]);
        }
    fprintf(report, "\n");
    }

return (0);
}

int Unmixing(FILE* fin, FILE* report, int rows, int cols, int bands, int nendmembers, char* enddir, double **a, int nvalue){

FILE     *fend;
int      i, j, k, l, m;
float    *buffer, *buffer2, *b,   c;
double   **N,**ata, **ata_inv, tmpd, *atal, *ata_invl;

/* Memory Allocation */

if ( (buffer=malloc(bands*cols*sizeof(float))) == NULL ){
    puts("Error in (buffer) Memory Allocation.");
    exit(-1);
    }

if ( (buffer2=malloc(nendmembers*cols*sizeof(float))) == NULL ){
    puts("Error in (buffer) Memory Allocation.");
    exit(-1);
    }

if ( (b=malloc(bands*sizeof(float))) == NULL ){
    puts("Error in (buffer) Memory Allocation.");
    exit(-1);
    }

if ( (ata=calloc(nendmembers,sizeof(double*))) == NULL ){
    puts("Error in (buffer) Memory Allocation.");
    exit(-1);
    }

if ( (ata_inv=calloc(nendmembers,sizeof(double*))) == NULL ){
    puts("Error in (buffer) Memory Allocation.");
    exit(-1);
    }

if ( (atal=calloc(nendmembers*nendmembers,sizeof(double))) == NULL ){
    puts("Error in (buffer) Memory Allocation.");
    exit(-1);
    }

if ( (ata_invl=calloc(nendmembers*nendmembers,sizeof(double))) == NULL ){
    puts("Error in (buffer) Memory Allocation.");
    exit(-1);
    }

if ( (N=calloc(nendmembers,sizeof(double*))) == NULL ){
    puts("Error in (buffer) Memory Allocation.");
    exit(-1);
    }

for (i=0; i<nendmembers; i++){
    if ( (ata[i]=calloc(nendmembers,sizeof(double))) == NULL ){
        puts("Error in (buffer) Memory Allocation.");
        exit(-1);
        }

    if ( (ata_inv[i]=calloc(nendmembers,sizeof(double))) == NULL ){
        puts("Error in (buffer) Memory Allocation.");
        exit(-1);
        }

    if ( (N[i]=calloc(bands,sizeof(double))) == NULL ){
        puts("Error in (buffer) Memory Allocation.");
        exit(-1);
        }
    }

fprintf(report, "\n\n %s\n", enddir);

if ( (fend=fopen(enddir,"wb+")) == NULL){
    printf("\nCannot create endmembers file.\n");
    exit (-1);
    }

for (i=0; i<bands; i++){
    for (j=0; j<nendmembers; j++){
        fprintf(report, " %f" ,(a[i][j]));
        }
    fprintf(report, "\n");
    }

/* aTa Calculation */

for (i=0; i<nendmembers ; i++){
    for (j=0; j<nendmembers ; j++){
        for (k=0; k<bands; k++){
            ata[i][j]+=a[k][i]*a[k][j];
            }
        }
    }

fprintf(report, "\nata\n");
for (i=0; i<nendmembers; i++){
    for (j=0; j<nendmembers; j++){
        fprintf(report, " %f" , ata[i][j]);
        }
    fprintf(report, "\n");
    }

/* aTa Inversion */

for (i=0; i< nendmembers; i++){
    for (j=0; j<nendmembers; j++){
        atal[i*nendmembers+j]=ata[i][j];
        }
    }

Matrix_inversion(atal,ata_invl,nendmembers);

for (i=0; i< nendmembers; i++){
    for (j=0; j<nendmembers; j++){
        ata_inv[i][j]=ata_invl[i*nendmembers+j];
        }
    }

fprintf(report, "\nata_inv\n");
for (i=0; i<nendmembers; i++){
    for (j=0; j<nendmembers; j++){
        fprintf(report, " %f" , ata_inv[i][j]);
        }
    fprintf(report, "\n");
    }
fflush(NULL);

fprintf(report, "\n");
fprintf(report, "\n");

/* N Matrix Calculation */

for (i=0; i<nendmembers; i++){
    for (j=0; j<bands; j++){
        tmpd=0.0;
        for (k=0; k<nendmembers; k++){
            tmpd+=ata_inv[i][k]*a[j][k];
            }
        N[i][j]=tmpd;
        }
    }

fprintf(report, "\nN :\n");
for(i=0; i<nendmembers; i++){
    for(j=0; j<bands; j++){
        fprintf(report, " %lf", N[i][j]);
        }
    fprintf(report, "\n");
    }

/* Abundance Fraction Calculation */

fprintf(report, " \nimage: \n");
rewind(fin);
for (i=0; i<rows; i++){
    fread(buffer, sizeof(float), cols*bands, fin);
    printf("Unmixing at row %d of %d\r", i+1, rows);
    for (j=0; j<cols; j++){
        for (k=0; k<bands; k++){
            b[k]=buffer[k*cols+j];
            }
        for (l=0; l<nendmembers; l++){
            c=0.0;
            if (b[0]==nvalue){
                buffer2[l*cols+j]=nvalue;
                }
            else {
                for (m=0; m<bands; m++){
                    c+=N[l][m]*b[m];
                    }
                buffer2[l*cols+j]=c;
                }
            }
        }
    fwrite(buffer2, sizeof(float), nendmembers*cols, fend);
    for (j=0; j<cols*nendmembers; j=j+30){
        fprintf(report, " %f", buffer2[j]);
        }
    }

/* Unmixing Closing */

for (i=0; i<nendmembers; i++){
    free(N[i]);
    free(ata_inv[i]);
    free(ata[i]);
    }
free(N);
free(ata_invl);
free(atal);
free(ata_inv);
free(ata);
free(b);
free(buffer2);
free(buffer);

return (0);
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
