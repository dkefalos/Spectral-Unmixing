/*  PCA - Performs Principal Components Transformation on Images
    Copyright (C) 2014  Kefalos Dimitris

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
int PCA(FILE*, FILE*, int, int, int, int, FILE*, int, int);

int compare (const void * x, const void * y){

double da = *(const double *)x;
double db = *(const double *)y;

if (da == db)
return 0;

return (da > db) ? -1 : 1;
}

int main(int argc, char *argv[]) {

FILE     *fin, *report, *fpca;
int      rows, cols, bands, i, j, k, m, nthreads, noc, nvalue, npixels;
char     ersfile[100], binfile[100], tiffile[100], sysstring[200], txtfile[100], pcadir[100];
float    *buffer;

if(argc!=7){
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * PCA:     Performs Principal Components Transformation on Images       *");
    puts(" *                                                                       *");
    puts(" *         -   Input parameters (at command line) :                      *");
    puts(" *                 -        Image filename                               *");
    puts(" *                 -        Number of Rows                               *");
    puts(" *                 -        Number of Columns                            *");
    puts(" *                 -        Number of Bands                              *");
    puts(" *                 -        Number of Components                         *");
    puts(" *                 -        Null Cell Value                              *");
    puts(" *                                                                       *");
    puts(" *         -   Works on the file types supported by GDAL                 *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * Output : - Report :    Imagename_PCA_report.txt                       *");
    puts(" *          - PCA Image : Imagename_PCA                                  *");
    puts(" *                                                                       *");
    puts(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
    puts(" *                                                                       *");
    puts(" * PCA      by Dimitris Kefalos  (dkefalos@gmail.com)                    *");
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
sprintf(txtfile, "%s_PCA_report.txt", binfile);
i=(strlen(binfile));
memcpy(pcadir, binfile, strlen(binfile));
pcadir[i]='\0';
sprintf(pcadir, "%s_PCA", binfile);

/* Argument Handling */

rows=atoi(argv[2]);
cols=atoi(argv[3]);
bands=atoi(argv[4]);
noc=atoi(argv[5]);
nvalue=atoi(argv[6]);

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

if (noc<0 || noc>bands){
    puts("Please give a valid number of components");
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

if ((fpca=fopen(pcadir , "wb+")) == NULL){
    printf("\nCannot open pca file.\n");
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


/* Run PCA */

PCA(fin, report, rows, cols, bands, noc, fpca, nvalue, npixels);
puts("PCA Executed\n");

/* Closing */

free(buffer);
fflush(NULL);
if (fclose(fpca)!=0){
   puts("Error in fpca4 closing");
   }

if (fclose(report)!=0){
   puts("Error in report closing");
   }

if (fclose(fin)!=0){
   puts("Error in fin closing");
   }

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
