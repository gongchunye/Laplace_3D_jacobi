#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#define  Max(a,b) ((a)>(b)?(a):(b))

/*
u(x,y,z) = x^3 + y^2 + z^2
f(x,y,z) = 6x + 4
x,y,z in [0,1]
*/

#define  N   (2*2*2*2*2*2+2)
double   maxeps = 0.1e-3;
int itmax = 5000;
int i,j,k;
double eps;
double A[N][N][N],  B[N][N][N],  F[N][N][N];
double dxyz = 1.0/(N-1), dx, dy, dz ;

void relax();
void resid();
void init();
void verify(); 
void wtime(double *t)
{
    *t = omp_get_wtime();
}


int main(int an, char **as)
{
	double time_start, time_fin;
	int it;
    init();
    wtime(&time_start);
	for(it=1; it<=itmax; it++)
	{
		eps = 0.;
		relax();
		resid();
		if(it%100==0) printf( "it=%4i   error=%f\n", it,eps);
		if (eps < maxeps) break;
	}
	wtime(&time_fin);
    printf("Time: %gs\t", time_fin - time_start);
	verify();
	return 0;
}

void init()
{
	for(i=0; i<=N-1; i++)
	for(j=0; j<=N-1; j++)
	for(k=0; k<=N-1; k++)
	{	
		dx = i*dxyz;
		dy = j*dxyz;
		dz = k*dxyz;
		A[i][j][k]=dx*dx*dx+dy*dy+dz*dz;
		F[i][j][k]= 6*dx + 4;
	}
	
	for(i=1; i<=N-2; i++)
	for(j=1; j<=N-2; j++)
	for(k=1; k<=N-2; k++)
	{
		A[i][j][k]=1;
		B[i][j][k]=1;
	}
	printf("init ok!\n");
} 

// jacobi iteration
void relax()
{
	for(i=1; i<=N-2; i++)
	for(j=1; j<=N-2; j++)
	for(k=1; k<=N-2; k++)
	{
		B[i][j][k]=(A[i-1][j][k]+A[i+1][j][k]+A[i][j-1][k]+A[i][j+1][k]+A[i][j][k-1]+A[i][j][k+1] -F[i][j][k]*dxyz*dxyz)/6.;
	}
}

void resid()
{ 
	for(i=1; i<=N-2; i++)
	for(j=1; j<=N-2; j++)
	for(k=1; k<=N-2; k++)
	{
		double e;
		e = fabs(A[i][j][k] - B[i][j][k]);         
		A[i][j][k] = B[i][j][k]; 
		eps = Max(eps,e);
	}
}

//
void verify()
{
	double err=0;
	for(i=0; i<=N-1; i++)
	for(j=0; j<=N-1; j++)
	for(k=0; k<=N-1; k++)
	{
		dx = dxyz*i;
		dy = dxyz*j;
		dz = dxyz*k;
		err = fmax(err, abs(dx*dx*dx+dy*dy+dz*dz - A[i][j][k]));
	}
	printf(" global max error = %f\n", err);
}
