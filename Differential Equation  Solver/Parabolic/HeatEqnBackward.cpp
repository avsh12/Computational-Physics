#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_permutation.h>

using namespace std;

//intitial profile
double phi(double x)
   {if(x>=0 && x<0.5)
       return x;
    else if(x>=0.5 && x<=1)
       return (1-x);
    else return 0;
    }

int main()
   {double xrange=1, trange=0.1, dx=0.01, s=0.51;
    double dt=s*dx*dx, a=2*s+1;

    int N = (int)ceil(xrange/dx)-1, Nt = (int)ceil(trange/dt);
    int signum;

    gsl_matrix *A = gsl_matrix_alloc(N,N);
    gsl_matrix *B = gsl_matrix_calloc(N,N);
    gsl_permutation *p = gsl_permutation_alloc(N);
    gsl_matrix *Xn1 = gsl_matrix_calloc(N,1);
    gsl_matrix *Xn = gsl_matrix_calloc(N,1);

    gsl_matrix_set(B, 0, 0, a);
    //set matrix B to be inverted
    for(int i=1; i<N; ++i)                 
       {gsl_matrix_set(B, i, i, a);
        gsl_matrix_set(B, i, i-1, -s);
        gsl_matrix_set(B, i-1, i, -s);
        }

    //calculate LU Decomposition
    gsl_linalg_LU_decomp(B, p, &signum);  
    //Compute inverse of B in A 
    gsl_linalg_LU_invert(B, p, A);         
    gsl_matrix_free(B);
    gsl_permutation_free(p);

    gsl_matrix *u = gsl_matrix_alloc(Nt,N+2);

    //set boundary values in matrix u
    for(int i=0; i<N+2; ++i)                    
       {gsl_matrix_set(u, 0, i, phi(i*dx));
        if(i>0 && i<N+1)
           //set initial values in Xn, at t=0
           gsl_matrix_set(Xn, i-1, 0, phi(i*dx));     
        }

    //set boundary values in matrix u
    for(int n=0; n<Nt; ++n)                
       {gsl_matrix_set(u, n, 0, 0);
        gsl_matrix_set(u, n, N+1, 0);
        }

    //compute Xn at succesive times
    for(int n=1; n<Nt; ++n)               
       {gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, Xn, 0.0, Xn1);
        gsl_matrix_memcpy(Xn, Xn1);
        for(int i=1; i<N; ++i)
            gsl_matrix_set(u, n, i, gsl_matrix_get(Xn, i-1, 0));
        }

    fstream fout;
    fout.open("HeatEqnBackwardData.txt", ios::out|ios::trunc);

    for(int i=0; i<N+2; ++i)
       {for(int n=0; n<Nt; n+=10)
           {fout<<fixed<<setprecision(2)<<i*dx<<"  ";
            fout<<fixed<<setprecision(5)<<n*dt<<"  "<<gsl_matrix_get(u, n, i)<<endl;
            }
        fout<<endl;
        }
    fout.close();

    gsl_matrix_free(A);
    gsl_matrix_free(Xn);
    gsl_matrix_free(Xn1);

    return 0;
    }
