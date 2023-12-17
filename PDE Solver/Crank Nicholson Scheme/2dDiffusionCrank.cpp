//Give initial profile at the origin at line 127 and 128
//To include source, go to 156
//Change time range at 93
//Change the number of cells in the grid at 90
//Vary value of s at 87
//Change which data to write at 154

#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_permutation.h>

using namespace std;

void matrixInitializer(gsl_matrix *Mat, int nx, double a, double b)
    {for(int i=0; i<nx; ++i)
       {for(int j=0; j<=i; ++j)
           {double row = i*(i+1)/2 + j;
            //(i,j)
            gsl_matrix_set(Mat, row, row, gsl_matrix_get(Mat, row, row)+2*b);
        
            //(i-1, j)
            if(i==0)
              {gsl_matrix_set(Mat, row, (i+1)*(i+2)/2 + j, gsl_matrix_get(Mat, row, (i+1)*(i+2)/2 + j)+a);
       	       }
      	    else if(i==j)
          	  {gsl_matrix_set(Mat, row, j*(j+1)/2 + i-1, gsl_matrix_get(Mat, row, j*(j+1)/2 + i-1)+a);
           	   }
        	else
          	  {gsl_matrix_set(Mat, row, i*(i-1)/2 + j, gsl_matrix_get(Mat, row, i*(i-1)/2 + j)+a);
               }

        	//(i+1, j)
        	if(i==nx-1)
          	  {gsl_matrix_set(Mat, row, i*(i+1)/2 + j, gsl_matrix_get(Mat, row, i*(i+1)/2 + j)+a);
           	   }
        	else
           	  {gsl_matrix_set(Mat, row, (i+1)*(i+2)/2 + j, gsl_matrix_get(Mat, row, (i+1)*(i+2)/2 + j)+a);
           	   }

        	//(i, j-1)
        	if(j==0)
          	  {if(i==0)
             	 {gsl_matrix_set(Mat, row, (j+1)*(j+2)/2 + i, gsl_matrix_get(Mat, row, (j+1)*(j+2)/2 + i)+a);
              	  }
           	   else
             	 {gsl_matrix_set(Mat, row, i*(i+1)/2 + j+1, gsl_matrix_get(Mat, row, i*(i+1)/2 + j+1)+a);
              	  }
           	   }
        	else
         	  {gsl_matrix_set(Mat, row, i*(i+1)/2 + j-1, gsl_matrix_get(Mat, row, i*(i+1)/2 + j-1)+a);
           	   }

       	    //(i, j+1)
        	if(i==j)
          	  {if(i==nx-1)
             	 {gsl_matrix_set(Mat, row, i*(i+1)/2 + j, gsl_matrix_get(Mat, row, i*(i+1)/2 + j)+a);
              	  }
           	   else
             	 {gsl_matrix_set(Mat, row, (j+1)*(j+2)/2 + i, gsl_matrix_get(Mat, row, (j+1)*(j+2)/2 + i)+a);
              	  }
           	   }
        	else
          	  {gsl_matrix_set(Mat, row, i*(i+1)/2 + j+1, (gsl_matrix_get(Mat, row, i*(i+1)/2 + j+1) + a));
           	   }
                }
       }
     }

double source(double x, double y)
    {return (x<1.0e-3 && y<1.0e-3) ? 1:0;}

void fileWrite(fstream &fout, gsl_matrix *Mat, double nx, double dx, int &count)
    {for(int i=0; i<nx; ++i)
        {for(int j=0; j<nx; ++j)
            {fout<<i*dx<<"  "<<j*dx<<"  "<<gsl_matrix_get(Mat, i, j)<<endl;
             ++count;
             }
         fout<<endl;
         }
     }

int main()
   {double s=0.1, Df =1;
    double a = s*Df;
    double b = -(1+2*a), c = 2*a-1;
    double x_min=0, x_max=1;
    int nx=40;
    int row = nx*(nx+1)/2;
    double dx = (x_max-x_min)/nx;
    double dt = s*dx*dx;         //for Source
    int t_range = 1300;
    int signum;

    gsl_matrix *A = gsl_matrix_calloc(row, row);
    gsl_matrix *B = gsl_matrix_calloc(row, row);
    gsl_matrix *C = gsl_matrix_calloc(row, row);
    gsl_matrix *D = gsl_matrix_calloc(row, 1);          //for Source
    gsl_matrix *S = gsl_matrix_calloc(row, 1);          //for Source
    gsl_matrix *Fn = gsl_matrix_calloc(row, 1);
    gsl_matrix *Fn1 = gsl_matrix_calloc(row, 1);
    gsl_permutation *p = gsl_permutation_alloc(row);
    
    matrixInitializer(A, nx, a, b);
    matrixInitializer(B, nx, -a, c);
    
    //for Source
    for(int i=0; i<nx; ++i)
       {for(int j=0; j<=i; ++j)
           {gsl_matrix_set(D, (i*(i+1)/2+j), 0, -2*dt*source(i*dx, j*dx));
            }
        }
    
    gsl_linalg_LU_decomp(A, p, &signum);    //compute LU decomposition in A
    gsl_linalg_LU_invert(A, p, C);             //compute inverse in C
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C, B, 0.0, A);   //compute CB and store in A
    //for Source
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, C, D, 0.0, S);   //compute CD and store in S
    
    gsl_permutation_free(p);
    gsl_matrix_free(B);
    gsl_matrix_free(C);
    //for Source
    gsl_matrix_free(D);

    gsl_matrix *u = gsl_matrix_alloc(nx, nx);
    //set initial profile at the origin
    gsl_matrix_set(Fn, 0,0, 1);
    gsl_matrix_set(Fn1,0,0, 1);
 
    fstream fout;
    fout.open("2dDiffusionCrankData.txt", ios::out|ios::trunc);
    int count=0;
    //fout<<endl<<endl;

    for(int n=0; n<=t_range; ++n)
       {for(int i=0; i<nx; ++i)
           {for(int j=0; j<=i; ++j)
               {double val = gsl_matrix_get(Fn1, i*(i+1)/2+j, 0);
                if(i==j)
                  {gsl_matrix_set(u, i, j, val);
                   }
                else
                  {gsl_matrix_set(u, i, j, val);
                   gsl_matrix_set(u, j, i, val);
                   }
                }
            }

        if(n==t_range)// || n==1000||n==500||n==100||n==50||n==0)
          {fileWrite(fout, u, nx, dx, count);
           fout<<endl;
           }

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, Fn, 0.0, Fn1);
        //uncomment to add source
        //gsl_matrix_add(Fn1, S);
        gsl_matrix_memcpy(Fn,Fn1);
        }
    cout<<count;
    return 0;
    }
