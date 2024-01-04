#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>

using namespace std;

double source(double x, double y)
    {return (x<1.0e-2 && y<1.0e-2) ? 1:0;}

void fileWrite(fstream &fout, double *Mat, double nx, double dx, int &count)
    {for(int i=0; i<nx; ++i)
        {for(int j=0; j<nx; ++j)
            {if(j<=i)
                fout<<i*dx<<"  "<<j*dx<<"  "<<Mat[i*(i+1)/2 + j]<<endl;
             else
                fout<<i*dx<<"  "<<j*dx<<"  "<<Mat[j*(j+1)/2 + i]<<endl;
             ++count;
             }
         fout<<endl;
         }
     }

void memcpy(double *A, double *B, int row)
    {for(int i=0; i<row; ++i)
        {A[i] = B[i];
         }
     }

int main()
   {double s=0.35, D=1;
    double x_min=0, x_max=1;
    int nx=40;
    int row = nx*(nx+1)/2;
    double dx = (x_max-x_min)/nx;
    double dt = s*dx*dx;         //for Source
    int t_range = 143;

    fstream fout;
    fout.open("2dDiffusionForwardData.txt", ios::out|ios::trunc);
    int count=0;
  
    double *Fn, *Fn1, *copy;
    Fn = new double[row];
    Fn1 = new double[row];
    
    //initial profile
    for(int i=0; i<row; ++i)
       {Fn[i]=0;
        Fn1[i]=0;
        }

    for(int n=1; n<=t_range; ++n)
       {if(n==t_range)
          {fileWrite(fout, Fn1, nx, dx, count);
           fout<<endl;
           }
        for(int i=0; i<nx; ++i)
           {for(int j=0; j<=i; ++j)
               {//(i,j)
                Fn1[i*(i+1)/2 + j] = (1-4*s*D)* Fn[i*(i+1)/2 + j];
            
                //(i-1, j)
                if(i==0)
                  {Fn1[i*(i+1)/2 + j] += s*D*Fn[(i+1)*(i+2)/2 + j];
       	           }
      	        else if(i==j)
                  {Fn1[i*(i+1)/2 + j] += s*D*Fn[j*(j+1)/2 + i-1];
                   }
                else
                  {Fn1[i*(i+1)/2 + j] += s*D*Fn[i*(i-1)/2 + j];
                   }
 
                //(i+1, j)
                if(i==nx-1)
                  {Fn1[i*(i+1)/2 + j] += s*D*Fn[i*(i+1)/2 + j];
                   }
                else
                  {Fn1[i*(i+1)/2 + j] += s*D*Fn[(i+1)*(i+2)/2 + j];
                   }

                //(i, j-1)
                if(j==0)
                  {if(i==0)
             	     {Fn1[i*(i+1)/2 + j] += s*D*Fn[(j+1)*(j+2)/2 + i];
              	      }
                   else
             	     {Fn1[i*(i+1)/2 + j] += s*D*Fn[i*(i+1)/2 + j+1];
              	      }
                   }
                else
                  {Fn1[i*(i+1)/2 + j] += s*D*Fn[i*(i+1)/2 + j-1];
                  }

       	        //(i, j+1)
                if(i==j)
                  {if(i==nx-1)
              	     {Fn1[i*(i+1)/2 + j] += s*D*Fn[i*(i+1)/2 + j];
              	      }
                   else
             	     {Fn1[i*(i+1)/2 + j] += s*D*Fn[(j+1)*(j+2)/2 + i];
              	      }
                   }
                else
                  {Fn1[i*(i+1)/2 + j] += s*D*Fn[i*(i+1)/2 + j+1];
                   }
                }
          }
        //add contribution from source
        Fn1[0] += dt;//*source(0,0);
        memcpy(Fn, Fn1, row);
        }

    delete[] Fn;
    delete[] Fn1;
    fout.close();
    return 0;
    }
