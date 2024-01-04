#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include"ODESolver.h"

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_blas.h>

#include<time.h>

using namespace std;

void acc(long double x, long double y, long double arr[])
    {arr[0] = -2*x*y*y;
     arr[1] = -2*x*x*y;
     }

class Trajectory
   {public:
    long double initial_val[5];
    long double final_val[5];
    
    void setInitialVal(RKSolver sol)
        {for(int i=0; i<5; ++i)
             initial_val[i] = sol.IC[i];
         }
 
    void setFinalVal(RKSolver sol)
        {for(int i=0; i<5; ++i)
             final_val[i] = sol.u[0][i];
         }
 
    void static diffInFinalVal(Trajectory traj1, Trajectory traj2, long double *diff)
        {for(int i=0; i<4; ++i)
             diff[i] = traj2.final_val[i+1] - traj1.final_val[i+1];
         }

    void static diffInInitialVal(Trajectory traj1, Trajectory traj2, long double *diff)
        {for(int i=0; i<5; ++i)
             diff[i] = traj2.initial_val[i] - traj1.initial_val[i];
         }
    };
//Trajectory class ends

class LyapunovSpectrum
    {public:
     gsl_matrix *Jacobian, *Lambda;
     gsl_vector *eigenvalues, *spectrum;
     gsl_eigen_symm_workspace *param;
     

     LyapunovSpectrum()
        {Jacobian = gsl_matrix_alloc(4,4);
         Lambda = gsl_matrix_alloc(4,4);
         eigenvalues = gsl_vector_alloc(4);
         spectrum = gsl_vector_alloc(4);
         param = gsl_eigen_symm_alloc(4); 
         }

     ~LyapunovSpectrum()
        {gsl_matrix_free(Jacobian);
         gsl_matrix_free(Lambda);
         gsl_vector_free(eigenvalues);
         gsl_vector_free(spectrum);
         gsl_eigen_symm_free(param);
         }

     void setJacobian(long double iDiff[5], long double fDiff[5])
        {for(int i=0; i<4; ++i)
            {for(int j=0; j<4; ++j)
                {long double val = fDiff[j]/iDiff[i];
                 gsl_matrix_set(Jacobian, i, j, val);
                 //CONSOLE OUTPUT
                 cout<<val<<"  ";
                 }
             //CONSOLE OUTPUT
             cout<<endl;
             }
         }
    
     //calculates Jacobian*Trans(Jacobian) matrix
     void setLambda()
         {gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Jacobian, Jacobian, 0.0, Lambda);

           //CONSOLE OUTPUT
           cout<<endl;
           for(int i=0; i<4; ++i)
            {for(int j=0; j<4; ++j)
                {cout<<gsl_matrix_get(Lambda, i, j)<<"  ";
                 }
             cout<<endl;
             }
          }

     //calculates eigenvalues of Jacobian*Trans(Jacobian) matrix
     void calcEigenVal()
         {gsl_eigen_symm(Lambda, eigenvalues, param);
          //CONSOLE OUTPUT
          cout<<endl;
          for(int i=0; i<4; ++i)
              cout<<gsl_vector_get(eigenvalues, i)<<"   ";
          cout<<endl;
          }

     //stores eigenvalues of Jacobian*Trans(Jacobian) matrix to an array
     void getEigenVal(long double *eigenval)
         {for(int i=0; i<4; ++i)
             {eigenval[i] = gsl_vector_get(eigenvalues, i);
              }
          }

     //calculates Lyapunov spectrum 
     void calcSpectrum(double time)
         {for(int i=0; i<4; ++i)
             {long double val = log(gsl_vector_get(eigenvalues, i))/(2*time);
              gsl_vector_set(spectrum, i, val);

              //CONSOLE OUTPUT
              cout<<val<<"  ";
              }
          cout<<endl;
          }

     //stores Lyapunov spectrum to an array
     void getSpectrum(long double *spec)
         {for(int i=0; i<4; ++i)
              spec[i] = gsl_vector_get(spectrum, i);
          }

     };
//LyapunovSpectrum class ends

int main()
   {clock_t start, end;
    start = clock();

    double time = 40000;
    RKSolver sol(time, 0.001);
    fstream fout;
    fout.open("LyapunovSpectrum.txt", ios::out|ios::trunc);

    Trajectory trajectory[2];
    //difference between the initial values of the two trajectories
    long double initial_diff[4] = {1.0e-20, 1.0e-20, 1.0e-20, 2.0e-20};
    long double final_diff[4];

    for(int traj=0; traj<=1; ++traj)
       {if(traj==0)
          {//initial condition for first trajectory
           sol.setInitialVal(0, 0.1, 0.5, 0, 0.1);
           }
        else
          {//initial condition for second trajectory
           sol.initialIncrement(initial_diff);
           }

        for(int i=0; i<5; ++i)
            trajectory[traj].initial_val[i] = sol.IC[i];
        
        for(int i=1; i<=sol.Nt; ++i)
           {//calculate position at time i*dt
            sol.nextPosition(acc);
            }

        for(int i=0; i<5; ++i)
            trajectory[traj].final_val[i] = sol.u[0][i];

        {//CONSOLE OUTPUT
         cout<<"coordinates at time "<<trajectory[traj].final_val[0]<<" of trajectory "<<traj+1<<endl;
         for(int i=0; i<5; ++i)
             cout<<trajectory[traj].final_val[i]<<"  ";
         cout<<endl;
         }
        }

    //calculate difference in final values of the two trajectories
    trajectory[0].diffInFinalVal(trajectory[0], trajectory[1], final_diff);

    {//CONSOLE OUTPUT
     cout<<"Diffrence in final coordinates of the two trajectories"<<endl;
     for(int i=0; i<4; ++i)
         cout<<final_diff[i]<<"  ";
     cout<<endl<<endl;
     }

    LyapunovSpectrum spectrum;
    //calculate Jacobian matrix
    {//CONSOLE OUTPUT
     cout<<"Jacobian matrix"<<endl;
     }
    spectrum.setJacobian(initial_diff, final_diff);
    //calculate Jacobian*Trans(Jacobian) matrix
    {//CONSOLE OUTPUT
     cout<<"Jacobian*Trans(Jacobian) matrix"<<endl;
     }
    spectrum.setLambda();
    //calculate eigenvalues of Jacobian*Trans(Jacobian) matrix
    {//CONSOLE OUTPUT
     cout<<"Eigenvalues of Jacobian*Trans(Jacobian) matrix"<<endl;
     }
    spectrum.calcEigenVal();
    //calculate Lyapunov spectrum
    {//CONSOLE OUTPUT
     cout<<"Lyapunov spectrum"<<endl;
     }
    spectrum.calcSpectrum(time);
    
    long double lspectrum[4];
    //store the Lyapunov spectrum in the array lspectrum
    spectrum.getSpectrum(lspectrum);

    fout.close();
    end = clock();
    //program run time
    cout<<endl<<double(end-start)/double(CLOCKS_PER_SEC);
    return 0;
    }
