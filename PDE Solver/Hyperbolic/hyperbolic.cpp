#include<iostream>
#include<math.h>
#include "PDESolver.h"

class BoundaryCondition:public Hyperbolic
     {public:
      static float xstep, xrange;
      static float phi(float x)
               {return exp(-(x-xrange/2)*(x-xrange/2));
               }
      static float u_t0(float x)
               {return phi(x);}

      static float u_t1(float x)                //initial distribution at t=1, i.e. u(x,t=1)
               {return s*(phi(x+xstep) + phi(x-xstep))/2 + (1-s)*phi(x);                //correct here
                }//t1_dist function ends

      static float xmin_dist(float x)              //initial distribution at x=xmin, i.e. u(xmin,t)
               {return sin(M_PI*x);
                }//xmin_dist function ends

      static float xmax_dist(float x)              //initial distribution at x=xmzx, i.e. u(xmax,t)
               {return sin(M_PI*x);
                }//xmax_dist function ends
      };

float BoundaryCondition::xstep;
float BoundaryCondition::xrange;

int main()
   {Matrix M(20,30,0.1);
    float s=0.9;
    M.tStep(sqrt(s)*M.xStep());
    M.updateMatrix();             //create an array

    BoundaryCondition BC;
    BoundaryCondition::xstep = M.xStep();
    BoundaryCondition::xrange = M.xRange();
    //cout<<BC.xrange;

    MatrixOperation MatOp(M.Nx, M.Nt, M.dx, M.dt);
    MatOp.tBoundary(&BC.u_t0, 0, M.u);
    MatOp.xBoundary(0.0, 0, M.u);
    MatOp.xBoundary(0.0, M.Nx-1, M.u);
    MatOp.tBoundary(&BC.u_t1, 1, M.u);

    MatOp.setMatrix(&BC.functionDef, 1, M.Nx-2, 2, M.Nt-1, M.u);
    
    MatOp.fileOutput("HyperbolicData.txt", M.u);
    
    return 0;
    }
