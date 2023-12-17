#include<iostream>
#include<math.h>
#include"PDESolver.h"

class BoundaryCondition: public Advection
     {public:
      BoundaryCondition(float s):Advection(s)
            {;}
      static float u_t0(float x)
            {if(x>=4.5 && x<=5.5)
                return 100*sin((x-4.5)*M_PI);
             else return 0;
             }
      /*static float u_t0(float x)
            {return exp(-x*x/(2*0.0001))/sqrt(0.0001);
             }*/
      };

int main()
   {Matrix M(10,2,0.02);
    float s=1;
    M.tStep(s*M.xStep());
    M.updateMatrix();
    
    BoundaryCondition BC(s);
    MatrixOperation MatOp(M.Nx, M.Nt, M.dx, M.dt);
    MatOp.tBoundary(&BC.u_t0, 0, M.u);
    
    MatOp.setMatrix(&BC.functionDef, 1, M.Nx-2, 1, M.Nt-1, M.u);
    
    MatOp.fileOutput("AdvectionData.txt", M.u);
    //system("/bin/gnuplot /home/adarshv/Documents/'Numerical Methods'/PDE/Advection/plot.txt");

    return 0;
    }
