#include<iostream>
#include<math.h>
#include"PDESolver.h"

//inherit Parabolic class within which the recurrence relation is defined
class BoundaryCondition: public Parabolic
     {public:
      BoundaryCondition(float s):Parabolic(s)
            {;}
      //initial profile, at t=0
      static float u_t0(float x)
            {if(x>=0 && x<=0.5)
                return x;
             else if(x>=0.5 && x<=1)
                return (1-x);
             else return 0;
             }
      };

int main()
   {float xmax =1, tmax=0.1, x_step=0.01;
    //create object that contains the matrix to store the solution
    Matrix M(xmax, tmax, x_step);
    float s=0.48;

    //set time-step size dt = s*dx*dx
    M.tStep(s*M.xStep()*M.xStep());

    //create matrix of appropriate size to store the data
    M.updateMatrix();

    BoundaryCondition BC(s);
    //Object to set boundary condition
    MatrixOperation MatOp(M.Nx, M.Nt, M.dx, M.dt);
    //set initial profile u_t0 at t=0 in the matrix u 
    MatOp.tBoundary(&BC.u_t0, 0, M.u);
    //set 0.0 at x=0 in the matrix u
    MatOp.xBoundary(0.0,0,M.u);
    //set 0.0 at x=xmax in the matrix u
    MatOp.xBoundary(0.0,M.Nx-1, M.u);

    //calculate the whole distribution by functionDef function implementing the recurrence relation
    MatOp.setMatrix(&BC.functionDef, 1, M.Nx-2, 1, M.Nt-1, M.u);
    
    //write the output to the file
    MatOp.fileOutput("HeatEqnForwardData.txt", M.u);

    return 0;
    }
