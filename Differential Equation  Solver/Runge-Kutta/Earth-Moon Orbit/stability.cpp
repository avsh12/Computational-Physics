#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include"Solver.h"

//Time is measured in days: 1d = 86400s
//Distance is measured in units of billion meters: 1bm = 10^9m
//Mass is measured in units of Earth mass: 1M = 5.9724e24kg

//G = 6.67408e-11m^3kg^(-1)s^(-2) = 2.9756e-3 bm^3M^(-1)d^(-2)
//Moon mass = 7.34767309e22kg = 1.2303e-2M
//Apogee distance = 405500000m = 0.4055bm
//Time period of moon = 393430

using namespace std;

int main()
   {//masses of Earth and Moon
    twoBody B(1,1.2303e-2);
    //file object
    fstream fout;
    fout.open("energystability.txt",ios::out|ios::app);

    fout<<endl<<endl;

    for(double dt=0.001; dt<=2.731; dt+=0.001)
       {//solver object. Contains solution at all times
        //RKSolver sol(27.31,dt);
        leapFrog sol(27.31,dt);
        
        //set initial value
        sol.setInitialVal(0, 0.4055, 0, 0, 8.3808e-2);

        //calculate energy from initial values
        double E0 = B.V(0.4055,0) + B.m*pow((sol.u[0][4]),2)/2;
        double sum_energy = 0, sum_energy_sq = 0;
    
        //RK method solver
        /*for(int i=1; i<sol.Nt; ++i)
           {sol.nextPosition(&B.a, i);
            }*/

        //leapFrog solver
        sol.solver(&B.a);
    
        double alpha = B.m2/(B.m1+B.m2), beta = B.m1/(B.m1+B.m2);
        double *Et = new double[sol.Nt];
             
        for(int i=0; i<sol.Nt; i+=1)
           {double x = sol.u[i][1], y = sol.u[i][2];
            double v1x = -alpha*sol.u[i][3], v1y = -alpha*sol.u[i][4], v2x = beta*sol.u[i][3], v2y = beta*sol.u[i][4];
            double KE = B.m1*(v1x*v1x+v1y*v1y)/2 + B.m2*(v2x*v2x+v2y*v2y)/2;
            double V = B.V(x,y);
        
            //sum numerical energies
            sum_energy += KE+V;
            Et[i] = KE+V;
            }

        //calculate average energy
        double avg_energy = sum_energy/sol.Nt;

        //calculate sum of squares of deviation of numerical energy from the numerical average
        for(int i=0; i<sol.Nt; ++i)
            sum_energy_sq += pow((avg_energy - Et[i]), 2);
        //standard deviation of numerical energy vs time distribution
        double st_dev = sqrt(sum_energy_sq/sol.Nt);
    
        fout<<dt<<"\t"<<E0<<"\t"<<avg_energy<<"\t"<<st_dev<<"\t"<<abs(st_dev/avg_energy)<<endl;
        delete[] Et;
        }
    fout.close();
    
    return 0;
    }
