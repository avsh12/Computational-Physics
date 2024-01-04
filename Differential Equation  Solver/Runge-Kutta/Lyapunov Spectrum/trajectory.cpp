#include<iostream>
#include<fstream>
#include<math.h>
#include"ODESolver.h"

using namespace std;

void acc(long double x, long double y, long double arr[])
    {arr[0] = -2*x*y*y;
     arr[1] = -2*x*x*y;
     }

int main()
   {RKSolver sol(95, 0.001);
    fstream fout;
    fout.open("trajectory.txt", ios::out|ios::trunc);

    //difference between the initial values of the two trajectories
    long double initial_diff[5] = {0, 1.0e-10, -1.0e-12, 0, 0};
    
    for(int traj=0; traj<=1; ++traj)
       {if(traj==0)
          {//initial condition for first trajectory
           sol.setInitialVal(0, 0.3, 0.1, 0, 1);
           }
        else
          {//initial condition for second trajectory
           sol.initialIncrement(initial_diff);
           }

        for(int i=1; i<=sol.Nt; ++i)
           {//calculate position at time i*dt
            sol.nextPosition(acc);

            //write position at time i*dt in the file
            long double x = sol.u[0][1], y = sol.u[0][2], vx = sol.u[0][3], vy = sol.u[0][4];
 	    fout<<sol.u[0][0]<<"\t"<<x<<"\t"<<y<<"\t"<<vx<<"\t"<<vy<<endl;
            }

        //double blank line to seperate datasets of the two trajectories
        fout<<endl<<endl;
        }

    fout.close();
    return 0;
    }
