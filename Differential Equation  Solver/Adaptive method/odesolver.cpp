#include<iostream>
#include<fstream>
#include<math.h>
#include<time.h>
#include<iomanip>

using namespace std;

double eulerSolver(double x0, double y0, double dt)
      {double time=x0;
       double y = y0;
       for(; time<=1; time += dt)
          {y = y*(1+dt) + dt*(1-4*time);
           }

       return y;
       }

int main()
   {clock_t start, end;
    start = clock();

    double dt1 = 0.5, dt2 = 0.6;
    double y1, y2;
    double tolerance = 1.0e-8;
    int index = 0, iteration = 0;

    y1 = eulerSolver(0, 1, dt1);
    y2 = eulerSolver(0, 1, dt2);
    double next_step;

    for(int i=1; ;++i, ++iteration)
       {next_step = 0.9*tolerance*abs(dt1 - dt2)/abs(y2 - y1);

        if(i%2==1)
          {dt1 = next_step;
           y1 = eulerSolver(0, 1, dt1);

           if((abs(y2-y1)/tolerance)<1)
             {index=1;
              break;
              }
           }
        else
          {dt2 = next_step;
           y2 = eulerSolver(0, 1, dt2);

           if((abs(y2-y1)/tolerance)<1)
             {index=2;
              break;
              }
           }
        }

    if(index==1)
       cout<<endl<<setprecision(10)<<y1;
    else
       cout<<endl<<setprecision(10)<<y2;

    cout<<endl<<"Iterations: "<<iteration;

    end = clock();
    //program run time
    cout<<endl<<"Run time: "<<double(end-start)/double(CLOCKS_PER_SEC)<<endl;

    return 0;
    }
