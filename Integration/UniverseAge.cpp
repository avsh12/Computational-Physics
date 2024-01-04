//Universe age calculated using Extended Trapezoid
//Calculated Universe age: 13.7976 by

#include<iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<stdlib.h>
#include "NInt.h"

using namespace std;

double omega_r = 2.47e-5/pow(0.6736, 2);
double omega_m = 0.3153;
double omega_l = 1 - omega_m - omega_r;
double omega_k = 1 - omega_m - omega_r - omega_l;
double H0 = 67.36;
double MPc = 3.085678e22;
double Km = 1000;
double Gyr = 3.1556926e16;

//function to be integrated in the range [0,1]
double func(double x)
      {return x/(H0*sqrt(omega_r + omega_m*x + omega_k*x*x + omega_l*x*x*x*x));
       }

int main()
    {double age=0;

     age = NIntExtendedRule1(&func, 0,1, 1000)*MPc/(Km*Gyr);
     cout<<"Age of the universe: "<<age<<endl;
     
     return 0;
     }
