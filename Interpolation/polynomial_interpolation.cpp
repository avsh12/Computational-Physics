#include<iostream>
#include<math.h>
#include<random>
#include<fstream>
#include<iomanip>
#include"interpolation.h"

using namespace std;

//function from which sample points are taken
double function(double x)
      {return 5*x*x*exp(-x*x)*cos(3*x);
       }

int main()
   {random_device dev;
    mt19937 rng(dev());

    //number of points to be sampled
    int data_len = 60;
    //generates uniform real random numbers in the specified range
    uniform_real_distribution<> dist6(-3.0,3.0);
    
    double x_points[data_len];
    //get data_len number of random sample x values
    for(int i=0; i<data_len; ++i)
       {x_points[i] = dist6(rng);
        }

    //create object for interpolation
    NevilleInterpolation poly_interp(data_len);
    //set the sampled x points
    poly_interp.setX(x_points);
    //set the values of function at the given x points using the given function
    poly_interp.setY(&function);

    fstream fout;
    fout.open("poly_interp_data.txt", ios::out|ios::trunc);

    //Interpolate the data
    for(double x=-3.0; x<=3.01; x+=0.05) {fout<<fixed<<setprecision(6)<<x<<"\t"<<function(x)<<"\t"<<poly_interp.InterpolatingPolynomial(x)<<endl;
        }

    fout<<endl<<endl;
    //write data that is interpolated 
    for(int i; i<=data_len; ++i)
       {fout<<poly_interp.X[i]<<"\t"<<poly_interp.Y[i]<<endl;
        }

    fout.close();

    return 0;
    }
