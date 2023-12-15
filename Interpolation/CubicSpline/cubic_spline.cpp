#include<iostream>
#include<fstream>
#include<math.h>
#include<random>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>

using namespace std;

//function from which sample points are taken
double function(double x)
      {return 5*x*x*exp(-x*x)*cos(3*x);
       }

int main()
  {int data_len = 60;
   double xrange = 6.0;
   double random_range = xrange/data_len;
   double *xdata, *ydata;
   xdata = new double[data_len];
   ydata = new double[data_len];

   random_device dev;
   mt19937 rng(dev());
   //generates uniform real random numbers in the specified range
   uniform_real_distribution<> dist6(0,random_range);

   //get data_len number of random sample x values and corresponding function values
   xdata[0] = -3.0;
   ydata[0] = function(xdata[0]);
   for(int i=1; i<data_len; ++i)
      {xdata[i] = xdata[0] + (i-1)*random_range + dist6(rng);
       ydata[i] = function(xdata[i]);
       }

   //create interpolation object of type cubic spline
   gsl_interp *workspace = gsl_interp_alloc(gsl_interp_cspline, data_len);
   gsl_interp_init(workspace, xdata, ydata, data_len);

   //Index Look-ip and Acceleration object
   gsl_interp_accel *accel_object = gsl_interp_accel_alloc();

   fstream fout;
   fout.open("cspline_data.txt", ios::out|ios::trunc);

   //write the interpolating values and exact values
   for(double x=xdata[0]; x<xdata[data_len-1]; x+=0.03)
      {function(x);
       fout<<x<<"\t"<<gsl_interp_eval(workspace, xdata, ydata, x, accel_object)<<"\t"<<function(x)<<endl;
       }

   fout<<endl<<endl;
   //write the sampled data
   for(int i=0; i<data_len; ++i)
      {fout<<xdata[i]<<"\t"<<ydata[i]<<endl;
       }

   fout.close();
   delete[] xdata, ydata;
   gsl_interp_free(workspace);
   gsl_interp_accel_free(accel_object);
   return 0;
   }
