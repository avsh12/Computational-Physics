#include<math.h>
#include <gsl/gsl_integration.h>
#include<functional>

//Takes function and range in which to be integrated
double NIntTrapezoid(double(*f)(double x), double a, double b)
   {double h=(b-a);
    return (f(a) + f(b))*h/2;
    }

double NIntOneThird(double(*f)(double x), double a, double b)
   {double h=(b-a)/2;
    return (f(a) + 4*f(a+h) + f(a+2*h))*h/3;
    }

double NIntThreeEighth(double(*f)(double x), double a, double b)
   {double h=(b-a)/3;
    return (3*f(a) + 9*(f(a+h) + f(a+2*h)) + 3*f(a+3*h))*h/8;
    }

double NIntBodeRule(double(*f)(double x), double a, double b)
   {double h=(b-a)/4;
    return (7*f(a) + 32*f(a+h)+ 12*f(a+2*h) + 32*f(a+3*h) + 7*f(a+4*h))*2*h/45;
    }

//Takes function, the range of integration, and the number of discretization
double NIntExtendedRule1(double(*f)(double x), double a, double b, int n)
   {//int n = int(ceil((b-a)/h));
    double h=(b-a)/n;
    double Int = (f(a) + f(a+h*n))*h/2;

    for(int i=1; i<n; ++i)
       {Int += h*f(a + h*i);
        }
    return Int;
    }

double NIntExtendedRule2(double(*f)(double x), double a, double b, int n)
   {//int n = int(ceil((b-a)/h));
    double h=(b-a)/n;
    double Int = (f(a) + f(a+n*h))*h/3;

    for(int i=1; i<n; ++i)
       {if(i%2 == 0)
           Int += f(a + h*i)*h*2/3;
        else 
           Int += f(a + h*i)*h*4/3;
        }
    return Int;
    }

double gaussianQuad(double(*f)(double x), double a, double b, int nodes)
      {double xi, wi, Int;

       gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(nodes);

       for(int i=0; i<nodes; ++i)
          {gsl_integration_glfixed_point(a, b, i, &xi, &wi, table);
           Int += f(xi)*wi;
           }

       return Int;
       }

double gaussianQuad(double a, double b, int nodes, std::function<double(double)> f)
      {double xi, wi, Int;

       gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(nodes);

       for(int i=0; i<nodes; ++i)
          {gsl_integration_glfixed_point(a, b, i, &xi, &wi, table);
           Int += f(xi)*wi;
           }

       return Int;
       }

