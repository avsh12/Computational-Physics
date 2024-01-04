//Take a function (1+ tanh((x-x0)/del))/2. Integrate it with a particular value of x0 and del. Now, taking the integral value find the value of x0 using bisection method.

//Output:
//Actual value of the abscissa of the midpoint of the function is 2.448
//Numerically obtained value using Bisection method is 2.447998


#include<iostream>
#include<math.h>
#include<iomanip>
#include "NInt.h"
#include "bisection.h"

using namespace std;

double f(double x, double x0, double del)
        {return (1+ tanh((x-x0)/del))/2;
         }

double intF(double a, double b, double x0, double del)
      {//gaussianQuad takes the range of integration, (a,b), number of nodes, and function to be integrated
       return gaussianQuad(a, b, 5, [x0, del](double x)->double{return f(x, x0, del);});
       }

int main()
   {double a, b, x0, del, Int, root, tolerance;

    tolerance = 1.0e-5;
    x0 = 2.448;
    del = 2;
    a = x0-6;
    b = x0 + 6;

    //Integrate with a particular value of x0 and del.
    Int = gaussianQuad(a, b, 5, [x0, del](double x)->double{ return f(x, x0, del); });

    //findRoot finds the zero of the function Int-intF(a, b, x, del) in the given range with given tolerance.
    //Here root is the value of x for which intF(a, b, x, del) = Int.
    root = findRoot(1, 4, tolerance, [Int, a, b, del](double x)->double{ return Int-intF(a, b, x, del); });

    cout<<"Actual value of the abscissa of the midpoint of the function is "<<setprecision(8)<<x0<<endl;
    cout<<"Numerically obtained value using Bisection method is "<<root<<endl;

    return 0;
    }
