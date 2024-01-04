#include<math.h>
#include<functional>

double findRoot(double(*func)(double x), double a, double b, double tolerance)
      {double midp = (a+b)/2;
       double val = func(midp);
       while(abs(val) > tolerance)
            {if(val < 0) a = midp;
             else b = midp;
               
             midp = (a+b)/2;
             val = func(midp);
             }

       return midp;
       }

double findRoot(double a, double b, double tolerance, std::function<double(double)> func)
      {double midp = (a+b)/2;
       double val = func(midp);
       while(abs(val) > tolerance)
            {if(val < 0) a = midp;
             else b = midp;
               
             midp = (a+b)/2;
             val = func(midp);
             }

       return midp;
       }
