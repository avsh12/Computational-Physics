#include<math.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>


class Interpolation
     {public:
      double *X, *Y;
      int data_len;

      Interpolation(int data_len)
          {this->data_len = data_len;
           X = new double[data_len];
           Y = new double[data_len];
           }
      Interpolation(double *x_arr, double *y_arr)
          {data_len = sizeof(x_arr)/sizeof(double);
           X = new double[data_len];
           Y = new double[data_len];
           }

      ~Interpolation()
          {delete[] X, Y;
           }

      void setY(double (*func)(double x))
          {for(int i=0; i<data_len; ++i)
              {Y[i] = func(X[i]);
               }
           }

      void setY(double *y_arr)
          {for(int i=0; i<data_len; ++i)
              {Y[i] = y_arr[i];
               }
           }
     
      void setX(double *x_arr)
          {for(int i=0; i<data_len; ++i)
              {X[i] = x_arr[i];
               }
           }
      };
//class Interpolation ends

class LagrangeInterpolation: public Interpolation
      {public:
       LagrangeInterpolation(int data_len): Interpolation(data_len){;}
       LagrangeInterpolation(double *x_arr, double *y_arr): Interpolation(x_arr, y_arr){;}

       double lagrangeBasisFunc(double x, int index)
          {double Product=1;
           for(int j=0; j<data_len; ++j)
              {//skip the index term
               if(j==index)
                  continue;
               Product *= (x-X[j])/(X[index] - X[j]);
               }

           return Product;
           }

      double InterpolatingPolynomial(double x)
          {double interpolating_value=0;
           for(int j=0; j<data_len; ++j)
              {interpolating_value += Y[j]*lagrangeBasisFunc(x, j);
               }

           return interpolating_value;
           }
      };
//class LagrangeInterpolation ends

class NevilleInterpolation: public Interpolation
     {public:
      NevilleInterpolation(int data_len): Interpolation(data_len){;}
      NevilleInterpolation(double *x_arr, double *y_arr): Interpolation(x_arr, y_arr){;}

      double recurrencePolynomial(double x, int index_i, int index_j)
            {if(index_i >=0 && index_j<=data_len-1 && index_j>=index_i)
               {if(index_i == index_j)
                  {return Y[index_i];
                   }
                else
                  {return ((x-X[index_i])*recurrencePolynomial(x,index_i+1, index_j) - (x-X[index_j])*recurrencePolynomial(x,index_i, index_j-1))/(X[index_j] - X[index_i]);
                   }
                }
             else
                return 0;
             }

      double InterpolatingPolynomial(double x)
            {double *poly_array1, *poly_array2;
             poly_array1 = new double[data_len];

             for(int i=0; i<data_len; ++i)
                {poly_array1[i] = Y[i];
                 }

             //traverse the triangle from left to right
             for(int j=1; j<data_len; ++j)
                {poly_array2 = new double[data_len-j];

                 //traverse the column of the triangle
                 for(int i=0; i<data_len-j; ++i)
                    {poly_array2[i] = ((x-X[i])*poly_array1[i+1] - (x-X[i+j])*poly_array1[i])/(X[i+j]-X[i]);
                     }

                 delete[] poly_array1;
                 poly_array1 = new double[data_len-j];
                 for(int k=0; k<data_len-j; ++k)
                    {poly_array1[k] = poly_array2[k];
                     }
                 delete[] poly_array2;
                 }
             
             double val = poly_array1[0];
             delete[] poly_array1;
             
             return val;
             }

      };
