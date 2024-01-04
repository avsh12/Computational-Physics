#include<iostream>
#include<math.h>
#include<iomanip>   //precision set
#include<fstream>

using namespace std;

class Derivative
    {public:
     float xmin, xmax, h;
     int arrsize;
     //Dynamical allocation for arrays
     float *X, *df, *Fdf, *Bdf, *Cdf; 

     //parameterized constructor. Takes minimum and maximum values of x and the step size
     Derivative(float a, float b, float c)
        {xmin = a;
         xmax = b;
         h = c;
         //calculate array size required
         arrsize  = (int)ceil((xmax-xmin)/h);  
         //allocate array to store x values
         X = new float[arrsize];        
         X[0]=xmin;
         for(int i=1; i!=arrsize; ++i)
             //assign x values in step of h
             X[i] = X[i-1]+h;           
        }
    
     //Default constructor with default xmin =0, xmax=1, h=0.01
     Derivative()
        {xmin = 0;
         xmax = 1;
         h = 0.01;
         //calculate array size required
         arrsize  = (int)ceil(xmax/h);  
         //allocate array to store x values 
         X = new float[arrsize];         
         X[0]=xmin;
         for(int i=1; i!=arrsize; ++i)
             //assign x values in step of h
             X[i] = X[i-1]+h;            
        }

     //Destructor, release memory acquired by pointers
     ~Derivative()
        {delete[] X, df, Fdf, Bdf, Cdf;
        }

    //forward derivative at a point
    float FDerx(float (*fx)(float x), float x)
         {return (fx(x+h)-fx(x))/h;
         }

    //Backward derivative at a point
    float BDerx(float (*fx)(float x), float x)
         {return (fx(x) - fx(x-h))/h;
         }

    //central derivative at a point
    float CDerx(float (*fx)(float x), float x)
         {return (fx(x+h)-fx(x-h))/(2*h);
         }

    //Analytic derivative values in a range
    void ADer(float (*ADerx)(float x))
         {df = new float[arrsize];       //allocate array to store analytic derivative values
          float x=xmin;
          for(int i=0; i!=arrsize; ++i)
             {df[i] = ADerx(x);          //assign derivative values
              x += h;
             }
         }

    //forward derivative in a range
    void FDer(float (*fx)(float x))
         {Fdf = new float[arrsize];      //allocate array to store forward derivative values
          float x=xmin;
          for(int i=0; i!=arrsize; ++i)
             {Fdf[i] = FDerx(fx, x);         //assign derivative values
              x += h;
             }
         }
  
    //backward derivative in a range
    void BDer(float (*fx)(float x))
         {Bdf = new float[arrsize];      //allocate array to store backward derivative values
          float x=xmin;
          for(int i=0; i!=arrsize; ++i)
             {Bdf[i] = BDerx(fx, x);         //assign derivative values
              x += h;
             }
         }

    //central derivative in a range
    void CDer(float (*fx)(float x))
         {//allocate array to store central derivative values
          Cdf = new float[arrsize];      
          float x=xmin;
          for(int i=0; i!=arrsize; ++i)
             {//assign derivative values
              Cdf[i] = CDerx(fx, x);         
              x += h;
             }
         }
    };

//function whose derivative is to be calculated
     float fx(float x)
        {float value = sin(x);
         return value - pow(value,3);
         }

//Analytic derivative at a point
    float ADerx(float x)
         {return cos(x)*(1-3*pow(sin(x),2));
         }

int main()
   {//object for with parameter minimum x , maximum x, step size
    Derivative Der(-M_PI/16,M_PI/8,0.01);      
    Der.ADer(&ADerx);
    Der.FDer(&fx);
    Der.BDer(&fx);
    Der.CDer(&fx);
    // fetch the size of arrsize from the object
    int size = Der.arrsize;         

    fstream fout;
    fout.open("derivativeData.txt",ios::out|ios::trunc);

    cout<<"Derivative of sin(x)-sin(x)^3 in the range x=0 to x="<<0.5<<" are:\n\n";
    cout<<"x\t df\t\t Fdf\t\t Bdf\t\t Cdf"<<endl<<endl;

    for(int i=0; i!=size; ++i)
       {fout<<fixed<<setprecision(10)<<Der.X[i]<<fixed<<setprecision(10)<<"\t"<<Der.df[i]<<"\t"<<Der.Fdf[i]<<"\t"<<Der.Bdf[i]<<"\t"<<Der.Cdf[i]<<endl;
        }

    return 0;
   }



