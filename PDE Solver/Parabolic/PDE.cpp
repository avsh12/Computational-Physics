#include<iostream>
#include<iomanip>
#include<math.h>
#include<fstream>

using namespace std;

float phi(float x)
     {if(x<=0.5)
         return 2*x;
      else if (x<=1)
         return 2*(1-x);
      else
         return 0;
      }

float u1(float x, float t)
     {float u;
      
      for(int m=1; m<=100; ++m)
         u += sin(m*M_PI/2) * sin(m*M_PI*x) * exp(-t*pow(m*M_PI,2))/(m*m);

      u *= 8/(M_PI*M_PI);
      return u;
      }

int main()
   {float xrange = 1, trange = 1;
    float dx = 0.01, s = 0.51, dt;
    dt = dx*dx*s;
  
    int Nx = floor(xrange/dx)+1, Nt = floor(trange/dt);
    cout<<Nx<<"\t"<<Nt<<typeid(Nt).name()<<endl;

    float **u;
    u = new float*[Nt];
    for(int i=0;i<Nt;++i)
       {u[i] = new float[Nx];
        u[i][0] = 0;          //boundary condition
        u[i][Nx-1] = 0;       //boundary condition
        }
    
    for(int i=1; i<Nx-1; ++i)
       {u[0][i] = phi(i*dx);  //boundary condition
        //cout<<"("<<i*dx<<","<<phi(i*dx)<<")"<<"\t";
        }

    for(int n=1; n<Nt; ++n)
       {for(int i=1; i<Nx-1; ++i)
           {u[n][i] = s*(u[n-1][i+1] + u[n-1][i-1]) + (1-2*s)*u[n-1][i];
            }
        }
    cout<<(Nx-2)*dx<<","<<u[0][Nx-2]<<","<<(Nx-1)*dx<<","<<u[0][Nx-1];

    fstream foutn, fouta;
    foutn.open("PDENData.txt",ios::out|ios::trunc);
    fouta.open("PDEAData.txt",ios::out|ios::trunc); 
    //fout<<Nx<<"\t"<<Nt<<typeid(Nt).name()<<endl;

    for(int i=0; i<Nx; ++i)
       {
        for(int n=0; n<Nt; n+=200)
           {foutn<<fixed<<setprecision(2)<<i*dx<<"\t";
            foutn<<fixed<<setprecision(5)<<n*dt<<"\t";
            fouta<<fixed<<setprecision(2)<<i*dx<<"\t";
            fouta<<fixed<<setprecision(5)<<n*dt<<"\t";
            foutn<<fixed<<setprecision(5)<<u[n][i]<<"\n";
            fouta<<fixed<<setprecision(5)<<u1(i*dx,n*dt)<<"\n";
            }
        foutn<<endl;
        fouta<<endl;
        }

    for(int n=0; n<Nt; ++n)
        delete[] u[n];
    delete[] u;
    return 0;
    }
