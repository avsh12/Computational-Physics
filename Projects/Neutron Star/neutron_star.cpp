//Masses of four Neutron stars:
//PRS J1614-2230 - M = 1.97 ± 0.04 M⊙
//PSR J0348+0432 - M = 2.01 ± 0.04 M⊙
//    J0740+6620 - M = 2.08 ± 0.07 M⊙
//PSR J1311-3430 - M = 2.55 ± 0.50 M⊙  black widow pulsar (BWP)


#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<time.h>
#include "bisection.h"

//constants
#define mass_solar 1.11544986e60; // MeV/c^2 : 1.98847e30 kg
#define mass_neutron 938.926; // MeV/c^2
#define speed_light 2.99792458e23 // fm s^-1
#define G 1.18976283e5; // (MeV/c^2)^-1 fm^3 s^-2
#define NSD 157.0684808 // MeV/c^2 fm^-3  nuclear saturation density (2.8e17 kg m^-3)

//conversion factors
#define kg_to_mevperc2 5.6095886e29; //MeV/c^2

//initial condition in dimensionless units
//r = 0, m = 0, rho = 560.95886

using namespace std;

struct conversion_factors
         {double R_0 = 1.0e19; //fm
          //(pow(speed_light, 2)*R_0/(G*mass_solar)) * C.mass_solar
          double M_0 = 6.7722184*mass_solar; //MeV/c^2
          //pow(speed_light, 2) / (4*pi*G* R_0*R_0)
          double rho_s = 601.133796; //MeV/c^2 fm^-3
          double P_s = rho_s;
          };

class massRadius
     {public:
      double R, P, rho, n;
      double m, r;
      //dimensionless initial values
      double P_0, m_0, r_0, rho_0;
      //conversion units
      conversion_factors factors;

      double dr;
      
      massRadius(double P_0, double m_0, double r_0)
                {P = this->P_0 = P_0/factors.P_s;
                 m = this->m_0 = m_0/factors.M_0;
                 r = this->r_0 = r_0/factors.R_0;
                 rho = rho_0 = density(this->P_0);
                 }

      massRadius()
                {}

      void updateCoreDensity(double P_0, double m_0, double r_0)
                 {P = this->P_0 = P_0/factors.P_s;
                  m = this->m_0 = m_0/factors.M_0;
                  r = this->r_0 = r_0/factors.R_0;
                  rho = rho_0 = density(this->P_0);
                  }

      //dimensionless
      double density(double P)
                 {//P/1.54 + P * mass_neutron/(pow(factors.rho_s, 1.54/2.54) * pow(363.44, 1/2.54))
                  return  P/1.54 + 1.904140406*pow(P, 0.3937);
                  }

      //dimensionless
      double Pressure(double rho)
                 {return findRoot(0, pow(0.52517*rho, 2.54), 1.0e-7, [rho](double x)->double{return x + 2.932376*pow(x, 0.3937) - 1.54*rho;});
                  }

      double mDer(double r)
                 {return rho*r*r;
                  }

      double PDer(double r)
                 {return -m*rho/(r*r);
                  }

      double PDerTOV(double r)
                 {return -(P + rho)*(pow(r,3)*P + m)/(r*r - 2*m*r);
                  }

      void nextStep()
              {double k1 = PDer(r);
               double k2 = PDer(r+dr/2);
               double k3 = PDer(r+dr/2);
               double k4 = PDer(r+dr);

               P += dr*(k1+2*k2 + 2*k3 + k4)/6;

               k1 = mDer(r);
               k2 = mDer(r+dr/2);
               k3 = mDer(r+dr/2);
               k4 = mDer(r+dr);

               m += dr*(k1+2*k2 + 2*k3 + k4)/6;

               rho = density(P);

               r += dr;
               }

      void fileWrite(fstream fout)
           {fout<<r<<"\t"<<m<<endl;
            }

      };

int count = 0;

int main()
     {clock_t start, end;
      start = clock();

      double initial_density_dim, initial_pressure_dim, initial_radius_dim, initial_mass_dim;
      double P_tolerance;
      int stride, index;
      massRadius nStar;

      fstream fout;
      fout.open("mass_radius.txt", ios::out|ios::trunc);

/*      
      cout<<endl<<"Initial values in dimensionless units: "<<endl;
      cout<<"radius\t mass\t pressure\t density"<<endl;
      cout<<nStar.r<<"\t"<<nStar.m<<"\t"<<nStar.P<<"\t"<<nStar.rho<<endl;
*/
      fout<<"Radius\t Mass\t pressure\tdensity"<<endl<<endl<<endl;

      initial_density_dim = 560.95886;
      initial_mass_dim = 0;
      initial_radius_dim = 1.0e-20; 
          
      P_tolerance = 1.0e-9;
      double step[] = {5.0e-8, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.2e-3};
      index = 1;

      for(int i=0; i<6; ++i)
         {stride = 1.5/(1000*step[i]);
          count=0;

          nStar.dr = step[i];
          initial_pressure_dim = findRoot(0, pow(0.01085*initial_density_dim, 2.54), 1.0e-7, 
                 [initial_density_dim](double x)
                        ->double{return x + 141.94095*pow(x, 0.3937) - 1.54*initial_density_dim;});
          
          nStar.updateCoreDensity(initial_pressure_dim, initial_mass_dim, initial_radius_dim);

          fout<<"#index"<<index++<<". dr = "<<nStar.dr<<endl;
          for(count=0; nStar.P >= P_tolerance;)
             {if(isnan(nStar.rho)) break;
              if(count%stride==0)
                  {fout<<setprecision(10)<<nStar.r<<"\t"<<nStar.m*6.7722184<<"\t";
                   //pressure in 10^32 Pa
                   fout<<(nStar.P*nStar.factors.rho_s*1.6021766)<<"\t";
                   fout<<(nStar.rho*nStar.factors.rho_s/NSD)<<endl;
                   }
              nStar.nextStep();
              
              count++;          
              }

          if(isnan(nStar.rho))
            {fout<<endl<<endl;
             continue;
             }
          fout<<setprecision(10)<<nStar.r<<"\t"<<nStar.m*6.7722184<<"\t";
          fout<<(nStar.P*nStar.factors.rho_s*1.6021766)<<"\t";
          fout<<(nStar.rho*nStar.factors.rho_s/NSD)<<endl;
          fout<<endl<<endl;

 //        fout<<nStar.r<<"\t"<<nStar.m*6.7722184<<"\t"<<nStar.P_0<<"\t";
 //        fout<<nStar.rho_0*nStar.factors.rho_s/NSD<<endl;
         }

      cout<<endl<<"Loop ended at: "<<endl;
      cout<<"radius\t\t mass\t\t pressure\t\t density"<<endl;
      cout<<setprecision(10)<<nStar.r<<"\t"<<nStar.m<<"\t"<<nStar.P<<"\t"<<nStar.rho<<endl;

      fout.close();
  
      end = clock();
      //program run time
      cout<<"Run time: "<<double(end-start)/double(CLOCKS_PER_SEC)<<endl; 
      return 0;
      }
