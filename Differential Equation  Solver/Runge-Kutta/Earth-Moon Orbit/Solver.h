#include<math.h>

using namespace std;

class twoBody
     {public:
      static double G;       //Gravitational constant
      double E;       //Total Energy
      double l;       //Total angular momentum
      static double m1;      //Mass of first body
      static double m2;      //Mass of second body
      static double m;       //Reduced mass
      
      twoBody(double m1, double m2)
             {this->m1 = m1;
              this->m2 = m2;
              m = m1*m2/(m1+m2);
              G = 2.9756e-3;//6.67408e-11;
              }
      
      //Potential of the two body
      double V(double x, double y)
              {return -G*m1*m2/sqrt(x*x+y*y);
               }

      //acceleration of the body
      static void a(double x, double y, double value[])
            {double k = G*m1*m2;
             double r = sqrt(x*x+y*y);
             value[0] = -(k/m)*x/(r*r*r);
             value[1] = -(k/m)*y/(r*r*r);
             }

      };
double twoBody::G;
double twoBody::m1;
double twoBody::m2;
double twoBody::m;

class RKSolver
     {public:
      double **u;
      double t_range;
      double dt;
      int Nt;
 
      RKSolver(double t_range, double dt)
            {this->t_range = t_range;
             this->dt = dt;
             
             Nt = ceil(t_range/dt)+ 1;
             u = new double*[Nt];
             for(int i=0; i<Nt; ++i)
                {u[i] = new double[5];
                 }             
             }

      ~RKSolver()
            {for(int i=0; i<Nt; ++i)
                 delete[] u[i];
             delete[] u;
             }

      void setInitialVal(long double t, long double x0, long double y0, long double vx0, long double vy0)
            {u[0][0] = t;
             u[0][1] = x0;
             u[0][2] = y0;
             u[0][3] = vx0;
             u[0][4] = vy0;
             }
      
      void nextPosition(void (*f)(double x, double y, double arr[]), int index)
            {double a[2];
             double x = u[index-1][1], y=u[index-1][2], vx=u[index-1][3], vy=u[index-1][4];
             
             f(x, y, a);
             double G1x = a[0];
             double G1y = a[1];
             double F1x = vx;
             double F1y = vy;

             f(x+ F1x*dt/2, y+ F1y*dt/2, a);
             double G2x = a[0];
             double G2y = a[1];
             double F2x = vx + G1x*dt/2;
             double F2y = vy + G1y*dt/2;

             f(x+ F2x*dt/2, y+ F2y*dt/2, a);
             double G3x = a[0];
             double G3y = a[1];
             double F3x = vx + G2x*dt/2;
             double F3y = vy + G2y*dt/2;

             f(x+ F3x*dt, y+ F3y*dt, a);
             double G4x = a[0];
             double G4y = a[1];
             double F4x = vx + G3x*dt;
             double F4y = vy + G3y*dt;
             
             u[index][0] = u[index-1][0] + dt;
             u[index][3] = vx + (G1x + 2*G2x + 2*G3x + G4x)*dt/6;
             u[index][4] = vy + (G1y + 2*G2y + 2*G3y + G4y)*dt/6;
             u[index][1] = x + (F1x + 2*F2x + 2*F3x + F4x)*dt/6;
             u[index][2] = y + (F1y + 2*F2y + 2*F3y + F4y)*dt/6;
             }
      };
//Numerical class ends

class leapFrog
     {public:
      double **u;
      double t_range;
      double dt;
      int Nt;
      
      leapFrog(double t_range, double dt)
            {this->t_range = t_range;
             this->dt = dt;
             
             Nt = ceil(t_range/dt)+ 1;
             u = new double*[Nt];
             for(int i=0; i<Nt; ++i)
                {u[i] = new double[5];
                 }             
             }

      ~leapFrog()
            {for(int i=0; i<Nt; ++i)
                 delete[] u[i];
             delete[] u;
             }

      void setInitialVal(double t, double x0, double y0, double vx0, double vy0)
            {u[0][0] = t;
             u[0][1] = x0;
             u[0][2] = y0;
             u[0][3] = vx0;
             u[0][4] = vy0;
             }
      
      void solver(void (*acc)(double x, double y, double arr[]))
            {double arr[2];
             u[0][0]=0;
             double vx_bhalf, vy_bhalf;
             
             for(int n=1; n<Nt; ++n)
                {u[n][0] = u[n-1][0] + dt;

                 acc(u[n-1][1], u[n-1][2], arr);
                 vx_bhalf = u[n-1][3] + arr[0]*dt/2;
                 vy_bhalf = u[n-1][4] + arr[1]*dt/2;

                 
                 u[n][1] = u[n-1][1] + vx_bhalf*dt;
                 u[n][2] = u[n-1][2] + vy_bhalf*dt;
 
                 acc(u[n][1], u[n][2], arr);
                 u[n][3] = vx_bhalf + arr[0]*dt/2;
                 u[n][4] = vy_bhalf + arr[1]*dt/2;
                 }
             }
      };
