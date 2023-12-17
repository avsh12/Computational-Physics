#include<iomanip>
#include<math.h>
#include<fstream>

using namespace std;

class Matrix
     {public:
      float xrange, trange;
      float dx, dt;
      int Nx, Nt;
      float **u;
      
      Matrix()                //default constructor
         {xrange = 20;
          trange = 50;
          dx = 0.1;
          u=NULL;
          }

      Matrix(float xmax, float tmax, float x_step)//, float s_param)    
         {xrange = xmax;
          trange = tmax;
          dx = x_step;
          u=NULL;
          }

      ~Matrix()
         {if(u!=NULL)                 //delete pointer if it's occupied
            {for(int n=0; n<Nt; ++n)
                 delete[] u[n];
             delete[] u;
             }
          }

      void xStep(float step)     //set x step
          {dx = step;}
      void tStep(float step)     //set t step
          {dt = step;}
      float xStep()              //return x step
          {return dx;}
      float tStep()              //return t step
          {return dt;}
      void xRange(float xrange)
          {this->xrange = xrange;}
      float xRange()
          {return xrange;}
      void tRange(float trange)
          {this->trange = trange;}
      float tRange()
          {return trange;}

      void updateMatrix()
           {Nx = floor(xrange/dx)+1;    //Number of x steps
            Nt = floor(trange/dt);      //Number of t steps 

            if(u!=NULL)                 //delete pointer if it's occupied
              {for(int n=0; n<Nt; ++n)
                   delete[] u[n];
               delete[] u;
               u=NULL;
               }
            u = new float*[Nt];
            for(int i=0;i<Nt;++i)
               {u[i] = new float[Nx];
                }
            }

      };

class MatrixOperation
     {public:
      static float Nx, Nt, dx, dt;


      MatrixOperation(float Nx, float Nt, float dx, float dt)
             {this->Nx = Nx;
              this->Nt = Nt;
              this->dx = dx;
              this->dt = dt;
              }
      
      //temporal boundary conditions as inputs
      static void tBoundary(float boundary, int t_index, float **u)    //at a given t_index boundary value is constant
          {for(int j=0; j<Nx; ++j)
              {u[t_index][j] = boundary;
               }
           }

      static void tBoundary(float (*boundary)(float x_value), int t_index, float **u)    //at a given t_index boundary value depends on x
          {for(int j=0; j<Nx; ++j)
              {u[t_index][j] = boundary(j*dx);
               }
           }

      static void tBoundary(float (*boundary)(float x_value, int t_value), int t_index, float **u)  //at a given t_index boundary value depends on x and t
          {for(int j=0; j<Nx; ++j)
              {u[t_index][j] = boundary(j*dx, t_index*dt);
               }
           }

      static void tBoundary(float boundary, int xstart_index, int xlast_index, int t_index, float **u)      //at a given t_index boundary value is constant. Boundary is set for given range of x
          {for(int j=xstart_index; j<=xlast_index; ++j)
              {u[t_index][j] = boundary;
               }
           }

      static void tBoundary(float (*boundary)(float x_value), int xstart_index, int xlast_index, int t_index, float **u)      //at a given t_index boundary value depends on x. Boundary is set for given range of x
          {for(int j=xstart_index; j<=xlast_index; ++j)
              {u[t_index][j] = boundary(j*dx);
               }
           }

      static void tBoundary(float (*boundary)(float x_value, int t_value), int xstart_index, int xlast_index, int t_index, float **u)   //at a given t_index boundary value depends on x and t. Boundary is set for given range of x
          {for(int j=xstart_index; j<=xlast_index; ++j)
              {u[t_index][j] = boundary(j*dx, t_index*dt);
               }
           }

      //spatial boundary conditions as inputs
      static void xBoundary(float boundary, int x_index, float **u)      //at a given x_index boundary value is constant
          {for(int j=0; j<Nt; ++j)
              {u[j][x_index] = boundary;
               }
           }

      static void xBoundary(float (*boundary)(float x_value, int t_value), int x_index, float **u)    //at a given x_index boundary value depends on x and t
          {for(int j=0; j<Nt; ++j)
              {u[j][x_index] = boundary(x_index*dx, j*dt);
               }
           }

      static void xBoundary(float (*boundary)(float t_value), int x_index, float **u)      //at a given x_index boundary value depends on t
          {for(int j=0; j<Nt; ++j)
              {u[j][x_index] = boundary(j*dt);
               }
           }

      static void xBoundary(float boundary, int tstart_index, int tlast_index, int x_index, float **u)      //at a given x_index boundary value is constant. Boundary is set for given range of t
          {for(int j=tstart_index; j<=tlast_index; ++j)
              {u[j][x_index] = boundary;
               }
           }

      static void xBoundary(float (*boundary)(int t_value), int tstart_index, int tlast_index, int x_index, float **u)      //at a given x_index boundary value depends on t. Boundary is set for given range of t
          {for(int j=tstart_index; j<=tlast_index; ++j)
              {u[j][x_index] = boundary(j*dt);
               }
           }

      static void xBoundary(float (*boundary)(float x_value, int t_value), int tstart_index, int tlast_index, int x_index, float **u)     //at a given x_index boundary value depends on x and t. Boundary is set for given range of t
          {for(int j=tstart_index; j<=tlast_index; ++j)
              {u[j][x_index] = boundary(x_index*dx, j*dt);
               }
           }
      //boundary condition as inputs end

      static void setMatrix(float (*function)(float **matrix, int x_index, int t_index, int row, int col), int xstart_index, int xlast_index, int tstart_index, int tlast_index, float **u)               //function
          {for(int n=tstart_index; n<=tlast_index; ++n)
              {for(int i=xstart_index; i<=xlast_index; ++i)
                  {u[n][i] = function(u, i, n, Nt, Nx);
                   }
               }
           }//setMatrix function ends

       static void fileOutput(string file_name, float **u)  //writes function values to a file 
           {fstream fout;
            fout.open(file_name,ios::out|ios::trunc);
             
            for(int i=0; i<Nx; ++i)
               {for(int n=0; n<Nt; n+=1)
                   {fout<<fixed<<setprecision(2)<<i*dx<<"\t";
                    fout<<fixed<<setprecision(5)<<n*dt<<"\t";
                    fout<<fixed<<setprecision(5)<<u[n][i]<<"\n";
                    }
                fout<<endl;
                }
            fout.close();
            }//file_output function ends
       };//PDE_Parabolic1 class ends

float MatrixOperation::Nx, MatrixOperation::Nt, MatrixOperation::dx, MatrixOperation::dt;

class Hyperbolic
     {public:
      static float c, s;
      
      Hyperbolic(float s, float c)
           {this->s = s;
            this->c = c;
            }      

      Hyperbolic()
           {s = 0.9;
            c = 1;
            }
       
      void setS(float s)                //set s parameter
               {this->s =s;}
      float getS()
               {return s;}
      void setC(float c)
               {this->c = c;}
      float getC()
               {return c;}

      static float functionDef(float **matrix, int x_index, int t_index, int row, int col)               //function
               {return s*(matrix[t_index-1][x_index+1] + matrix[t_index-1][x_index-1]) + 2*(1-s)*matrix[t_index-1][x_index]-matrix[t_index-2][x_index];
                }//function ends

      };//Hyperbolic class ends

float Hyperbolic::c;
float Hyperbolic::s;

class Parabolic
     {public:
      static float s;

      Parabolic(float s)
            {this->s = s;}
      Parabolic()
            {s = 0.48;}

      float getS()
            {return s;}

      void setS(float s)
            {this->s = s;}

      static float functionDef(float **matrix, int x_index, int t_index, int row, int col)
            {return s*(matrix[t_index-1][x_index+1] + matrix[t_index-1][x_index-1]) + (1-2*s)*matrix[t_index-1][x_index];
             }
      };

float Parabolic::s;

class Advection
     {public:
      static float s;

      Advection(float s)
            {this->s = s;}
      Advection()
            {s = 0.9;}

      float getS()
            {return s;}

      void setS(float s)
            {this->s = s;}

      static float functionDef(float **matrix, int x_index, int t_index, int row, int col)
            {return -s*(matrix[t_index-1][x_index] - matrix[t_index-1][x_index-1]) + matrix[t_index-1][x_index];
             }
      };

float Advection::s;
