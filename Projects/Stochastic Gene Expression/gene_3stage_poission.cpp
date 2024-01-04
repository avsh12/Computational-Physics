#include<iostream>
#include<fstream>
#include<math.h>
#include<time.h>

using namespace std;

/*
poisson distribution
k0 = 6, k1 = 2
a = 1, b = 40
gamma = 1
dt = 0.001
iter = 20000
*/

/*
near gaussian distribution
k0 = 6, k1 = 2
a = 40, b = 2
gamma = 1
dt = 0.001
iter = 20000
*/

/*
bimodal distribution
k0 = 0.6, k1 = 0.2
a = 4, b = 10
gamma = 1
dt = 0.001
iter = 20000
*/

//parameters used. Default for bimodal distribution
struct param
      {float k0 = 6, k1 = 2, d1 = 0.0005;
       float a = 1, b = 40;
       float gama = 1;
       };

class Distribution
     {public:
      int row, col;
      //probabilities, inactive case
      float **P;
      //probabilities, active case
      float **Q;
      //probabilities, inactive case
      float **P1;
      //probabilities, active case
      float **Q1;

      Distribution()
          {row = 100;
           col = 100;

           P = new float*[row];
           Q = new float*[row];
           P1 = new float*[row];
           Q1 = new float*[row];

           for(int i=0; i<row; ++i)
              {P[i] = new float[col];
               Q[i] = new float[col];
               P1[i] = new float[col];
               Q1[i] = new float[col];
               }
           }
      Distribution(int row, int col)
          {this->row = row;
           this->col = col;
           
           P = new float*[row];
           Q = new float*[row];
           P1 = new float*[row];
           Q1 = new float*[row];

           for(int i=0; i<row; ++i)
              {P[i] = new float[col];
               Q[i] = new float[col];
               P1[i] = new float[col];
               Q1[i] = new float[col];
               }
           }

       ~Distribution()
          {//free pointers
          for(int i=0; i<row; ++i)
             {delete[] P[i], Q[i], P1[i], Q1[i];
              }
    
          delete[] P, Q, P1, Q1;
          //pointers freed
           }
       };
      
class RungeKuttaO4:public Distribution
     {public:
      float k0, k1, d1;
      float a, b;
      float gama;
      float dt;

      RungeKuttaO4():Distribution()
          {;}
      RungeKuttaO4(param p):Distribution()
          {k0 = p.k0;
           k1 = p.k1;
           d1 = p.d1;
           a = p.a;
           b = p.b;
           gama = p.gama;
           }
      RungeKuttaO4(int row, int col, param p):Distribution(row, col)
          {k0 = p.k0;
           k1 = p.k1;
           d1 = p.d1;
           a = p.a;
           b = p.b;
           gama = p.gama;
           }
          
      void solver(float **P, float **Q, float **P1, float **Q1)
            {float f1=0, f2=0;
             float result[2];
             //set alias for row and col
             int &rmax = this->row;
             int &cmax = this->col;

             for(int row=0; row<rmax; ++row)
                {for(int col=0; col<cmax; ++col)
                    {   f1 = k1*Q[row][col];
                        f1 += (row == rmax-1) ? 0 : gama*(row+1) * P[row+1][col];
                        f1 += (col == cmax-1) ? 0 : (col+1) * P[row][col+1];
                        f1 += (col == 0) ? 0 : b*gama*row * P[row][col-1];
 
                        //f1 *= d1;                

                        float coeff = -(k0 + col + (b+1)*row*gama);
                        float y = P[row][col];
                        float F1 = dt * (f1 + coeff*y);
                        float F2 = dt * (f1 + coeff*(y + F1/2));
                        float F3 = dt * (f1 + coeff*(y + F2/2));
                        float F4 = dt * (f1 + coeff*(y + F3));

                        result[0] =  y + (F1 + 2*F2 + 2*F3 + F4)/6;

                        f2 = k0*P[row][col];
                        f2 += (row == rmax-1) ? 0 : gama*(row+1) * Q[row+1][col];
                        f2 += (col == cmax-1) ? 0 : (col+1) * Q[row][col+1];
      	                f2 += (row == 0) ? 0 : a * Q[row-1][col];
                        f2 += (col == 0) ? 0 : b*gama*row* Q[row][col-1];

                        //f2 *= d1;
                
            	        coeff = -(k1 + a + col + (b+1)*row*gama);
              	        y = Q[row][col];
                        F1 = dt * (f2 + coeff*y);
                        F2 = dt * (f2 + coeff*(y + F1/2));
                        F3 = dt * (f2 + coeff*(y + F2/2));
                        F4 = dt * (f2 + coeff*(y + F3));

                       result[1] =  y + (F1 + 2*F2 + 2*F3 + F4)/6;
                       
                    P1[row][col] = result[0];
                    Q1[row][col] = result[1];
                    }
                }
             }//solver ends
      };

int main()
   {clock_t start, end;
    start = clock();

    int row = 50, col = 200;
    param parameters;
    RungeKuttaO4 RK4(row, col, parameters);
    
    //Initialization
     for(int i=0; i<row; ++i)
        {for(int j=0; j<col; ++j)
            {RK4.P[i][j] = 0;
             RK4.Q[i][j] = 0;
             }
         }
     RK4.P[3][5] = 0.5;
     RK4.Q[3][5] = 0.5;
     //Initialized
    
    //set time step size
    RK4.dt = 0.0001;

    //solve for steady state solution
    for(int iter=1; iter<2e5; ++iter)
       {if(iter%2 == 1)
          {RK4.solver(RK4.P, RK4.Q, RK4.P1, RK4.Q1);
           }
        else
          {RK4.solver(RK4.P1, RK4.Q1, RK4.P, RK4.Q);
           }
        }
    
    fstream fout;
    fout.open("poission_data.txt", ios::out|ios::app);

    //stores protein distribution
    float protein_dist[col];
    float total_unnormalized_prob=0;

    //run through protein numbers 0 to col-1
    for(int i=0; i<col; ++i)
       {float unnormalized_prob=0;
        for(int j=0; j<row; ++j)
           {//sum over all mRNAs both for active and inactive cases
            unnormalized_prob += RK4.P1[j][i] + RK4.Q1[j][i];
            }
        //unnormalized probability of i proteins
        protein_dist[i] = unnormalized_prob;
        //sum unnormalized probabilities for each protein numbers
        total_unnormalized_prob += unnormalized_prob;
        }

    for(int i=0; i<col; ++i)
       {//write the probability after normalizing
        fout<<i<<"\t"<<protein_dist[i]/total_unnormalized_prob<<endl;
        }
fout<<endl<<endl;
cout<<total_unnormalized_prob;
    
    fout.close();
    end = clock();
    //program run time
    cout<<endl<<double(end-start)/double(CLOCKS_PER_SEC);
    return 0;
    }
