//add the p_offset(m_offset) to the column(row) index to obtain number of proteins(mRNA).

#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<time.h>

using namespace std;
int count =0;
/*
poisson distribution
a = 0.5, b = 100
gama = 10

dt = 0.7
iter = 2.3e4
file: gamma10_poission_data.txt
*/

/*
near gaussian distribution
a = 20, b = 2.5
gama = 1

dt = 0.001
iter = 20000
*/

//parameters used. Default for bimodal distribution
struct param
      {float d1 = 0.0005;
       float a = 0.5, b = 100;
       float gama = 10;
       float d0 = gama*d1;
       float nu0 = a*d1, nu1 = b*d0;
       };

class Distribution
     {public:
      int row, col;
      int p_offset, m_offset;
      vector<vector<float>> P, P1;
      
      Distribution()
          {row = 100;
           col = 100;

           P.resize(row, vector<float>(col));
           P1.resize(row, vector<float>(col));
           }
      Distribution(int row, int col)
          {this->row = row;
           this->col = col;

           P.resize(row, vector<float>(col));
           P1.resize(row, vector<float>(col));
           }

       Distribution(int row, int col, int m_offset, int p_offset)
          {this->row = row;
           this->col = col;
           this->p_offset = p_offset;
           this->m_offset = m_offset;

           P.resize(row, vector<float>(col));
           P1.resize(row, vector<float>(col));
           }
       };
      
class RungeKuttaO4:public Distribution
     {public:
      float d0, d1;
      float nu0, nu1;
      float a, b;
      float gama;
      float dt;
      //indicators for resizing different rows and columns
      int trow_push, lrow_push;
      int lcol_push, rcol_push;
      //minimum probability at the boundary at which the matrix need to be resized
      float local_tolerance;

      RungeKuttaO4():Distribution()
          {;}
      RungeKuttaO4(param p):Distribution()
          {d0 = p.d0;
           d1 = p.d1;
           a = p.a;
           b = p.b;
           gama = p.gama;
           nu0 = p.nu0;
           nu1 = p.nu1;

           trow_push = 0;
           lrow_push = 0;
           lcol_push = 0;
           rcol_push = 0;

           local_tolerance = 1.0e-4;
           }
      RungeKuttaO4(int row, int col, param p):Distribution(row, col)
          {d0 = p.d0;
           d1 = p.d1;
           a = p.a;
           b = p.b;
           gama = p.gama;
           nu0 = p.nu0;
           nu1 = p.nu1;

           trow_push = 0;
           lrow_push = 0;
           lcol_push = 0;
           rcol_push = 0;
 
           local_tolerance = 1.0e-4;
           }

      RungeKuttaO4(int row, int col, int m_offset, int p_offset, param p):Distribution(row, col, m_offset, p_offset)
          {d0 = p.d0;
           d1 = p.d1;
           a = p.a;
           b = p.b;
           gama = p.gama;
           nu0 = p.nu0;
           nu1 = p.nu1;

           trow_push = 0;
           lrow_push = 0;
           lcol_push = 0;
           rcol_push = 0;

           local_tolerance = 5.0e-5;
           }

        //normalizes the matrix passed to it
        void normalize(vector<vector<float>> &P)
            {float total_prob = 0;

             for(vector<float> vect : P)
                {for(float prob : vect)
                     total_prob += prob;
                 }

             for(int i=0; i<P.size(); ++i)
                {for(int j=0; j<P[i].size(); ++j)
                    {P[i][j] = P[i][j] / total_prob;
                     }
                 }
             }
          
        void solver(vector<vector<float>> &P, vector<vector<float>> &P1)
            {float f1=0;
             //set alias for row and col
             int &rmax = this->row;
             int &cmax = this->col;
             //store the present dimension
             int temp_row = rmax;
             int temp_col = cmax;

             //decide if matrix is to be resized
             for(int col=0; col<temp_col; ++col)
                {if(m_offset != 0 && P[0][col]>local_tolerance)
                   {trow_push = 1;
                    }
                 if(P[temp_row-1][col]>local_tolerance)
                   {lrow_push = 1;
                    }
                 }
             for(int row=0; row<temp_row; ++row)
                {if(p_offset != 0 && P[row][0]>local_tolerance)
                   {lcol_push = 1;
                    }
                 if(P[row][temp_col-1]>local_tolerance)
                   {rcol_push = 1;
                    }
                 }
 
             //resize P1 matrix if needed
             if(trow_push||lrow_push||lcol_push||rcol_push)
               {//free P1
                P1.clear();

                //expand row from the top
                if(trow_push==1)
                  {int incr = (m_offset ==1) ? 1 : 2;
 
                   rmax += incr;
                   m_offset -= incr;
                   }

                //expand column from the left
                if(lcol_push==1)
                  {int incr = (p_offset ==1) ? 1 : 2;
 
                   cmax += incr;
                   p_offset -= incr;
                   }

                //expand row from the bottom
                if(lrow_push==1)
                  {//increment row by two
                   rmax += 2;
                   }

                //expand column from the right
                if(rcol_push==1)
                  {//increment col by two
                   cmax += 2;
                   }

                //create P1 with new dimension
                P1.resize(rmax, vector<float>(cmax));
                }

             for(int row=0, i, j; row<rmax; ++row)
                {for(int col=0; col<cmax; ++col)
                    {i = trow_push == 1 ? row -(rmax-temp_row) : row;
                     j = lcol_push == 1 ? col -(cmax-temp_col) : col;

                     f1 = (i >= temp_row-1)||(i<-1)||(j<0)||(j>temp_col-1) ? 0 : d0*(m_offset+i+1) * P[i+1][j];
                     f1 += (j >= temp_col-1)||(j<-1)||(i<0)||(i>temp_row-1) ? 0 : d1*(p_offset+j+1) * P[i][j+1];
                     f1 += (i <= 0)||(i > temp_row)||(j<0)||(j>temp_col-1) ? 0 : nu0 * P[i-1][j];
                     f1 += (j <= 0)||(j > temp_col)||(i<0)||(i>temp_row-1) ? 0 : nu1*(m_offset+i) * P[i][j-1];
                
                     float coeff = -(nu0 + (m_offset+i)*(nu1+d0) + (p_offset+j)*d1);
                     float y = (0<=i&&i<temp_row)&&(0<=j&&j<temp_col) ? P[i][j] : 0;

                     float F1 = dt * (f1 + coeff*y);
                     float F2 = dt * (f1 + coeff*(y + F1/2));
                     float F3 = dt * (f1 + coeff*(y + F2/2));
                     float F4 = dt * (f1 + coeff*(y + F3));

                     P1[row][col] =  y + (F1 + 2*F2 + 2*F3 + F4)/6;
                    }
                }

             //resize P if needed.
             if(temp_row != rmax || temp_col != cmax)
               {//free P
                P.clear();
  
                trow_push = 0;
                lrow_push = 0;
                lcol_push = 0;
                rcol_push = 0;

                P.resize(rmax, vector<float>(cmax));
                }
             }//solver ends
      };

float fileWrite(fstream &fout, RungeKuttaO4 RK4)
        {float total_prob = 0;
    	 float prob_protein = 0;
    	 float prob_mRNA = 0;

    	 fout<<"#index "<<count++<<". Protein Distribution"<<endl;
    	 //run through protein numbers 0 to col-1
    	 for(int i=0; i<RK4.col; ++i)
      	    {prob_protein = 0;
     	     for(int j=0; j<RK4.row; ++j)
     	        {//sum over all mRNAs
     	         prob_protein += RK4.P1[j][i];
     	         }
    	     fout<<i<<"\t"<<prob_protein<<endl;
    	     total_prob += prob_protein;
     	     }

   	 fout<<endl<<endl;
         fout<<"#index "<<count++<<". mRNA Distribution"<<endl;

   	 for(int i=0; i<RK4.row; ++i)
    	    {prob_mRNA = 0;
     	     for(int j=0; j<RK4.col; ++j)
      	        {//sum over all proteins
     	         prob_mRNA += RK4.P[i][j];
     	         }
     	     fout<<i<<"\t"<<prob_mRNA<<endl;
    	     }
         fout<<endl<<endl;
    	 return total_prob;
    	 }

int main()
   {clock_t start, end;
    start = clock();

    //initial distribution
    int protein = 3, mrna = 5;
    int p_offset = protein-2, m_offset = mrna-2;

    int row = 5, col = 5;
    param parameters;
    RungeKuttaO4 RK4(row, col, m_offset, p_offset, parameters);
    
    fstream fout;
    fout.open("gamma10_poission_data.txt", ios::out|ios::trunc);
   
    //Initialization
     for(int i=0; i<row; ++i)
        {for(int j=0; j<col; ++j)
            {RK4.P[i][j] = 0;
             }
         }
     RK4.P[mrna-m_offset][protein-p_offset] = 1;
     //Initialized
    
    //set time step size
    RK4.dt = 0.7;

    //solve for steady state solution
    for(int iter=1; iter<5e5; iter += 1)
       {if(iter%2 == 1)
          {RK4.solver(RK4.P, RK4.P1);
           }
        else
          {RK4.solver(RK4.P1, RK4.P);
           }
        }

    cout<<"Unnormalized probability: "<<fileWrite(fout, RK4)<<endl;
    RK4.normalize(RK4.P1);
    cout<<"Normalized probability: "<<fileWrite(fout, RK4)<<endl;

    fout.close();

    end = clock();
    //program run time
    cout<<"Run time: "<<double(end-start)/double(CLOCKS_PER_SEC)<<endl;
    return 0;
    }
