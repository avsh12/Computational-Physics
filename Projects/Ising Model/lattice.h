#include<fstream>
#include<vector>
#include<random>
#include<algorithm>

using namespace std;

random_device dev;
mt19937 rng(dev());
uniform_real_distribution<> prob(0, 1);

//Prescription for 2d to 1d index conversion
//(x, y) --> y*stride + x

struct param
   {fstream mTout;
    fstream thermalization_out;
    fstream lattice_config_out;
    float temp_start = 4;
    float temp_end = 0;
    float temp_step = 0.1;
    };

class Lattice
     {public:
      vector<int> lattice, shuffled_lattice;
      int xcount, ycount, N;
      int stride, MC_step;

      double temperature;
      float lattice_energy;
      
      Lattice(int xdim, int ydim)
         {xcount = xdim;
          ycount = ydim;
          N = xdim*ydim;
          stride = xdim;
          MC_step = N;

          lattice.resize(N);
          shuffled_lattice.resize(N);
          iota(shuffled_lattice.begin(), shuffled_lattice.end(), 0);
          }

      Lattice()
         {xcount = ycount = 0;
          N = stride = MC_step = 0;
          }

      void ResizeLattice(int xdim, int ydim);
      void RandomAssignment();
      void OrderedAssignment();
      void OrderedAssignment(int s);
      void GetNeighbours(int x, int y, int* neighbours);
      void WriteLattice(fstream &fout);
      int EnergyChangeAtFlippingSite(int x, int y);
      };
      
void Lattice::ResizeLattice(int xdim, int ydim)
    {xcount = xdim;
     ycount = ydim;
     N = xdim*ydim;
     stride = xdim;
     MC_step = N;

     lattice.resize(N);
     shuffled_lattice.resize(N);
     iota(shuffled_lattice.begin(), shuffled_lattice.end(), 0);
     }

void Lattice::RandomAssignment()
    {for(int i=0; i<N; ++i)
        {lattice[i] = prob(rng)>0.5 ? 1 : -1;
         }
     }

void Lattice::OrderedAssignment()
    {fill(lattice.begin(), lattice.end(), 1);
     }

void Lattice::OrderedAssignment(int s)
    {fill(lattice.begin(), lattice.end(), s);
     }

//spin over boundaries
//boundary 1: y=0, x=0 to N-1
//boundary 2: x=N-1, y=0 to N-1
//boundary 3: y=N-1, x=0 to N-1
//boundary 4: x=0, y=0 to N-1
void Lattice::GetNeighbours(int x, int y, int* neighbours)
    {neighbours[0] = x + (y-1)*stride;
     neighbours[1] = x+1 + y*stride;
     neighbours[2] = x + (y+1)*stride;
     neighbours[3] = x-1 + y*stride;

     //apply periodic boundary condition
     if(y==0)
        neighbours[0] = (ycount-1)*stride + x;
     if(x==xcount-1)
        neighbours[1] = y*stride;
     if(y==ycount-1)
        neighbours[2] = x;
     if(x==0)
        neighbours[3] = y*stride + xcount-1;
     }

void Lattice::WriteLattice(fstream &fout)
    {for(int i=0; i<ycount; ++i)
        {for(int j=0; j<xcount; ++j)
            {fout<<j<<"\t"<<i<<"\t"<<lattice[i*stride + j]<<endl;
             }
         }
     }

int Lattice::EnergyChangeAtFlippingSite(int x, int y)
   {int neighbours[4];
    int neighbour_spin_sum=0;

    //get neighbours
    GetNeighbours(x, y, neighbours);
  
    for(int i=0; i<4; ++i)
       {neighbour_spin_sum += lattice[neighbours[i]];
        }
          
    int spin_at_xy = lattice[x + y*stride];
    //int flipped_spin_at_xy = -spin_at_xy;
          
    //return -(flipped_spin_at_xy - spin_at_xy)*neighbour_spin_sum;
    return 2*spin_at_xy*neighbour_spin_sum;
    }
      
