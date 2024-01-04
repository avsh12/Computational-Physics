#include"lattice.h"

using namespace std;

class IsingModel: public Lattice
     {public:
      long double magnetization, m1, m2, m3, m4, E2;
      long double variance, susceptibility, binder_cumulant, specific_heat;  
      random_device rd;

      IsingModel(): Lattice()
         {magnetization = 0;
          }   

      IsingModel(int xdim, int ydim): Lattice(xdim, ydim)
         {magnetization = 0;
          }

      void UpdateLatticeEnergyAndMagnetization();
      void MinimizationForOneMCStep();
      void ThermalizeLattice(long double fluctuation, int MC_in_block, int comparison_attempt);
      void ThermalizeLattice(fstream &fout, long double fluctuation, int MC_in_block, int comparison_attempt);
      void WriteMomentsVsTemperature(param &parameters);
      };


void IsingModel::UpdateLatticeEnergyAndMagnetization()
     {int neighbour_spin_sum, neighbours[4];

      lattice_energy = 0;
      magnetization = 0;

      for(int y = 0; y < ycount; ++y)
         {for(int x=0; x < xcount; ++x)
             {GetNeighbours(x, y, neighbours);

              neighbour_spin_sum = lattice[neighbours[0]] + lattice[neighbours[1]] + lattice[neighbours[2]] + lattice[neighbours[3]];
 
              lattice_energy += (-lattice[x + y*stride]) * neighbour_spin_sum;
              magnetization += lattice[x + y*stride];
              }
           }

      lattice_energy = lattice_energy/(2*N);
      magnetization /= N;
      }

void IsingModel::MinimizationForOneMCStep()
    {long double deltabyT;

     shuffle(shuffled_lattice.begin(), shuffled_lattice.end(), rd);

     for(int i=0; i<N; ++i)
        {int y = shuffled_lattice[i]/stride;
         int x = shuffled_lattice[i] - y*stride;

         deltabyT = (long double)(EnergyChangeAtFlippingSite(x, y))/temperature;
                 
         if(deltabyT <= 0 || (deltabyT > 0 && prob(rng) < exp(-deltabyT)))
           {lattice[stride*y + x] = -lattice[stride*y + x];
            }
         }
      }

//Thermalize lattice 
//Tries to thermalize until change in magnetization does not exceed "fluctuation"
//after "MC_in_block" MC steps or "comparison_attempt" is reached
void IsingModel::ThermalizeLattice(long double fluctuation, int MC_in_block, int comparison_attempt)
     {long double mag1MC = magnetization;
      long double avg_mag, avg_mag2;
      int avg_count = 100;

      start:
      for(int i=1; i<=MC_in_block; ++i)
         {MinimizationForOneMCStep();
          UpdateLatticeEnergyAndMagnetization();
          }

      //calculate average magnetization using avg_count points
      avg_mag = avg_mag2 = 0;
      for(int i=1; i<=avg_count; ++i)
         {avg_mag += magnetization;
          avg_mag2 += magnetization*magnetization;
               
          MinimizationForOneMCStep();
          UpdateLatticeEnergyAndMagnetization();
          }
       avg_mag /= avg_count;
       avg_mag2 /= avg_count;
           
       //thermalize again if fluctuation is larger than threshold
       if(sqrt(avg_mag2 - avg_mag*avg_mag) > fluctuation && comparison_attempt--)
       //if((abs(mag1MC - avg_mag) > fluctuation) || comparison_attempt--)
         {//mag1MC = avg_mag;
          goto start;
          }
               
       //extra thermalization to ensure thermalized state
       for(int i=1; i<=15*MC_in_block; ++i)
          {MinimizationForOneMCStep();
           UpdateLatticeEnergyAndMagnetization();
           }
       }

//also writes thermalization data to the given file
void IsingModel::ThermalizeLattice(fstream &fout, long double fluctuation, int MC_in_block, int comparison_attempt)
    {long double mag1MC = magnetization;
     long double avg_mag, avg_mag2;
     int avg_count = 100;
     int MC_count = 0;

     start:
     for(int i=1; i<=MC_in_block; ++i)
        {MinimizationForOneMCStep();
         UpdateLatticeEnergyAndMagnetization();

         fout<<++MC_count<<"\t"<<magnetization<<endl;
         }

      //calculate average magnetization using avg_count points
      avg_mag = avg_mag2 = 0;
      for(int i=1; i<=avg_count; ++i)
         {avg_mag += magnetization;
          avg_mag2 += magnetization*magnetization;

          MinimizationForOneMCStep();
          UpdateLatticeEnergyAndMagnetization();

          fout<<++MC_count<<"\t"<<magnetization<<endl;
          }
       avg_mag /= avg_count;
       avg_mag2 /= avg_count;
           
       if(sqrt(avg_mag2 - avg_mag*avg_mag) > fluctuation && comparison_attempt--)
       //if((abs(mag1MC - avg_mag) > fluctuation) || comparison_attempt--)
         {//mag1MC = avg_mag;
          goto start;
          }
               
       for(int i=1; i<=3*MC_in_block; ++i)
          {MinimizationForOneMCStep();
           UpdateLatticeEnergyAndMagnetization();

           fout<<++MC_count<<"\t"<<magnetization<<endl;
           }
        }

//void WriteMomentsVsTemperature()(fstream &mTout, fstream &thermalization_out, fstream &lattice_config, long double temp_start, long double temp_end, long double temp_step)
void IsingModel::WriteMomentsVsTemperature(param &parameters)
    {long double fluctuation = 0.01;
     int no_of_MC_blocks = 100;
     int MC_in_block = 10000;
     int total_MC = no_of_MC_blocks*MC_in_block;

     long double avg_magnetization, avg_abs_magt, avg_energy;
     int comparison_attempt = 8;

     //write lattice configuration before thermalization
   //  WriteLattice(parameters.lattice_config_out);
   //  parameters.lattice_config_out<<endl<<endl;

     for(temperature = parameters.temp_start; temperature >= parameters.temp_end; temperature -= parameters.temp_step)
        {avg_magnetization = avg_abs_magt = avg_energy = 0;
         m1 = m2 = m3 = m4 = E2 = 0;

         //parameters.temp_step = (temperature < 3)&&(temperature>1) ? 0.01 : 0.05;
         /*if(temperature < 2.34 && temperature > 2.24)
           {if(temperature <= 2.322 && temperature >= 2.258)
               parameters.temp_step = 0.001;
            else
               parameters.temp_step = 0.005;
            }*/
         //16*16
         /*if(temperature < 2.34 && temperature > 1.8)
           {parameters.temp_step = 0.001;
            }*/

         /*if(temperature < 2.4 && temperature > 2.1)
           {parameters.temp_step = 0.001;
            }*/        

       //  MC_in_block = (temperature < 3.5)&&(temperature>1.5) ? 1000 : 500;
       //  no_of_MC_blocks = (temperature < 3.5)&&(temperature>1.5) ? ((temperature < 2.4)&&(temperature>2.1) ? 150 : 100) : 60;
       //  total_MC = no_of_MC_blocks*MC_in_block;

         //thermalize the lattice
         if(temperature == parameters.temp_start)
           {ThermalizeLattice(parameters.thermalization_out, fluctuation, MC_in_block, comparison_attempt);
            }
         else
            ThermalizeLattice(fluctuation, MC_in_block, comparison_attempt);
            //lattice thermalized

            //First four magnetization moments about zero
            //and first two moments of lattice energy about zero
            for(int i=0; i<total_MC; ++i)
               {MinimizationForOneMCStep();
                UpdateLatticeEnergyAndMagnetization();

                avg_magnetization += magnetization;
                m1 = abs(magnetization);
                avg_abs_magt += m1;
                m2 += (m1*m1);
                m3 += (m1*m1*m1);
                m4 += (m1*m1*m1*m1);
                avg_energy += lattice_energy;
                E2 += lattice_energy*lattice_energy;
                }

          avg_energy /= total_MC;
          E2 /= total_MC;
          avg_magnetization /= total_MC;
          avg_abs_magt /= total_MC;
          m1 = avg_abs_magt;
          m2 /= total_MC;
          m3 /= total_MC;
          m4 /= total_MC;
          variance = (m2 - m1*m1);
          susceptibility = variance/temperature;
          binder_cumulant = 1 - m4/(3*m2*m2);
          specific_heat = (E2 - avg_energy*avg_energy)/(temperature*temperature);
          //Moments calculated

          parameters.mTout<<temperature<<"\t"<<avg_magnetization<<"\t"<<avg_abs_magt<<"\t";
          parameters.mTout<<m2<<"\t"<<m3<<"\t"<<m4<<"\t";
          parameters.mTout<<variance<<"\t"<<susceptibility<<"\t"<<binder_cumulant<<"\t";
          parameters.mTout<<avg_energy<<"\t"<<E2<<"\t"<<specific_heat<<"\t"<<endl;

          //write thermalized lattice configuration at temperature temp_start
       //   WriteLattice(parameters.lattice_config_out);
       //   parameters.lattice_config_out<<endl<<endl;               
          }
       //write thermalized lattice configuration at temperature temp_end
     //  WriteLattice(parameters.lattice_config_out);
       }

