#include<iostream>
#include<time.h>
#include"metropolis.h"

using namespace std;

int main()
   {clock_t start, end;
    start = clock();
    
    int dims[] = {2, 4, 8, 16, 32, 64};
    param parameters[7];
    IsingModel model[7];

    for(int i=0; i<=0; ++i)
       {//Create lattice of the desired dimension
        model[i].ResizeLattice(dims[i], dims[i]);
        
        //create various files for the data
        //file for moments and other observables
        parameters[i].mTout.open("mtcurve" + to_string(dims[i]) + ".txt", ios::out|ios::trunc);
        //file for thermalization of the lattice with MC steps
        parameters[i].thermalization_out.open("thermalization" + to_string(dims[i]) + ".txt", ios::out|ios::trunc);
        //file for the configuration of the lattice with temperature
        parameters[i].lattice_config_out.open("lattice_config" + to_string(dims[i]) + ".txt", ios::out|ios::trunc);

        //set initial temperature
        model[i].temperature = 4;
        //assign a random configuration
        model[i].RandomAssignment();
        //update the lattice for its energy and magnetizxation at the present configuration
        model[i].UpdateLatticeEnergyAndMagnetization();

        parameters[i].mTout<<"Temperature\t <m>\t <|m|>\t <m^2>\t <|m|^3>\t <m^4>\t";
        parameters[i].mTout<<"Variance\t Susceptibility\t Binder Cumulant"<<"\t";
        parameters[i].mTout<<"<Energy>\t <Energy^2>\t Specific Heat\t"<<endl;

        //write the various quantities in the files present in the parameters structure
        model[i].WriteMomentsVsTemperature(parameters[i]);
        }

    end = clock();
    cout<<"\nRun time: "<<double(end-start)/double(CLOCKS_PER_SEC)<<endl;

    return 0;
    }

