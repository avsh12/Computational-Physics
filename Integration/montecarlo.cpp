//Problem: Given a torus z*z + pow((sqrt(x*x + y*y) - 3), 2) and a box {1<x<4, -3<y<4, -1<z<1}, find the mass and centre of mass of the torus enclosed in the box

//Output:
//Fraction of the points inside the torus: 0.526223.
//Volume of the torus enclosed in the box = 22.1014.
//Total mass of the torus lying in the box is 438.066.
//The centre of mass of the region: X_cm = 2.36484, Y_cm = 0.100719, Z_cm = 0.722681.

#include<iostream>
#include<math.h>
#include<random>
#include<time.h>

//tells if the point is inside the torus
#define inTorus(x, y, z) (z*z + pow((sqrt(x*x + y*y) - 3), 2) < 1)

using namespace std;

//mass density of the torus
double density(double x, double y, double z)
      {return 2*exp(5*z);
       }

int main()
   {clock_t start, end;
    start = clock();

    random_device dev;
    mt19937 rng(dev());
    uniform_real_distribution<> dist6(0,1.0);

    double x_cm, y_cm, z_cm;
    double mass, volume_box, volume;
    double x_box, y_box, z_box;
    double x, y, z;
    int N, N_in;

    x_cm = y_cm = z_cm = 0;
    mass = 0;
    x_box = 3;
    y_box = 7;
    z_box = 2;
    volume_box = x_box*y_box*z_box;
    N = 10000000;
    N_in = 0;
    
    for(int i = 0; i < N; ++i)
       {x = 1 + 3*dist6(rng);
        y = -3 + 7*dist6(rng);
        z = -1 + 2*dist6(rng);

        if(inTorus(x,y,z))
          {mass += density(x,y,z);
           x_cm += x*density(x,y,z);
           y_cm += y*density(x,y,z);
           z_cm += z*density(x,y,z);

           N_in++;
           }
        }

    volume = N_in*volume_box/N;

    mass *= volume/N_in;
    x_cm *= volume/(mass*N_in);
    y_cm *= volume/(mass*N_in);
    z_cm *= volume/(mass*N_in);

    cout<<"Fraction of the points inside the torus: "<<double(N_in)/N<<"."<<endl;
    cout<<"Volume of the torus enclosed in the box = "<<volume<<"."<<endl;
    cout<<"Total mass of the torus lying in the box is "<<mass<<"."<<endl;
    cout<<"The centre of mass of the region: X_cm = "<<x_cm<<", Y_cm = "<<y_cm<<", Z_cm = "<<z_cm<<"."<<endl;
    
    end = clock();
    //program run time
    cout<<"Run time: "<<double(end-start)/double(CLOCKS_PER_SEC)<<endl;
    return 0;
    }
