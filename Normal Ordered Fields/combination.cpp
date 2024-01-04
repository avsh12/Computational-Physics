#include<iostream>

using namespace std;

void cprint(int *a, int n)
    {cout<<"(";
     for(int i=0; i<n; ++i)
        {cout<<a[i]<<" ";
         //if(i!=n-1) cout<<",";
         }
     cout<<")"<<"    ";
     }
     
bool incrementAtPos(int *a, int pos, int r, int n)
   {if(pos==r)
      {if(a[pos-1] < n)
         {++a[pos-1];
          return true;
         }
       else return false;
       }
    else 
      {if(r-pos+a[pos-1] < n)
         {++a[pos-1];
          for(int i=pos+1; i<=r; ++i)
             {a[i-1] = a[i-2]+1;
              }
          return true;
          }
        else return false;
       }
    }

int main()
   {int n = 6, r=4;
    int ncomb=0;
    int *a = new int[r];
    
    for(int i=1; i<=r; ++i)
       {a[i-1]=i;
        }
    cprint(a,r);
    ncomb++;
    
    for(int i=r; i>=1;)
       {if(incrementAtPos(a, i, r, n))
          {cprint(a,r);
           ncomb++;
           if(i<r)
              i=r;
           }
        else
          {if(i>1)
             {--i;
              }
           else 
              break;
           }
        }
        
    cout<<endl<<endl<<"Nunber of possible values at each entry is "<<n<<endl;
    cout<<"Number of entries is "<<r<<endl;
    cout<<"Number of distinct combinations generated = "<<ncomb;
    
    delete[] a;
    }
