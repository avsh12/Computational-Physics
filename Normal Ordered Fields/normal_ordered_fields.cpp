#include<iostream>
#include<fstream>

using namespace std;

int tcount=0;

void cprint(int *a, int r, int n)
    {
     
     cout<<"(";
     for(int i=1; i<=r; ++i)
        {cout<<a[i-1]<<" ";
         //if(i!=n-1) cout<<",";
         }
         
     bool exists=false;
     int val, m=1;
         
     for(int i=r+1; i<=n; ++i)
        {for(; m<=n; ++m)
            {for(int j=1; j<=r; ++j)
        	    {if(a[j-1]==m)
        	       {exists = true;
        	       	break;
        	       	}
        	     }
        	     
        	 if(exists==false)
        	   {val = m++;
        	   	break;
        	   	}
        	   
        	 exists = false;
        	 }
        	 
         cout<<val<<" ";
         }
         
     cout<<")"<<"    ";
     }
     
void fprint(fstream &fout, int *a, int r, int n)
    {for(int i=1; i<=r; ++i)
        {fout<<"\\phi_{"<<a[i-1]<<"}^-";
         //if(i!=n-1) cout<<",";
         }
         
     bool exists=false;
     int val, m=1;
         
     for(int i=r+1; i<=n; ++i)
        {for(; m<=n; ++m)
            {for(int j=1; j<=r; ++j)
        	    {if(a[j-1]==m)
        	       {exists = true;
        	       	break;
        	       	}
        	     }
        	     
        	 if(exists==false)
        	   {val = m++;
        	   	break;
        	   	}
        	   
        	 exists = false;
        	 }
        	 
         fout<<"\\phi_{"<<val<<"}^+";
         }
     
     ++tcount;
     if(tcount==5)
       {fout<<" \\\\ ";
       	fout<<"\n& & ";
       	tcount=0;
       	}
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

void distinctFields(fstream &fout, int n, int r)
   {//int n = 10, r=4;
   	int ncomb=0;
    int *a = new int[r];
    
    for(int i=1; i<=r; ++i)
       {a[i-1]=i;
        }
        
    cout<<"Number of creation operators on the left is "<<r<<endl<<endl;
    fprint(fout, a,r, n);
    fout<<" + ";
    ncomb++;
    
    for(int i=r; i>=1;)
       {if(incrementAtPos(a, i, r, n))
          {fprint(fout, a,r, n);
           fout<<" + ";
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
       
    delete[] a;
    }

int main()
   {int nfields = 8, nminus;
   	
    //LaTeX code starts
    fstream fout;
    fout.open("normalorderedproduct.tex", ios::out|ios::trunc);
     
    fout<<"\\documentclass{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage[a3paper, left=5mm,right=5mm,top=1cm,bottom=1cm]{geometry}\n\\begin{document}\n";
     
    fout<<"\\begin{eqnarray*}\n N(";
    for(int i=1; i<=nfields; ++i)
       {fout<<"\\phi_{"<<i<<"}";
        }
    fout<<") &=& ";
        
    //LaTeX code ends
   	
    for(nminus=0; nminus<=nfields; ++nminus)
       {distinctFields(fout, nfields, nminus);
   	       	
   	if(nminus!=nfields)
   	  {fout<<" \\\\ ";
   	   fout<<"\n & & ";
   	   fout<<" \\\\ ";
   	   fout<<"\n& & ";
   	   }
   	 }
   	      	
     fout<<"\\end{eqnarray*}\n\\end{document}";

     //system("/bin/pdflatex /home/adarshv/Documents/Numerical/QFT/normalorderedproduct.tex");
     }
