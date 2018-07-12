#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>
using namespace std;
using namespace arma;
#define pi 3.14159265
int main ()
{
char const* ch1="Ar";
char const* ch2="Ai"; 
char const* ch3="Br";
char const* ch4="Bi";
int    N=121;
double R=10000.0;
double K=1.0;
double  yC;
double   U;
double  dU;
double d2U;


mat cheb    = zeros<mat>(5,N)   ;// contains Chebyschev polynomials and derivs
cx_mat A    = zeros<cx_mat>(N,N);
cx_mat B    = zeros<cx_mat>(N,N);
cx_mat iota = zeros<cx_mat>(1,1);
iota=cx_double(0,1);

// Spectral collocation points yC=cos(m*pi/(N-1)) for m=0..N-1
for (int m=N-2;m>1;m--)
{
 yC = cos(m*pi/(N-1));

 
 // Cheb Poly @ yC

 cheb(0,0)=1.0;
 cheb(0,1)=yC;

 for (int jj=1; jj <N-1;jj++)
     {
      cheb(0,jj+1)=2.0*yC*cheb(0,jj)-cheb(0,jj-1);
     }

// Cheb derivs @ yC
 for (int ii=1; ii<5; ii++)
     {
      cheb(ii,0)=0.0;
      cheb(ii,1)=cheb(ii-1,0);
      cheb(ii,2)=4.0*cheb(ii-1,1);
        
      for (int jj=3;jj<N;jj++)
          {
           cheb(ii,jj)=2.0*jj*cheb(ii-1,jj-1)+(jj)/(jj-2.0)*cheb(ii,jj-2);
          }
     }

// Base flow @ yC
U   = 1.0-pow(yC,2);
dU  =-2.0*yC;
d2U =-2.0;
// A and B of Aa=cBa
for (int jj=0;jj<N;jj++)
    {
      A(N-m-1,jj)=U*(cheb(2,jj)-(K*K)*cheb(0,jj))-d2U*cheb(0,jj)+(iota(0,0)/(K*R))*(cheb(4,jj)-2.0*K*K*cheb(2,jj)+(pow(K,4))*cheb(0,jj));
      B(N-m-1,jj)=   cheb(2,jj)-(K*K)*cheb(0,jj);
    }
}
// BCs
for (int j=0;j<N;j++ )
{
A(0,j)  = pow(+1,j);
A(1,j)  = pow(+1,j-1)*pow(j,2);
A(N-2,j)= pow(-1,j-1)*pow(j,2);
A(N-1,j)= pow(-1,j);
B(0,j)  = 0.0;
B(1,j)  = 0.0;
B(N-2,j)= 0.0;
B(N-1,j)= 0.0;
}
/*
cx_vec eigval = eig_pair( A, B );
cout<<eigval<<endl;
*/
char fname1[50];
sprintf(fname1, "/home/vikas/Desktop/hw_5430/%s",ch1);
ofstream fout1;
fout1.open(fname1);

char fname2[50];
sprintf(fname2, "/home/vikas/Desktop/hw_5430/%s",ch2);
ofstream fout2;
fout2.open(fname2);

char fname3[50];
sprintf(fname3, "/home/vikas/Desktop/hw_5430/%s",ch3);
ofstream fout3;
fout3.open(fname3);

char fname4[50];
sprintf(fname4, "/home/vikas/Desktop/hw_5430/%s",ch4);
ofstream fout4;
fout4.open(fname4);


 for(int i = 0; i < N; i++)//row wise
    {
     for(int j = 0; j < N; j++)
         {

          fout1.precision(18); 
          fout1 << real(A(i,j)); 
          fout1 << "\n";  

          fout2.precision(18); 
          fout2 << imag(A(i,j)); 
          fout2 << "\n";      

          fout3.precision(18); 
          fout3 << real(B(i,j));; 
          fout3 << "\n";      

          fout4.precision(18); 
          fout4 << imag(B(i,j)); 
          fout4 << "\n";      

         }
       
    }

fout1.close();
fout2.close();
fout3.close();
fout4.close();


return 0;

}
