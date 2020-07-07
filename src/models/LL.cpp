#include "LL.h"

LL::LL(){

strcpy(this->name,"LL");

int dim_sys = 7;
// List of state variables
symbol x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5"), x6("x6"), x7("x7");

lst vars;
vars = {x1, x2, x3, x4, x5, x6, x7};

float delta = 0.05;

ex dx1 =  x1 + (1.4*x3 - 0.9*x1)*delta;
ex dx2 = x2 + (2.5*x5 - 1.5*x2)*delta;
ex dx3 = x3 + (0.6*x7 - 0.8*x2*x3)*delta;
ex dx4 = x4 + (2 - 1.3*x3*x4)*delta;
ex dx5 = x5 + (0.7*x1 - x4*x5)*delta;
ex dx6 = x6 + (0.3*x1 - 3.1*x6)*delta;
ex dx7 = x7 + (1.8*x6 - 1.5*x2*x7)*delta;

lst dyns;
dyns = {dx1, dx2, dx3, dx4, dx5, dx6, dx7};

this->vars = vars;
this->dyns = dyns;

int num_dirs = 8;
int num_temp = 2;

vector< double > Li (dim_sys,0);
vector< vector< double > > L (num_dirs,Li);

for(int i = 0; i < dim_sys; i++) {
      L[i][i] = 1;
}

L[7][2] = 1; L[7][3] = 1;
//L[8][3] = 1; L[8][4] = 1;

// Template matrix
vector< int > Ti (dim_sys,0);
vector< vector< int > > T (num_temp,Ti);

for(int i=0; i<dim_sys; i++) {
  T[0][i] = i;
}

for(int i=0; i<dim_sys; i++) {
  T[1][i] = i+1;
}

vector< double > offu (num_dirs,0);
vector< double > offl (num_dirs,0);

offu[0] = 1.21; offl[0] = -1.19;
offu[1] = 1.06; offl[1] = -1.04;
offu[2] = 1.51; offl[2] = -1.49;
offu[3] = 2.41; offl[3] = -2.39;
offu[4] = 1.01; offl[4] = -0.99;
offu[5] = 0.11; offl[5] = -0.09;
offu[6] = 0.46; offl[6] = -0.44;
offu[7] = 10; offl[7] = 10;
//offu[8] = 10; offl[8] = 10;

Bundle *B = new Bundle(L,offu,offl,T);
this->reachSet = B;

}
