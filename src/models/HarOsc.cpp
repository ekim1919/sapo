#include "HarOsc.h"

 HarOsc::HarOsc(){

   // System dimension (number of variables)
   strcpy(this->name,"HarOsc");

   int dim_sys = 2;
   symbol x("x"), y("y");
   lst vars, dyns;
   vars = {x,y};

   ex dx =  y;
   ex dy = -x;
   dyns = {dx, dy};

   this->vars = vars;
   this->dyns = dyns;

   ///// Parallelotope bundle for reachable set representation /////
   int num_dirs = 2;
   int num_temp = 1;

   // Directions matrix
   vector< double > Li (dim_sys,0);
   vector< vector< double > > L (num_dirs,Li);
   //initialize box 1
   for(int i=0; i<dim_sys; i++){
     L[i][i] = 1;
   }

   // Template matrix
   vector< int > Ti (dim_sys,0);
   vector< vector< int > > T (num_temp,Ti);
   for(int i=0; i<dim_sys; i++){
     T[0][i] = i;
   }
   vector< double > offp (num_dirs,0);
   vector< double > offm (num_dirs,0);

   // Offsets for the set of initial conditions
   offp[0] = 3; offm[0] = -2;		// h
   offp[1] = 1; offm[1] = 1;				// init quaternion		// init quaternion

   Bundle *B = new Bundle(L,offp,offm,T);
   this->reachSet = B;

 }
