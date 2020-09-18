
#include "OscPart.h"

 OscPart::OscPart(){

   // System dimension (number of variables)
   strcpy(this->name,"OscPart");

   int dim_sys = 3;
   symbol x("x"), y("y"), z("z");
   lst vars, dyns;
   vars = {x,y,z};

   ex dx =  -0.1 * x - y;
   ex dy =  x - 0.1 * y;
   ex dz = -0.15 * z;

   dyns = {dx, dy, dz};

   this->vars = vars;
   this->dyns = dyns;

   ///// Parallelotope bundle for reachable set representation /////
   int num_dirs = 3;
   int num_temp = 1;

   // Directions matrix
   vector< double > Li (dim_sys,0);
   vector< vector< double > > L (num_dirs,Li);
   //initialize box
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
   offp[0] = 0.1; offm[0] = 0.1;		// h
   offp[1] = 1; offm[1] = -0.8;				// init quaternion
   offp[2] = 1; offm[2] = -0.9;		// init quaternion

   Bundle *B = new Bundle(L,offp,offm,T);
   this->reachSet = B;

 }
