#include "BestFit.h"

 BestFit::BestFit(bool box){


   // Initialize model
   //Supplies GiNaC expressions
   strcpy(this->name,"BestFit");
   int dim_sys = 3;
   // List of state variables
   symbol s("s"), i("i"), r("r");
   lst vars;
   vars = {s, i, r};

   // System's dynamics
   ex ds = s - (0.34*s*i)*0.1;				// susceptible
   ex di = i + (0.34*s*i - 0.05*i)*0.1;	// infected
   ex dr = r + 0.05*i*0.1;					// removed
   lst dyns;
   dyns = {ds,di,dr};

   this->vars = vars;
   this->dyns = dyns;


   // Init reach set
   Bundle *B;
   if(box){
     int num_dirs = 3;		// number of bundle directions
     int num_temps = 1;		// number of bundle templates

     // Directions matrix
     vector< double > Li (dim_sys,0); //Initialze with 0
     vector< vector< double > > L (num_dirs,Li);

     #if DEFAULT
     L[0][0] = 1;  // [1 0 0 ]^T
     L[1][1] = 1; // [0 1 0 ]^T
     L[2][2] = 1; // [0 0 1 ]^T
     #else
     L[0][0] = 1;
     L[1][1] = 1;
     L[2][2] = 1; L[2][0] = 1;
     #endif

     // Template matrix
     vector< int > Ti (dim_sys,0);
     vector< vector< int > > T (num_temps,Ti);
     T[0][0] = 0; T[0][1] = 1; T[0][2] = 2;


     // Offsets for the set of initial conditions
     vector< double > offp (num_dirs,0);
     vector< double > offm (num_dirs,0);
     offp[0] = 0.8; offm[0] = -0.79;
     offp[1] = 0.2; offm[1] = -0.19;
     offp[2] = 0.0001; offm[2] = -0.000099;

     B = new Bundle(L,offp,offm,T);
 }else{
   int num_dirs = 5;		// number of bundle directions
   int num_temps = 3;		// number of bundle templates

   // Directions matrix
   vector< double > Li (dim_sys,0);
   vector< vector< double > > L (num_dirs,Li);

   #if DEFAULT
   L[0][0] = 1;
   L[1][1] = 1;
   L[2][2] = 1;
   L[3][0] = 1; L[3][1] = 1;
   L[4][0] = 1; L[4][2] = 1;
   #else
   L[0][0] = 1;
   L[1][1] = 1;
   L[2][2] = 1;
   L[3][0] = 1; L[3][1] = -1;  L[3][2] = 1;
   L[4][0] = 1; L[4][2] = 1;   L[4][1] = 1;
   #endif

   // Template matrix
   vector< int > Ti (dim_sys,0);
   vector< vector< int > > T (num_temps,Ti);
   T[0][0] = 0; T[0][1] = 1; T[0][2] = 2;
   T[1][0] = 1; T[1][1] = 2; T[1][2] = 3;
   T[2][0] = 2; T[2][1] = 3; T[2][2] = 4;

   // Offsets for the set of initial conditions
   vector< double > offp (num_dirs,0);
   vector< double > offm (num_dirs,0);
   offp[0] = 0.8; offm[0] = -0.79;
   offp[1] = 0.2; offm[1] = -0.19;
   offp[2] = 0.0001; offm[2] = -0.000099;
   offp[3] = 1; offm[3] = 0;
   offp[4] = 1; offm[4] = 0;

   B = new Bundle(L,offp,offm,T);
 }

 this->reachSet = B;

}
