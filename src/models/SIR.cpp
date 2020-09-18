/**
 * @file SIR.cpp
 * SIR epidemic model
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#define DEFAULT 1
#define TWODTEST 0
#define THREEDTEST 0

#include "SIR.h"

 SIR::SIR(bool box){


   // Initialize model
   //Supplies GiNaC expressions
   strcpy(this->name,"SIR");
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
     offp[0] = 0.8; offm[0] = -0.75;
     offp[1] = 0.2; offm[1] = -0.15;
     offp[2] = 0.0001; offm[2] = -0.000099;

     B = new Bundle(L,offp,offm,T);
 }else{

   #if DEFAULT

   int num_dirs = 5;		// number of bundle directions
   int num_temps = 3;		// number of bundle templates

   // Directions matrix
   vector< double > Li (dim_sys,0);
   vector< vector< double > > L (num_dirs,Li);

   L[0][0] = 1;
   L[1][1] = 1;
   L[2][2] = 1;
   L[3][0] = 1; L[3][1] = 0.5;
   L[4][0] = 0.5; L[4][2] = 0.5;

   cout << "Multiple temps" << endl;

   #elif TWODTEST

   int num_dirs = 7;		// number of bundle directions
   int num_temps = 4;		// number of bundle templates

   // Directions matrix
   vector< double > Li (dim_sys,0);
   vector< vector< double > > L (num_dirs,Li);

   L[0][0] = 1;
   L[1][1] = 1;
   L[2][2] = 1;
   L[3][0] = 1; L[3][1] = 1;
   L[4][0] = 1; L[4][2] = 1;
   L[5][0] = 1; L[5][1] = -1;
   L[6][0] = 1; L[6][2] = -1; //2D-Rotations

   #elif THREEDTEST

   int num_dirs = 14;		// number of bundle directions
   int num_temps = 4;		// number of bundle templates

   // Directions matrix
   vector< double > Li (dim_sys,0);
   vector< vector< double > > L (num_dirs,Li);

   L[0][0] = 1;
   L[1][1] = 1;
   L[2][2] = 1;
   L[3][0] = 1; L[3][1] = 1;
   L[4][0] = 1; L[4][2] = 1;

   L[5][0] = 1; L[5][1] = 1;  L[5][2] = 1;
   L[6][0] = 1; L[6][1] = -1;   L[6][2] = 1;
   L[7][0] = 1; L[7][1] = 1;   L[7][2] = -1; //+x-invariant 3D rotations

   L[8][0] = -1; L[8][1] = 1;  L[8][2] = 1;
   L[9][0] = 1; L[9][1] = 1;   L[9][2] = 1;
   L[10][0] = 1; L[10][1] = 1;   L[10][2] = -1; //+y-invariant 3D rotations

   L[11][0] = -1; L[11][1] = 1;  L[11][2] = 1;
   L[12][0] = 1; L[12][1] = -1;   L[12][2] = 1;
   L[13][0] = 1; L[13][1] = 1;   L[13][2] = 1; //+z-invariant 3D rotations

   #endif

   #if DEFAULT

   // Template matrix
   vector< int > Ti (dim_sys,0);
   vector< vector< int > > T (num_temps,Ti);
   T[0][0] = 0; T[0][1] = 1; T[0][2] = 2;
   T[1][0] = 1; T[1][1] = 2; T[1][2] = 3;
   T[2][0] = 2; T[2][1] = 3; T[2][2] = 4;

   // Offsets for the set of initial conditions
   vector< double > offp (num_dirs,0);
   vector< double > offm (num_dirs,0);
   offp[0] = 0.8; offm[0] = -0.75;
   offp[1] = 0.2; offm[1] = -0.15;
   offp[2] = 0.0001; offm[2] = -0.000099;
   offp[3] = 1; offm[3] = 0;
   offp[4] = 1; offm[4] = 0;

   #elif TWODTEST

      // Template matrix
      vector< int > Ti (dim_sys,0);
      vector< vector< int > > T (num_temps,Ti);
      T[0][0] = 0; T[0][1] = 1; T[0][2] = 2;
      T[1][0] = 1; T[1][1] = 2; T[1][2] = 3;
      T[2][0] = 2; T[2][1] = 3; T[2][2] = 4;
      T[3][0] = 4; T[3][1] = 5; T[3][2] = 6;

      // Offsets for the set of initial conditions
      vector< double > offp (num_dirs,0);
      vector< double > offm (num_dirs,0);
      offp[0] = 0.8; offm[0] = -0.79;
      offp[1] = 0.2; offm[1] = -0.19;
      offp[2] = 0.0001; offm[2] = -0.000099;
      offp[3] = 1; offm[3] = 1;
      offp[4] = 1; offm[4] = 1;
      offp[5] = 1; offm[5] = 1;
      offp[6] = 1; offm[6] = 1;

   #elif THREEDTEST
   // Template matrix
   vector< int > Ti (dim_sys,0);
   vector< vector< int > > T (num_temps,Ti);
   T[0][0] = 0; T[0][1] = 1; T[0][2] = 2;
   //T[1][0] = 2; T[1][1] = 3; T[1][2] = 4;
   //T[2][0] = 1; T[2][1] = 2; T[2][2] = 3;

   T[1][0] = 5; T[1][1] = 6; T[1][2] = 7;

   T[2][0] = 8; T[2][1] = 9; T[2][2] = 10;

   T[3][0] = 11 ; T[3][1] = 12; T[3][2] = 13;
   // Offsets for the set of initial conditions
   vector< double > offp (num_dirs,0);
   vector< double > offm (num_dirs,0);
   offp[0] = 0.8; offm[0] = -0.79;
   offp[1] = 0.2; offm[1] = -0.19;
   offp[2] = 0.0001; offm[2] = -0.000099;
   offp[3] = 1; offm[3] = 1;
   offp[4] = 1; offm[4] = 1;

   offp[5] = 2.8; offm[5] = 2.8;
   offp[6] = 2.8; offm[6] = 2.8;
   offp[7] = 2.8; offm[7] = 2.8;

   offp[8] = 2.8; offm[8] = 2.8;
   offp[9] = 2.8; offm[9] = 2.8;
   offp[10] = 2.8; offm[10] = 2.8;

   offp[11] = 2.8; offm[11] = 2.8;
   offp[12] = 2.8; offm[12] = 2.8;
   offp[13] = 2.8; offm[13] = 2.8;
   #endif
   B = new Bundle(L,offp,offm,T);
 }

 this->reachSet = B;

}
