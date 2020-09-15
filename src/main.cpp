/**
 *
 * @file main.cpp
 * main: This main file reproduces the experiments reported in "Sapo: Reachability Computation and Parameter Synthesis of Polynomial Dynamical Systems"
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 *
 */

#define PHOS 0
#define SIRF 1
#define LKF 0
#define QUADF 0
#define ROSSLERF 0

#include <stdio.h>
#include <iostream>

#include "Common.h"
#include "Bundle.h"
#include "Sapo.h"

#include "VanDerPol.h"
#include "Rossler.h"
#include "SIR.h"
#include "LotkaVolterra.h"
#include "Phosphorelay.h"
#include "Quadcopter.h"
#include "LL.h"

#include "SIRp.h"
#include "Influenza.h"
#include "Ebola.h"

using namespace std;

int main(int argc, char** argv) {

  // Sapo's options
  sapo_opt options;
  options.trans = 0;			 // Set transformation (0=OFO, 1=AFO)
  options.decomp = 0;			  // Template decomposition (0=no, 1=yes)
  //options.alpha = 0.5;		// Weight for bundle size/orthgonal proximity
  options.verbose = false;


  cout<<"TABLE 1"<<endl;
  // Load modles
  vector< Model* > reach_models;
  vector< int > reach_steps;
  //reach_models.push_back(new VanDerPol());      reach_steps.push_back(300);
  //reach_models.push_back(new Rossler());        reach_steps.push_back(300);
  //reach_models.push_back(new SIR(true));        reach_steps.push_back(300);
  //reach_models.push_back(new LotkaVolterra());  reach_steps.push_back(300);
  //reach_models.push_back(new Phosphorelay());   reach_steps.push_back(300);
  //reach_models.push_back(new Quadcopter());     reach_steps.push_back(300);
  //reach_models.push_back(new LL());               reach_steps.push_back(300);

  // Compute reach sets
  for(int i=0; i<reach_models.size(); i++){

    cout<<"Model: "<<reach_models[i]->getName()<<"\tReach steps: "<<reach_steps[i]<<"\t";

    Sapo *sapo = new Sapo(reach_models[i],options);
    Flowpipe* flowpipe = sapo->reach(reach_models[i]->getReachSet(),reach_steps[i]);
  }	// reachability analysis
  cout<<"\n";

/*
  cout<<"TABLE 2"<<endl;
  // Load models
  vector< Model* > synth_models;
  synth_models.push_back(new SIRp());
  synth_models.push_back(new Influenza());
  synth_models.push_back(new Ebola());

  // Synthesize parameters
  for(int i=0; i<synth_models.size(); i++){

    cout<<"Model: "<<synth_models[i]->getName()<<"\t";

    Sapo *sapo = new Sapo(synth_models[i], options);
    LinearSystemSet *synth_parameter_set = sapo->synthesize(synth_models[i]->getReachSet(),synth_models[i]->getParaSet(),synth_models[i]->getSpec());	// parameter synthesis
  }
  cout<<"\n";
*/
/*
  cout<<"FIGURE 5"<<endl;
  Model* ll = new LL();
  Sapo *sapo6 = new Sapo(ll, options);

  // Compute reach set with box template
  cout <<"Model: "<< ll->getName() <<"\tReach steps: 300\t";
  Flowpipe* flowpipeC = sapo6->reach(ll->getReachSet(), 300);

  // Generate matlab script to plot flowpipe
  char figC[] = "plotFigureLL.m";
  flowpipeC->plotProjToFile(3, 1, figC, 'b');
  // Set picture appearence
  ofstream matlab_script;
  matlab_script.open (figC, ios_base::app);
  matlab_script<<"xlabel('t');\n";
  matlab_script<<"ylabel('h');\n";
  //matlab_script<<"zlabel('r');\n";
  //matlab_script<<"axis([0 1 0 0.7 0 0.8]);\n";
  //matlab_script<<"view([74 23]);\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<figC<<" generated\n"<<endl;
*/
#if PHOS
  cout<<"FIGURE 5"<<endl;
  Model* ll = new Phosphorelay();
  Sapo *sapo6 = new Sapo(ll, options);

  // Compute reach set with box template
  cout <<"Model: "<< ll->getName() <<"\tReach steps: 300\t";
  Flowpipe* flowpipeC = sapo6->reach(ll->getReachSet(), 200);

  // Generate matlab script to plot flowpipe
  char figC[] = "plotPhos.m";
  flowpipeC->plotProjToFile(2, 1, figC, 'b');
  // Set picture appearence
  ofstream matlab_script;
  matlab_script.open (figC, ios_base::app);
  matlab_script<<"xlabel('t');\n";
  matlab_script<<"ylabel('x3');\n";
  //matlab_script<<"zlabel('r');\n";
  //matlab_script<<"axis([0 1 0 0.7 0 0.8]);\n";
  //matlab_script<<"view([74 23]);\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<figC<<" generated\n"<<endl;
#endif

  #if QUADF
    cout<<"FIGURE 5"<<endl;
    Model* ll = new Quadcopter();
    Sapo *sapo6 = new Sapo(ll, options);

    // Compute reach set with box template
    cout <<"Model: "<< ll->getName() <<"\tReach steps: 300\t";
    Flowpipe* flowpipeC = sapo6->reach(ll->getReachSet(), 3);

    // Generate matlab script to plot flowpipe 
    char figC[] = "plotQuad.m";
    flowpipeC->plotProjToFile(13, 1, figC, 'b');
    // Set picture appearence
    ofstream matlab_script;
    matlab_script.open (figC, ios_base::app);
    matlab_script<<"xlabel('t');\n";
    matlab_script<<"ylabel('hI');\n";
    //matlab_script<<"zlabel('r');\n";
    //matlab_script<<"axis([0 1 0 0.7 0 0.8]);\n";
    //matlab_script<<"view([74 23]);\n";
    matlab_script<<"grid on;";
    matlab_script.close();
    cout<<figC<<" generated\n"<<endl;
  #endif


#if LKF
  cout<<"FIGURE 5"<<endl;
  Model* lk = new LotkaVolterra();
  Sapo *sapo5 = new Sapo(lk,options);

  // Compute reach set with box template
  cout <<"Model: "<< lk->getName() <<"\tReach steps: 300\t";
  Flowpipe* flowpipeLK = sapo5->reach(lk->getReachSet(),300);

  // Generate matlab script to plot flowpipe
  char figLK[] = "plotFigureLK.m";
  flowpipeLK->plotProjToFile(0,1,figLK,'b');
  // Set picture appearence
  ofstream matlab_script;
  matlab_script.open (figLK, ios_base::app);
  matlab_script<<"xlabel('t');\n";
  matlab_script<<"ylabel('x3');\n";
  //matlab_script<<"zlabel('r');\n";
  //matlab_script<<"axis([0 1 0 0.7 0 0.8]);\n";
  //matlab_script<<"view([74 23]);\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<figLK<<" generated\n"<<endl;
#endif

#if ROSSLERF
  cout<<"FIGURE 5"<<endl;
  Model* lk = new Rossler();
  Sapo *sapo5 = new Sapo(lk,options);

  // Compute reach set with box template
  cout <<"Model: "<< lk->getName() <<"\tReach steps: 300\t";
  Flowpipe* flowpipeLK = sapo5->reach(lk->getReachSet(),300);

  // Generate matlab script to plot flowpipe
  char figLK[] = "plotFigureRoss.m";
  flowpipeLK->plotProjToFile(2,1,figLK,'b');
  // Set picture appearence
  ofstream matlab_script;
  matlab_script.open (figLK, ios_base::app);
  matlab_script<<"xlabel('t');\n";
  matlab_script<<"ylabel('z');\n";
  //matlab_script<<"axis([0 1 0 0.7 0 0.8]);\n";
  //matlab_script<<"view([74 23]);\n";
  //matlab_script<<"grid on;";
  matlab_script.close();
  cout<<figLK<<" generated\n"<<endl;
#endif

#if SIRF
  /*
  cout<<"FIGURE 3a"<<endl;
  Model* sir3a = new SIR(true);
  Sapo *sapo3a = new Sapo(sir3a,options);

  // Compute reach set with box template
  cout<<"Model: "<<sir3a->getName()<<"\tReach steps: 300\t";
  Flowpipe* flowpipe3a = sapo3a->reach(sir3a->getReachSet(),300);

  // Generate matlab script to plot flowpipe
  char fig3a[] = "plotSIR.m";
  flowpipe3a->plotProjToFile(1,1,fig3a,'b');
  ofstream matlab_script;
  // Set picture appearence
  matlab_script.open (fig3a, ios_base::app);
  matlab_script<<"xlabel('t');\n";
  matlab_script<<"ylabel('r');\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<fig3a<<" generated\n"<<endl;
  */

  cout<<"FIGURE 3b"<<endl;
  Model* sir3b = new SIR(false);
  Sapo *sapo3b = new Sapo(sir3b,options);

  // Compute reach set with bundle template
  cout<<"Model: "<<sir3b->getName()<<"\tReach steps: 300\t";
  Flowpipe* flowpipe3b = sapo3b->reach(sir3b->getReachSet(),3);

  // Generate matlab script to plot flowpipe
  char fig3b[] = "plotFigure3b.m";
  flowpipe3b->plotProjToFile(1,1,fig3b,'b');
  // Set picture appearence
  ofstream matlab_script;
  matlab_script.open (fig3b, ios_base::app);
  matlab_script<<"xlabel('t');\n";
  matlab_script<<"ylabel('i');\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<fig3b<<" generated\n"<<endl;
#endif

/*
  cout<<"FIGURE 4a"<<endl;
  Model *sirp = new SIRp();
  Sapo *sapo_synth = new Sapo(sirp,options);

  // Synthesize parameters
  cout<<"Model: "<<sirp->getName()<<"\t";
  LinearSystemSet *synth_parameter_set = sapo_synth->synthesize(sirp->getReachSet(),sirp->getParaSet(),sirp->getSpec());

  // Generate matlab script to plot the parameter set
  char fig4a[] = "plotFigure4a.m";
  sirp->getParaSet()->at(0)->plotRegionToFile(fig4a,'w');
  synth_parameter_set->at(0)->plotRegionToFile(fig4a,'k');
  // Set picture appearence
  matlab_script.open (fig4a, ios_base::app);
  matlab_script<<"xlabel('\\beta');\n";
  matlab_script<<"ylabel('\\gamma');\n";
  matlab_script<<"axis([0.17 0.21 0.045 0.065]);\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<fig4a<<" generated\n"<<endl;


  cout<<"FIGURE 4b"<<endl;
  cout<<"Model: "<<sirp->getName()<<"\t";

  // Compute reach set with contrained parameters
  Sapo *sapo_reach = new Sapo(sirp,options);
  Flowpipe* flowpipe = sapo_reach->reach(sirp->getReachSet(),synth_parameter_set->at(0),300);

  // Generate matlab script to plot the reach set
  char fig4b[] = "plotFigure4b.m";
  flowpipe->plotRegionToFile(fig4b,'w');
  // Set picture appearence
  matlab_script.open (fig4b, ios_base::app);
  matlab_script<<"xlabel('s');\n";
  matlab_script<<"ylabel('i');\n";
  matlab_script<<"zlabel('r');\n";
  matlab_script<<"axis([0 1 0 1 0 1]);\n";
  matlab_script<<"view([126 35]);\n";
  matlab_script<<"grid on;";
  matlab_script.close();
  cout<<fig4b<<" generated"<<endl;
*/
  exit(EXIT_SUCCESS);
}
