#include <iostream>
#include <fstream>
#include <map>
#include <iomanip>
#include <cstring>


#include "Functional.h"
#include "Common.h"
#include "System.h"
#include "Setup.h"
#include "Potential.h"


using namespace std;




int main(int argc, char *argv[]) {
  //  feenableexcept(FE_INVALID | FE_OVERFLOW);  
  string input_file ="input_file";
  string run_type = "train";

  for (int i=1; i<argc; i++) {
    
    if (!strcmp(argv[i],"-in") || !strcmp(argv[i], "-input")) {
      if (argc > i+1) {
	input_file = argv[++i];
      }
    } else if (!strcmp(argv[i],"-init")) {
      
      run_type = "initialize";
    } else if (!strcmp(argv[i],"-run")) {
      run_type = "run";
    } else if (!strcmp(argv[i],"-train")) {
      run_type = "train";
    } else if (!strcmp(argv[i],"-validate")) {
      run_type = "validate";
    } else if (!strcmp(argv[i],"-optimize")) {
      run_type = "optimize";
    } else if (!strcmp(argv[i],"-forces")) {
        run_type = "forces";
    } else {
      ERROR("Command line argument \""+(string)argv[i]+"\" not recognized");
    }

  }
  
  Setup run_params(input_file);
  
  vector<System*> systems(run_params.Nsystems());
  
  for (int i=0; i<run_params.Nsystems(); i++) {
    systems.at(i) = new System(run_params.system(i), &run_params.F);
  }
  
  
  if (run_params.is_potential()) {
    
    Potential pot(systems, run_params.F);
    
    if (run_type == "train") {
      pot.train();
    } else if (run_type == "run") {
      pot.evaluate();
    } else if (run_type == "validate") {
      pot.validate();
    } else if (run_type == "optimize") {
      pot.optimize_Gs();
    } else if (run_type == "forces") {
        pot.forces();
    }

  } else {
  
    Functional F(systems, run_params.F);
  
    if (run_type == "train") {
      F.train();
    } else if (run_type == "run") {
      F.evaluate();
    } else if (run_type == "validate") {
      F.validate();
    } 
    
  }
  
  
  for (int i=0; i<run_params.Nsystems(); i++) {
    delete systems.at(i);
  }

  
}



