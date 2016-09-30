//     _____________________________________      _____   |    
//     ___/ __ \__/ __ \_/ __ \__/ __ \__/ /________/ /   |
//     __/ /_/ /_/ /_/ // / / /_/ /_/ /_/ __ \/ _ \/ __/  |
//     _/ ____/_/ _, _// /_/ /_/ ____/_/ / / /  __/ /_    |
//     /_/     /_/ |_| \____/ /_/     /_/ /_/\___/\__/    |
//---------------------------------------------------------

/*
  This file is part of the PROPhet code, which was written
  by Brian Kolb and Levi Lentz in the group of Alexie
  Kolpak at MIT.

  PROPhet is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.
  
  PROPhet is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with PROPhet.  If not, see <http://www.gnu.org/licenses/>.
*/


// ####################################################################
//                         CLASS DESCRIPTION
// ####################################################################
// This is the main class for network functions. As it inherits from 
// Network.h, it uses that interface.
//  <<<< THIS CLASS IS NOT FULLY IMPLEMENTED. >>>>
// ####################################################################






#ifndef __NETWORK_FUNCTION
#define __NETWORK_FUNCTION

#include <iostream>
#include <vector>
#include <complex>

#include "Parallel.h"
#include "Common.h"
#include "Network.h"
#include "Functional_params.h"
#include "System.h"
#include "Optimizer.h"
#include "NF_node.h"

using namespace std;


class Network_function : public Network {
  
 public:
  
  Network_function(const vector<System*> &systems, Functional_params network_details);
  ~Network_function();
  
  virtual bool train();
  
  virtual void load(ifstream &input_file);

  virtual void print(ostream &output = std::cout);

  vector<REAL> evaluate();


  // MD-specific routines
  virtual REAL evaluate_MD(int i_sys);
  virtual REAL evaluate_MD(vector<REAL> &dE_dG);
  virtual REAL train_MD(int i_sys);
  virtual void scale_gradient(REAL scale_factor);

  vector<complex<REAL> > get_parameters();
  void set_parameters(vector<complex<REAL> > new_params);

  virtual void attach_optimizer(Optimizer *opt);

  virtual vector<REAL> get_state();
  virtual void set_state(const vector<REAL>& new_state);

  virtual vector<REAL> count();
  virtual vector<REAL> mean();
  inline virtual void set_means(vector<REAL> means) { params.input_mean = means; }
  inline virtual void set_variances(vector<REAL> variances) { params.input_variance = variances; }
  inline virtual bool is_preconditioned() { return params.input_mean.size(); }
  
  virtual vector<REAL> variance(vector<REAL> means);
  virtual void Normalize(vector<REAL> means, vector<REAL> variances);
  inline virtual vector<REAL> get_targets() {
    vector<REAL> targets;
    for (int i=0; i<systems.size(); i++) {
      targets.push_back(systems[i]->properties.target());
    }
    return targets;
  }
  
  virtual void center_outputs(REAL bias_correction);

  virtual inline int Nparameters() { return Nparams; }
  virtual inline const vector<REAL>& get_gradient() { return my_dOutput_dParameters; }
  virtual void init_optimizer(Optimizer* opt);

  inline void set_output_mean(REAL new_mean) { output_mean = new_mean; }
  inline void set_output_variance(REAL new_variance) { output_variance = new_variance; }

 private:
  
  void set_parameters();
  
  REAL SSE;
  REAL SSE_mod;
  vector<REAL> phi; 

  REAL output_mean;
  REAL output_variance;

  vector<vector<complex<REAL> > > values;
  vector<vector<Network_function_node*> > nodes;
  
  vector<REAL> my_dOutput_dParameters;
  vector<complex<REAL> > my_dOutput_dInputs;

  vector<int> NNodes;
  vector<int> deriv_offsets;
  
  const vector<System*> systems;

  vector<REAL> my_means;
  vector<REAL> my_variances;
  
  Functional_params params;

  int Nparams;
  int NNodes_max;

  Optimizer *optimizer;
 
  inline vector<complex<REAL> > vector_sum(vector<complex<REAL> > in1, vector<complex<REAL> > in2) {
    for (int i=0; i<in1.size(); i++) {
      in1.at(i) += in2.at(i);
    }
    return in1;
  }

  REAL numeric_correction;
  vector<complex<REAL> > numeric_corr;

  inline int sign(REAL value) {
    return ((0 < value) - (value < 0));
  }

};


#endif

