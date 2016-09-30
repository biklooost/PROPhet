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
// This is the main neural network class. It inherits from
// Network.h so it follows that interface. 
// ####################################################################




#ifndef __NEURAL_NETWORK
#define __NEURAL_NETWORK

#include <iostream>
#include <vector>

#include "Parallel.h"
#include "Common.h"
#include "Network.h"
#include "NN_node.h"
#include "Functional_params.h"
#include "System.h"
#include "Optimizer.h"

using namespace std;


class Neural_network : public Network {
  
 public:
  
  Neural_network(const vector<System*> &systems, Functional_params network_details);
  ~Neural_network();
  
  virtual bool train();
  
  virtual void load(ifstream &input_file);

  virtual void print(ostream &output = std::cout);

  vector<REAL> evaluate();

  // MD-specific routines
  virtual REAL evaluate_MD(int i_sys);
  virtual REAL evaluate_MD(vector<REAL> &dE_dG);
  virtual REAL train_MD(int i_sys);
  virtual void scale_gradient(REAL scale_factor);

  virtual vector<REAL> gradient() { return my_dOutput_dParameters; }
  
  vector<REAL> get_parameters();
  void set_parameters(vector<REAL> new_params);

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

  virtual inline void set_output_mean(REAL new_mean) { 
    output_mean = new_mean; 
    params.output_mean = new_mean;
  }
  virtual inline void set_output_variance(REAL new_variance) { 
    output_variance = new_variance; 
    params.output_variance = new_variance;
  }
  
  
 private:
  
  void set_parameters();
  
  REAL SSE, SSE_mod;
  int last_i_sys;
  
  REAL lambda, lambda_max;
  REAL t;

  vector<REAL> phi; 
  REAL d_phi;

  vector<vector<REAL> > values;
  vector<vector<Neural_network_node*> > nodes;
  
  vector<REAL> my_dOutput_dParameters;
  vector<REAL> my_dOutput_dInputs;

  vector<int> NNodes;
  vector<int> deriv_offsets;
  
  const vector<System*> systems;
  
  REAL output_mean;
  REAL output_variance;

  int Nparams;
  int NNodes_max;

  Optimizer *optimizer;
  vector<REAL> results;
 
  inline vector<REAL> vector_sum(vector<REAL> in1, vector<REAL> in2) {
    for (int i=0; i<in1.size(); i++) {
      in1.at(i) += in2.at(i);
    }
    return in1;
  }
  

  inline int sign(REAL value) {
    return ((0 < value) - (value < 0));
  }


  Functional_params  params;

};


#endif
