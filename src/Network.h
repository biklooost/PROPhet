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
//                           CLASS DESCRIPTION
// ####################################################################
// This is an interface class to machine-learning algorithms
// (currently neural networks and network functions). It is intended
// to abstract away details of a particular machine-learning approach.
// ####################################################################







#ifndef __NETWORK
#define __NETWORK

#include <iostream>
#include <fstream>

#include "Network_node.h"
#include "Functional_params.h"
#include "Common.h"
#include "System.h"
#include "Optimizer.h"

using namespace std;


class Network {
  
  
 public:
  
  Network() {
    
    

  }
    
  virtual ~Network() {

  }

  virtual bool train() = 0;
  virtual vector<REAL> evaluate() = 0;


  // MD specific routines
  virtual REAL evaluate_MD(int i_sys) = 0;
  virtual REAL evaluate_MD(vector<REAL> &dE_dG) = 0;
  virtual REAL train_MD(int i_sys) = 0;
  virtual void scale_gradient(REAL scale_factor) = 0;

  virtual void load(ifstream &input_file) = 0;

  virtual int Nparameters() = 0;
  virtual const vector<REAL>& get_gradient() = 0;
  
  virtual void print(ostream &output=std::cout) = 0;
  virtual void attach_optimizer(Optimizer *opt) = 0;
  

  virtual vector<REAL> get_state() = 0;
  virtual void set_state(const vector<REAL>& new_state) = 0;
  
  virtual vector<REAL> count() = 0;
  virtual vector<REAL> mean() = 0;
  virtual vector<REAL> variance(vector<REAL> means) = 0;
  virtual void Normalize(vector<REAL> means,vector<REAL> variances) = 0;
  virtual vector<REAL> get_targets() = 0;
  virtual void center_outputs(REAL bias_correction) = 0;

  virtual void init_optimizer(Optimizer* opt) = 0;

  virtual void set_means(vector<REAL> means) = 0;
  virtual void set_variances(vector<REAL> variances) = 0;

  virtual bool is_preconditioned() = 0;

  virtual void set_output_mean(REAL new_mean) = 0;
  virtual void set_output_variance(REAL new_variance) = 0;
  
 private:
    
    
  

  

};

#endif
