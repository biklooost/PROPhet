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
// This is the main driver for using machine-learning for fitting
// properties to other properties.  It handles all inputs and outputs
// except structure, which is handled by the Potential class.
// ####################################################################







#ifndef __FUNCTIONAL
#define __FUNCTIONAL


#include <vector>
#include <iostream>
#include <map>

#include "Error.h"
#include "Parallel.h"
#include "Network.h"
#include "Functional_params.h"
#include "System.h"

using namespace std;

class Functional {

 public:
  
  Functional(const vector<System*> &systems, Functional_params F);
  
  ~Functional();
  
  void train();
  void evaluate();
  void validate();
  
 private:
  
  Parallel *mpi;
  Network *net;
  
  int Nsystems;
  
  Functional_params params;
  REAL output_mean;
  REAL output_variance;
  
  vector<REAL> dE_dParameters;
  
  void save();
  void backup(int niter);
  void load(string file);
  
  void syncronize();

  void Normalize_data(vector<REAL> in_mean, vector<REAL> in_variance);

  void create_system_map();
  vector<int> system_map;
  void bernoulli_sample(REAL p,bool update);

};




#endif
