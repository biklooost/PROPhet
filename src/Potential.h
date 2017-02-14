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
//                          CLASS DESCRIPTION
// ####################################################################
// This class handles the creation and use of any mapping involving
// structure (e.g. analytical potentials: structure -> energy ). It is
// akin to the Functional class which handles all non-structure
// mappings. Mappings of structure require a special class because of
// the specifics of the algorithm to deal with structure as an input.
// ####################################################################





#ifndef __POTENTIAL__
#define __POTENTIAL__

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>

#include "Parallel.h"
#include "Optimizer.h"
#include "Network.h"
#include "Functional_params.h"
#include "System.h"


using namespace std;

class Potential {
  
 public:
  Potential(vector<System*> systems_in, Functional_params F_in);
  Potential(Functional_params F_in);
  Potential();
  ~Potential();
  
  REAL train();
  vector<REAL> evaluate();
  void validate();
  void forces();
  
  void insert_atom_type(int atom_number, char* filename, System* system);
  
  //double evaluate_MD(int index, int type, vector<REAL> &dE_dG);
  double evaluate_MD(int index, int type, vector<REAL> &dE_dG);
  void add_system(System* new_system);
  inline Functional_params ret_params() { return params; };
  
  void set_F_params(Functional_params F);
  REAL energy_shift() { return output_mean; }

  void optimize_Gs();

 private:
  
  Parallel* mpi;
  map<int, Network*> nets;
  Optimizer* opt;
  
  int Nsystems;

  REAL output_mean;
  REAL output_variance;

  void init_optimizer();
  vector<string> atom_names;
  vector<int> atom_types;

  Functional_params params;
  vector<REAL> gradient;

  map<int,vector<REAL>::iterator> indices;
  void save();
  void backup(int niter);
  void load(string filename);
  void Precondition();
  
  void syncronize();

  vector<int> get_all_atom_types();
  
  vector<System*> systems;

  int Ntotal_params;
  

  void create_system_map();
  vector<int> system_map;






  
};


#endif

