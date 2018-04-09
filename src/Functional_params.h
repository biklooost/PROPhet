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
//                        CLASS DESCRIPTION
// ####################################################################
// This is effectively a container class that holds information about
// what the user wants to do. It allows objects to have access to run
// parameters set by the user or by other objects. Most details of a
// network are held in here.
// ####################################################################







#ifndef __FUNCTIONAL_PARAMS
#define __FUNCTIONAL_PARAMS

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <vector>
#include <map>

#include "Common.h"
#include "Atom.h"

using namespace std;

class Functional_params
{

  friend class Setup;

public:

  Functional_params();
  ~Functional_params();

  inline int is_key(string key)
  {
    return my_keys.count(key);
  }


  inline int hidden(int i)
  {
    return my_hidden.at(i);
  }
  inline void set_hidden(vector<int> Node_counts)
  {
    my_hidden.clear();
    for (int i=0; i<Node_counts.size(); i++) {
      my_hidden.push_back(Node_counts.at(i));
    }
  }

  inline string inputs(int i)
  {
    return my_inputs.at(i);
  }
  inline vector<string> inputs()
  {
    return my_inputs;
  }
  inline void set_inputs(vector<string> new_inputs)
  {
    my_inputs.clear();
    for (int i=0; i<new_inputs.size(); i++) {
      my_inputs.push_back(new_inputs.at(i));
    }
  }



  inline int Ninputs()
  {
    return my_inputs.size();
  }
  inline int Ninput_nodes()
  {
    if (my_Ninput_nodes) {
      return my_Ninput_nodes;
    } else {
      my_Ninput_nodes = my_inputs.size();
      return my_Ninput_nodes;
    }
  }

  inline void Ninput_nodes(int new_NNodes)
  {
    my_Ninput_nodes = new_NNodes;
  }

  inline string output()
  {
    return my_output;
  }
  inline void set_output(string new_out)
  {
    my_output = new_out;
  }

  inline int Nlayers()
  {
    return my_hidden.size();
  }
  inline int NNodes(int i)
  {
    if (i<0) {
      return 0;
    }
    return my_hidden.at(i);
  }

  inline REAL threshold()
  {
    return my_threshold;
  }

  inline string transfer_function()
  {
    return my_transfer_function;
  }

  inline string input_file(string filename = "")
  {
    if (filename.empty()) {
      return my_input_file;
    } else {
      my_input_file = filename;
      return filename;
    }
  }

  inline string output_file(string filename = "")
  {
    if (filename.empty()) {
      return my_output_file;
    } else {
      my_output_file = filename;
      return filename;
    }
  }

  inline string Network_type()
  {
    return my_network_type;
  }

  void print(ostream &out);
  void read(istream &in);

  inline int Niterations()
  {
    return my_Niterations;
  }
  inline int Ncheckpoint()
  {
    return my_Ncheckpoint;
  }
  inline int Nprint()
  {
    return my_Nprint;
  }
  inline int Ncycles()
  {
    return my_Ncycles;
  }
  inline int lambda_max()
  {
    return my_lambda_max;
  }

  inline REAL line_min_epsilon()
  {
    return my_line_min_epsilon;
  }
  inline REAL sd_momentum()
  {
    return my_sd_momentum;
  }

  inline string training_algorithm()
  {
    return my_training_algorithm;
  }

  inline REAL debug()
  {
    return my_debug;
  }

  vector<REAL> input_mean;
  vector<REAL> input_variance;


  inline bool output_precondition()
  {
    return my_output_precondition;
  }
  inline bool input_precondition()
  {
      return my_input_precondition;
  }
  REAL output_mean;
  REAL output_variance;

  inline int sample_step()
  {
    return my_downsample;
  }

  inline bool output_is_intensive()
  {
    if (my_output == "dft_gap" || my_output=="gw_gap") {
      return true;
    } else {
      return false;
    }
  }

  inline bool NormCD()
  {
    return my_norm_cd;
  }
  inline REAL NormCD_val()
  {
    return my_norm_cd_val;
  }
  inline bool Tbackup()
  {
    return my_Tbackup;
  }
  inline int Nbackup()
  {
    return my_Nbackup;
  }
  inline int Nradial()
  {
    return my_Nradial;
  }
  inline int Nangular()
  {
    return my_Nangular;
  }
  inline REAL Rcut()
  {
    return my_Rcut;
  }
  inline bool SGD()
  {
    return my_SGD;
  }
  inline int SGD_cnt()
  {
    return my_SGD_cnt;
  }
  inline void Rcut(REAL cutoff)
  {
    my_Rcut = cutoff;
  }
  inline map<string,REAL> FE()
  {
    return my_FE;
  }
  inline vector <int> var_bounds()
  {
    return my_bounds;
  }
  inline REAL dropout()
  {
    return my_dropoutP;
  }
  inline int Nconv()
  {
    return my_conv;
  }
  inline int Nconv_stride()
  {
    return my_conv_stride;
  }
  inline bool PGvectors()
  {
    return my_printgvectors;
  }

  vector<vector<REAL> >  G1;
  vector<vector<REAL> >  G2;
  vector<vector<REAL> >  G3;
  vector<vector<REAL> >  G4;

  bool is_MD;

  inline REAL regularization()
  {
    return my_regularization;
  }

  inline vector<int> atomic_numbers()
  {
    return my_atomic_numbers;
  }
  inline vector<string> atomic_symbols()
  {
    return my_atomic_symbols;
  }
  inline void atomic_numbers(vector<int> new_types)
  {
    my_atomic_numbers = new_types;
    my_atomic_symbols.clear();
    for (int i=0; i<my_atomic_numbers.size(); i++) {
      Atom atom(my_atomic_numbers[i]);
      my_atomic_symbols.push_back(atom.atomic_symbol());
    }
  }

  string current_atom_type;

  // Parameters relating to the Monte Carlo optimization of the G
  // functions of the Behler method
  REAL T0;
  REAL Tf;
  REAL dT;
  int Nsteps_per_T;
  REAL step_size;



private:

  std::set<string> my_keys;
  std::set<string> my_valid_params;
  bool my_SGD;
  int my_SGD_cnt;
  bool my_output_precondition;
  bool my_input_precondition;
  REAL my_regularization;
  bool my_norm_cd;
  REAL my_norm_cd_val;

  int my_Ninput_nodes;

  string my_input_file;
  string my_output_file;
  string my_transfer_function;
  string my_network_type;

  REAL my_threshold;
  REAL my_lambda_max;

  REAL my_line_min_epsilon;
  REAL my_sd_momentum;

  vector<int> my_hidden;
  vector<string> my_inputs;
  string my_output;
  int my_Ncycles;
  int my_Niterations;
  int my_Ncheckpoint;
  int my_Nprint;
  int my_Nbackup; //How often to save the functional file
  bool my_Tbackup;
  REAL my_debug;
  vector <int> my_bounds; //vector of lbound and rbound used in variance in Charge Density

  REAL my_dropoutP; //dropout probability

  string my_training_algorithm;

  int my_downsample;
  vector<int> my_atomic_numbers;
  vector<string> my_atomic_symbols;

  // Parameters relating to the potential surface method
  // of Behler and Parrinello
  REAL my_Rcut;
  int my_Nradial;
  int my_Nangular;
  std::map<string,REAL> my_FE;
  //Printing the g-vectors to *.gvector files
  bool my_printgvectors;

  int my_conv;
  int my_conv_stride;


};

#endif

