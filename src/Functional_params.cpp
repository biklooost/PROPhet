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



#include "Functional_params.h"

#include <sstream>
#include <limits>
#include <iomanip>



// ########################################################
//                       Constructor
// ########################################################
//

Functional_params::Functional_params()
{

  // default values
  this->my_input_file = "";
  this->my_output_file = "functional.save";
  this->my_network_type = "neural_network";
  this->my_transfer_function = "tanh_spline";
  this->my_training_algorithm = "default";

  this->my_Niterations = std::numeric_limits<int>::max();
  this->my_Ncheckpoint = 100;
  this->my_Nprint = 100;
  this->my_Ncycles = 1;
  this->my_lambda_max = 0;
  this->my_downsample = 1;
  this->my_Ninput_nodes = 0;
  this->my_line_min_epsilon = 0;
  this->my_sd_momentum = 0.0;
  this->my_SGD = false;
  this->my_SGD_cnt = 0;

  this->my_Rcut = 6.0;
  this->my_Nradial = 8;
  this->my_Nangular = 8;

  this->T0 = 1;
  this->Tf = 0.1;
  this->dT = 0.98;
  this->Nsteps_per_T = 50;
  this->step_size = 0.2;

  this->is_MD = false;

  this->my_regularization = 0.0;

  this->my_output_precondition = false;
  this->my_input_precondition = true;
  this->output_mean = 0;
  this->output_variance = 1;

  this->my_Nbackup = 0;
  this->my_Tbackup = false;
  this->my_norm_cd = false;
  this->my_bounds = vector <int> (2,0);

  this->my_dropoutP = 1.0;

  this->my_conv = 250; //This needs to just be some large number
  this->my_conv_stride = 1;

  this->my_printgvectors = false;
}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

Functional_params::~Functional_params()
{

}

// ########################################################
// ########################################################




// ########################################################
//                       PRINT
// ########################################################
// Write state to checkpoint file

void Functional_params::print(ostream &out)
{

  out << this->my_network_type << endl;

  for (int i=0; i<my_inputs.size(); i++) {
    out << my_inputs.at(i) <<"  ";
  }
  out << endl;
  if (!current_atom_type.empty()) {
    out << current_atom_type<<":  ";
    for (int i=0; i<my_atomic_symbols.size(); i++) {
      out << my_atomic_symbols[i] << " ";
    }
    out << endl;
    if (!this->my_FE.empty()) {
      out.precision(9);
      out << "FE:  " << this->my_FE[current_atom_type] << endl;
      out.precision(6);
    }
    out << G1.size()+G2.size()+G3.size()+G4.size()<<endl;
    for  (int i=0; i<G1.size(); i++) {
      out << setprecision(12)<<"G1 "<<G1.at(i).at(0) << " "<<G1.at(i).at(1)
          << " "<<G1.at(i).at(2)<<endl;
    }
    for (int i=0; i<G2.size(); i++) {
      out << setprecision(12)<<"G2 "<<G2.at(i).at(0) << " "<<G2.at(i).at(1)
          <<" "<<G2.at(i).at(2)<< " "<<G2.at(i).at(3)<<endl;
    }
    for (int i=0; i<G3.size(); i++) {
      out << setprecision(12)<<"G3 "<<G3.at(i).at(0)<<" "<<G3.at(i).at(1)
          <<" "<<G3.at(i).at(2)<<" "<< G3.at(i).at(3)<<" "<<G3.at(i).at(4)
          <<endl;
    }
    for (int i=0; i<G4.size(); i++) {
      out << setprecision(12)<<"G4 "<<G4.at(i).at(0)<<" "<<G4.at(i).at(1)
          <<" "<<G4.at(i).at(2)<<" "<< G4.at(i).at(3)<<" "<<G4.at(i).at(4)
          <<endl;
    }
  }
  for (int i=0; i<input_mean.size(); i++) {
    out <<input_mean.at(i) << "  ";
  }
  out << endl;
  for (int i=0; i<input_variance.size(); i++) {
    out << input_variance.at(i) << "  ";
  }
  out << endl;

  out << my_output << endl;
  out << output_mean << endl;
  out << output_variance << endl;
  if (my_hidden.size()) {
    for (int i=0; i<my_hidden.size(); i++) {
      out << my_hidden.at(i) << "  ";
    }
  } else {
    out << 0;
  }
  out << endl;
}

// ########################################################
// ########################################################





// ########################################################
//                       READ
// ########################################################
// Read state from checkpoint file

void Functional_params::read(istream &in)
{
  string line, input;
  istringstream Line;
  int hidden;
  bool is_potential = false;
  REAL values;


  getline(in,line);
  this->my_network_type = line;

  getline(in,line);
  Line.str(line);

  vector<string> input_names;
  Line >> input;
  if (input.empty()) {
    ERROR("Could not read network inputs from file");
  }
  input_names.push_back(input);
  if (input == "structure") {
    is_potential = true;
  }
  while (Line >> input) {
    input_names.push_back(input);
    if (input == "structure") {
      is_potential = true;
    }
  }
  this->set_inputs(input_names);
  Line.clear();

  if (is_potential) {
    my_atomic_symbols.clear();
    my_atomic_numbers.clear();
    getline(in,line);
    Line.str(line);
    string symbol;
    Line >> symbol;
    symbol.erase(symbol.end()-1,symbol.end());
    this->current_atom_type = symbol;
    while (Line >> symbol) {
      Atom atom(symbol);
      my_atomic_symbols.push_back(atom.atomic_symbol());
      my_atomic_numbers.push_back(atom.atomic_number());
    }
    Line.clear();

    getline(in, line);
    if (line.find("FE") != std::string::npos) {
      Line.str(line);
      REAL FE;
      string temp;
      Line.precision(9);
      Line >> temp >> FE;
      this->my_FE[this->current_atom_type] = FE;
      Line.precision(6);
      Line.clear();
      getline(in,line);
    }
    Line.str(line);
    int Nfunctions, function_count=1;
    Line >> Nfunctions;
    G1.clear();
    G2.clear();
    G3.clear();
    G4.clear();
    Line.clear();

    REAL Rcut_max = 0;
    string function_type;
    while (function_count <= Nfunctions) {
      getline(in,line);
      Line.str(line);
      vector<REAL> temp;
      Line >> function_type;
      while (Line >> values) {
        temp.push_back(values);
      }
      Line.clear();
      if (function_type == "G1") {
        G1.push_back(temp);
      } else if (function_type == "G2") {
        G2.push_back(temp);
      } else if (function_type == "G3") {
        G3.push_back(temp);
      } else if (function_type == "G4") {
        G4.push_back(temp);
      } else {
        ERROR("Symmetry functions not read correctly from checkpoint file");
      }
      function_count++;
      if (temp[0] > Rcut_max) {
        Rcut_max = temp[0];
      }
    }
    if (G1.size()+G2.size()+G3.size()+G4.size()) {
      this->my_Ninput_nodes = my_atomic_numbers.size()*(G1.size()+G2.size())
                              +0.5*my_atomic_numbers.size()*(my_atomic_numbers.size()+1)*(G3.size()+G4.size());
    } else {
      this->my_Ninput_nodes = this->inputs().size();
    }
    if (Rcut_max > 0) {
      this->Rcut(Rcut_max);
    }
  }

  input_mean.clear();
  input_variance.clear();
  getline(in,line);
  string numbers = " \t-+0123456789eEdD.";
  if (line.find_first_not_of(numbers) != string::npos) {
    ERROR("Input means read incorrectly");
  }
  Line.str(line);
  while (Line >> values) {
    input_mean.push_back(values);
  }
  Line.clear();

  getline(in,line);
  if (line.find_first_not_of(numbers) != string::npos) {
    ERROR("Input variances read incorrectly");
  }
  Line.str(line);
  while (Line >> values) {
    input_variance.push_back(values);
  }
  Line.clear();

  if (input_mean.size() != input_variance.size() || input_mean.empty()) {
    ERROR("Input means/variances read incorrectly");
  }

  getline(in,line);
  Line.str(line);
  Line >> this->my_output;
  Line.clear();

  getline(in, line);
  Line.str(line);
  Line >> this->output_mean;
  Line.clear();
  getline(in,line);
  Line.str(line);
  Line >> this->output_variance;
  Line.clear();

  vector<int> layer_config;

  getline(in,line);
  Line.str(line);
  if (line.empty()) {
    ERROR("Could not read node configuration from save file");
  }
  Line >> hidden;

  layer_config.push_back(hidden);
  while (Line >> hidden) {
    layer_config.push_back(hidden);
  }
  this->my_hidden = layer_config;
  Line.clear();

}

// ########################################################
// ########################################################




