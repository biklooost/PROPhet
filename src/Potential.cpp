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



#include "Potential.h"

#include <set>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "NN.h"
#include "NF.h"
#include "Atom.h"
#include "Analysis.h"


// ########################################################
//                       Constructor
// ########################################################
// 

Potential::Potential(vector<System*> systems_in, Functional_params F_in) : systems(systems_in),params(F_in) {

  mpi = new Parallel;

  this->atom_types = this->get_all_atom_types();

  params.atomic_numbers(atom_types);
  
  this->Nsystems = mpi->Reduce(systems.size(), MPI_SUM);
  Nsystems = mpi->Bcast(Nsystems);
  
  Ntotal_params = 0;
  
  for (int i=0; i<atom_types.size(); i++) {
    string filename;
    Atom atom(atom_types[i]);
    if (!F_in.input_file().empty()) {
      ifstream file;
      filename = params.input_file() + "_" + atom.atomic_symbol();
      file.open(filename.c_str());
      if (!file.is_open()) { ERROR("Found atom type "+atom.atomic_symbol()+" but could not open restart file '"+filename+"'"); }
      if (mpi->io_node()) { 
	cout << "Reading parameters for "+atom.atomic_symbol()+" from file '"+filename+"'"<<endl;
      }
      
      params.read(file);

      if (mpi->io_node()) {
	cout << "Hidden layers :  ";
	for (int layer=0; layer<params.Nlayers(); layer++) {
	  cout << params.NNodes(layer)<<"  ";
	}
        cout << endl;
        map<string,REAL> FE = params.FE();
        if(!FE.empty()) {
            cout << params.current_atom_type << " Free Energy : " << FE[params.current_atom_type] << endl;
            
        }
	cout << endl << endl;
      }
    }
    
    if (!systems[0]->structure.is_initialized()) {
      for (int i_sys=0; i_sys<systems.size(); i_sys++) {
	systems[i_sys]->properties.set_inputs(systems[i_sys]->structure.init_G(&params));
      }
    }
    
    
    params.current_atom_type = atom.atomic_symbol();
    
    if (params.Network_type() == "neural_network" || params.Network_type() == "nn") {
      nets.insert(pair<int, Network*>(atom_types[i], new Neural_network(systems, params)));
    } else if (params.Network_type() == "network_functions" || params.Network_type() == "nf") {
      nets.insert(pair<int, Network*>(atom_types[i], new Network_function(systems, params)));
    }
    if (!F_in.input_file().empty()) {
      ifstream file;
      file.open(filename.c_str());
      if (!file.is_open()) { ERROR("Could not open network restart file '"+filename+"'"); }
      nets[atom_types[i]]->load(file);
    }
    Ntotal_params += nets[atom_types[i]]->Nparameters();
  }
  
  gradient.reserve(Ntotal_params);
  
  Ntotal_params = 0;
  for (map<int,Network*>::iterator it = nets.begin(); it != nets.end(); ++it) {
    gradient.insert(gradient.end(), it->second->Nparameters(), 0.0);
    indices.insert(pair<int,vector<REAL>::iterator>(it->first, gradient.begin()+Ntotal_params));
    Ntotal_params += it->second->Nparameters();
  }
  
  if (mpi->io_node()) {

    double ratio = (double)Nsystems/Ntotal_params;
    if (ratio <= 1) {
      cout << "<< WARNING >> ";
    }
    cout << "Ratio of training data to parameters =  "<<ratio<<endl<<endl<<endl;
  
  }
  
  this->output_mean = params.output_mean;

}

// ########################################################
// ########################################################




// ########################################################
//                       Constructor
// ########################################################
//

Potential::Potential(Functional_params F_in):params(F_in) {
  
}

// ########################################################
// ########################################################




// ########################################################
//                       Constructor
// ########################################################
//

Potential::Potential() {

}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

Potential::~Potential() {

  if (mpi) { delete mpi; }
  for (map<int,Network*>::iterator it=nets.begin(); it != nets.end(); ++it) {
    if (it->second) { delete it->second; }
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       SET_F_PARAMS
// ########################################################
// Sets the internal version of the Functional_params object

void Potential::set_F_params(Functional_params F_in) {
  params = F_in;
}

// ########################################################
// ########################################################



// ########################################################
//                       SAVE
// ########################################################
// Write the checkpoint file

void Potential::save() {
  
  if (!mpi->io_node()) { return; }

  Atom atom_type;

  for (map<int, Network*>::iterator it=nets.begin(); it != nets.end(); ++it) {
    
    atom_type.set_type(it->first);

    ofstream output;
    string filename = params.output_file()+"_"+atom_type.atomic_symbol();
    output.open(filename.c_str());
    it->second->print(output);
    output.close();
  
  }

}

// ########################################################
//                       BACKUP
// ########################################################
// Write the Backup Potential file
void Potential::backup(int niter) {
  
  if (!mpi->io_node()) { return; }
  ostringstream convert; 
  convert << niter - 1;
  Atom atom_type;
  
  for (map<int, Network*>::iterator it=nets.begin(); it != nets.end(); ++it) {
    
    atom_type.set_type(it->first);

    ofstream output;
    string filename = params.output_file()+"_"+atom_type.atomic_symbol()+"_"+convert.str();
    output.open(filename.c_str());
    it->second->print(output);
    output.close();
  
  }

}




// ########################################################
//                       LOAD
// ########################################################
// Read the checkpoint file

void Potential::load(string in) {

  Atom atom_type;
  
  for (map<int, Network*>::iterator it = nets.begin(); it != nets.end(); ++it) {
    ifstream input_file;
    atom_type.set_type(it->first);
    string filename = in + "_" + atom_type.atomic_symbol();
    input_file.open(filename.c_str());
    if (!input_file.is_open()) {ERROR("Could not open input file"+filename+"\"");}
    it->second->load(input_file);
    input_file.close();
  }
  
}

// ########################################################
// ########################################################





// ########################################################
//                       PRECONDITION
// ########################################################
// Perform input and (optionally) output preconditioning.
// This helps prevent node saturation and eases training.

void Potential::Precondition() {
			      
  if (mpi->io_node()) {
    cout << "Preconditioning ... "<<endl<<endl;
  }
  
  vector<REAL>  means;
  vector<REAL> variances;

  for (int atom=0; atom<atom_types.size(); atom++) {
    int type = atom_types[atom];
    
    if (nets.at(type)->is_preconditioned()) {
      continue;
    }

    int count = 0;
    means.assign(systems[0]->structure.NG(),0.0);
    
    
    for (int system=0; system<systems.size(); system++) {
      vector<REAL> temp_means = systems[system]->structure.mean(type);
      count += systems[system]->structure.count(type);
      for (int i=0; i<temp_means.size(); i++) {
	means[i] += temp_means[i];
      }
    }
    
    
    means = mpi->Reduce(means, MPI_SUM);
    count = mpi->Reduce(count, MPI_SUM);
    if (mpi->rank() == 0 && count != 0) {
      for (int i=0; i<means.size(); i++) {
	means[i] /= count;
      }
    }
    
    means = mpi->Bcast(means,means.size());
    nets.at(type)->set_means(means);
    
    variances.assign(systems[0]->structure.NG(),0.0);
    for (int system=0; system<systems.size(); system++) {
      vector<REAL> temp_variances = systems[system]->structure.variance(type, means);
      for (int i=0; i<temp_variances.size(); i++) {
	variances[i] += temp_variances[i];
      }
    }
    variances = mpi->Reduce(variances, MPI_SUM);
    if (mpi->rank() == 0 && count != 0) {
      for (int i=0; i<variances.size(); i++) {
	variances[i] = sqrt(variances[i]/count);
	if (variances[i] < 1e-5) { variances[i] = 1e-5; }
      }
    }
    variances = mpi->Bcast(variances, variances.size());
    nets.at(type)->set_variances(variances);
  }
  
  
  this->syncronize();
  
  if ( params.output_precondition() ) {
    
    // Precondition the network outputs
    vector<REAL> targets = nets.begin()->second->get_targets();
    targets = mpi->Gatherv(targets);
    
    this->output_mean = 0;
    this->output_variance = 0;
    
    if (mpi->io_node()) {
      int N = targets.size();
      for (int i=0; i<N; i++) {
	output_mean += targets.at(i);
      }
      output_mean /= N;
      for (int i=0; i<N; i++) {
	output_variance += (targets.at(i) - output_mean)*(targets.at(i) - output_mean);
      }
      output_variance /= N;
    }
    
    output_mean = mpi->Bcast(output_mean);
    output_variance = mpi->Bcast(output_variance);
   
  } else {
    output_mean = 0;
    output_variance = 1;
  }
 
  for (map<int,Network*>::iterator it=nets.begin(); it != nets.end(); ++it) {
    it->second->set_output_mean(output_mean);
    it->second->set_output_variance(output_variance);
  }
  
}

// ########################################################
// ########################################################




// ########################################################
//                       EVALUATE
// ########################################################
// Fire the potential to make predictions

vector<REAL> Potential::evaluate() {
  
  if (params.input_file().empty()) {
    ERROR("Cannot do predictions without checkpoint file (checkpoint_in flag)");
  }

  vector<REAL> my_outputs(systems.size(), 0.0);
  
  this->syncronize();

  for (int i_sys=0; i_sys<systems.size(); i_sys++) {
    for (int atom=0; atom<systems[i_sys]->structure.Natom; atom++) {
      my_outputs[i_sys] += nets[systems[i_sys]->structure.types[atom].atomic_number()]->evaluate_MD(i_sys);
    }
  }

  vector<REAL> output = mpi->Gatherv(my_outputs);
  
  if (mpi->io_node()) {
    cout << endl;
    cout << "System         Prediction"<<endl;
    cout << "-----------------------------"<<endl;
    int Nroot = my_outputs.size();
    int count = 1;

    create_system_map();
    for (int i=0; i<system_map.size(); i++) {
      int index = system_map[i];
	cout << count++ << "            "<<output.at(index)+output_mean<<endl;
    }
    cout << endl;
  }
  
  return output;
  
}

// ########################################################
// ########################################################




// ########################################################
//                       SYNCRONIZE
// ########################################################
// Make sure all cores have the same network parameters.

void Potential::syncronize() {

  for (map<int,Network*>::iterator it = nets.begin(); it != nets.end(); ++it) {
    vector<REAL> current_params = it->second->get_state();
    it->second->set_state(mpi->Bcast(current_params,current_params.size()));
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       INIT_OPTIMIZER
// ########################################################
// Initialize the optimizer we're going to use to train the
// potential.

void Potential::init_optimizer() {

  for (map<int,Network*>::iterator it=nets.begin(); it != nets.end(); ++it) {
    it->second->init_optimizer(opt);
  }
  opt->internal_init();

}

// ########################################################
// ########################################################




// ########################################################
//                       TRAIN
// ########################################################
// Train the potential

REAL Potential::train() {
  
  this->Precondition();
  opt = new Optimizer(mpi,params,Nsystems);
  init_optimizer();  
  
  REAL Error;
  REAL output;
  
  this->syncronize();

  while (!opt->is_converged()) {
    
    REAL SSE = 0;
    gradient.assign(Ntotal_params, 0.0);
    
    for (int i_sys=0; i_sys<systems.size(); i_sys++) {
      
      output = 0;
      for (int atom=0; atom<systems[i_sys]->structure.Natom; atom++) {
	output += nets[systems[i_sys]->structure.types[atom].atomic_number()]->train_MD(i_sys);
      }
      
      Error = (output + output_mean) - systems[i_sys]->properties.target();
      SSE += Error*Error;

      vector<REAL>::iterator vect_it;
      
      for (map<int,Network*>::iterator it=nets.begin(); it != nets.end(); ++it) {
	const vector<REAL> temp_grad = it->second->get_gradient();
	vect_it = indices[it->first];
	for (int i=0; i<temp_grad.size(); i++) {
	  (*vect_it) += 2*Error*temp_grad[i];
	  ++vect_it;
	}
	
      }
      
    }
    
    opt->update_network(SSE, gradient);
    
    if (opt->is_checkpoint()) {
      this->save();
    }
    
    if (opt->is_backup()){
        this->backup(opt->iteration());
    }
    
  }
  
  this->save();

  REAL final_error = sqrt(opt->get_error());

  delete opt;
  
  return final_error;

}

// ########################################################
// ########################################################




// ######################################################## 
//                      OPTIMIZE_GS
// ########################################################
// This is a routine to use a Monte Carlo algorithm to 
// find the set of G functions for a potential.  The 
// algorithm is described in :
// J. Chem. Phys. 144, 224103 (2016)

void Potential::optimize_Gs() {
  
  int counter = 0;
  vector<vector<REAL> > G1, G2, G3, G4;
  
  REAL T = params.T0;
  int Naccepted;
  
  REAL last_E = 9e9;

  ofstream log_file, parameter_file;

  log_file.open("logfile.out",std::fstream::out);

  if (mpi->rank() == 0) {
    
    log_file << "  T    accepted    Error"<<endl;
    log_file << " -------------------------"<<endl;
    
  }
  
  G1 = systems[0]->structure.G1p;
  G2 = systems[0]->structure.G2p;
  G3 = systems[0]->structure.G3p;
  G4 = systems[0]->structure.G4p;
  
  
  REAL max_eta = 65;
  
  while ( T >= params.Tf ) {
    
    if (mpi->rank() == 0) {
      Naccepted = 0;
      log_file << T << " : ";
    }    
    
    for (int Step = 1; Step <= params.Nsteps_per_T; Step++) {
      if (mpi->rank() == 0) {
	cout << "Step : " << Step << endl;
      }
      if (mpi->rank() == 0 && counter) {
	
	for (int i=0; i<params.G2.size(); i++) {
	  REAL parameter = RAND::Normal(0,params.step_size) + G2[i][2];
	  while (parameter <= 0 || parameter > max_eta) {
	    parameter = RAND::Normal(0,params.step_size) + G2[i][2];
	  }
	  params.G2[i][2] = parameter;
	  
	  parameter = RAND::Normal(0,params.step_size) + G2[i][3];
	  while (parameter < 0 || parameter > params.Rcut()) {
	    parameter = RAND::Normal(0,params.step_size) + G2[i][3];
	  }
	  params.G2[i][3] = parameter;
	  
	}
	
	for (int i=0; i<params.G3.size(); i++) {
	  REAL parameter = RAND::Normal(0,params.step_size) + G3[i][2];
	  while (parameter <= 0 || parameter > max_eta) {
	    parameter = RAND::Normal(0,params.step_size) + G3[i][2];
	  }
	  params.G3[i][2] = parameter;
	  
	  parameter = RAND::Normal(0,params.step_size) + G3[i][3];
	  while (parameter < 0 || parameter > params.Rcut()) {
	    parameter = RAND::Normal(0,params.step_size) + G3[i][3];
	  }
	  params.G3[i][3] = parameter;
	}

	for (int i=0; i<params.G4.size(); i++) {
	  REAL parameter = RAND::Normal(0,params.step_size) + G4[i][2];
	  while (parameter <= 0 || parameter > max_eta) {
	    parameter = RAND::Normal(0,params.step_size) + G4[i][2];
	  }
	  params.G4[i][2] = parameter;
	  
	  parameter = RAND::Normal(0,params.step_size) + G4[i][3];
	  while (parameter <= 0 || parameter > max_eta) {
	    parameter = RAND::Normal(0,params.step_size) + G4[i][3];
	  }
	  params.G4[i][3] = parameter;
	}

      }
      
      if (counter) {
	for (int i=0; i<params.G2.size(); i++) {
	  params.G2[i] = mpi->Bcast(params.G2[i],params.G2[i].size());
	}
	for (int i=0; i<params.G3.size(); i++) {
	  params.G3[i] = mpi->Bcast(params.G3[i],params.G3[i].size());
	}
	for (int i=0; i<params.G4.size(); i++) {
	  params.G4[i] = mpi->Bcast(params.G4[i],params.G4[i].size());
	}
	
	for (int i=0; i<systems.size(); i++) {
	  systems[i]->properties.set_inputs(systems[i]->structure.init_G(&params));
	}
      }
      
      REAL Error = this->train();
      
      if (mpi->rank() == 0) {
	cout << "Error =  "<<Error << endl;
	if (RAND::Uniform(0,1) < exp(-(Error - last_E)/T)) {
	  G2 = params.G2;
	  G3 = params.G3;
	  G4 = params.G4;
	  last_E = Error;
	  Naccepted++;

	  parameter_file.open("G_parameters.out",std::fstream::out);

	  for (int i=0; i<G1.size(); i++) {
	    parameter_file << setprecision(10)<<"G1  ";
	    for (int j=0; j<G1[i].size(); j++) {
	      parameter_file << G1[i][j] <<" ";
	    }
	    parameter_file << endl;
	  }
	  
	  for (int i=0; i<G2.size(); i++) {
	    parameter_file << "G2  ";
	    for (int j=0; j<G2[i].size(); j++) {
	      parameter_file << G2[i][j] <<" ";
	    }
	    parameter_file << endl;
	  }  
	  
	  for (int i=0; i<G3.size(); i++) {
	    parameter_file << "G3  ";
	    for (int j=0; j<G3[i].size(); j++) {
		  parameter_file << G3[i][j] <<" ";
	    }
	    parameter_file << endl;
	  }

	  for (int i=0; i<G4.size(); i++) {
	    parameter_file << "G4  ";
	    for (int j=0; j<G4[i].size(); j++) {
	      parameter_file << G4[i][j] <<" ";
	    }
	    parameter_file << endl;
	  }
	    
	  parameter_file.close();

	}

      }
     
      counter++;
 
    }
    
    if (mpi->rank() == 0) {
      log_file << (double)(Naccepted)/(double)(params.Nsteps_per_T)<<"  "<<last_E <<endl;
    }
    
    T *= params.dT;
    
  }
  
  parameter_file.close();
  
}

// ######################################################## 
// ######################################################## 






// ########################################################
//                       GET_ALL_ATOM_TYPES
// ########################################################
// Figures out what atomic species are present in the 
// union of all systems and broadcasts the result.

vector<int> Potential::get_all_atom_types() {

  set<int> atomic_numbers;
  
  for (int sys=0; sys<systems.size(); sys++) {
    for (int type=0; type<systems[sys]->structure.types.size(); type++) {
      atomic_numbers.insert(systems[sys]->structure.types.at(type).atomic_number());
    }
  }
  
  vector<int> numbers;
  numbers.assign(atomic_numbers.begin(), atomic_numbers.end());
  
  numbers = mpi->Gatherv(numbers);
  
  // This filters out repeated atom types
  atomic_numbers.clear();
  for (int i=0; i<numbers.size(); i++) {
    atomic_numbers.insert(numbers[i]);
  }
  
  // Broadcast the unique atomic numbers to all cores
  numbers.assign(atomic_numbers.begin(), atomic_numbers.end());
  std::sort(numbers.begin(),numbers.end());
  int N = mpi->Bcast(numbers.size());
  numbers = mpi->Bcast(numbers, N);
  
  return numbers;

}

// ########################################################
// ########################################################



// ########################################################
//                       VALIDATE
// ########################################################
// This is an analysis routine to test the fidelity of the
// potential. 

void Potential::validate() {

  if (params.input_file().empty()) {
    ERROR("Cannot do validation without checkpoint file (checkpoint_in flag)");
  }

  vector<REAL> my_outputs(systems.size(), 0.0);
  vector<REAL> unraveled(systems.size(),0.0);
  map<string,REAL> FE = params.FE();
  bool unrav = !FE.empty();
  this->syncronize();
  
  for (int i_sys=0; i_sys<systems.size(); i_sys++) {
    for (int atom=0; atom<systems[i_sys]->structure.Natom; atom++) {
      my_outputs[i_sys] += nets[systems[i_sys]->structure.types[atom].atomic_number()]->evaluate_MD(i_sys);
    }
    if (unrav) {
        unraveled[i_sys] = systems[i_sys]->structure.unravel_Energy(&params,my_outputs[i_sys]);
    }
  }
  
  vector<REAL> output = mpi->Gatherv(my_outputs);
  unraveled = mpi->Gatherv(unraveled);
  
  vector<REAL> targets = nets.begin()->second->get_targets();
  targets = mpi->Gatherv(targets);
  vector<REAL> unTargets(unraveled.size(),0.0);
  if (unrav) {
        for (int i_sys=0; i_sys<systems.size(); i_sys++) {
            unTargets[i_sys] = systems[i_sys]->structure.unravel_Energy(&params,targets.at(i_sys));
        }
  }

  double SSE = 0.0;
  int count = 0;

  if (mpi->io_node()) {
    cout << endl;
    if (unrav){
        cout << "System         Prediction(per Atom)       Target(per Atom)"<<endl;
        cout << "----------------------------------------------------------"<<endl;  
    }else {
        cout << "System         Prediction       Target"<<endl;
        cout << "--------------------------------------"<<endl;
    }
    REAL Error, min_E = 9e9, max_E = 0.0;
    int Nroot = my_outputs.size(), count=1;  
    
    create_system_map();

    cout << setprecision(10);
    
    for (int i=0; i<system_map.size(); i++) {
      int index = system_map[i];
      cout <<count++ << "              "<<output.at(index)+output_mean << "              "<<targets.at(index) <<  endl;
      Error = output.at(index)+output_mean - targets.at(index);
      SSE += pow(Error,2);
      if (abs(Error)<abs(min_E)) {
	min_E = Error;
      } 
      if (abs(Error)>abs(max_E)) {
	max_E = Error;
      }
    }
    if (unrav){
        cout << endl;
        cout << "Predicted Total Energies" << endl;
        cout << "System         Prediction(total)       Target(total)"<<endl;
        cout << "----------------------------------------------------"<<endl; 
        for (int i=0; i<system_map.size(); i++) {
            int index = system_map[i];
            cout <<count++ << "              "<<unraveled.at(index) << "              "<<unTargets.at(index) <<  endl;
        }
    }
    cout << setprecision(6);
    cout << endl;
    cout << "Minimum Error =  "<<min_E<<endl;
    cout << "Maximum Error =  "<<max_E<<endl;
    cout << "RMS Error =  "<< sqrt(SSE/(double)(Nsystems)) << endl;
    
    Analysis A;
    REAL D = A.KS(output, targets);
    cout << endl<<endl;
    cout << "-- Approximate Kolmogorov-Smirnov test --"<<endl;
    cout << "  D       =  "<<D<<endl;
    cout << "  p value = "<<A.p <<endl<<endl;
    
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       ADD_SYSTEM
// ########################################################
// Push a new system to the back of the stack. 

void Potential::add_system(System* new_system) {
  this->systems.push_back(new_system);
}

// ########################################################
// ########################################################



// ########################################################
//                       EVALUATE_MD
// ########################################################
// Fire the potential in the event that we also need forces.
// Note that forces are not computed here, but dE_dG is.

double Potential::evaluate_MD(int index, int type, vector<REAL> &dE_dG) {

  if (!systems[0]->structure.is_initialized()) {
    systems[0]->properties.set_inputs(systems[0]->structure.init_G(&params));
  }
  systems[0]->structure.Calc_G(index);
  return (double)(nets.at(type)->evaluate_MD(dE_dG));
}

// ########################################################
// ########################################################




// ########################################################
//                       CREATE_SYSTEM_MAP
// ########################################################
// Basically sorts given examples so they are printed in 
// the order that they were specified during a "run" or
// "validate" run.

void Potential::create_system_map() {

  if (!this->system_map.empty()) { return; }

  vector<vector<int> > by_node(mpi->Nprocs(),vector<int>(0));
  this->system_map.assign(Nsystems,0);

  int proc = 0;
  for (int i=0; i<Nsystems; i++) {
    by_node[proc].push_back(i);
    proc++;
    proc = (proc < mpi->Nprocs() ? proc : 0);
  }

  int count = 0;
  for (int proc=0; proc<by_node.size(); proc++) {
    for (int element=0; element<by_node[proc].size(); element++) {
      system_map[by_node[proc][element]] = count++;
    }
  }

}

// ########################################################
// ########################################################

// ########################################################
//                       INSERT_ATOM_TYPE
// ########################################################
// Adds a new atom type to the potential. Necessary for the
// LAMMPS interface.

void Potential::insert_atom_type(int atom_number, char* filename, System* system) {
  
  ifstream file;
  
  if (nets.count(atom_number)) {
    return;
  }
  
  file.open(filename);
  if (!file.is_open()) { ERROR("Could not open parameter file \""+string(filename)+"\""); }
  
  params.read(file);
  this->output_mean = params.output_mean;

  nets.insert(pair<int,Network*>(atom_number,new Neural_network(vector<System*>(1,system), params)));
  
  nets[atom_number]->load(file);
  file.close();

  atom_names.push_back(params.current_atom_type);
  atom_types.push_back(atom_number);
  
  if (nets.size() == params.atomic_numbers().size()) {
    vector<string> all_types = params.atomic_symbols();
    vector<int> temp(params.atomic_numbers().size(),0);
    for (int i=0; i<atom_names.size(); i++) {
      for (int j=0; j<temp.size(); j++) {
	if (atom_names[i] == all_types[j]) {
	  temp[j] = atom_types[i];
	}
      }
    }
    atom_types = temp;
    params.atomic_numbers(atom_types);
  }
	  
}

// ########################################################
// ########################################################



