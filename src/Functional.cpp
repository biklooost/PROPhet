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



#include "Functional.h"
#include "Optimizer.h"
#include "NN.h"
#include "NF.h"
#include "Parallel.h"
#include "Analysis.h"



// ########################################################
//                       Constructor
// ########################################################
//

Functional::Functional(const vector<System*> &systems, Functional_params F_in) : params(F_in) { 
  
  if (params.Network_type() == "neural_network" || params.Network_type() == "nn") {
    
    this->net = new Neural_network(systems, F_in);
    
  } else if (params.Network_type() == "network_function" || params.Network_type() == "nf") {

    this->net = new Network_function(systems, F_in);	     
    
  } else {
    
    ERROR("Unsupported network type '"+params.Network_type()+"'");
    
  }
  
  mpi = new Parallel;
  
  if (!F_in.input_file().empty()) {
    this->load(F_in.input_file());
  }
  
  this->Nsystems = mpi->Reduce(systems.size(), MPI_SUM);
  int Ntotal_params = this->net->Nparameters();
  REAL ratio = (REAL)this->Nsystems/(REAL)Ntotal_params;
  if (mpi->io_node()) {
    if (ratio <= 1) {
        cout << "<< WARNING >> ";
    }
    cout << "Ratio of training data to parameters =  "<< (REAL)this->Nsystems/Ntotal_params << endl << endl;
  }

  
}

// ########################################################
// ########################################################





// ########################################################
//                       Destructor
// ########################################################
//

Functional::~Functional() { 

  delete this->net;
  delete this->mpi;
  
}

// ########################################################
// ########################################################





// ########################################################
//                       NORMALIZE_DATA
// ########################################################
// Routine to carry output input/output preconditioning

void Functional::Normalize_data(vector<REAL> in_mean, vector<REAL> in_variance) {

  if (this->params.Network_type() == "nf" || params.Network_type() == "network_functions") {
    return;
  }

  if (in_mean.size() == 0) {

    if (mpi->io_node()) {
      cout << "Preconditioning ... "<<endl;
    }

    vector<REAL> counts;
    counts = net->count();

    vector<REAL> means;
    means = net->mean();

    counts = mpi->Reduce(counts, MPI_SUM);
    means = mpi->Reduce(means, MPI_SUM);
    
    if (mpi->rank() == 0) {
      for (int i=0; i<counts.size(); i++) {
	means.at(i) /= counts.at(i);
      }
    }

    means = mpi->Bcast(means,means.size());

    vector<REAL> variances;
    variances = net->variance(means);
    variances = mpi->Reduce(variances, MPI_SUM);

    if (mpi->rank() == 0) {
      for (int i=0; i<variances.size(); i++) {
	variances.at(i) /= counts.at(i);
	//	cout << means.at(i) << "  " << variances.at(i)<<endl;
      }
    }
    
    variances = mpi->Bcast(variances, variances.size());

    net->Normalize(means,variances);
    
  } else {
    net->Normalize(in_mean, in_variance);

  }

  if (mpi->io_node()) {
    cout << endl;
  }

  this->syncronize();



  if (params.output_precondition() ) {
    // Precondition the network outputs                                                                              
    vector<REAL> targets = net->get_targets();
    targets = mpi->Gatherv(targets);
    
    output_mean = 0;
    output_variance = 0;
    if (mpi->io_node()) {
      this->output_mean = 0;
      this->output_variance = 0;
      
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
  net->set_output_mean(output_mean);
  net->set_output_variance(output_variance);
  
}

// ########################################################
// ########################################################




// ########################################################
//                       TRAIN
// ########################################################
// Trains the given functional

void Functional::train() { 

  
  Optimizer *opt = new Optimizer(mpi);
  
  opt->set_params(params);
  opt->set_Nsystems(Nsystems);
  opt->set_Niterations(params.Niterations());
  opt->set_Ncheckpoint(params.Ncheckpoint());
  opt->set_training_algorithm(params.training_algorithm());
  opt->set_threshold(params.threshold());
  opt->set_debug(params.debug());
  opt->set_Nbackup(params.Nbackup());
  
  
  net->attach_optimizer(opt);

  this->syncronize();

  this->Normalize_data(params.input_mean, params.input_variance);
  
  while (net->train()) {
      
    this->save();
    if (opt->is_backup()){
        this->backup(opt->iteration());
    }
    if (opt->iteration() %100 == 0){
        this->bernoulli_sample(params.dropout(),true);
    } else {
        this->bernoulli_sample(params.dropout(),false);
    }
  }
  
  delete opt;

  this->save();
  
}

// ########################################################
// ########################################################




// ########################################################
//                       EVALUATE
// ########################################################
// Fires the functional to get a prediction

void Functional::evaluate() { 

  if (params.input_file().empty()) {
    ERROR("Cannot do prediction without checkpoint file (checkpoint_in flag)");
  }

  this->Normalize_data(params.input_mean, params.input_variance);
  
  this->syncronize();

  vector<REAL> my_output = net->evaluate();
  int Nroot=my_output.size();
  vector<REAL> output = mpi->Gatherv(my_output);

  if (mpi->io_node()) {
    cout << endl;
    cout << "System         Prediction"<<endl;
    cout << "-----------------------------"<<endl;
    int count = 1;
    for (int i=0; i<Nroot; i++) {
      for (int index=i; index<output.size(); index+=(double)(output.size())/(double)(mpi->Nprocs())) {
	cout << count++ << "             " << output.at(i)<<endl;
      }
    }
    cout << endl;
  }  
  
}

// ########################################################
// ########################################################



// ########################################################
//                       VALIDATE
// ########################################################
// Performs statistics on a fit to help determine quality

void Functional::validate() {
 
  if (params.input_file().empty()) {
    ERROR("Cannot validate without checkpoint file (checkpoint_in flag)");
  }
 
  this->Normalize_data(params.input_mean, params.input_variance);
  this->syncronize();
  vector<REAL> output = net->evaluate();
  int Nroot = output.size();
  output = mpi->Gatherv(output);

  vector<REAL> my_targets = net->get_targets();
  vector<REAL> targets = mpi->Gatherv(my_targets);

  double SSE = 0.0;
  int count = 1;
  
  if (mpi->io_node()) {
    cout << endl;
    cout << "System         Prediction       Target"<<endl;
    cout << "-----------------------------------------"<<endl;
    
    REAL min_E = 9e9, max_E = -9e9, Error;
    create_system_map();

    for (int i=0; i<system_map.size(); i++) {
      int index = system_map[i];
      
      cout << setprecision(10);
      cout << count++ << "          " << output.at(index) << "  "<<targets.at(index)<<endl;
	Error = output.at(index) - targets.at(index);
	if (Error < min_E) {
	  min_E = Error;
	}
	if (Error > max_E) {
	  max_E = Error;
	}
	SSE += pow(Error,2);
      }
      
    
    cout << setprecision(6);
    cout << endl;
    cout << "Minimum error =  "<<min_E<<endl;
    cout << "Maximum error =  "<<max_E<<endl;
    cout << "RMS Error =  "<< sqrt(SSE/(double)(Nsystems))<<endl;
    
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
//                       SAVE
// ########################################################
// Write functional to checkpoint file

void Functional::save() { 
  
  if (mpi->io_node()) {
    ofstream output;
    output.open(params.output_file().c_str());

    net->print(output);
    
    output.close();
  }

}

// ########################################################
// ########################################################

// ########################################################
//                       BACKUP
// ########################################################
// Backup the functional every nsave iterations

void Functional::backup(int niter) { 
  
  if (mpi->io_node()) {
    ofstream output;
    ostringstream convert; 
    convert << niter - 1;
    string filename = params.output_file()+"_"+convert.str();
    //output.open(params.output_file().c_str());
    output.open(filename.c_str());
    net->print(output);
    
    output.close();
  }

}

// ########################################################
//                       LOAD
// ########################################################
// Read functional from checkpoint file

void Functional::load(string in) { 
  
  if (mpi->io_node()) {
      ifstream input_file;
      input_file.open(in.c_str());
      if (!input_file.is_open()) {ERROR("Could not open input file"+in+"\"");}
      net->load(input_file);
    }
    syncronize();
  
}

// ########################################################
// ########################################################





// ########################################################
//                       SYNCRONIZE
// ########################################################
// Sync network parameters across all cores during parallel
// run

void Functional::syncronize() {
  
  vector<REAL> params = net->get_state();
  net->set_state(mpi->Bcast(params,params.size()));
  
}

// ########################################################
// ########################################################




// ########################################################
//                       CREATE_SYSTEM_MAP
// ########################################################
// Basically sorts given examples so they are printed in 
// the order that they were specified during a "run" or
// "validate" run.

void Functional::create_system_map() {

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

void Functional::bernoulli_sample(REAL p,bool update) {
    vector <double> dropout;
    for (int layer=0; layer<params.Nlayers(); layer++) {
        for (int node=0; node<params.NNodes(layer); node++) {
                dropout.push_back(1.0);
        }
    }
    int cnt = 0;
    if (this->mpi->io_node() && update) {
        for (int layer=0; layer<params.Nlayers(); layer++) {
            for (int node=0; node<params.NNodes(layer); node++) {
                REAL r = ((REAL) rand() / (RAND_MAX));
                if (r >= p) {
                    dropout[cnt] = 0.0;
                } else {
                    dropout[cnt] = 1.0;;
                } 
                cnt++;
            }
        }
    }
    dropout = this->mpi->Bcast(dropout,dropout.size());
    this->net->set_dropout(dropout);
    //this->net->    
}

// ########################################################
// ########################################################




