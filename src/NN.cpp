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



#include "NN.h"
#include <iomanip>
#include <sstream>





// ########################################################
//                       Constructor
// ########################################################
//

Neural_network::Neural_network(const vector<System*> &in_systems, Functional_params network_details) : systems(in_systems), params(network_details) { 
  
  int current_node_index = 0;

  vector<REAL> input_layer(network_details.Ninput_nodes());
  values.push_back(input_layer);
  NNodes.push_back(network_details.Ninput_nodes());
  vector <REAL> temp(network_details.Ninput_nodes(),1.0);
  //dropout.push_back(temp);
  
  for (int i=0; i<network_details.Nlayers(); i++) {
    vector<REAL> layer_values(network_details.NNodes(i),0.0);
    vector<Neural_network_node *> layer_nodes;
    for (int j=0; j<network_details.NNodes(i); j++) {
      layer_nodes.push_back(new Neural_network_node(values.at(i).size(), network_details.transfer_function()));
      layer_nodes.back()->set_index(current_node_index++);
    }
    values.push_back(layer_values);
    nodes.push_back(layer_nodes);
    NNodes.push_back(network_details.NNodes(i));
    //vector <REAL> layer_dropout(network_details.NNodes(i),1.0);
    //dropout.push_back(layer_dropout);
  }
  vector<Neural_network_node*> output_layer;
  output_layer.push_back(new Neural_network_node(values.back().size(), "linear"));
  output_layer.back()->set_index(current_node_index++);
  nodes.push_back(output_layer);
  NNodes.push_back(1);


  Nparams = 0;
  for (int layer=0; layer<network_details.Nlayers(); layer++) {
    for (int node=0; node<network_details.NNodes(layer); node++) {
      deriv_offsets.push_back(Nparams);
      Nparams = Nparams + values.at(layer).size() + 1;
    }
  }
  deriv_offsets.push_back(Nparams);
  if (network_details.Nlayers()) {
    Nparams = Nparams + network_details.NNodes(network_details.Nlayers()-1)+1;
  } else {
    Nparams = Nparams + NNodes[0]+1;
  }
  
  for (int i=0; i<Nparams; i++) {
    my_dOutput_dParameters.push_back(0.0);
  }

  this->SSE = 0;

  this->NNodes_max = 0;
  for (int i=0; i<NNodes.size(); i++) {
    if (NNodes.at(i) > NNodes_max) {
      NNodes_max = NNodes.at(i);
    }
  }

  for (int i=0; i<in_systems.size(); i++) {
    phi.push_back(0.0);
  }
  
  last_i_sys = -1;
  
  /*if (!network_details.SGD()){
      for(int i = 0; i < in_systems.size(); i++){
          this->training_set.push_back(i);
      }
  }*/

  
}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

Neural_network::~Neural_network() { 
  
  for (int i=0; i<nodes.size(); i++) {
    for (int j=0; j<nodes.at(i).size(); j++) {
      delete nodes[i][j];
    }
  }

}

// ########################################################
// ########################################################





// ########################################################
//                       LOAD
// ########################################################
// Reads network details from checkpoint file

void Neural_network::load(ifstream &input_file) { 
  
  string line,transfer_function,junk;
  REAL param;
  
  for (int layer=0; layer<nodes.size(); layer++) {
    if (nodes.at(layer).size()==0) { continue; }
    getline(input_file,line);
    while (!input_file.eof() && line.find("[[ layer ") == std::string::npos) {
      getline(input_file, line);
    }

    for (int node=0; node<nodes.at(layer).size(); node++) {
      vector<REAL> parameters;
      istringstream Line;
      getline(input_file,line);
      while (!input_file.eof() && line.find("[ node ") == std::string::npos) {
	getline(input_file, line);
      }

      Line.str(line);
      Line >> junk >> junk >> junk >> junk >> transfer_function;
      nodes.at(layer).at(node)->set_transfer_function(transfer_function);
      Line.clear();
      getline(input_file,line);
      Line.str(line);
      while (Line >> param) {
      	parameters.push_back(param);
      }
      Line.clear();
      getline(input_file,line);
      Line.str(line);
      Line >> param;
      Line.clear();
      parameters.push_back(param);
      nodes[layer][node]->set_parameters(parameters);
    }
  }
  
}

// ########################################################
// ########################################################


// ########################################################
//                       COUNT
// ########################################################
// Determines how many values for each network input there 
// are.  Needed for input preconditioning.

vector<REAL> Neural_network::count() {
  
  vector<REAL> counts(values[0].size(), 0);
  for (int i=0; i<systems.size(); i++) {
      if (systems.at(i)->train == "train") {
        counts = vector_sum(counts, systems.at(i)->properties.count());
      }
  }
  return counts;

}

// ########################################################
// ########################################################




// ########################################################
//                       MEAN
// ########################################################
// Calculates the sum (not actually mean) of each input 
// over its training examples.  Needed for preconditioning.

vector<REAL> Neural_network::mean() {
  
  vector<REAL> means(values[0].size(), 0.0);
  for (int i=0; i<systems.size(); i++) {
    if (systems.at(i)->train == "train") {
        means = vector_sum(means, systems.at(i)->properties.mean());
    }
  }
  return means;

}

// ########################################################
// ########################################################




// ########################################################
//                       VARIANCE
// ########################################################
// Computes the sum of the squared differences between 
// the mean and each input value over this cores training
// examples.  Needed for preconditioning.

vector<REAL> Neural_network::variance(vector<REAL> means) {
  
  vector<REAL> variances(values[0].size(), 0.0);
  for (int i=0; i<systems.size(); i++) {
    if(systems.at(i)->train == "train") {
        variances = vector_sum(variances, systems.at(i)->properties.variance(means));
    }
  }
  return variances;
  
}

// ########################################################
// ########################################################




// ########################################################
//                       NORMALIZE
// ########################################################
// Preconditions inputs

void Neural_network::Normalize(vector<REAL> means, vector<REAL> variances) { 
  
  params.input_mean = means;
  params.input_variance = variances;
  
  for (int i=0; i<systems.size(); i++) {
    systems.at(i)->properties.Normalize(means, variances);
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       ATTACH_OPTIMIZER
// ########################################################
// Attach an optimizer object (needed for training) 
// to this network.

void Neural_network::attach_optimizer(Optimizer *opt) {
  
  this->optimizer = opt;
  
}

// ########################################################
// ########################################################




// ########################################################
//                       TRAIN
// ########################################################
// Actually train the network.

bool Neural_network::train() {
  
  optimizer->init((vector<vector<Network_node*> >*)(&nodes));
  
  int debug_count = 0;
  bool early_stop = false; 
  int nval = 0;
  
  while (!optimizer->is_converged() && !optimizer->is_checkpoint()) {
    REAL debug_error = 0.0;
    this->SSE = 0.0;
    this->SSE_mod = 0.0;
    REAL output = 0.0;
    REAL vSSE = 0.0;
    nval = 0.0;
    my_dOutput_dParameters.assign(my_dOutput_dParameters.size(), 0);
    //this->bernoulli_sample(this->params.dropout());
    for (int i_sys=0; i_sys<systems.size(); i_sys++) {
        if(systems.at(i_sys)->train == "train"){
            output = 0.0; 
            vector<REAL> temp_dOutput_dParameters(Nparams, 0.0);
            vector<REAL> dOut_dIn; 
            dOut_dIn.reserve(NNodes_max);
            vector<REAL> temp_dOut_dIn;
            temp_dOut_dIn.reserve(this->NNodes_max);

            while (systems[i_sys]->properties.iterate(values[0])) {
              for (int layer=0; layer<nodes.size()-1; layer++) {
                for (int node=0; node<nodes[layer].size(); node++) {
                  values[layer+1][node] = nodes[layer][node]->evaluate(values[layer]);
                }
              }

              output += nodes.back().at(0)->evaluate(values.back());

              //Back propagate for network derivative with respect to parameters
              ///////////////////////////////////////////////////////////////////

              // Start with the output node
              dOut_dIn.assign(NNodes.at(NNodes.size()-2),0.0);

              for (int i=deriv_offsets.back(); i<Nparams; i++) {
                temp_dOutput_dParameters.at(i) += nodes.back().at(0)->dOutput_dParameters(i-deriv_offsets.back());
              }
              for (int i=0; i<dOut_dIn.size(); i++) {
                dOut_dIn.at(i) = nodes.back().at(0)->dOutput_dInputs(i);
              }

              for (int layer=nodes.size()-1; layer>=1; layer--) {
                temp_dOut_dIn.assign(NNodes[layer-1],0.0);
                for (int node=0; node<NNodes[layer]; node++) {
                  for (int i=0; i<NNodes[layer-1]+1; i++) {
                    temp_dOutput_dParameters.at(deriv_offsets.at(nodes.at(layer-1).at(node)->index())+i)
                      += dOut_dIn.at(node)*nodes.at(layer-1).at(node)->dOutput_dParameters(i);
                  }
                  for (int i=0; i<NNodes[layer-1]; i++) {
                    temp_dOut_dIn.at(i) += dOut_dIn.at(node)*nodes.at(layer-1).at(node)->dOutput_dInputs(i);
                  }
                }
                dOut_dIn = temp_dOut_dIn;
              }

            }

            output *= systems[i_sys]->Prefactor*params.output_variance;
            output += params.output_mean;

            REAL dSSE_dOutput;
            REAL error = output - systems[i_sys]->properties.target();
            REAL lambda = optimizer->alpha();

            SSE += error*error;

            if (lambda == 0) {

              dSSE_dOutput = 2*error*params.output_variance;

            } else {

              REAL sign_phi = sign(error)*phi.at(i_sys);
              REAL cos_error = cos(error+sign_phi);

              SSE_mod += systems[i_sys]->Prefactor*error*error*(1+lambda*cos_error*cos_error);
              dSSE_dOutput = params.output_variance*systems[i_sys]->Prefactor*(2*lambda*(error*cos_error*cos_error - error*error*cos_error*sin(error+sign_phi)) + 2*error);
              phi.at(i_sys) += optimizer->dphi_dt();
            }

            for (int i=0; i<Nparams; i++) {
              my_dOutput_dParameters.at(i) += systems[i_sys]->Prefactor*dSSE_dOutput*temp_dOutput_dParameters.at(i);
            }
        } else if(systems.at(i_sys)->train == "val") {
            early_stop = true;
            output = 0.0; 
            while (systems[i_sys]->properties.iterate(values[0])) {
              for (int layer=0; layer<nodes.size()-1; layer++) {
                for (int node=0; node<nodes[layer].size(); node++) {
                  values[layer+1][node] = nodes[layer][node]->evaluate(values[layer]);
                }
              }
              output += nodes.back().at(0)->evaluate(values.back()); 
            }
            output *= systems[i_sys]->Prefactor*params.output_variance;
            output += params.output_mean;
            REAL error = output - systems[i_sys]->properties.target();
            //cout << output << " " << systems[i_sys]->properties.target() << endl;
            vSSE += error*error;
            nval++;
        }
    }
    optimizer->set_val_sse(vSSE,nval);
    if (SSE_mod) {
      optimizer->update_network(SSE_mod, my_dOutput_dParameters, SSE);
    } else {
      optimizer->update_network(SSE, my_dOutput_dParameters);
    }
    
  }

  return !optimizer->is_converged();
  
}

// ########################################################
// ########################################################


  
  
// ########################################################
//                       EVALUATE
// ########################################################
// Use the network to make predictions.

vector<REAL> Neural_network::evaluate() {
  
  results.clear();
  
  for (int i_sys=0; i_sys<systems.size(); i_sys++) {
    int count = 0;
    REAL output = 0.0;
    vector<REAL> in;

    REAL cheap_sum = 0.0;

    while (systems[i_sys]->properties.iterate(values[0])) {

      for (int layer=0; layer<nodes.size()-1; layer++) {
	for (int node=0; node<nodes[layer].size(); node++) {
	  values[layer+1][node] = nodes[layer][node]->evaluate(values[layer]);
	}
      }

      count++;
      output += params.output_variance*nodes.back().at(0)->evaluate(values.back());
      
    }
    
    results.push_back(output*systems[i_sys]->Prefactor+params.output_mean);
    
  }
  
  return results;
  
}

// ########################################################
// ########################################################


// ########################################################
//                       EVALUATE_MD
// ########################################################
// Use the network to make predictions for a potential
// (Structure-Energy relationship). 

REAL Neural_network::evaluate_MD(vector<REAL> &dE_dG) {
  
  REAL output = 0.0;
  vector<REAL> dOut_dIn;
  dOut_dIn.reserve(NNodes_max);
  vector<REAL> temp_dOut_dIn;
  temp_dOut_dIn.reserve(this->NNodes_max);
  
  if (!systems[0]->properties.iterate(values[0])) {
    systems[0]->properties.iterate(values[0]);
  }
  
  for (int i=0; i<values[0].size(); i++) {
    if (!(values[0][i]-params.input_mean[i])) {
      values[0][i] = 0;
    } else {
      values[0][i] = (values[0][i]-params.input_mean[i])/params.input_variance[i];
    }
  }

  for (int layer=0; layer<nodes.size()-1; layer++) {
    for (int node=0; node<nodes[layer].size(); node++) {
      values[layer+1][node] = nodes[layer][node]->evaluate(values[layer]);
    }
  }

  output += params.output_variance*nodes.back().at(0)->evaluate(values.back());
  
  
  dOut_dIn.assign(NNodes.at(NNodes.size()-2),0.0);
  
  
  for (int i=0; i<dOut_dIn.size(); i++) {
    dOut_dIn.at(i) = nodes.back().at(0)->dOutput_dInputs(i);
  }
  
  for (int layer=nodes.size()-1; layer>0; layer--) {
    temp_dOut_dIn.assign(NNodes[layer-1],0.0);
    for (int node=0; node<NNodes[layer]; node++) {
      for (int i=0; i<NNodes[layer-1]; i++) {
	temp_dOut_dIn.at(i) += dOut_dIn.at(node)*nodes.at(layer-1).at(node)->dOutput_dInputs(i);
      }
    }
    dOut_dIn = temp_dOut_dIn;
  }
  dE_dG = dOut_dIn;
  
  for (int i=0; i<dE_dG.size(); i++) {
    dE_dG[i] *= params.output_variance/params.input_variance[i];
  }
  
  return output;
  
}

// ########################################################
// ########################################################





// ########################################################
//                       EVALUATE_MD
// ########################################################
// Use network to make predictions for a potential
// (struture-energy relationship).

REAL Neural_network::evaluate_MD(int i_sys) {
  
  REAL output;
  
  if (!systems[i_sys]->properties.iterate(values[0])) {
    systems[i_sys]->properties.iterate(values[0]);
  }
  
  for (int i=0; i<values[0].size(); i++) {
    if (!(values[0][i]-params.input_mean[i])) {
      values[0][i] = 0;
    } else {
      values[0][i] = (values[0][i] - params.input_mean[i])/params.input_variance[i];
    }
  }

  for (int layer=0; layer<nodes.size()-1; layer++) {
    for (int node=0; node<nodes[layer].size(); node++) {
      values[layer+1][node] = nodes[layer][node]->evaluate(values[layer]);
    }
  }
  
  output = params.output_variance*nodes.back().at(0)->evaluate(values.back());
  
  return output;
  
}

// ########################################################
// ########################################################




// ########################################################
//                       TRAIN_MD
// ########################################################
// Actually train the network for a potential (structure-
// energy relationship).

REAL Neural_network::train_MD(int i_sys) {
  
  if (i_sys != last_i_sys) {
    my_dOutput_dParameters.assign(my_dOutput_dParameters.size(), 0);
    last_i_sys = i_sys;
  }
  
  
  REAL output = 0.0; 
  
  vector<REAL> temp_dOutput_dParameters(Nparams, 0.0);
  vector<REAL> dOut_dIn; 
  dOut_dIn.reserve(NNodes_max);
  vector<REAL> temp_dOut_dIn;
  temp_dOut_dIn.reserve(this->NNodes_max);
  
  if (!systems[i_sys]->properties.iterate(values[0])) {
    systems[i_sys]->properties.iterate(values[0]);
    my_dOutput_dParameters.assign(my_dOutput_dParameters.size(), 0);
  }
  
  for (int i=0; i<values[0].size(); i++) {
    if (!(values[0][i]-params.input_mean[i])) {
      values[0][i] = 0;
    } else {
      values[0][i] = (values[0][i] - params.input_mean[i])/params.input_variance[i];
    }
  }
  
  for (int layer=0; layer<nodes.size()-1; layer++) {
    for (int node=0; node<nodes[layer].size(); node++) {
      values[layer+1][node] = nodes[layer][node]->evaluate(values[layer]);
    }
  }
  
  output = params.output_variance*nodes.back().at(0)->evaluate(values.back())*systems[i_sys]->Prefactor;
  
  //Back propagate for network derivative with respect to parameters
  ///////////////////////////////////////////////////////////////////

  // Start with the output node
  dOut_dIn.assign(NNodes.at(NNodes.size()-2),0.0);
  
  for (int i=deriv_offsets.back(); i<Nparams; i++) {
    temp_dOutput_dParameters[i] += nodes.back()[0]->dOutput_dParameters(i-deriv_offsets.back());
  }
  for (int i=0; i<dOut_dIn.size(); i++) {
    dOut_dIn[i] = nodes.back()[0]->dOutput_dInputs(i);
  }
  
  
  for (int layer=nodes.size()-1; layer>=1; layer--) {
    temp_dOut_dIn.assign(NNodes[layer-1],0.0);
    for (int node=0; node<NNodes[layer]; node++) {
      for (int i=0; i<NNodes[layer-1]+1; i++) {
	temp_dOutput_dParameters[deriv_offsets[nodes[layer-1][node]->index()]+i] += dOut_dIn[node]*nodes[layer-1][node]->dOutput_dParameters(i);
      }
      for (int i=0; i<NNodes[layer-1]; i++) {
	temp_dOut_dIn[i] += dOut_dIn[node]*nodes[layer-1][node]->dOutput_dInputs(i);
      }
    }
    dOut_dIn = temp_dOut_dIn;
  }
  
  for (int i=0; i<Nparams; i++) {
    my_dOutput_dParameters[i] += params.output_variance*systems[i_sys]->Prefactor*temp_dOutput_dParameters[i];
  }
  
  return output;
  
}

// ########################################################
// ########################################################




// ########################################################
//                       INIT_OPTIMIZER
// ########################################################
// Initialize the attached optimizer.  Needed for training.

void Neural_network::init_optimizer(Optimizer* opt) {

  vector<Network_node*> temp_nodes;
  for (int i=0; i<nodes.size(); i++) {
    for (int j=0; j<nodes[i].size(); j++) {
      temp_nodes.push_back(nodes[i][j]);
    }
  }
  opt->init(temp_nodes);
}
  
// ########################################################
// ########################################################




// ########################################################
//                       SCALE_GRADIENT
// ########################################################
// Multiply the gradient by a scalar.

void Neural_network::scale_gradient(REAL scale_factor) {

  for (int i=0; i<my_dOutput_dParameters.size(); ++i) {
    my_dOutput_dParameters[i] *= scale_factor;
  }

}

// ########################################################
// ########################################################





// ########################################################
//                       GET_STATE
// ########################################################
// Return the parameters of the network. Needed for parallel
// syncronization.

vector<REAL> Neural_network::get_state() {
  
  vector<REAL> output;

  for (int layer=0; layer<nodes.size(); layer++) {
    for (int node=0; node<nodes.at(layer).size(); node++) {
      vector<REAL> node_params = nodes.at(layer).at(node)->get_parameters();
      for (vector<REAL>::iterator it=node_params.begin(); it!=node_params.end(); ++it) {
     	output.push_back(*it);
      }
    }
  }
  return output;
}

// ########################################################
// ########################################################




// ########################################################
//                       CENTER_OUTPUTS
// ########################################################
//

void Neural_network::center_outputs(REAL bias_correction) {
  vector<REAL> params = nodes.back().at(0)->get_parameters();
  params.back() += bias_correction;
  nodes.back().at(0)->set_parameters(params);
}

// ########################################################
// ########################################################




// ########################################################
//                       SET_STATE
// ########################################################
// Set the network parameters. Needed to syncronize across
// cores in a parallel run.

void Neural_network::set_state(const vector<REAL> &new_state) {

  if (new_state.size() != Nparams) { ERROR("Wrong number of parameters in set_state");}
  vector<REAL>::const_iterator begin = new_state.begin(), end = new_state.begin();
  for (int layer=0; layer<nodes.size(); layer++) {
    for (int node=0; node<nodes.at(layer).size(); node++) {
      end = begin+NNodes.at(layer)+1;
      nodes.at(layer).at(node)->set_parameters(vector<REAL>(begin,end));
      begin = end;      
    }
  }
}

// ########################################################
// ########################################################




// ########################################################
//                       PRINT
// ########################################################
// Write network details to checkpoint file

void Neural_network::print(ostream &output) {
  
  this->params.print(output);
  
  for (int layer=0; layer<nodes.size(); layer++) {
    output << "[[ layer "<<layer<<" ]] " << endl;
    for (int node=0; node<nodes.at(layer).size(); node++) {
        /*nodes[layer][node]->print(output);*/
        
        if (layer == 0 || layer == nodes.size() - 1) {
            nodes[layer][node]->print(output,1.0);
        }  else {
            nodes[layer][node]->print(output,this->params.dropout());
        }
    }
        
  }
  
}

void Neural_network::get_training_set() {
    if (this->params.SGD()) {
        this->training_set.clear();
        for(int i = 0; i < this->params.SGD_cnt(); i++ ) {
            this->training_set.push_back(rand() % systems.size());
        }
    }
}

void Neural_network::set_dropout(vector <double> dropout) {
    int cnt = 0;
    for (int layer=0; layer<nodes.size()-1; layer++) {
        for (int node=0; node<nodes[layer].size(); node++) {
                nodes[layer][node]->dropout = dropout[cnt];
                cnt++;
        }
    }
    
}
// ########################################################
// ########################################################



