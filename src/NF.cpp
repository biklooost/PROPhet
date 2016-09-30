#include "NF.h"








// ########################################################
//                       Constructor
// ########################################################
//

Network_function::Network_function(const vector<System*> &in_systems, Functional_params network_details) : systems(in_systems) { 
  
  vector<complex<REAL> > input_layer(network_details.Ninput_nodes());
  values.push_back(input_layer);
  NNodes.push_back(network_details.Ninput_nodes());
  
  for (int i=0; i<network_details.Nlayers(); i++) {
    vector<complex<REAL> > layer_values(network_details.NNodes(i),0.0);
    vector<Network_function_node *> layer_nodes;
    for (int j=0; j<network_details.NNodes(i); j++) {
      layer_nodes.push_back(new Network_function_node(values.at(i).size()));
    }
    values.push_back(layer_values);
    nodes.push_back(layer_nodes);
    NNodes.push_back(network_details.NNodes(i));
  }
  vector<Network_function_node*> output_layer;
  output_layer.push_back(new Network_function_node(values.back().size()));
  nodes.push_back(output_layer);
  NNodes.push_back(1);


  Nparams = 0;
  for (int layer=0; layer<network_details.Nlayers(); layer++) {
    for (int node=0; node<network_details.NNodes(layer); node++) {
      deriv_offsets.push_back(Nparams);
      Nparams = Nparams + 3*values.at(layer).size() + 1;
    }
  }
  deriv_offsets.push_back(Nparams);
  Nparams = Nparams + 3*network_details.NNodes(network_details.Nlayers()-1)+1;

  my_dOutput_dParameters.assign(2*Nparams, 0.0);


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
    
}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

Network_function::~Network_function() { 
  
  for (int i=0; i<nodes.size(); i++) {
    for (int j=0; j<nodes.at(i).size(); j++) {
      delete nodes[i][j];
    }
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       INIT_OPTIMIZER
// ########################################################
// Initilizes the attached optimizer.

void Network_function::init_optimizer(Optimizer* opt) {
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
//                       CENTER_OUTPUTS
// ########################################################
//

void Network_function::center_outputs(REAL bias_correction) {
  vector<REAL> params = nodes.back().at(0)->get_parameters();
  params.back() += bias_correction;
  nodes.back().at(0)->set_parameters(params);
}

// ########################################################
// ########################################################



  

// ########################################################
//                       LOAD
// ########################################################
// Reads the current network from a checkpoint file.

void Network_function::load(ifstream &input_file) { 
  
  string line,junk;
  complex<REAL> param;
  vector<REAL>  parameters;

  for (int layer=0; layer<nodes.size(); layer++) {
    
    getline(input_file,line);
    while (!input_file.eof() && line.find("layer "+layer) == std::string::npos) {
      getline(input_file, line);
    }
    
    for (int node=0; node<nodes.at(layer).size(); node++) {
      istringstream Line;
      getline(input_file,line);
      while (!input_file.eof() && line.find("[ node ") == std::string::npos) {
	getline(input_file, line);
      }
      Line.clear();
      getline(input_file,line);
      Line.str(line);
      while (Line >> param) {
      	parameters.push_back(param.real());
	parameters.push_back(param.imag());
      }
      Line.clear();
    }
  }
  
  if (parameters.size() != 2*Nparams) { ERROR("Wrong number of parameters found in input file"); }

  this->set_state(parameters);

}

// ########################################################
// ########################################################





void Network_function::set_parameters() { }

vector<REAL> Network_function::count() {

}

vector<REAL> Network_function::mean() { }
vector<REAL> Network_function::variance(vector<REAL> means) { }
void Network_function::Normalize(vector<REAL> means, vector<REAL> variances) { }

void Network_function::attach_optimizer(Optimizer *opt) {
  
  this->optimizer = opt;
  
}


// ########################################################
//                       TRAIN
// ########################################################
// Trains the network.

bool Network_function::train() {

  optimizer->init( (vector<vector<Network_node*> >*)(&nodes));

  int debug_count = 0;

  while (!optimizer->is_converged() && !optimizer->is_checkpoint()) {

    this->SSE = 0.0;
    this->SSE_mod = 0.0;
    complex<REAL> output;
    my_dOutput_dParameters.assign(my_dOutput_dParameters.size(), 0.0);
    
    for (int i_sys=0; i_sys<systems.size(); i_sys++) {
      
      output = complex<REAL>(0.0,0.0); 
      vector<complex<REAL> > temp_dOutput_dParameters(Nparams, complex<REAL>(0.0,0.0));
      vector<complex<REAL> > dOut_dIn; 
      dOut_dIn.reserve(NNodes_max);
      vector<complex<REAL> > temp_dOut_dIn;
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
	dOut_dIn.assign(NNodes.at(NNodes.size()-2),complex<REAL>(0.0,0.0));

	for (int i=deriv_offsets.back(); i<Nparams; i++) {
	  temp_dOutput_dParameters.at(i) += nodes.back().at(0)->dOutput_dParameters(i-deriv_offsets.back());
	}
	for (int i=0; i<dOut_dIn.size(); i++) {
	  dOut_dIn.at(i) = nodes.back().at(0)->dOutput_dInputs(i);
	}
	
	for (int layer=nodes.size()-1; layer>=1; layer--) {
	  temp_dOut_dIn.assign(NNodes[layer-1],complex<REAL>(0.0,0.0));
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
      
      complex<REAL> dSSE_dOutput;
      
      complex<REAL> error = output - complex<REAL>(systems[i_sys]->properties.target(), 0.0);

      REAL lambda = optimizer->alpha();
      
      SSE += error.real()*error.real() + error.imag()*error.imag();
      
      if (lambda == 0) {

	dSSE_dOutput = conj(error);
	
      } else {
	
	REAL error1 = sqrt(error.real()*error.real() + error.imag()*error.imag());
	REAL error2 = error1*error1;
	REAL sign_phi = phi.at(i_sys);
	REAL cos_error = cos(error1 + phi.at(i_sys));
	
	SSE_mod += error2*(1+lambda*cos_error*cos_error);
	dSSE_dOutput = 2*lambda*(error1*cos_error*cos_error - error2*cos_error*sin(error1+sign_phi) + 2*error2);
	
	phi.at(i_sys) += optimizer->dphi_dt();
	  
      }

      for (int i=0; i<Nparams; i++) {
	complex<REAL> deriv = dSSE_dOutput*temp_dOutput_dParameters.at(i);
	my_dOutput_dParameters.at(2*i) += 2*deriv.real();
	my_dOutput_dParameters.at(2*i+1) -= 2*deriv.imag();
      }
      
    }
    
    optimizer->update_network(SSE, my_dOutput_dParameters);
    
  }
  
  return !optimizer->is_converged();
  
}

// ########################################################
// ########################################################




// ########################################################
//                       EVALUATE
// ########################################################
// Fires the network to make predictions.

vector<REAL> Network_function::evaluate() {
  
  vector<REAL> results;

  for (int i_sys=0; i_sys<systems.size(); i_sys++) {
    int count = 0;
    complex<REAL> output = complex<REAL>(0.0,0.0);
    vector<complex<REAL> > in;

    while (systems[i_sys]->properties.iterate(values[0])) {

      for (int layer=0; layer<nodes.size()-1; layer++) {
	for (int node=0; node<nodes[layer].size(); node++) {
	  values[layer+1][node] = nodes[layer][node]->evaluate(values[layer]);
	}
      }
      count++;
      output += nodes.back().at(0)->evaluate(values.back());
      
    }

    results.push_back(output.real());
    
  }
  
  return results;
  
}

// ########################################################
// ########################################################





REAL Network_function::evaluate_MD(int i_sys) {

}

REAL Network_function::evaluate_MD(vector<REAL> &dE_dG) {

}

REAL Network_function::train_MD(int i_sys) {

}


void Network_function::scale_gradient(REAL scale_factor) {

}


vector<REAL> Network_function::get_state() {
  
  vector<REAL> output;
  output.reserve(Nparams);

  for (int layer=0; layer<nodes.size(); layer++) {
    for (int node=0; node<nodes.at(layer).size(); node++) {
      vector<REAL> node_params = nodes.at(layer).at(node)->get_parameters();
      for (vector<REAL>::iterator it=node_params.begin(); it<node_params.end(); it++) {
     	output.push_back(*it);
      }
    }
  }
  return output;
}

void Network_function::set_state(const vector<REAL>& new_state) {

  vector<REAL>::const_iterator begin = new_state.begin(), end;
  for (int layer=0; layer<nodes.size(); layer++) {
    for (int node=0; node<nodes.at(layer).size(); node++) {
      end = begin + nodes.at(layer).at(node)->Nparameters();
      nodes.at(layer).at(node)->set_parameters(vector<REAL>(begin,end));
      begin = end;
    }
  }

}



void Network_function::print(ostream &output) {
  
  for (int layer=0; layer<nodes.size(); layer++) {
    output << "[[ layer "<<layer<<" ]] " << endl;
    for (int node=0; node<nodes.at(layer).size(); node++) {
      nodes[layer][node]->print(output);
    }
  }
  
}

