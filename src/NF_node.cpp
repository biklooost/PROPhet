#include "NF_node.h"
#include <cstdlib>
#include <time.h>
#include <iomanip>
#include <sstream>

Network_function_node::Network_function_node(int Ninputs) : Ninputs(Ninputs), my_Nparameters(3*Ninputs+1) { 

  for (int i=0; i<my_Nparameters; i++) {
    Parameters.push_back(complex<REAL>(RAND::Normal(0.0,1.0),RAND::Normal(0.0,1.0)));
  }
  
  for (int i=0; i<my_Nparameters; i++) {
    my_dOutput_dParameters.push_back(complex<REAL>(0.0,0.0));
  }
  for (int i=0; i<Ninputs; i++) {
    my_dOutput_dInputs.push_back(complex<REAL>(0.0,0.0));
  }
  
  static int current_index = 0;
  my_index = current_index++;

  for (int i=0; i<Ninputs; i++) {
    last_imag_ln.push_back(0.0);
  }
  
}



Network_function_node::~Network_function_node() { }

void Network_function_node::set_parameters(vector<REAL> new_parameters) {
  if (new_parameters.size() != 2*my_Nparameters) { ERROR("Parameters set wrong"); }
  for (int i=0; i<my_Nparameters; i++) {
    Parameters.at(i).real() = new_parameters.at(2*i);
    Parameters.at(i).imag() = new_parameters.at(2*i+1);
  }

}

vector<REAL> Network_function_node::get_parameters() {
  
  vector<REAL> params;
  params.reserve(2*my_Nparameters);
  for (int i=0; i<my_Nparameters; i++) {
    params.push_back(Parameters.at(i).real());
    params.push_back(Parameters.at(i).imag());
  }
  return params;

}

complex<REAL> Network_function_node::evaluate(const vector<complex<REAL> > &in) {
  
  if (in.size() != Ninputs) { ERROR("Node evaluation went wrong");}
  
  complex<REAL> output;
  complex<REAL> expn, expn1, ln;
  my_dOutput_dInputs.assign(Ninputs,complex<REAL>(0.0,0.0));
  
  for (int i=0; i<Ninputs; i++) {

    expn = exp(in.at(i));
    expn1 = exp(-1.0e-3*pow(in.at(i).real(),4));
	
    ln = log(in.at(i));
    
    // C1*exp(x)
    output += Parameters.at(3*i)*expn1*expn;
    my_dOutput_dParameters.at(3*i) = expn1*expn;
    my_dOutput_dInputs.at(i) += Parameters.at(3*i)*expn1*expn*(1 - 4.0e-3*pow(in.at(i).real(),3));
    
    // C2*ln(x)
    output += Parameters.at(3*i+1)*ln;
    my_dOutput_dParameters.at(3*i+1) = ln;
    my_dOutput_dInputs.at(i) += Parameters.at(3*i+1)/in.at(i);
    
    // C3*x
    output += Parameters.at(3*i+2)*in.at(i);
    my_dOutput_dParameters.at(3*i+2) = in.at(i);
    my_dOutput_dInputs.at(i) += Parameters.at(3*i+2);
    
  }  
  
  output += Parameters.back();
  my_dOutput_dParameters.back() = complex<REAL>(1.0,0.0);
  
  output = complex<REAL>(10.0,0.0)*tanh(output/complex<REAL>(10.0,0.0));

  return output;

}
  
  
void Network_function_node::read(istream &input) {
  string line, skip, test;
  istringstream Line;
  complex<REAL> param;
  
  this->Parameters.clear();
  this->my_dOutput_dParameters.clear();
  this->my_dOutput_dInputs.clear();
  
  getline(input, line);
  Line.str(line);
  Line >> skip >> test;
  if (!test.compare("layer")) {
    getline(input, line);
    Line.str(line);
    Line >> skip >> skip >> this->my_index;
  } else {
    Line >> this->my_index;
  }
  
  getline(input, line);
  Line.str(line);
  while (Line >> param) {
    Parameters.push_back(param);
    my_dOutput_dParameters.push_back(0);
    my_dOutput_dInputs.push_back(0);
  }
  this->Ninputs = Parameters.size();
}

void Network_function_node::print(ostream &output) {
  
  output << "  [ node "<<this->index()<<" ]  " << endl;
  
  for (int i=0; i<Parameters.size(); i++) {
    output << "   " << setprecision(12) << Parameters.at(i) << "  ";
  }
  output << endl;

}

