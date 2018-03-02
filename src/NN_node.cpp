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
// A class defining the behavior of a node in a neural network.
// ####################################################################



#include "Network_node.h"
#include "NN_node.h"
#include <cstdlib>
#include <time.h>
#include <iomanip>
#include <sstream>

Tables Neural_network_node::Table = Tables("tanh");



// ########################################################
//                       Constructor
// ########################################################
//

Neural_network_node::Neural_network_node(int Ninputs, string which_transfer_function) : my_transfer_function(which_transfer_function),my_Nparameters(Ninputs+1), Ninputs(Ninputs)
{

  this->set_transfer_function(which_transfer_function);

  for (int i=0; i<Ninputs; i++) {
    Parameters.push_back(RAND::Normal(0,10.5/sqrt(Ninputs)));
  }

  // add a bias
  Parameters.push_back(RAND::Normal(0,1));

  for (int i=0; i<my_Nparameters; i++) {
    my_dOutput_dParameters.push_back(0);
  }
  for (int i=0; i<Ninputs; i++) {
    my_dOutput_dInputs.push_back(0);
  }
  this->dropout = 1.0;
}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

Neural_network_node::~Neural_network_node() { }

// ########################################################
// ########################################################



// ########################################################
//                       SET_PARAMETERS
// ########################################################
// Sets the node's parameter values.  Needed for parallel
// syncronization.

void Neural_network_node::set_parameters(vector<REAL> new_parameters)
{

  if (Parameters.size() != new_parameters.size()) {
    ERROR("Parameters not set correctly");
  }
  for (int i=0; i<Parameters.size(); i++) {
    Parameters.at(i) = new_parameters.at(i);
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       SET_TRANSFER_FUNCTION
// ########################################################
// Sets a node's transfer function.

void Neural_network_node::set_transfer_function(string new_function)
{

  if (new_function == "tanh_spline") {

    transfer_function = &Neural_network_node::table_function;

  } else if (new_function == "tanh") {

    transfer_function = &Neural_network_node::tanh;

  } else if (new_function == "tanh_crazy") {

    transfer_function = &Neural_network_node::tanh_crazy;

  } else if (new_function == "cubic") {

    transfer_function = &Neural_network_node::cubic;

  } else if (new_function == "pulse") {

    transfer_function = &Neural_network_node::pulse;

  } else if (new_function == "linear") {

    transfer_function = &Neural_network_node::linear;

  } else if (new_function == "sin") {

    transfer_function = &Neural_network_node::sin;

  } else if (new_function == "periodic_tanh") {

    transfer_function = &Neural_network_node::periodic_tanh;

  } else if (new_function == "random") {

    transfer_function = &Neural_network_node::random;

  } else if (new_function == "relu") {

      transfer_function = &Neural_network_node::relu;

  } else {

    ERROR("Unrecognized transfer function '"+new_function+"'");

  }

  my_transfer_function = new_function;

}

// ########################################################
// ########################################################




// ########################################################
//                       READ
// ########################################################
// Reads node's parameters from checkpoint file.

void Neural_network_node::read(istream &input)
{
  string line, skip, test;
  istringstream Line;
  REAL param;

  this->Parameters.clear();
  this->my_dOutput_dParameters.clear();
  this->my_dOutput_dInputs.clear();

  getline(input, line);
  Line.str(line);
  Line >> skip >> test;
  if (!test.compare("layer")) {
    getline(input, line);
    Line.str(line);
    Line >> skip >> skip >> this->my_index >> skip >> this->my_transfer_function;
  } else {
    Line >> this->my_index >> skip >> this->my_transfer_function;
  }

  getline(input, line);
  Line.str(line);
  while (Line >> param) {
    Parameters.push_back(param);
    my_dOutput_dParameters.push_back(0);
    my_dOutput_dInputs.push_back(0);
  }
  this->Ninputs = Parameters.size();

  getline(input, line);
  Line.str(line);
  Line >> param;
  this->Parameters.push_back(param);

}

// ########################################################
// ########################################################




// ########################################################
//                       PRINT
// ########################################################
// Writes node's parametrs to checkpoint file.

void Neural_network_node::print(ostream &output,REAL p)
{

  output << "  [ node "<<this->index()<<" ]  " << my_transfer_function << endl;

  for (int i=0; i<Parameters.size()-1; i++) {
    output << "   " << setprecision(18) << p*Parameters.at(i) << "  ";
  }
  output << endl;
  output << "   " << setprecision(18) << p*Parameters.back() << endl;

}

/*void Neural_network_node::print(ostream &output) {

  output << "  [ node "<<this->index()<<" ]  " << my_transfer_function << endl;

  for (int i=0; i<Parameters.size()-1; i++) {
    output << "   " << setprecision(18) << Parameters.at(i) << "  ";
  }
  output << endl;
  output << "   " << setprecision(18) << Parameters.back() << endl;

}*/
// ########################################################
// ########################################################


