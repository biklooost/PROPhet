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
// 
// ####################################################################


// A class defining the behavior of a node in a network function.


#ifndef __NETWORK_FUNCTION_NODE
#define __NETWORK_FUNCTION_NODE

#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <complex>

#include "Tables.h"
#include "Common.h"
#include "Network_node.h"

using namespace std;


class Network_function_node : public Network_node {

public:
  
  Network_function_node(int Ninputs);
  ~Network_function_node();
  
  
  complex<REAL> evaluate(const vector<complex<REAL> > &in);

  inline complex<REAL> dout_din(const vector<REAL> &in, REAL &deriv) { }
  
  inline vector<complex<REAL> > dOutput_dParameters() { return my_dOutput_dParameters; }
  inline complex<REAL> dOutput_dParameters(int i) { return my_dOutput_dParameters[i]; }
  
  inline vector<complex<REAL> > dOutput_dInputs() { return my_dOutput_dInputs; }
  inline complex<REAL> dOutput_dInputs(int i) { return my_dOutput_dInputs[i]; }
  
  inline int index() { return my_index; }
  
  virtual void set_parameters(vector<REAL> new_params);
  virtual vector<REAL> get_parameters();
  
  virtual void print(ostream &stream = std::cout, REAL p = 1.0);
  virtual void read(istream &input);
  
  virtual inline int Nparameters() { return 2*my_Nparameters; }
    
 private:
  
  int my_Nparameters;
  int Ninputs;
  vector<complex<REAL> > my_dOutput_dInputs;
  vector<complex<REAL> > my_dOutput_dParameters;
  
  int my_index;
  vector<REAL> last_imag_ln;
  
  
  vector<complex<REAL> > Parameters;
    
};



#endif
