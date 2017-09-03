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



// This is an interface to a node in a network type machine-learning approach. Currently it is used for neural networks and network functions. It is intended to abstract away details of the particular behavior of a node in the given machine-learning approach.

#ifndef __NETWORK_NODE
#define __NETWORK_NODE

#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include "Common.h"


using namespace std;


class Network_node {
  
public:
  
  Network_node() { }
  ~Network_node() { }
  
  virtual void set_parameters(vector<REAL> new_params) = 0;
  virtual vector<REAL> get_parameters() = 0;
  
//  virtual void print(ostream &stream = std::cout) = 0;
  virtual void print(ostream &steam = std::cout, REAL p=1.0) = 0;
  virtual void read(istream &input) = 0;
  
  virtual int Nparameters() = 0;
  
};



#endif
