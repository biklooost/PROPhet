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





#ifndef __NEURAL_NETWORK_NODE
#define __NEURAL_NETWORK_NODE

#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include "Tables.h"
#include "Common.h"
#include "Network_node.h"


using namespace std;


class Neural_network_node : public Network_node
{

public:

  Neural_network_node(int Ninputs, string which_transfer_function);
  ~Neural_network_node();


  inline REAL evaluate(const vector<REAL> &in)
  {
    if (Ninputs != in.size()) {
      ERROR("Wrong number of inputs passed to node");
    }
    REAL X = 0.0, dT_dX, out;
    for (int i=0; i<Ninputs; i++) {
      X += in[i]*Parameters[i];
    }
    X += Parameters[Ninputs];
    out = (this->*transfer_function)(X, &dT_dX);
    for (int i=0; i<Ninputs; i++) {
      my_dOutput_dParameters[i] = this->dropout*dT_dX*in[i];
      my_dOutput_dInputs[i] = this->dropout*dT_dX*Parameters[i];
    }
    my_dOutput_dParameters[Ninputs] = dT_dX;
    return this->dropout*out;
  }

  inline REAL dout_din(const vector<REAL> &in, REAL &deriv) { }

  inline vector<REAL> dOutput_dParameters()
  {
    return my_dOutput_dParameters;
  }
  inline REAL dOutput_dParameters(int i)
  {
    return my_dOutput_dParameters[i];
  }

  inline vector<REAL> dOutput_dInputs()
  {
    return my_dOutput_dInputs;
  }
  inline REAL dOutput_dInputs(int i)
  {
    return my_dOutput_dInputs[i];
  }

  inline int set_index(int new_index)
  {
    my_index = new_index;
  }
  inline int index()
  {
    return my_index;
  }

  virtual void set_parameters(vector<REAL> new_params);
  virtual inline vector<REAL> get_parameters()
  {
    return Parameters;
  }

  virtual void print(ostream &stream = std::cout,REAL p=1.0);
  //virtual void print(ostream &stream = std::cout);
  virtual void read(istream &input);

  void set_transfer_function(string new_transfer_function);

  virtual inline int Nparameters()
  {
    return my_Nparameters;
  }

  REAL dropout;

private:

  int my_Nparameters;
  int Ninputs;
  vector<REAL> my_dOutput_dInputs;
  vector<REAL> my_dOutput_dParameters;


  REAL (Neural_network_node::*transfer_function)(REAL,REAL*);


  int my_index;

  inline REAL table_function(REAL in, REAL *deriv)
  {

    if (deriv != NULL) {
      (*deriv) = Table.derivative(in);//3*Table(row,0)*x2+2*Table(row,1)*x+Table(row,2);
      return Table(in);//Table(row,0)*x3+Table(row,1)*x2+Table(row,2)*x+Table(row,3);
    } else {
      return Table(in);//Table(row,0)*x3+Table(row,1)*x2+Table(row,2)*x+Table(row,3);
    }
  }

  inline REAL periodic_tanh(REAL in, REAL *deriv)
  {
    REAL inner = 4.0*std::sin(in);
    if (deriv != NULL) {
      (*deriv) = 4.0*std::cos(in)/pow(std::cosh(inner),2);
    }
    return std::tanh(inner);
  }


  inline REAL random(REAL in, REAL *deriv)
  {
    if (deriv != NULL) {
      (*deriv) = RAND::Uniform(10);
    }
    return RAND::Uniform(10);
  }

  inline REAL sin(REAL in, REAL *deriv)
  {
    if (deriv != NULL) {
      (*deriv) = std::cos(in);
    }
    return std::sin(in);
  }

  inline REAL pulse(REAL in, REAL *deriv)
  {
    REAL a=1;
    REAL b=1;
    REAL c=1;
    if (deriv != NULL) {
      (*deriv) = c*exp(-b*in*in)*(-2*b*in*std::sin(a*in)+a*std::cos(a*in));
    }
    return c*std::sin(a*in)*exp(-b*in*in);
  }

  // Standard library tanh transfer function
  inline REAL tanh(REAL in,REAL *deriv)
  {
    if (abs(in) > 6) {
      int sign = ((REAL)(0) < in) - (in < (REAL)(0));
      REAL exptl = exp(-2*sign*in);
      if (deriv != NULL) {
        (*deriv) = 4.0*exptl*(1-2*exptl);
      }
      return sign*(1-2*exptl*(1-exptl));
    }
    if (deriv != NULL) {
      (*deriv) = 1.0/pow(std::cosh(in),2);
    }
    return std::tanh(in);
  }

  // Crazy tanh
  inline REAL tanh_crazy(REAL in,REAL *deriv)
  {
    if (deriv != NULL) {
      (*deriv) = (1-pow(0.1*in,4))/pow(std::cosh(in),2) - std::tanh(in)*0.1*4*pow(0.1*in,3);
    }
    return (1-pow(0.1*in,4))*std::tanh(in);
  }

  //Cubic function
  inline REAL cubic(REAL in, REAL *deriv)
  {
    REAL x2 = in*in;
    if (deriv != NULL) {
      (*deriv) = 3*x2;
    }
    return x2*in;
  }



  // Linear transfer function
  inline REAL linear(REAL in,REAL *deriv)
  {
    if (deriv != NULL) {
      (*deriv) = 1.0;
      return in;
    } else {
      return in;
    }
  }

  inline REAL relu(REAL in, REAL *deriv)
  {
    if ( deriv != NULL) {
      if (in <= 0.0) {
        (*deriv) = 0.0;
      } else {
        (*deriv) = 1.0;
      }
    }
    if (in <= 0.0) {
      return 0.0;
    } else {
      return in;
    }
  }



  vector<REAL> Parameters;
  REAL bias;

  string my_transfer_function;

  static Tables Table;



};



#endif
