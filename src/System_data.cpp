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
// This is a low-level class that handles passing the requested
// properties of a system to a machine-learning algorithm as input. As
// iterating over system properties is the very inner loop of training
// a mapping, it needs to be as efficient as possible.
// ####################################################################


#include "System_data.h"




// ########################################################
//                       Constructor
// ########################################################
//

System_data::System_data()
{

  is_locked = false;

}

// ########################################################
// ########################################################






// ########################################################
//                       Destructor
// ########################################################
//

System_data::~System_data()
{

}

// ########################################################
// ########################################################




// ########################################################
//                       ITERATE
// ########################################################
// Performs the inner loop over the properties used in the
// fit.

bool System_data::iterate(vector<complex<REAL> > &output_data)
{

  if (indices.empty()) {
    my_N = data.size();
    for (int i=0; i<my_N; i++) {
      indices.push_back(0);
      N_max.push_back(data.at(i)->size());
      #if __cplusplus >= 201103L
        output_data.at(i).real(data.at(i)->at(0));// = data.at(i)->at(0);
        output_data.at(i).imag(0.0);// = 0.0;
      #else
        output_data.at(i).real() = data.at(i)->at(0);
        output_data.at(i).imag() = 0.0;
      #endif
    }
    return true;
  }

  indices.at(my_N-1)++;
  for (int i=my_N-1; (i>0)&&(indices[i]==N_max[i]); i--) {
    indices[i] = 0;
    indices[i-1]++;
  }

  if (indices[0] == N_max[0]) {
    for (int i=0; i<my_N; i++) {
      indices[i] = 0;
    }
    indices.back() = -1;
    return false;
  }

  for (int i=0; i<my_N; i++) {
    #if __cplusplus >= 201103L
      output_data[i].real(data[i]->at(indices[i]));
      output_data[i].imag(0.0);
    #else
      output_data[i].real() = data[i]->at(indices[i]);
      output_data[i].imag() = 0.0;
    #endif
  }

  return true;
}

// ########################################################
// ########################################################




// ########################################################
//                       ITERATE
// ########################################################
// Performs the inner loop over the properties used in the
// fit.

bool System_data::iterate(vector<REAL>& output_data)
{


  if (indices.empty()) {
    my_N = data.size();
    for (int i=0; i<my_N; i++) {
      indices.push_back(0);
      N_max.push_back(data.at(i)->size());
      output_data.at(i) = data.at(i)->at(0);
    }
    return true;
  }

  if (this->is_locked) {

    for (int i=0; i<my_N; i++) {
      indices[i]++;
    }

  } else {

    indices.at(my_N-1)++;
    for (int i=my_N-1; (i>0)&&(indices[i]==N_max[i]); i--) {
      indices[i] = 0;
      indices[i-1]++;
    }

  }


  if (indices[0] == N_max[0]) {
    for (int i=0; i<my_N; i++) {
      if (is_locked) {
        indices[i] = -1;
      } else {
        indices[i] = 0;
      }
    }
    indices.back() = -1;
    return false;
  }

  for (int i=0; i<my_N; i++) {
    output_data[i] = data[i]->at(indices[i]);
  }

  return true;
}

// ########################################################
// ########################################################



