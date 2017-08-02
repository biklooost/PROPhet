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
// Class to handle data that are on a grid. Mainly this is intended
// for charge denisities and (possibly) wave functions, but any 3D 
// data can use this class.
// ####################################################################



#include "Grid_data.h"
#include <fstream>
#include <sstream>





// ########################################################
//                       Constructor
// ########################################################
//

Grid_data::Grid_data() { 

  this->N1 = 0;
  this->N2 = 0;
  this->N3 = 0;

  for (int i=0; i<3; i++) {
    v1.push_back(0.0);
    v2.push_back(0.0);
    v3.push_back(0.0);
  }
  

}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

Grid_data::~Grid_data() {  }

// ########################################################
// ########################################################





// ########################################################
//                       SET_DV
// ########################################################
// Sets the volume element within a grid-based property.
// Used for proper scaling.

void Grid_data::set_dV() {
  
  volume = this->v3[0]*(this->v1[1]*this->v2[2] - this->v1[2]*this->v2[1])
    -this->v3[1]*(this->v1[0]*this->v2[2] - this->v1[2]*this->v2[0])
    +this->v3[2]*(this->v1[0]*this->v2[1] - this->v1[1]*this->v2[0]);
  
  volume *= pow(this->a,3);
  dV = volume/double(N1*N2*N3);
  
}

// ########################################################
// ########################################################




// ########################################################
//                       INTEGRATE
// ########################################################
// Integrates the values on the grid. 

double Grid_data::integrate() {
  
  double integral = 0.0;
  
  for (int i=0; i<data.size();i++) {
    integral += data.at(i);
  }
  
  return integral;
  
}

// ########################################################
// ########################################################




// ########################################################
//                       NORMALIZE
// ########################################################
// Precondition the values on the grid.

void Grid_data::normalize(REAL A) {

  double integral = this->integrate();
  
  double N = A/integral;
  
  for (int i=0; i<this->data.size(); i++) {
    this->data.at(i) *= N;
  }

}

// ########################################################
// ########################################################


