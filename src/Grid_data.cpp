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
#include <algorithm>
#include <math.h>
#include <valarray>




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
//                       VARIANCE
// ########################################################
// This is a quick and dirty function to remove the tails of the 
// charge density vector. May have memory leaks

void Grid_data::variance(vector <int> bounds) {
    if ((bounds[0] == 0 && bounds[1] == 0 ) || bounds[0] == bounds[1]) { return; }
    vector <REAL> old_data = this->data;
    sort(old_data.begin(),old_data.end());
    REAL lower = old_data[floor(old_data.size()*bounds[0]/100.)];
    REAL upper = old_data[ceil(old_data.size()*bounds[1]/100.)];
    vector <REAL> new_data;
    for(int i = 0; i<this->data.size(); i++) {
        if ((this->data.at(i) < upper) && (this->data.at(i) > lower)) {
            new_data.push_back(this->data.at(i));
        }
    }
    this->data = new_data; 
}

// ########################################################
// ########################################################

void Grid_data::conv_matrix(int n) {
    vector <REAL> tmp = this->data;
    REAL *conv; 
    int size;
    REAL sigma;
    switch (n) {
        case 1: //Edge detection
            size = 3;
            conv = (REAL*)malloc(size*size*size*sizeof(REAL));
            for (int i = 0; i < size*size*size; i++) { conv[i] = -1.0; }
            conv[size*size*1 + size*1 + 1] = 8.0;   
        case 2: //Normalization/average
            conv = (REAL*)malloc(size*size*size*sizeof(REAL));
            for (int i = 0; i < size*size*size; i++) { conv[i] = 1.0/(size*size*size); }
        case 3: //5x5 gaussian blur
            size = 5;
            sigma = 1.132;
            conv = (REAL*)malloc(size*size*size*sizeof(REAL));
            for (int i = 0; i < size; i++) { 
                for (int j = 0; j <size; j++) {
                    for (int k = 0; k < size; k++) {
                        conv[size*size*k + size*j * i] = 1/(2*3.14*sigma)*exp(-(pow((k-2),2) + pow((j-2),2) + pow((i-2),2))/(2*sigma));
                    }
                }
            }
        case 4: //unsharpen mask
            size = 5;
            sigma = 1.132;
            conv = (REAL*)malloc(size*size*size*sizeof(REAL));
            for (int i = 0; i < size; i++) { 
                for (int j = 0; j <size; j++) {
                    for (int k = 0; k < size; k++) {
                        conv[size*size*k + size*j * i] = -1/(2*3.14*sigma)*exp(-(pow((k-2),2) + pow((j-2),2) + pow((i-2),2))/(2*sigma));
                    }
                }
            }
            conv[size*size*2 + size*2 + 2 ] *= -13.5; 
        default:
            size = 3;
            conv = (REAL*)malloc(size*size*size*sizeof(REAL));
            for (int i = 0; i < size*size*size; i++) { conv[i] = 0.0; }
            conv[size*size*1 + size*1 + 1] = 1.0;
    }
    int shift = (size-1)/2;
    for (int i = 0; i<this->N1; i++) {
        for (int j = 0; j <this->N2; j++) {
            for (int k = 0; k<this->N3; k++) {
                REAL sum = 0.0;
                for (int ii =0; ii < size; ii ++){
                    for (int jj = 0; jj < size; jj++) {
                        for (int kk = 0; kk < size; kk++) {
                            sum += conv[size*size*kk + size*jj + ii]*(*this)(i+ii-shift,j+jj-shift,k+kk-shift); //shift is used to shift convolutional matrix  to center on (i,j,k) in tmp
                        }
                    }
                }
                tmp[N1*N2*k + N1*j + i] = sum;
            }
        }
    }
    this->data = tmp;
}
