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






#ifndef __GRID_DATA
#define __GRID_DATA

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "Common.h"

using namespace std;



class Grid_data {
  
  friend class VASP;
  friend class QE;
  friend class FHIAIMS;

 public:
  
  Grid_data();
  ~Grid_data();
  
  double integrate();
  
  void variance(vector <int> bounds);
  
  
  
  inline REAL& operator()(int x, int y, int z) {
    
    if (x<0) { x += N1;}
    if (y<0) { y += N2;}
    if (z<0) { z += N3;}
    if (x>=N1) { x -= N1;}
    if (y>=N2) { y -= N2;}
    if (z>=N3) { z -= N3;}
    
    return data.at(N1*N2*z + N1*y + x);
    
  }
  
  inline int& Nx() { return N1; }
  inline int& Ny() { return N2; }
  inline int& Nz() { return N3; }
  
  inline int N() { return N1*N2*N3; }
  
  REAL& operator[](int i) { return data.at(i); }
  
  inline vector<REAL>* as_vector_ptr() { return &this->data; }
  
  inline void push_back(REAL value) { data.push_back(value); }

  inline void reserve(int N) {
    data.reserve(N);
  }
  
  inline void resize(int x, int y, int z) { 
    data.resize(x*y*z); 
    N1 = x; N2 = y; N3 = z;
  }

  inline REAL get_dV() { return dV;}

  void normalize(REAL A = 1.0);

  void set_dV();
  inline void set_dV(REAL new_dV) { dV = new_dV; }
  double volume;

 private:
  
  double a;
  vector<double> v1, v2, v3;
  

  
  REAL dV;
  
  
  int N1, N2, N3;
  
  vector<REAL> data;
  
};

#endif
