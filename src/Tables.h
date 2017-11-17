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
// A class to define functions as look-up tables for efficiency. It
// was originally created to facilitate a tanh_spline lookup-table to
// efficiently implement a tanh-like transfer function. Note that
// this does not need to be an excellent approximation
// to tanh, as all that is really required for a neural network
// transfer function is an odd function that saturates at large
// inputs.  Other table functions to 1) efficiently implement
// expensive-to-compute functions or 2) implement non-analytical or
// extremely complex functions can easily be added. Functions can be
// entered into a table of coefficients to a cubic spline
// interpolation. This allows rapid evaluate of any desired
// function and its derivative.
// ####################################################################






#ifndef __TABLES
#define __TABLES

#include <vector>
#include <string>
#include <cmath>

#include "Common.h"



class Tables
{

public:

  Tables();
  Tables(string function);
  ~Tables();

  void init(string type);

  int slow_row(REAL value);

  int fast_row(REAL value)
  {

    if (value < _X[0]) {
      return -1;
    } else if (value >= _X[my_Nvalues]) {
      return my_Nvalues + 1;
    }

    int bin = (int)(m*value+b);
    if (value < _X[bin] || value > _X[bin+1]) {
      // Just in case something goes wrong
      int B = slow_row(value);
      cout << bin << "  "<<B<<endl;
      bin = slow_row(value);
      my_shift = -_X[bin];
      return bin;
    }

    my_shift = -_X[bin];
    return bin;
  }



  inline REAL operator()(int row, int power)
  {
    return coefficients[row][power];
  }

  inline REAL operator()(REAL in_x)
  {
    int irow = fast_row(in_x);
    if (irow < 0) {
      return small_x();
    }
    if (irow > my_Nvalues) {
      return large_x();
    }
    x = in_x + my_shift;
    return ((coefficients[irow][0]*x + coefficients[irow][1])*x + coefficients[irow][2])*x + coefficients[irow][3];
  }

  inline REAL derivative(REAL in_x)
  {
    int irow = fast_row(in_x);
    if (irow < 0 ) {
      return 0;
    }
    if (irow > my_Nvalues) {
      return 0;
    }
    x = in_x + my_shift;
    return (3*coefficients[irow][0]*x + 2*coefficients[irow][1])*x + coefficients[irow][2];
  }


  inline REAL shift()
  {
    return my_shift;
  }

  inline int Nvalues()
  {
    return my_Nvalues;
  }

  inline REAL small_x()
  {
    if (small_x_saturation != NOT_SET) {
      return small_x_saturation;
    } else {
      ERROR("value out of bounds for table");
    }
  }

  inline REAL large_x()
  {
    if (large_x_saturation != NOT_SET) {
      return large_x_saturation;
    } else {
      ERROR("value out of bounds for table");
    }
  }


private:

  REAL my_shift;
  REAL **coefficients;
  REAL* _X;

  // Terms for linear interpolation of the bin
  REAL m, b;

  REAL x,x2,x3;

  REAL small_x_saturation;
  REAL large_x_saturation;

  REAL min_value;
  REAL max_value;
  int my_Nvalues;

  void init_tanh();
  void init_cos();
  void init_cutoff_1();

  void set_coefficients(REAL X[], REAL in[][4]);


};

#endif


