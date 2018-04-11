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


#include "Tables.h"
#include <iostream>
#include <iomanip>




// ########################################################
//                       Constructor
// ########################################################
//

Tables::Tables()
{

  _X = NULL;
  coefficients = NULL;

  small_x_saturation = NOT_SET;
  large_x_saturation = NOT_SET;

}

// ########################################################
// ########################################################




// ########################################################
//                       Constructor
// ########################################################
//

Tables::Tables(string function)
{
  this->init(function);
}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

Tables::~Tables()
{

  if (coefficients != NULL) {
    for (int i=0; i<my_Nvalues-1; i++) {
      delete [] coefficients[i];
    }
    delete [] coefficients;
  }
  coefficients = NULL;

  if (_X != NULL) {
    delete [] _X;
  }
  _X = NULL;
}

// ########################################################
// ########################################################



// ########################################################
//                       INIT
// ########################################################
// Initialize the table funtion to the function of choice.

void Tables::init(string type)
{

  if (!type.compare("tanh")) {

    init_tanh();

  } else if (type == "cos") {

    init_cos();

  } else if (type == "cutoff_1") {

    init_cutoff_1();

  } else {

    ERROR("Unrecognized table type '"+type+"'");

  }

}

// ########################################################
// ########################################################



// ########################################################
//                       SET_COEFFICIENTS
// ########################################################
// Set the internal coefficients for the function after reading.

void Tables::set_coefficients(REAL X[], REAL a[][4])
{

  this->coefficients = new REAL* [my_Nvalues];
  this->_X = new REAL [my_Nvalues+1];

  for (int i=0; i<my_Nvalues; i++) {
    this->_X[i] = X[i];
    this->coefficients[i] = new REAL [4];
    for (int j=0; j<4; j++) {
      this->coefficients[i][j] = a[i][j];
    }
  }

  this->_X[my_Nvalues] = X[my_Nvalues];

  this->m = (REAL)(my_Nvalues)/(_X[my_Nvalues]-_X[0]);
  this->b = -m*_X[0];

}

// ########################################################
// ########################################################






// ########################################################
//                       INIT_COS
// ########################################################
// Make a cosine function.

void Tables::init_cos()
{

#include "Table_cos.data"

  this->my_Nvalues = 999;
  set_coefficients(X,a);

}

// ########################################################
// ########################################################





// ########################################################
//                       INIT_TANH
// ########################################################
// Make a tanh function.

void Tables::init_tanh()
{

#include "Table_tanh_accurate.data"
  this->my_Nvalues = 1000;
  this->small_x_saturation = a[0][3];
  this->large_x_saturation = a[my_Nvalues-1][0]+a[my_Nvalues-1][1]+a[my_Nvalues-1][2]+a[my_Nvalues-1][3];
  set_coefficients(X,a);

}

// ########################################################
// ########################################################




// ########################################################
//                       INIT_CUTOFF_1
// ########################################################
// Make the function 0.5*(cos(pi*r/Rcut)+1);

void Tables::init_cutoff_1()
{

#include "Table_cutoff_1.data"
  this->my_Nvalues = 3;
  set_coefficients(X,a);
  this->small_x_saturation = 0;
  this->large_x_saturation = 1;
}

// ########################################################
// ########################################################



// ########################################################
//                       SLOW_ROW
// ########################################################
// A slow but robust method to find the row in the table
// that corresponds to a given input. Used as a fallback
// if fast_row fails.

int Tables::slow_row(REAL x)
{

  int min = 0;
  int max = my_Nvalues;
  int delta = max - min;
  int irow = (int)(0.5*(double)(max+min));

  if (x < _X[0]) {
    return min-1;
  } else if (x > _X[max]) {
    return max+1;
  }

  while (max-min > 15) {
    if (_X[irow] > x) {
      max = irow;
    } else if (_X[irow] <= x) {
      min = irow;
    }
    irow = (int)(0.5*(double)(max+min));
  }



  for (irow=min; irow<max; irow++) {
    if (_X[irow] <= x && _X[irow+1] > x) {
      this->my_shift = -_X[irow];
      return irow;
    }
  }

  ERROR("Could not find x value in Table function");

}

// ########################################################
// ########################################################













