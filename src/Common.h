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
// This class is included in all other *.h files. It defines or
// includes a few common routines such as warnings and errors and
// some random number functions. It also defines the precision we're
// using, as set by the USE_FLOAT or USE_LONG_DOUBLE defines.
// ####################################################################







#ifndef __COMMON
#define __COMMON

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <fenv.h>

#include "Error.h"

#ifdef USE_FLOAT
#define REAL float
#else
#ifdef USE_LONG_DOUBLE
#define REAL  long double
#else
#define REAL double
#endif
#endif

#define NOT_SET -23482.249283420

namespace RAND
{

inline REAL Uniform(REAL lower=0, REAL upper=1)
{
  return (upper-lower)*((REAL)rand()/RAND_MAX) + lower;
}
inline REAL Normal(REAL mean=0, REAL std_dev=1)
{
  return mean+std_dev*std::sqrt(-2.0*log(RAND::Uniform(0,1)))*cos(2.0*3.141592653589793*RAND::Uniform(0,1));
}

};


#endif

