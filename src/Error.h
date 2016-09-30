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
// Class to define a general warning and error interface. The functions
// here automatically provide location information as to where they 
// were called.
// ####################################################################







#ifndef __ERROR
#define __ERROR

using namespace std;

#define INFO __FILE__,__FUNCTION__,__LINE__
#define ERROR(a) fatal_error(a,INFO)
#define WARNING(a) warning(a,INFO)

#include <iostream>
#include <string>
#include <cstdlib>

void fatal_error(string message="", string in_class="", string in_function="",int on_line=0);
void warning(string message="",string in_class="",string in_function="",int on_line=0);

  
#endif
