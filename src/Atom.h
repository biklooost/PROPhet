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
// This class is mainly to identify atoms by either their atomic number
// or chemical symbol.  It uses a look-up table to do this.
// ####################################################################



#ifndef __ATOM__
#define __ATOM__

#include "Common.h"

using namespace std;

class Atom {
  
 public:
  
  Atom(int number);
  Atom(string symbol);
  Atom();

  ~Atom();
  
  void set_type(int number);
  void set_type(string symbol);
  
  inline int atomic_number() { return my_atomic_number; }
  inline string atomic_symbol() { return my_atomic_symbol; }
    
  
 private:
  
  int my_atomic_number;
  string my_atomic_symbol;

  // An atomic symbol table
  static char _symbols[118][3];
  
};

#endif
