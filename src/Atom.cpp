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




#include "Atom.h"

#include <iostream>
#include <algorithm>
#include <ctype.h>



// ########################################################
//                       Constructor
// ########################################################
Atom::Atom()
{

  my_atomic_number = 0;
  my_atomic_symbol.clear();

}

// ########################################################
// ########################################################




// ########################################################
//                       Constructor
// ########################################################

Atom::Atom(int number)
{
  this->set_type(number);
}

// ########################################################
// ########################################################





// ########################################################
//                       Constructor
// ########################################################

Atom::Atom(string symbol)
{
  this->set_type(symbol);
}

// ########################################################
// ########################################################






// ########################################################
//                     SET_TYPE
// ########################################################
// Set the type of atom given the atomic number as input

void Atom::set_type(int number)
{

  my_atomic_number = number;
  my_atomic_symbol = _symbols[number-1];

  if (my_atomic_symbol.empty()) {
    ERROR("Could not find atomic symbol for atomic number "+number);
  }

}

// ########################################################
// ########################################################





// ########################################################
//                     SET_TYPE
// ########################################################
// Set the type of atom given the atomic symbol as input

void Atom::set_type(string symbol)
{

  my_atomic_symbol = symbol;

  string whitespace = " \t\n";
  symbol = symbol.substr(symbol.find_first_not_of(whitespace), symbol.find_last_not_of(whitespace) - symbol.find_first_not_of(whitespace)+1);
  transform(symbol.begin()+1, symbol.end(), symbol.begin()+1, ::tolower);
  transform(symbol.begin(), symbol.begin()+1, symbol.begin(), ::toupper);

  for (int i=0; i<118; i++) {
    if (_symbols[i] == symbol) {
      my_atomic_number = i+1;
      break;
    }
  }

  if (!my_atomic_number) {
    ERROR("Could not find atomic number corresponding to the element "+symbol);
  }

}

// ########################################################
// ########################################################





// ########################################################
//                       Destructor
// ########################################################

Atom::~Atom()
{

}

// ########################################################
// ########################################################





// Atomic chart to convert between atomic symbol and atomic number
char Atom::_symbols[118][3] = {
  {'H'}, {'H','e'}, {'L','i'}, {'B','e'}, {'B'}, {'C'}, {'N'}, {'O'}, {'F'}, {'N','e'},
  {'N','a'}, {'M','g'}, {'A','l'}, {'S','i'}, {'P'}, {'S'}, {'C','l'}, {'A','r'},
  {'K'}, {'C','a'}, {'S','c'}, {'T','i'}, {'V'}, {'C','r'}, {'M','n'}, {'F','e'},
  {'C','o'}, {'N','i'}, {'C','u'}, {'Z','n'}, {'G','a'}, {'G','e'}, {'A','s'}, {'S','e'},
  {'B','r'}, {'K','r'}, {'R','b'}, {'S','r'}, {'Y'}, {'Z','r'}, {'N','b'}, {'M','o'},
  {'T','c'}, {'R','u'}, {'R','h'}, {'P','d'}, {'A','g'}, {'C','d'}, {'I','n'}, {'S','n'},
  {'S','b'}, {'T','e'}, {'I'}, {'X','e'}, {'C','s'}, {'B','a'}, {'L','a'}, {'C','e'},
  {'P','r'}, {'N','d'}, {'P','m'}, {'S','m'}, {'E','u'}, {'G','d'}, {'T','b'}, {'D','y'},
  {'H','o'}, {'E','r'}, {'T','m'}, {'Y','b'}, {'L','u'}, {'H','f'}, {'T','a'}, {'W'},
  {'R','e'}, {'O','s'}, {'I','r'}, {'P','t'}, {'A','u'}, {'H','g'}, {'T','l'}, {'P','b'},
  {'B','i'}, {'P','o'}, {'A','t'}, {'R','n'}, {'F','r'}, {'R','a'}, {'A','c'}, {'T','h'},
  {'P','a'}, {'U'}, {'N','p'}, {'P','u'}, {'A','m'}, {'C','m'}, {'B','k'}, {'C','f'},
  {'E','s'}, {'F','m'}, {'M','d'}, {'N','o'}, {'L','r'}, {'R','f'}, {'D','b'}, {'S','g'},
  {'B','h'}, {'H','s'}, {'M','t'}, {'D','s'}, {'R','g'}, {'C','n'}, {'U','u','t'}, {'F','l'},
  {'U','u','p'}, {'L','v'}, {'U','u','s'}, {'U','u','o'}
};
