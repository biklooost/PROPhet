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
// A general class to hold information about system properties.
// ####################################################################


// This class reads and stores information about system properties for any training or testing example system. It uses DFT_IO to get properties from DFT runs, though it can also read user_input files to obtain user-defined properties.

#ifndef __SYSTEM
#define __SYSTEM

#include <iostream>
#include <string>
#include <map>

#include "Common.h"
#include "Grid_data.h"
#include "DFT_IO.h"
#include "Functional_params.h"
#include "System_data.h"
#include "Structure.h"

using namespace std;

class System {
  
public:
  
  System(map<string,string> files, Functional_params* F);
  System();
  ~System();
  
  System_data properties;

  REAL Prefactor;
  
  Structure structure;

 private:
  
  Grid_data Density;
  map<string,vector<REAL> > data;

  vector<REAL> GW_gap;
  vector<REAL> DFT_gap;
  vector<REAL> volume;

  string directory;
  
  
};

#endif
