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
// An interface to the VASP code. This class follows the DFT_IO 
// interface.
// ####################################################################




#ifndef __VASP
#define __VASP


#include "DFT_IO.h"
#include "Common.h"
#include "xml_reader.h"
#include "Structure.h"


class VASP : public DFT_IO {

 public:

  VASP();
  ~VASP();
  
  virtual Grid_data get_density(string directory, int step = 1);
  virtual Structure read_structure(string filename);
  virtual REAL read_band_gap(string filename);
  virtual REAL read_Nelectrons(string filename);

  virtual REAL get_property(string property, string directory);

 private:
  
  xml_reader xml;
  

};

#endif

