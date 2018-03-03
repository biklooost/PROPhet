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
// Interface to the FHIAims code. It inherits from DFT_IO.  This class
// reads properties from an FHIAims run.
// ####################################################################






#ifndef __FHIAIMS
#define __FHIAIMS

#include "DFT_IO.h"
#include "Structure.h"
#include "xml_reader.h"
#include <iostream>

class FHIAIMS : public DFT_IO
{
public:
  FHIAIMS();
  ~FHIAIMS();
  virtual Grid_data get_density(string prefix, int step=1);
  virtual Structure read_structure(string prefix);
  virtual REAL read_band_gap(string prefix);
  virtual REAL read_Nelectrons(string prefix);
  virtual REAL get_property(string property,string directory);
private:
  xml_reader xml;
};
#endif
