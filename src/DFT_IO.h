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
// This is an interface class that handles IO from DFT codes 
// (currently QE, VASP, and FHIAims). Any new DFT interface 
// should use this class.
// ####################################################################






#ifndef __DFT_IO
#define __DFT_IO

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Common.h"
#include "Grid_data.h"
#include "Structure.h"

using namespace std;

class DFT_IO {
 
 public:
  
  DFT_IO() {};
  virtual ~DFT_IO() {};
  
  
  virtual Grid_data get_density(string directory, int step=1) = 0;
  virtual Structure read_structure(string filename) = 0;
  
  virtual REAL read_band_gap(string filename) = 0;
  virtual REAL read_Nelectrons(string filename) = 0;

  virtual REAL get_property(string property,string directory) = 0;
  
  inline vector<REAL> get_user_property(int property, string directory) {
    
    vector<REAL> data;

    this->open(directory+"/user_input");
    this->skip_lines(property-1);
    this->get_line();
    if (line.empty()) { 
      ERROR("User property does not exist"); 
    }
    
    REAL value;
    Line.clear();
    Line >> noskipws >> value;
    if (!Line.eof() && Line.fail()) {
      ERROR("Invalid value for " + directory + ", check user_input file");
    } else {
      data.push_back(value);
    }
    file.close();
    
    return data;
  }
 
  
 protected:
  
  ifstream file;
  
  string line;
  istringstream Line;
  
  REAL Nelect;

  string directory;

  inline void skip_lines(int N) {
    for (int i=0; i<N; i++) {
      getline(this->file,line);
    }
  }
  
  inline void get_line() {
    Line.clear();
    getline(this->file,line);
    string whitespace = " \t\n";
    std::size_t position = line.find_first_not_of(whitespace);
    if (position != string::npos) {
      line = line.substr(position,line.find_last_not_of(whitespace)-position+1);
    }
    Line.str(line);
  }
  
  
  inline void open(string filename) {
    if (file.is_open()) {
      file.close();
    }
    file.open(filename.c_str());
    if (!file.is_open()) {
      ERROR("Could not open file "+filename);
    }
  }
  


  
  
};




#endif
