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
//                          CLASS DESCRIPTION
// ####################################################################
// This class handles user input from the main input file and sets up
// the run, including locating files, parallelization, etc.
// ####################################################################




#ifndef __SETUP
#define __SETUP

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include "Parallel.h"
#include "Error.h"
#include "Functional_params.h"

using namespace std;

class Setup {
  
public:
 
  Setup(string input_file);
  ~Setup();
  
  
  
  inline int Nsystems() { return this->systems.size(); }
  inline map<string,string> system(int i) { 
    if (this->systems.at(i).empty()) {
      ERROR("No file specifiec for input field"); 
    } 
    return this->systems.at(i); 
  }

  inline bool is_potential() { return my_is_potential; }
  

  Functional_params F;  
  
 private:
  
  vector<map<string,string> > systems;
  
  map<string,string> property_map;
  
  Parallel *mpi;

  vector<ifstream*> input_files;
  
  void read_input(string filename);
  
  void split_systems();
  
  void print_details();

  string clean_line(string line);
  
  enum {none, functional, input};
  
  bool my_is_potential;
  
  
  inline string plural(string in, int N) {
    if (N == 1) { return in; }
    else if (N != 1) { return in+"s"; }
  }
  
  inline string to_lower(const string in) {
    string output = in;
    transform(output.begin(), output.end(), output.begin(), ::tolower);
    return output;
  }


 
  
};













#endif
