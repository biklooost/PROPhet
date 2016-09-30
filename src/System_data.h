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
//                        CLASS DESCRIPTION
// ####################################################################
// This is a low-level class that handles passing the requested
// properties of a system to a machine-learning algorithm as input. As
// iterating over system properties is the very inner loop of training
// a mapping, it needs to be as efficient as possible.
// ####################################################################




#ifndef __SYSTEM_DATA
#define __SYSTEM_DATA

#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <complex>

#include "Common.h"

using namespace std;

class System_data {
  
 public:
  
  System_data();
  ~System_data();
  
  inline void push_back(vector<REAL> *to_add) { this->data.push_back(to_add); }
  inline void set_inputs(vector<vector<REAL>* > to_add) { 
    data.clear();
    for (int i=0; i<to_add.size(); i++) {
      this->push_back(to_add[i]);
    }
  }

  bool iterate(vector<REAL> &output_data);
  bool iterate(vector<complex<REAL> > &output_data);
  
  inline REAL& at(string key) { 
    for (int i=0; i<keys.size(); i++) {
      if (!keys.at(i).compare("key")) {
	return data[i]->at(0);
      }
    }
    ERROR("Could not find key '"+key+" in system data");
  }
  
  inline REAL target() { return my_target.at(0); }
  inline void target(REAL target_value) { my_target.push_back(target_value); }

  inline vector<REAL> count() {
    vector<REAL> counts(data.size(), 0);
    for (int i=0; i<data.size(); i++) {
      counts.at(i) = data.at(i)->size();
    }
    return counts;
  }
  
  inline vector<REAL> mean() {
    vector<REAL> means(data.size(), 0.0);
    for (int i=0; i<data.size(); i++) {
      for (int j=0; j<(*data.at(i)).size(); j++) {
	means.at(i) += (*data.at(i)).at(j);
      }
    }
    return means;
  }
  
  inline vector<REAL> variance(vector<REAL> means) {
    vector<REAL> variances(data.size(), 0.0);
    for (int i=0; i<data.size(); i++) {
      for (int j=0; j<(*data.at(i)).size(); j++) {
	variances.at(i) += pow( (*data.at(i)).at(j)-means.at(i) ,2);
      }
    }
    return variances;
  }

  inline void Normalize(vector<REAL> means, vector<REAL> variances) {
    
    for (int i=0; i<data.size(); i++) {
      for (int j=0; j<(*data.at(i)).size(); j++) {
	(*data.at(i)).at(j) -= means.at(i);
	if (variances.at(i) >= 0.01) {
	  (*data.at(i)).at(j) /= sqrt(variances.at(i));
	}
      }
    }
  
  }


  inline void lock(bool new_is_locked=true) { this->is_locked = new_is_locked; }

 private:

  bool is_locked;
  
  vector<vector<REAL>* > data;
  vector<int> indices;
  vector<REAL> my_target;
  int my_N;
  vector<int> N_max;
  vector<REAL> empty_vector;

  vector<string> keys;
  
};

#endif

