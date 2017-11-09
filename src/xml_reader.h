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
// A wrapper to the pugixml library (included with this package). This
// is just for convenience.
// ####################################################################



#ifndef __XML_READER
#define __XML_READER

#include <iostream>
#include "pugixml.hpp"
#include <string>
#include <utility>
#include <vector>

using namespace std;

class xml_reader {
  
 public:
  
  xml_reader(string filename="");
  ~xml_reader();

  pugi::xml_node get_node_by_name(string node_name, string attribute_name="", string attribute_value="");
  pugi::xml_node get_node_by_name(string node_name, string attribute_name, string attribute_value, pugi::xml_node in_node);

  vector<pugi::xml_node> get_all_nodes_by_name(string node_name, string attribute_name="", string attribute_value="");
  vector<pugi::xml_node> get_all_nodes_by_name(string node_name, string attribute_name, string attribute_value, pugi::xml_node in_node);
  
  void read(string filename);
  void read_string(string parse);
  inline void print() {
      this->doc.save(std::cout, "", pugi::format_raw);
  }
  
 private:
  
  pugi::xml_document doc;
  pugi::xml_node node;
  pugi::xml_attribute attribute;


  bool is_correct(pugi::xml_node node, string node_name, string attribute_name, string attribute_value);
  
};


#endif
