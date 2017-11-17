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



#include <cstdlib>

#include "xml_reader.h"
#include "Error.h"
#include <sstream>



// ########################################################
//                       Constructor
// ########################################################
//

xml_reader::xml_reader(string filename)
{


  if (!filename.empty()) {
    this->read(filename);
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

xml_reader::~xml_reader()
{

}

// ########################################################
// ########################################################





// ########################################################
//                       READ
// ########################################################
// Reads an XML file for parsing.

void xml_reader::read(string filename)
{

  pugi::xml_parse_result result = this->doc.load_file(filename.c_str());

  if (!result) {
    ERROR("Could not load \""+filename+"\"");
  }

}

void xml_reader::read_string(string parse)
{

  pugi::xml_parse_result result = this->doc.load_string(parse.c_str());


  if (!result) {
    std::stringstream msg;
    msg << "XML [" << parse << "] parsed with errors, attr value: [" << doc.child("node").attribute("attr").value() << "]\n";
    msg << "Error description: " << result.description() << "\n";
    ERROR(msg.str());
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       GET_NODE_BY_NAME
// ########################################################
// Finds a node with a given name in the XML tree.

pugi::xml_node xml_reader::get_node_by_name(string node_name, string attribute_name, string attribute_value)
{

  return get_node_by_name(node_name, attribute_name, attribute_value,this->doc);

}

// ########################################################
// ########################################################




// ########################################################
//                       GET_NODE_BY_NAME
// ########################################################
// Finds a node with a given name in the XML tree.

pugi::xml_node xml_reader::get_node_by_name(string node_name, string attribute_name, string attribute_value, pugi::xml_node in_node)
{


  pugi::xml_node child = in_node.first_child();

  while(child) {
    if (is_correct(child, node_name, attribute_name, attribute_value)) {
      return child;
    }
    pugi::xml_node grandchild = get_node_by_name(node_name, attribute_name, attribute_value, child);
    if (grandchild) {
      return grandchild;
    }
    child = child.next_sibling();
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       GET_ALL_NODES_BY_NAME
// ########################################################
// Finds all nodes with a given name in the XML tree.

vector<pugi::xml_node> xml_reader::get_all_nodes_by_name(string node_name, string attribute_name, string attribute_value)
{

  return get_all_nodes_by_name(node_name, attribute_name, attribute_value, doc);

}

// ########################################################
// ########################################################





// ########################################################
//                       GET_ALL_NODES_BY_NAME
// ########################################################
// Finds all nodes with a given name in the XML tree.

vector<pugi::xml_node> xml_reader::get_all_nodes_by_name(string node_name, string attribute_name, string attribute_value, pugi::xml_node in_node)
{

  vector<pugi::xml_node> nodes;

  pugi::xml_node child = in_node.first_child();
  while (child) {
    if (is_correct(child,node_name,attribute_name,attribute_value)) {
      nodes.push_back(child);
    }
    vector<pugi::xml_node> grandchildren = get_all_nodes_by_name(node_name,attribute_name,attribute_value,child);
    for (int i=0; i<grandchildren.size(); i++) {
      nodes.push_back(grandchildren[i]);
    }
    child = child.next_sibling();
  }

  return nodes;

}

// ########################################################
// ########################################################




// ########################################################
//                       IS_CORRECT
// ########################################################
// Checks that the node is indeed what was asked for.

bool xml_reader::is_correct(pugi::xml_node node, string node_name, string attribute_name, string attribute_value)
{

  if (!node_name.empty() && node.name() != node_name) {
    return false;
  }

  if (!attribute_name.empty()) {
    for (pugi::xml_attribute node_attribute = node.first_attribute(); node_attribute; node_attribute = node_attribute.next_attribute()) {
      if (node_attribute.name() == attribute_name) {
        if (!attribute_value.empty()) {
          if (node_attribute.value() == attribute_value) {
            return true;
          } else {
            return false;
          }
        } else {
          return true;
        }

      }

    }

    return false;
  } else {

    return true;

  }

}

// ########################################################
// ########################################################


