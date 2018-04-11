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
// An interface to the PROPhet XML file. This package facilitates
// reading directly from CUSTOM files.  The package follows the DFT_IO
// interface.
// ####################################################################



#include "CUSTOM.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>

#include "Atom.h"

// ########################################################
//                       Constructor
// ########################################################
//

CUSTOM::CUSTOM()
{
  this->xml_process = false;

}

// ########################################################
// ########################################################



// ########################################################
//                       Destructor
// ########################################################
//

CUSTOM::~CUSTOM()
{

}

// ########################################################
// ########################################################

// ########################################################
//                       READ_STRUCTURE
// ########################################################
// Reads the structure from a PROPhet XML run.

Structure CUSTOM::read_structure(string prefix)
{
  if (this->xml_process == false) {
    Structure temps;
    temps.CART = true;
    temps.periodic = true;
    ifstream fname("PROPhet.xml",std::ifstream::in);
    string ltemp;
    stringstream cstring;
    int id =0;
    stringstream temp(prefix);
    temp >> id;
    temp.str("");
    temp.clear();
    cstring << "<system id=\"" << id + 1<< "\">";
    while(std::getline(fname, ltemp)) {
      if (ltemp.find(cstring.str()) != std::string::npos) {
        temp << ltemp;
        do {
          std::getline(fname,ltemp);
          temp << ltemp;
        } while(ltemp.find("</system>") == std::string::npos);
        break;
      }
    }
    fname.close();
    this->xml.read_string(temp.str());
    stringstream tvector,tcord;
    REAL tnumber;
    pugi::xml_node lattice = xml.get_node_by_name("lattice");
    vector <pugi::xml_node> train = xml.get_all_nodes_by_name("train");
    if (train.size() > 0) {
      temps.train = string(train.back().child_value());
    } else {
      temps.train = "train";
    }
    //cout << temps.train << endl;

    tvector << string(lattice.child("a").child_value());
    while (tvector >> tnumber) {
      temps.a.push_back(tnumber);
    }
    tvector.str("");
    tvector.clear();
    tvector << string(lattice.child("b").child_value());
    while (tvector >> tnumber) {
      temps.b.push_back(tnumber);
    }
    tvector.str("");
    tvector.clear();
    tvector << string(lattice.child("c").child_value());
    while (tvector >> tnumber) {
      temps.c.push_back(tnumber);
    }
    tvector.str("");
    tvector.clear();
    vector <pugi::xml_node> node = xml.get_all_nodes_by_name("atoms");
    pugi::xml_node atoms = node.back();
    REAL tmp;
    temps.Natom = atoi(atoms.child("natoms").child_value());
    for (pugi::xml_node t_node = atoms.child("atom"); t_node; t_node = t_node.next_sibling("atom")) {
      string specie = string(t_node.attribute("specie").value());
      temps.types.push_back(Atom(specie));
      tcord<< t_node.child_value();
      vector<REAL> row;
      for (int i = 0; i < 3; i++ ) {
        tcord >> tmp;
        row.push_back(tmp);
      }
      temps.pos.push_back(row);
      tcord.str("");
      tcord.clear();
    }
    this->xml_process = true;
    this->xstruct = temps;
    return temps;
  } else {
    return this->xstruct;
  }
}
REAL CUSTOM::get_property(string property,string directory)
{

  vector<pugi::xml_node> nodes = xml.get_all_nodes_by_name("target");
  if (nodes.size() == 0) {
    ERROR("For system " + directory + " no target key");
  }
  REAL target = atof(nodes.back().child_value());
  return target;

}
REAL CUSTOM::read_band_gap(string prefix)
{
  return 0.0;
}

// ########################################################
// ########################################################




// ########################################################
//                       GET_DENSTIY
// ########################################################
// Reads the charge density from a cube file.

Grid_data CUSTOM::get_density(string prefix, int step)
{
  if (this->xml_process == false) {
    Structure x = this->read_structure(prefix);
    this->xstruct = x;
    this->xml_process = true;
  }
  vector<pugi::xml_node> nodes = xml.get_all_nodes_by_name("charge-density");
  if (nodes.size() == 0) {
    ERROR("For system " + prefix + " charge density tag");
  }
  string chd = nodes.back().child_value();
  Grid_data grid, downsample;
  int natom;
  string token;
  grid.a = 1;
  this->open(chd);
  this->skip_lines(2);
  this->get_line();
  Line >> natom;
  this->get_line();
  Line >> grid.Nx() >> grid.v1[0] >> grid.v1[1] >> grid.v1[2];
  this->get_line();
  Line >> grid.Ny() >> grid.v2[0] >> grid.v2[1] >> grid.v2[2];
  this->get_line();
  Line >> grid.Nz() >> grid.v3[0] >> grid.v3[1] >> grid.v3[2];
  for(int i = 0; i<3; i++) {
    grid.v1[i] *= grid.Nx();
    grid.v2[i] *= grid.Ny();
    grid.v3[i] *= grid.Nz();
  }
  this->skip_lines(natom);
  do {
    this->get_line();
    while(getline(this->Line,token,' ')) {
      if(token == "") {
        continue;
      } else {
        grid.push_back(atof(token.c_str()));
      }
    }
  } while(!this->file.eof());
  grid.set_dV();
  grid.train = this->xstruct.train;
  return grid;
}

// ########################################################
// ########################################################



// ########################################################
//                       READ_NELECTRONS
// ########################################################
// Reads the number of electrons from a PROPhet XML run.

REAL CUSTOM::read_Nelectrons(string prefix)
{
  string datafile = prefix + "/data-file.xml";
  this->xml.read(datafile);
  vector<pugi::xml_node> nodes = xml.get_all_nodes_by_name("NUMBER_OF_ELECTRONS");
  this->Nelect = atof(nodes.back().child_value());
  return this->Nelect;
}

// ########################################################
// ########################################################


