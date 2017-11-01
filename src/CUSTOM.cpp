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
// An interface to the Quantum Espresso package. This package facilitates
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

CUSTOM::CUSTOM() { }

// ########################################################
// ########################################################



// ########################################################
//                       Destructor
// ########################################################
//

CUSTOM::~CUSTOM() { }

// ########################################################
// ########################################################




// ########################################################
//                       READ_ENERGIES
// ########################################################
// Read energies from 

void CUSTOM::read_energies(string prefix,int nbnd, REAL  *cbm, REAL *vbm){
	vector<REAL> energy;
	vector<int> occu;
	int i,indexVBM;
	this->xml.read(prefix + "/data-file.xml");
	vector<pugi::xml_node> nodes1 = xml.get_all_nodes_by_name("EIGENVALUES");
	string temp = nodes1.front().child_value();
	stringstream ss(temp.c_str());
	while(getline(ss,temp,'\n')) {
		if(std::string::npos != temp.find_first_of("0123456789")){
		  energy.push_back(atof(temp.c_str()));
		}
	}
	nodes1 = xml.get_all_nodes_by_name("OCCUPATIONS");
	temp = nodes1.front().child_value();
	stringstream ss1(temp);
	while(getline(ss1,temp,'\n')) {
		if(std::string::npos != temp.find_first_of("0123456789")) {
		  occu.push_back(atoi(temp.c_str()));
		}
	}
	for (i = 0; i < occu.size(); i++) {
		if (occu[i] == 1) {
			indexVBM = i;
		}
	}
	file.close();
	*vbm = energy[indexVBM];
	*cbm = energy[indexVBM + 1];
}

// ########################################################
// ########################################################



// ########################################################
//                       READ_STRUCTURE
// ########################################################
// Reads the structure from a Quantum Espresso run.

Structure CUSTOM::read_structure(string prefix) { 
	Structure temps;

	temps.CART = true;
	temps.periodic = true;
	
	stringstream tvector;
	vector< vector<REAL> > coordinates;
	string datafile = prefix + "/data-file.xml";
	this->xml.read(datafile);
	vector<pugi::xml_node> node = xml.get_all_nodes_by_name("NUMBER_OF_ATOMS");
	temps.Natom = atoi(node.back().child_value());
	stringstream tcord;
	stringstream atom;
	node = this->xml.get_all_nodes_by_name("a1");
	tvector << node.back().child_value();
	REAL tnumber;
	while (tvector >> tnumber) { temps.a.push_back(tnumber*0.529177249); }
	tvector.str("");
	tvector.clear();	
	node = this->xml.get_all_nodes_by_name("a2");
	tvector << node.back().child_value();
	while (tvector >> tnumber) { temps.b.push_back(tnumber*0.529177249); }
	tvector.str("");
	tvector.clear();	
	node = this->xml.get_all_nodes_by_name("a3");
	tvector << node.back().child_value();
	while (tvector >> tnumber) { temps.c.push_back(tnumber*0.529177249); }
	tvector.str("");
	tvector.clear();	
	for (int i = 1; i<=temps.Natom; i++) {
		vector<REAL> row;
		REAL tmp;
		atom << "ATOM." << i;
		node = xml.get_all_nodes_by_name(atom.str());
		temps.types.push_back(Atom(node.back().attribute("SPECIES").value()));
		tcord<< node.back().attribute("tau").value();		
		for (int i = 0; i < 3; i++ ){
			tcord >> tmp;
			row.push_back(tmp*0.529177249);
		}
		temps.pos.push_back(row);
		coordinates.push_back(row);
		atom.str("");
		atom.clear();
		tcord.str("");
		tcord.clear();
	}	
	return temps;
}
REAL CUSTOM::get_property(string property,string directory) { 
  string filename = directory+"/data-file.xml";
  this->xml.read(filename);

  if (property == "energy") {
      try {
        vector<pugi::xml_node> nodes = xml.get_all_nodes_by_name("TOTAL_ENERGY");
        if (nodes.size() == 0) {
          throw 3;
        }
        REAL energy = atof(nodes.back().child_value());
        return energy;
      } catch (int e) {
        ERROR("Unable to parse " +directory + "/data-file.xml, it does not contain \"TOTAL ENERGY\" key");
      }
  } else if ( (property == "dft_gap") || (property == "gw_gap") ) {
    REAL bndgap = read_band_gap(directory);
    return bndgap;
  }
  

}
REAL CUSTOM::read_band_gap(string prefix) { 
		string datafile = prefix + "/data-file.xml";
		int nKpoints, i, nbnd;
		REAL CBM=9.0e9;
		REAL VBM=-9.0e-9;
		REAL tCBM, tVBM, Bandgap;
		stringstream save_file;
		this->xml.read(datafile);
		vector<pugi::xml_node> nodes = xml.get_all_nodes_by_name("NUMBER_OF_K-POINTS");	
		nKpoints = atoi(nodes.back().child_value());
		if(nodes.empty()) { cout << "Unable to read K-Points from " << datafile << "\n"; }
		nodes = xml.get_all_nodes_by_name("NUMBER_OF_BANDS");	
		if(nodes.empty()) { cout << "Unable to read Number of bands from " << datafile << "\n"; }
		nbnd = atoi(nodes.back().child_value());	
		for (i = 1; i <= nKpoints; i++) {
			if (i < 10) {
				save_file << prefix << "/K0000" << nKpoints << "/eigenval.xml";		
			}else if (i < 100) {
				save_file << prefix << "/K000" << nKpoints << "/eigenval.xml";		
			}else if (i < 1000) {
				save_file << prefix << "/K00" << nKpoints << "/eigenval.xml";		
			}else if (i < 10000) {
				save_file << prefix << "/K0" << nKpoints << "/eigenval.xml";		
			}else if (i < 100000) {
				save_file << prefix << "/K" << nKpoints << "/eigenval.xml";		
			}
			read_energies(save_file.str(),nbnd, &tCBM, &tVBM);
			if (tCBM < CBM) {CBM = tCBM;}
			if (tVBM > VBM) {VBM = tVBM;}
			save_file.str("");
		}	
		return (CBM - VBM)*27.211396132;
}

// ########################################################
// ########################################################




// ########################################################
//                       GET_DENSTIY
// ########################################################
// Reads the charge density from a Quantum Espresso run.

Grid_data CUSTOM::get_density(string prefix, int step) { 
    Grid_data x;
    return x;
}

// ########################################################
// ########################################################



// ########################################################
//                       READ_NELECTRONS
// ########################################################
// Reads the number of electrons from a Quantum Espresso run.

REAL CUSTOM::read_Nelectrons(string prefix) {
	string datafile = prefix + "/data-file.xml";
	this->xml.read(datafile);
	vector<pugi::xml_node> nodes = xml.get_all_nodes_by_name("NUMBER_OF_ELECTRONS");	
	this->Nelect = atof(nodes.back().child_value());
	return this->Nelect;
}

// ########################################################
// ########################################################


