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



#include <cstdlib>

#include "VASP.h"
#include "Atom.h"



// ########################################################
//                       Constructor
// ########################################################
//

VASP::VASP() { 

}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

VASP::~VASP() { }

// ########################################################
// ########################################################




// ########################################################
//                       READ_STRUCTURE
// ########################################################
// Reads the atomic configuration from a VASP run.

Structure VASP::read_structure(string filename) { 

  Structure new_struct;
  new_struct.CART = false;
  new_struct.periodic = true;
  
  string line;
  istringstream Line;
  double temp;

  filename = filename + "/vasprun.xml";
  this->xml.read(filename);
  
  pugi::xml_node basis = xml.get_all_nodes_by_name("varray", "name", "basis").back();
  vector<pugi::xml_node> vectors = xml.get_all_nodes_by_name("v","","",basis);


  Line.str(vectors.at(0).text().get());
  while (Line >> temp) {
    new_struct.a.push_back(temp);
  }						
  Line.clear();
  
  Line.str(vectors.at(1).text().get());
  while (Line >> temp) {
    new_struct.b.push_back(temp);
  }
  Line.clear();

  Line.str(vectors.at(2).text().get());
  while (Line >> temp) {
    new_struct.c.push_back(temp);
  }
  Line.clear();
  
  new_struct.Natom = atoi(xml.get_node_by_name("atoms").child_value());
  
  
  pugi::xml_node positions = xml.get_all_nodes_by_name("varray","name","positions").back();
  vectors = xml.get_all_nodes_by_name("v","","",positions);
  
  if (vectors.size() != new_struct.Natom) { ERROR("Atoms were not read correctly"); }
  
  for (int i=0; i<new_struct.Natom; i++) {
    vector<REAL> atom;
    Line.str(vectors.at(i).text().get());
    while(Line >> temp) { atom.push_back(temp); }
    new_struct.pos.push_back(atom);
    Line.clear();
  }
  
  pugi::xml_node species = xml.get_node_by_name("atominfo");
  species = xml.get_node_by_name("array","name","atoms",species);
  vectors = xml.get_all_nodes_by_name("c","","",species);

  
  for (int i=1; i<vectors.size(); i+=2) {
    for (int j=0; j<atoi(vectors.at(i).child_value()); j++) {
      new_struct.types.push_back(Atom(vectors.at(i-1).child_value()));
    }
  }


  return new_struct;
}

// ########################################################
// ########################################################





// ########################################################
//                       GET_PROPERTY
// ########################################################
// Reads various properties from VASP run.

REAL VASP::get_property(string property, string directory) {
  
  string filename = directory+"/vasprun.xml";
  this->xml.read(filename);

  if (property == "energy") {
    
    vector<pugi::xml_node> nodes = xml.get_all_nodes_by_name("i","name","e_wo_entrp");
    if (nodes.empty()) {ERROR("Could not get energy from "+directory+"/vasprun.xml");}
    return atof(nodes.back().child_value());
    
  } else if ( (property == "dft_gap") || (property == "gw_gap") ) {
        
    int VB = (int)(0.5*atof(xml.get_node_by_name("i","name","NELECT").text().get()))-1;
    pugi::xml_node eigenvalues = xml.get_node_by_name("eigenvalues");
    if (!eigenvalues) { ERROR("Could not read eigenvalues in "+filename);}
    vector<pugi::xml_node> kpoints = xml.get_all_nodes_by_name("set","","",xml.get_node_by_name("set","comment","spin 1",eigenvalues));
    if (!kpoints.front()) {ERROR("Could not read kpoints from "+filename);}
    
    int Nkpoints = kpoints.size();
    
    double VBM = -9e9;
    double CBM = 9e9;
    
    for (int i=0; i<Nkpoints; i++) {
      vector<pugi::xml_node> Ek = xml.get_all_nodes_by_name("r","","",kpoints[i]);
      
      if (Ek.size() < VB+1) {ERROR("Problem reading k-points from "+filename);}
      double E_vb = atof(Ek[VB].text().get());
      double E_cb = atof(Ek[VB+1].text().get());
      if (E_vb > VBM) {
    	VBM = E_vb;
      }
      if (E_cb < CBM) {
    	CBM = E_cb;
      }
    }
    
    return CBM - VBM;
  
  } else if (property == "fermi") {

    return (atof(xml.get_node_by_name("i","name","efermi").text().get()));
    
  } else if (property == "volume") {
    
    return (atof(xml.get_node_by_name("i","name","volume").text().get()));
    
  } else if (property == "nelectrons") {
    
    return atoi(xml.get_node_by_name("i","name","NELECT").text().get());

  } else {
    ERROR("Property \'"+property+"\' not implemented");
  }
  
}

// ########################################################
// ########################################################



// ########################################################
//                       GET_DENSITY
// ########################################################
// Reads the charge density from the VASP CHGCAR file.

Grid_data VASP::get_density(string directory, int step) {
  
  int temp_int, Natoms=0;
  
  Grid_data density, downsample;
  
  this->open(directory+"/CHGCAR");
  
  this->get_line();
  this->get_line();
  Line >> density.a;

  this->get_line();
  Line >> density.v1[0] >> density.v1[1] >> density.v1[2];
  this->get_line();
  Line >> density.v2[0] >> density.v2[1] >> density.v2[2];
  this->get_line();
  Line >> density.v3[0] >> density.v3[1] >> density.v3[2];

  this->get_line();
  this->get_line();

  while (!Line.eof()) {
    Line >> temp_int; 
    Natoms += temp_int;
  }

  this->get_line();

  this->skip_lines(Natoms+1);

  this->get_line();
  Line >> density.Nx() >> density.Ny() >> density.Nz();
  
  density.set_dV();

  int Ntotal = density.N();
  if (Ntotal <= 0) { ERROR("Did not read density correctly from file \""+directory+"/CHGCAR"+"\"");}

  int count = 0;
  double temp_double;
  density.reserve(Ntotal);
  this->get_line();
  while(count < Ntotal) {
    while(!Line.eof()) {
      Line >> temp_double;
      density.push_back(temp_double);
      count++;
    }
    this->get_line();
  }
  
  if (!file.eof() && !line.empty() && line.substr(0,12).compare("augmentation")) {
    ERROR("Charge density was not read correctly");
  }
  
  file.close();
  
  REAL Nelectrons = density.integrate()*density.dV/density.volume;
  
  REAL average;
  int count1;

  downsample.Nx() = 0;
  downsample.Ny() = 0;
  downsample.Nz() = 0;
  for (int i=0; i<density.Nx(); i+=step) {
    downsample.Nx() = downsample.Nx() + 1; 
    for (int j=0; j<density.Ny(); j+=step) {
      if (i == 0) { downsample.Ny() = downsample.Ny() + 1; }
      for (int k=0; k<density.Nz(); k+=step) {
	if (i == 0 && j == 0) {downsample.Nz() = downsample.Nz() + 1; }
	count1 = 0;
	average = 0.0;
	for (int ii=i; (ii<i+step && ii<density.Nx()); ii++) {
	  for (int jj=j; (jj<j+step && jj<density.Ny()); jj++) {
	    for (int kk=k; (kk<k+step && kk<density.Nz()); kk++) {
	      ++count1;
	      average += density(ii,jj,kk);
	    }
	  }
	}
	downsample.push_back(average/(density.volume*(double)(count1)));
      }
    }
  }
  
  
  downsample.volume = density.volume;
  double old_Ncells = density.Nx()*density.Ny()*density.Nz();
  double new_Ncells = downsample.Nx()*downsample.Ny()*downsample.Nz();
  downsample.set_dV(density.get_dV()*old_Ncells/new_Ncells);

  return downsample;
  
}

// ########################################################
// ########################################################




// ########################################################
//                       READ_BAND_GAP
// ########################################################
// Determines the band gap from a VASP run.

REAL VASP::read_band_gap(string filename) {
  
  this->open(filename);

  this->skip_lines(5);
  
  this->get_line();
  Line >> this->Nelect;
  
  int valence = int(Nelect/2);
  
  for (int i=0; i<2; i++) {
    this->get_line();
  }

  int band;
  double CBM = 9.0e9;
  double VBM = -9.0e9;
  double E;

  this->get_line();
  while(!file.eof()) {
    Line >> band;
    
    if (band == valence) {
      Line >> E;
      if (E > VBM) { VBM = E; }
    } else if (band == valence+1) {
      Line >> E;
      if (E < CBM) { CBM = E; }
    }

    this->get_line();
  }

  file.close();

  return CBM-VBM;

}
    
// ########################################################
// ########################################################



   
  

// ########################################################
//                       READ_NELECTRONS
// ########################################################
// Reads the number of electrons from a VASP run.

REAL VASP::read_Nelectrons(string filename) {
  
  if (this->Nelect) {
    
    return this->Nelect;
    
  }
  
  this->open(filename);
  
  for (int i=0; i<5; i++) {
    this->get_line();
  }
  
  this->get_line();
  Line >> this->Nelect;
  return Nelect;

}

// ########################################################
// ########################################################


