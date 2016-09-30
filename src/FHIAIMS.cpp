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



#include "FHIAIMS.h"
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

FHIAIMS::FHIAIMS() { }

// ########################################################
// ########################################################





// ########################################################
//                       Destructor
// ########################################################
//

FHIAIMS::~FHIAIMS() { }

// ########################################################
// ########################################################




// ########################################################
//                       READ_STRUCTURE
// ########################################################
// Routine to read the atomic configuration from an 
// FHI-Aims run

Structure FHIAIMS::read_structure(string file) { 
  int Natoms;
  Structure temps;
  temps.CART = true;
  stringstream tvector;
  for (int i = 0; i<3; i++) {temps.a.push_back(0.0);temps.b.push_back(0.0);temps.c.push_back(0.0);}
  vector< vector<REAL> > coordinates;
  
  string species;
  this->open(file);
  stringstream temp;
  string trash;
  REAL tmp1;
  REAL tmp2;
  REAL tmp3;
  string tmp4;
  vector <REAL> row;
  stringstream HOMO;
  vector < vector <REAL> > atom_array;
  temps.periodic = false;
  for (int i=0; i<3; i++) { row.push_back(0.0);}
  do {
    this->get_line();
    if(this->Line.str().find("Number of atoms") != std::string::npos) {
      Line >> trash >> trash >> trash >> trash >> trash >> temps.Natom;
      if (temps.pos.size() == 0 ) {
        for (int i=0; i < temps.Natom; i++) { 
          temps.pos.push_back(row); 
          temps.types.push_back(1);
        } 
      } 
    }
    if(this->Line.str().find("lattice_vector") != std::string::npos && this->Line.str().find("#") == std::string::npos) {
      Line >> trash >> temps.a[0] >> temps.a[1] >> temps.a[2];
      this->get_line();
      Line >> trash >> temps.b[0] >> temps.b[1] >> temps.b[2];
      this->get_line();
      Line >> trash >> temps.c[0] >> temps.c[1] >> temps.c[2];
      temps.periodic = true;
    }
  } while (!this->file.eof());
  if (temps.pos.size() == 0) {
    ERROR("There is somethine wrong with " + file + "\n");
  }else{
    this->open(file);
    do {
      this->get_line();
      if(this->Line.str().find("x [A]") != std::string::npos && this->Line.str().find("#") == std::string::npos) {
        for(int i=0; i< temps.Natom; i++) {
          this->get_line();
          if(this->Line.str().find(":") != std::string::npos) {
            Line >> trash >> trash >> trash >> species >> temps.pos[i][0] >> temps.pos[i][1] >> temps.pos[i][2];
            temps.types[i] = Atom(species);
          }else if(this->Line.str().find("atom") != std::string::npos && this->Line.str().find("#") == std::string::npos) {
            Line >> trash >> temps.pos[i][0] >> temps.pos[i][1] >> temps.pos[i][2] >> species;
            temps.types[i] = Atom(species);
          }
        }
      }
    } while(!this->file.eof());
  }
  return temps;
}

// ########################################################
// ########################################################




// ########################################################
//                       GET_PROPERTY
// ########################################################
// Reads various properties from FHI-Aims run

REAL FHIAIMS::get_property(string property,string file) { 
  stringstream ENERGY;
  string trash;
  REAL energy = 0.0;
  if (property == "energy") {
    string energyline = "| Total energy of the DFT / Hartree-Fock s.c.f. calculation";
    this->open(file);
    do {
      this->get_line();
      if(this->Line.str().find(energyline) != std::string::npos) {
        ENERGY.str(Line.str());
      }
    } while(!this->file.eof());
    ENERGY >> trash >> trash >> trash >> trash >> trash >> trash >> trash >> trash >> trash >> trash >> trash >> energy;
    return energy;
  } else if ( (property == "dft_gap") || (property == "gw_gap") ) {
    REAL bndgap = read_band_gap(file);
    return bndgap;
  }

}

// ########################################################
// ########################################################





// ########################################################
//                       READ_BAND_GAP
// ########################################################
//

REAL FHIAIMS::read_band_gap(string file) { 
	this->open(file);
	stringstream temp,HOMO;
	string trash;
	REAL bandgap;
	do {
		this->get_line();
		if(this->Line.str().find("Overall HOMO-LUMO") != std::string::npos) {
			HOMO.str(Line.str());
		}
	} while (!this->file.eof());
	HOMO >>	trash >> trash >> trash >> bandgap;
	return bandgap;
}

// ########################################################
// ########################################################




// ########################################################
//                       GET_DENSITY
// ########################################################
// Reads the charge density from an FHI-Aims run

Grid_data FHIAIMS::get_density(string CubeFile, int step) { //This is applicable to general cube file
  Grid_data grid, downsample;
	int natom;
	string token;
	grid.a = 1;
	this->open(CubeFile);
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
	do{
	  this->get_line();
	  while(getline(this->Line,token,' ')) {
	    if(token == "") { 
	      continue;
	      	}else{
       		  grid.push_back(atof(token.c_str()));
	    }
	  }
	}while(!this->file.eof());
	grid.set_dV();
	REAL Nelectrons = grid.integrate()*grid.dV;
  
	REAL average;
	int count1;

	downsample.Nx() = 0;
	downsample.Ny() = 0;
	downsample.Nz() = 0;
	for (int i=0; i<grid.Nx(); i+=step) {
	  downsample.Nx() = downsample.Nx() + 1; 
	  for (int j=0; j<grid.Ny(); j+=step) {
	    if (i == 0) { downsample.Ny() = downsample.Ny() + 1; }
	    for (int k=0; k<grid.Nz(); k+=step) {
	      if (i == 0 && j == 0) {downsample.Nz() = downsample.Nz() + 1; }
	      count1 = 0;
	      average = 0.0;
	      for (int ii=i; (ii<i+step && ii<grid.Nx()); ii++) {
          for (int jj=j; (jj<j+step && jj<grid.Ny()); jj++) {
            for (int kk=k; (kk<k+step && kk<grid.Nz()); kk++) {
              ++count1;
              average += grid(ii,jj,kk);
            }
          }
	      }
	      downsample.push_back(average/((double)(count1)));
	    }
	  }
	}
  
	downsample.volume = grid.volume;
	double old_Ncells = grid.Nx()*grid.Ny()*grid.Nz();
	double new_Ncells = downsample.Nx()*downsample.Ny()*downsample.Nz();
	downsample.set_dV(grid.get_dV()*old_Ncells/new_Ncells);

	return downsample;

}

// ########################################################
// ########################################################




// ########################################################
//                       READ_NELECTRONS
// ########################################################
// Gets number of electrons from an FHI-Aims run

REAL FHIAIMS::read_Nelectrons(string file) {
	this->open(file);
	stringstream temp;
	REAL nelect;
	string trash;
	int inttrash;
	do {
		this->get_line();
		if(this->Line.str().find("electrons") != std::string::npos) {
			Line >> trash >> trash >> trash >> inttrash >> trash >> trash >> trash >> trash >> trash >> this->Nelect;
		}
	} while (!this->file.eof());
	return this->Nelect;
}

// ########################################################
// ########################################################


