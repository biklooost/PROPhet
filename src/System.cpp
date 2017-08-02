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
// A general class to hold information about system properties.
// ####################################################################



#include "System.h"

#include "VASP.h"
#include "QE.h"
#include "FHIAIMS.h"




// ########################################################
//                       Constructor
// ########################################################
//

System::System(map<string,string> files, Functional_params *F) {
  
  DFT_IO *DFT;
  
  this->Prefactor = 1.0;
  
  if (files["code"] == "vasp") {
    
    DFT = new VASP();
    
  } else if (files["code"] == "qe") {

    DFT = new QE();

  } else if (files["code"] == "fhiaims") {

    DFT = new FHIAIMS();
    
  } else {
    
    ERROR("Interface to code '"+files["code"]+"' has not been implemented");
    
  }
  
  
  // Get input values
  for (int i=0; i<F->Ninputs(); i++) {
    string input = F->inputs(i);
    
    if (!input.compare("density")) {
      if (!Density.N()) { this->Density = DFT->get_density(files["density"],F->sample_step()); }
      if(F->NormCD()) { this->Density.normalize(F->NormCD_val());}
      properties.push_back(this->Density.as_vector_ptr());
      Prefactor *= Density.get_dV();
      if (F->output_is_intensive()) {
	Prefactor /= Density.volume;
      } 
    } else if (input == "density^2") {
      
      //Not implemented 
      
    } else if (input.find("user") == 0) {
      int prop = atoi(input.substr(4).c_str());
      data.insert(pair<string,vector<REAL> >(input,DFT->get_user_property(prop, files["user"])));
      properties.push_back(&data[input]);
    } else if (input == "structure") {
      this->structure = DFT->read_structure(files["structure"]);
      properties.lock(true);
    } else if (input == "random") {
      data.insert(pair<string,vector<REAL> >(input, vector<REAL>(1,RAND::Uniform())));
      properties.push_back(&data[input]);
    } else {
      data.insert(pair<string,vector<REAL> >(input,vector<REAL>(1,DFT->get_property(input,files[input]))));
      properties.push_back(&data[input]);
    }
  }


  // Get output values
  if (F->output() == "gw_gap") {
      REAL GW_gap = DFT->get_property("gw_gap",files["gw_gap"]);
      properties.target(GW_gap);
  } else if (F->output().find("user")==0) {
    int prop = atoi(F->output().substr(4).c_str());
    properties.target(DFT->get_user_property(prop, files["user"]).at(0));
  } else if (F->output() == "energy") {
      REAL energy = DFT->get_property(F->output(), files[F->output()]);;
      if (!F->FE().empty() || !this->structure.FE.empty()){
          energy = this->structure.train_Local(F,energy);
      } 
      //cout << "What we are training to: " << energy << endl;
      properties.target(energy);
    } else {
    properties.target(DFT->get_property(F->output(), files[F->output()]));
  }
  delete DFT;

}

// ########################################################
// ########################################################




// ########################################################
//                       Constructor
// ########################################################
//

System::System() {
  
}

// ########################################################
// ########################################################



// ########################################################
//                       Destructor
// ########################################################
//

System::~System() {
  
}

// ########################################################
// ########################################################








