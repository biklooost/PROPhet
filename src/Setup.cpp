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


#include "Setup.h"
#include "Parallel.h"
#include <time.h>
#include "Art.h"
#include <sys/stat.h>
#include <limits>
#include "xml_reader.h"


// ########################################################
//                       Constructor
// ########################################################
//

Setup::Setup(string input_file)
{


  mpi = new Parallel;

  if (mpi->io_node()) {
    cout << endl;
    cout << ART::PROPhet <<endl;
    cout << "run on:  ";
    time_t rawtime;
    time (&rawtime);
    cout << asctime(localtime(&rawtime))<<endl;
  }

  unsigned int seed = time(NULL);
  if (mpi->io_node()) {
    cout << "random seed =  "<<seed << endl;
  }
  srand(seed);

  my_is_potential = false;

  read_input(input_file);

  if (F.Ninputs() == 0) {
    ERROR("No functional inputs given in input file");
  }
  if (F.output().empty()) {
    ERROR("No functional output given in input file");
  }

  if (!F.input_file().empty()) {

    if (!my_is_potential) {
      ifstream input_file;
      input_file.open(F.input_file().c_str());
      if (!input_file.is_open()) {
        ERROR("Could not open network restart file \""+F.input_file()+"\"");
      }

      if (mpi->io_node()) {
        cout << endl << endl << "Loading state from file \""+F.input_file()+"\""<<endl<<endl;
      }
      F.read(input_file);
      input_file.close();
    }

  }

  split_systems();

  print_details();

}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

Setup::~Setup()
{

  delete mpi;

}

// ########################################################
// ########################################################




// ########################################################
//                       READ_INPUT
// ########################################################
// Reads the input file and sets up the run.

void Setup::read_input (string filename)
{

  string line;
  istringstream Line;
  string KEY, key;
  string string_value;
  int int_value;
  REAL real_value;
  REAL Rcut_max = 0;

  ifstream* file = new ifstream();

  file->open(filename.c_str());
  if (!file->is_open()) {
    ERROR("Could not open input file \""+filename+"\"");
  }

  getline(*file, line);

  input_files.push_back(file);

  while (1) {

    line = clean_line(line);

    if (input_files.front()->eof() && line.empty()) {
      break;
    }

    if (input_files.back()->eof()) {
      input_files.back()->close();
      delete input_files.back();
      input_files.pop_back();
    }

    if (line.empty()) {
      getline(*input_files.back(), line);
      continue;
    }
    Line.clear();
    istringstream Line;
    Line.str(line);
    Line >> KEY;

    key = to_lower(KEY);

    if (key == "data") {

      if (mpi->io_node()) {
        cout << endl << "Found data block"<<endl<<endl;
      }
      property_map.clear();
      while (Line >> KEY >> string_value) {
        key = to_lower(KEY);
        if (property_map.count("key") != 0) {
          ERROR("Key \""+key+"\" specified more than once in data block");
        }
        if (to_lower(KEY) == "code") {
          property_map.insert(pair<string,string>(key,to_lower(string_value)));
          if (to_lower(string_value) == "prophet") {
            ifstream fname("PROPhet.xml",std::ifstream::in);
            string ltemp;
            do {
              std::getline(fname,ltemp);
            } while(ltemp.find("<nsystem>") == std::string::npos);
            xml_reader xml;
            xml.read_string(ltemp);
            int nsys = atoi(xml.get_node_by_name("nsystem").child_value());
            for (int isys = 0; isys < nsys; isys++) {
              map<string,string> new_system;
              ostringstream t;
              t << isys;
              new_system.insert(pair<string,string>("base",t.str()));
              for (map<string,string>::iterator it=property_map.begin(); it!=property_map.end(); ++it) {
                if (!it->first.compare("code")) {
                  new_system.insert(pair<string,string>(it->first,it->second));
                } else {
                  new_system.insert(pair<string,string>(it->first,new_system["base"]+"/"+it->second));
                }
              }
              this->systems.push_back(new_system);
            }
            fname.close();
          }

        } else {
          property_map.insert(pair<string,string>(key, string_value));
        }
        if (mpi->io_node()) {
          cout << "   " << key << " --> " << string_value << endl;
        }
      }

    } else if (key == "input") {

      while (Line >> string_value) {
        F.my_inputs.push_back(to_lower(string_value));
        if (string_value == "structure") {
          my_is_potential = true;
        }
      }
    } else if ( key == "variance_bounds" ) {
      for (int cnt = 0; cnt < 2; cnt++) {
        Line >> int_value;
        F.my_bounds[cnt] = int_value;
      }
    } else if (key == "conv_cd") {
      Line >> int_value;
      F.my_conv = int_value;
    } else if (key == "conv_cd_stride") {
      Line >> int_value;
      F.my_conv_stride = int_value;
    } else if (key == "norm_cd") {
      Line >> real_value;
      F.my_norm_cd = true;
      F.my_norm_cd_val = real_value;
    } else if (key == "sgd") {
      Line >> string_value;
      if (string_value == "true" || string_value == "t" || string_value == "1") {
        F.my_SGD = true;
      } else {
          F.my_SGD = false;
      }
      //F.my_SGD = true;
      //F.my_SGD_cnt = int_value;
    } else if (key == "nsave") {
      Line >> int_value;
      F.my_Nbackup = int_value;
      F.my_Tbackup = true;
    } else if (key == "free_energy" ) {
      //}else if (key.compare("free_energy" ) ||key.compare("formation_energy")||key.compare("cohesive_energy")||key.compare("coh_energy")){
      //ERROR("Free Energy Training currently under development");
      Line.precision(9);
      while (Line >> KEY >> real_value) {
        F.my_FE[KEY] = real_value;
      }
      Line.precision(6);
    } else if (key == "output") {

      Line >> string_value;
      F.my_output = to_lower(string_value);

    } else if (key == "hidden" ) {

      while (Line >> int_value) {
        if (int_value == 0) {
          continue;
        }
        F.my_hidden.push_back(int_value);
      }

    } else if (key == "transfer") {

      Line >> string_value;
      F.my_transfer_function = to_lower(string_value);

    } else if (key == "checkpoint_in") {

      Line >> string_value;
      F.my_input_file = string_value;

    } else if (key == "checkpoint_out") {

      Line >> string_value;
      F.my_output_file = string_value;

    } else if (key == "niterations") {

      Line >> int_value;
      F.my_Niterations = int_value;

    } else if (key == "ncheckpoint") {

      Line >> int_value;
      F.my_Ncheckpoint = int_value;

    } else if (key == "precondition_output") {

      Line >> string_value;
      if (string_value == "true" || string_value == "t" || string_value == "1") {
        F.my_output_precondition = true;
      }
      
    }else if (key == "precondition_input") {
      Line >> string_value;
      if (string_value == "true" || string_value == "t" || string_value == "1") {
        F.my_input_precondition = true;
      } else {
          F.my_input_precondition = false;
      }  

    } else if (key == "print_gvector" || key == "print_gvectors") {
      Line >> string_value;
      if (string_value == "true" || string_value == "t" || string_value == "1") {
        F.my_printgvectors = true;
      }
    } else if (key == "network_type") {

      Line >> string_value;
      F.my_network_type = to_lower(string_value);

    } else if (key == "algorithm") {

      Line >> string_value;
      F.my_training_algorithm = to_lower(string_value);

    } else if (key == "threshold") {

      Line >> real_value;
      F.my_threshold = real_value;
      if (F.my_line_min_epsilon == 0) {
        F.my_line_min_epsilon = (0.1*real_value > 100*std::numeric_limits<REAL>::epsilon() ? 0.1*real_value : 100*std::numeric_limits<REAL>::epsilon());
      }

    } else if (key == "sd_momentum") {

      Line >> real_value;
      F.my_sd_momentum = real_value;

    } else if (key == "debug") {

      Line >> real_value;
      F.my_debug = real_value;

    } else if (key == "nprint") {

      Line >> int_value;
      F.my_Nprint = int_value;

    } else if (key == "lambda_max") {

      Line >> real_value;
      F.my_lambda_max = real_value;

    } else if (key == "ncycles") {

      Line >> int_value;
      F.my_Ncycles = int_value;

    } else if (key == "rcut") {

      Line >> real_value;
      F.my_Rcut = real_value;

    } else if (key == "nradial") {

      Line >> int_value;
      F.my_Nradial = int_value;

    } else if (key == "nangular") {

      Line >> int_value;
      F.my_Nangular = int_value;

    } else if (key == "g1") {
      vector<REAL> temp;
      while (Line >> real_value) {
        temp.push_back(real_value);
      }
      if (temp[0] > Rcut_max) {
        Rcut_max = temp[0];
      }
      if (temp.size() != 2) {
        ERROR("A G1 function requires 2 parameters");
      }
      Line.clear();
      F.G1.push_back(temp);

    } else if (key == "g2") {

      vector<REAL> temp;
      while (Line >> real_value) {
        temp.push_back(real_value);
      }
      if (temp[0] > Rcut_max) {
        Rcut_max = temp[0];
      }
      if (temp.size() != 4) {
        ERROR("A G2 function requires 4 parameters");
      }
      Line.clear();
      F.G2.push_back(temp);

    } else if (key == "g3") {

      vector<REAL> temp;
      while (Line >> real_value) {
        temp.push_back(real_value);
      }
      if (temp[0] > Rcut_max) {
        Rcut_max = temp[0];
      }
      if (temp.size() != 5) {
        ERROR("A G3 function requires 5 parameters");
      }
      Line.clear();
      F.G3.push_back(temp);

    } else if (key == "g4") {

      vector<REAL> temp;
      while (Line >> real_value) {
        temp.push_back(real_value);
      }
      if (temp[0] > Rcut_max) {
        Rcut_max = temp[0];
      }
      if (temp.size() != 5) {
        ERROR("A G4 function requires 5 parameters");
      }
      Line.clear();
      F.G4.push_back(temp);

    } else if (key == "include") {

      Line >> string_value;
      filename = string_value;
      input_files.push_back(new ifstream);
      input_files.back()->open(string_value.c_str());
      if (!input_files.back()->is_open()) {
        ERROR("Could not open nested input file '"+string_value+"'");
      }

    } else if (key == "downsample") {

      Line >> int_value;
      F.my_downsample = int_value;

    } else if (key == "seed") {

      Line >> int_value;
      srand(int_value);

    } else if (key == "regularization") {

      Line >> real_value;
      F.my_regularization = real_value;

    } else if (key == "line_min_epsilon") {

      Line >> real_value;
      F.my_line_min_epsilon = real_value;

    } else if (key == "t0") {

      Line >> real_value;
      F.T0 = real_value;

    } else if (key == "tf") {

      Line >> real_value;
      F.Tf = real_value;

    } else if (key == "nsteps_per_t") {

      Line >> int_value;
      F.Nsteps_per_T = int_value;

    } else if (key == "step_size") {

      Line >> real_value;
      F.step_size = real_value;

    } else if (key == "dt") {

      Line >> real_value;
      F.dT = real_value;


    } else if (key == "dropout") {

      Line >> real_value;
      F.my_dropoutP = real_value;

    } else {

      struct stat info;
      if (stat(KEY.c_str(), &info) != 0) {
        ERROR("Unrecognized flag \""+KEY+"\" in file \""+filename+"\"");
      }

      if (Rcut_max) {
        F.Rcut(Rcut_max);
      }

      map<string,string> new_system;
      while (Line >> string_value) { }
      string train;
      if (string_value == "train" || string_value == "val" || string_value == "test") {
        train = string_value;
      } else {
        train = "train";
      }


      new_system.insert(pair<string,string>("base",KEY));
      new_system.insert(pair<string,string>("train",train));
      for (map<string,string>::iterator it=property_map.begin(); it!=property_map.end(); ++it) {
        if (!it->first.compare("code")) {
          new_system.insert(pair<string,string>(it->first,it->second));
        } else {
          new_system.insert(pair<string,string>(it->first,new_system["base"]+"/"+it->second));
        }
      }

      this->systems.push_back(new_system);

    }

    getline(*input_files.back(), line);

  }

  input_files.front()->close();
  delete input_files.front();

}

// ########################################################
// ########################################################





// ########################################################
//                       SPLIT_SYSTEMS
// ########################################################
// For a parallel run this splits up the work over cores.

void Setup::split_systems()
{

  vector<map<string,string> > my_systems;

  if (mpi->Nprocs() > systems.size()) {
    if (mpi->io_node()) {
      cout << endl;
      cout << "---Warning---"<<endl;
      cout << "There are more processors (" << mpi->Nprocs() << ")"
           << " than systems (" << systems.size() << ")"<< endl;
      cout << "There will be unused processors" << endl;
      cout << "and the calculation might be unstable"<<endl;
    }
  }

  for (int i=mpi->rank(); i<systems.size(); i+=mpi->Nprocs()) {
    for (int j=0; j<F.Ninputs(); j++) {
      if (systems.at(i).count(F.inputs(j)) == 0) {
        systems.at(i).insert(pair<string,string>(F.inputs(j),systems.at(i)["base"]));
      }
    }
    if (systems.at(i).count(F.output())==0) {
      systems.at(i).insert(pair<string,string>(F.output(),systems.at(i)["base"]));
    }
    systems.at(i).insert(pair<string,string>("user",systems.at(i)["base"]));
    systems.at(i).erase("base");
    my_systems.push_back(systems.at(i));
  }



  if (systems.size() != my_systems.size()) {
    systems.assign(my_systems.begin(),my_systems.end());
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       PRINT_DETAILS
// ########################################################
// Tells the user what it found in the input file.

void Setup::print_details()
{

  int Ntotal_systems = mpi->Reduce(systems.size(), MPI_SUM);

  if (mpi->io_node()) {
    cout << endl << "Found " << Ntotal_systems << plural(" system",Nsystems()) << endl  << endl;
  }

  if (mpi->io_node()) {
    if (!F.my_FE.empty()) {
      map<string,REAL>::iterator it;
      cout << "############" << endl << "TRAINING TO FORMATION ENERGIES" << endl;
      cout << "Species Energies: " << endl;
      cout.precision(9);
      for (it = F.my_FE.begin(); it != F.my_FE.end(); it++) {
        cout << it->first << " : " << it->second << endl;
      }
      cout.precision(6);
      cout <<"############" << endl << endl;
    }
    cout << "Functional parameters" <<endl;
    cout << F.Ninputs() << plural(" input",F.Ninputs()) << " :  ";
    for (int i=0; i<F.Ninputs(); i++) {
      cout << F.inputs(i) << "   ";
    }
    cout << endl;
    cout << F.Nlayers() << plural(" hidden layer",F.Nlayers()) << " :  ";
    for (int i=0; i<F.Nlayers(); i++) {
      cout << F.hidden(i) << "   ";
    }
    //cout << endl << "Dropout Probability: " << 1 - F.dropout() << endl;
    cout << endl << "output :  " << F.output() << endl;

    //cout << endl;
    /*
    if (F.step_size == 0.2) {
        cout << "Default step size of " << F.step_size << endl;
    } else {
        cout << "step size : " << F.step_size << endl;
    }*/
    //cout << endl;
  }

  int min,max,temp;

  min = mpi->Reduce(systems.size(), MPI_MIN);
  max = mpi->Reduce(systems.size(), MPI_MAX);

  if (mpi->io_node()) {
    cout << endl;
    cout << "Running on " << mpi->Nprocs() << " core";
    if (mpi->Nprocs()==1) {
      cout << endl << endl;
    } else {
      cout << "s..." << endl;
      cout << "Minimum core load = " << min <<endl;
      cout << "Maximum core load = " << max << endl <<endl;
    }
  }

}

// ########################################################
// ########################################################





// ########################################################
//                       CLEAN_LINE
// ########################################################
// Removes unnecessary characters from a line for easier
// parsing.

string Setup::clean_line(string line_in)
{

  if (line_in.empty()) {
    return line_in;
  }

  for (int i=0; i<line_in.size(); i++) {
    if (line_in[i] == '=') {
      line_in.replace(i,1," ");
    }
    if (line_in[i] == '#') {
      line_in.resize(i);
    }
  }

  if (!line_in.empty()) {
    string whitespace=" \t\n";
    if (line_in.find_first_not_of(whitespace) != string::npos) {
      line_in = line_in.substr(line_in.find_first_not_of(whitespace),line_in.find_last_not_of(whitespace)-line_in.find_first_not_of(whitespace)+1);
    } else {
      line_in.clear();
    }
  }

  return line_in;
}

// ########################################################
// ########################################################




