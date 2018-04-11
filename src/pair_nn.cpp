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
// This is an interface to the LAMMPS MD code from Sandia National
// Labs. This allows one to use PROPhet potentials in MD simulations
// within LAMMPS.
// ####################################################################


#ifdef PAIR_CLASS
PairStyle(nn,PairNN)
#else

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_nn.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include "Common.h"
#include "Setup.h"
#include "System.h"
#include "Network.h"
#include "NN.h"


using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */


// ########################################################
//                       Constructor
// ########################################################
//

PairNN::PairNN(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  feenableexcept(FE_INVALID | FE_OVERFLOW);

}

// ########################################################
// ########################################################


// ########################################################
//                       Destructor
// ########################################################
//

PairNN::~PairNN()
{
  //  if (potential) { delete potential; }
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

  }
}

// ########################################################
// ########################################################


// ########################################################
//                       COMPUTE
// ########################################################
// Determine the energy and forces for the current structure.

void PairNN::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;

  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else {
    evflag = vflag_fdotr = 0;
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double Energy;
  vector<vector<REAL> > dE_dG(nlocal,vector<REAL>(1,0.0));

  /* ----------------------------------------------------------------------
      set eflag,vflag for current iteration
      invoke matchstep() on all timestep-dependent computes to clear their arrays
      eflag/vflag based on computes that need info on this ntimestep
      eflag = 0 = no energy computation
      eflag = 1 = global energy only
      eflag = 2 = per-atom energy only
      eflag = 3 = both global and per-atom energy
      vflag = 0 = no virial computation (pressure)
      vflag = 1 = global virial with pair portion via sum of pairwise interactions
      vflag = 2 = global virial with pair portion via F dot r including ghosts
      vflag = 4 = per-atom virial only
      vflag = 5 or 6 = both global and per-atom virial
  */
  sysdata->structure.inum = list->inum;
  sysdata->structure.ilist.assign(list->ilist, list->ilist + list->inum );
  sysdata->structure.numneigh.assign(list->numneigh, list->numneigh + list->inum );
  sysdata->structure.firstneigh.assign(list->inum, vector<int>(1,0));

  sysdata->structure.set_structure(x,type,atom->nlocal,atom->nghost);

  for (int ii=0; ii<list->inum; ii++) {
    sysdata->structure.firstneigh.at(ii).assign(list->firstneigh[ii], list->firstneigh[ii] + list->numneigh[ii]);
  }
  sysdata->structure.periodic = 0;

  for(int ii = 0; ii< sysdata->structure.ilist.size(); ii++) {
    Energy = potential->evaluate_MD(ii, type[ii], dE_dG[ii]);
    if (eflag_global) {
      eng_vdwl += Energy;
    }
    if (eflag_atom) {
      eatom[ii] += Energy;
    }

  }
  /*
  for(int ii=0; ii<dE_dG.size();ii++){
      //cout << ii << endl;
      for (int jj=0; jj<dE_dG[ii].size(); jj++) {
          cout << dE_dG[ii][jj] << '\t';
      }
      cout << endl;
  }*/

  if (comm->me == 0) {
    if (eflag_global) {
      eng_vdwl += potential->energy_shift();
    }
  }


  sysdata->structure.Get_Forces(dE_dG, f);

  if (vflag_fdotr) {
    virial_fdotr_compute();
  }

}

// ########################################################
// ########################################################

// ########################################################
//                       ALLOCATE
// ########################################################
// Allocates all necessary arrays.

void PairNN::allocate()
{

  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 0;
    }
  }

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

}

// ########################################################
// ########################################################

// ########################################################
//                       COEFF
// ########################################################
// Assigns the proper network to an atom type.

void PairNN::coeff(int narg, char **arg)
{

  if (narg != 2 && narg != 3) {
    error->all(FLERR,"Format of pair_coeff command is\npair_coeff atom_type parameter_file\n");
  }

  int this_type = force->numeric(FLERR, arg[0]);
  double this_cut;
  if (narg == 3) {
    this_cut = force->numeric(FLERR,arg[2]);
  }

  if (!allocated) {
    allocate();
  }

  if (comm->me==0) {
    cout << "Initializing atom type "<<this_type<<" with file "<<arg[1]<<endl;
  }


  potential->insert_atom_type(this_type, arg[1], sysdata);
  sysdata->structure.lammps_conv[this_type] = potential->ret_params().current_atom_type; //testing
  stringstream ss;
  ss << arg[1];
  potential_names.push_back(ss.str());
  potential_numbers.push_back(this_type);

  for (int i=1; i<=atom->ntypes; i++) {
    setflag[i][this_type] = 1;
    if (narg==3) {
      cutsq[i][this_type] = this_cut*this_cut;
    }
  }

}

// ########################################################
// ########################################################

// ########################################################
//                       INIT_STYLE
// ########################################################
// Set up the pair style to be a NN potential.

void PairNN::init_style()
{

  int irequest;
  neighbor->cutneighmin = 1.0;
  neighbor->cutneighmax = params.Rcut();
  neighbor->delay = 0;
  neighbor->every = 10;
  neighbor->skin = 1.0;
  irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 1;
  neighbor->requests[irequest]->id=1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 0;

}

// ########################################################
// ########################################################

// ########################################################
//                       INIT_LIST
// ########################################################
//

void PairNN::init_list(int id, NeighList *ptr)
{
  if(id == 1) {
    list = ptr;
  }
}

// ########################################################
// ########################################################

// ########################################################
//                       init_one
// ########################################################
// Initilize 1 pair interaction.  Needed by LAMMPS but not
// used in this style.

double PairNN::init_one(int i, int j)
{
  return sqrt(cutsq[i][j]);//params.Rcut();
}

// ########################################################
// ########################################################



// ########################################################
//                       WRITE_RESTART
// ########################################################
// Writes restart file but not yet implemented.

void PairNN::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  double Rcut = params.Rcut();
  int rank;
  rank = comm->me;
  if (rank == 0 ) {
    fwrite(&Rcut,sizeof(double),1,fp);
    for (int i = 0; i < potential_numbers.size(); i++) {
      int size = potential_names[i].size();
      int type = potential_numbers[i];
      fwrite(&type,sizeof(int),1,fp);
      fwrite(&size,sizeof(int),1,fp);
      for (int j =0; j < potential_names[i].size(); j++) {
        fwrite(&potential_names[i][j],sizeof(char),1,fp);
      }
    }
  }

}

// ########################################################
// ########################################################


// ########################################################
//                       READ_RESTART
// ########################################################
// Reads from restart file.  Tested to work, but use with caution

void PairNN::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
  int rank;
  rank = comm->me;
  char* array[2];
  double Rcut;
  int type;
  int size;
  stringstream Line;
  if (rank == 0) {
    fread(&Rcut,sizeof(double),1,fp);
  }
  MPI_Bcast(&Rcut,1,MPI_DOUBLE, 0, world);
  Line << Rcut;
  array[0] = const_cast<char*>(Line.str().c_str());
  this->settings(1, array);
  params.Rcut(Rcut);
  for (int i = 1; i <= atom->ntypes; i++) {
    if (rank == 0) {
      fread(&type, sizeof(int),1,fp);
      fread(&size, sizeof(int),1,fp);
    }
    MPI_Bcast(&type,1,MPI_INT, 0, world);
    MPI_Bcast(&size,1,MPI_INT, 0, world);
    char tmp[size];
    if (rank == 0) {
      for (int j =0; j < size; j++) {
        fread(&tmp[j],sizeof(char),1,fp);
      }
    }
    MPI_Bcast(&tmp,size,MPI_CHAR, 0, world);
    Line.str(std::string());
    Line << type;
    array[0] = const_cast<char*>(Line.str().c_str());
    array[1] = tmp;
    this->coeff(2, array);
  }
}

// ########################################################
// ########################################################





// ########################################################
//                       WRITE_RESTART_SETTINGS
// ########################################################
// Writes settings to restart file.

void PairNN::write_restart_settings(FILE *fp)
{
  //fwrite(&cut_global,sizeof(double),1,fp);
  //fwrite(&offset_flag,sizeof(int),1,fp);
  //fwrite(&mix_flag,sizeof(int),1,fp);
  //fwrite(&tail_flag,sizeof(int),1,fp);
}

// ########################################################
// ########################################################



/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

// ########################################################
//                       READ_RESTART_SETTINGS
// ########################################################
// Reads settings from restart file.

void PairNN::read_restart_settings(FILE *fp)
{
  /*
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
  */


}

// ########################################################
// ########################################################




/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairNN::write_data(FILE *fp)
{
  /*
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon[i][i],sigma[i][i]);
  */
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairNN::write_data_all(FILE *fp)
{
  /*
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j]);
  */
}

/* ---------------------------------------------------------------------- */

double PairNN::single(int i, int j, int itype, int jtype, double rsq,
                      double factor_coul, double factor_lj,
                      double &fforce)
{
  return 1;
}

/* ---------------------------------------------------------------------- */



// ########################################################
//                       Settings
// ########################################################
// Initializes settings.

void PairNN::settings(int narg, char* argv[])
{
  //cout << comm->me<<" starting"<<endl;
  if (narg != 1) {
    error->all(FLERR,"");
  }
  if (narg == 1) {
    params.Rcut(force->numeric(FLERR,argv[0]));
  } else {

  }
  //Functional_params F;
  params.is_MD = true;
  vector<string> inputs(1,"structure");
  params.set_inputs(inputs);
  params.set_output("energy");

  // Need to set up the system here
  this->sysdata = new System();

  sysdata->structure.set_structure(atom->x, atom->type, atom->nlocal, atom->nghost);

  sysdata->properties.lock(true);
  sysdata->structure.periodic = 0;

  vector<System*> system(1,sysdata);
  this->potential = new Potential(params);
  potential->add_system(sysdata);

  if (!allocated) {
    allocate();
  }

  for (int i=1; i<=atom->ntypes; i++) {
    for (int j=i; j<=atom->ntypes; j++) {
      cutsq[i][j] = params.Rcut()*params.Rcut();
    }
  }
  //  this->params = F;

}

// ########################################################
// ########################################################





#endif
