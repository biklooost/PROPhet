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
// This is a wrapper to MPI. It handles the case where the user
// requests not to use MPI by mimicking the interface of MPI.
// ####################################################################



#include "Parallel.h"





// ########################################################
//                       Constructor
// ########################################################
//

Parallel::Parallel()
{

  my_rank = 0;
  my_Nprocs = 1;

#ifdef USE_MPI
  int init;
  MPI_Initialized(&init);
  if (!init) {
    MPI_Init(NULL, NULL);
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &my_Nprocs);
#endif
}

// ########################################################
// ########################################################




#ifndef USE_MPI
int Parallel::io_node()
{
#ifndef NOT_MD
  return 0;
#else
  return 1;
#endif
}
#else
int Parallel::io_node()
{
  return !my_rank;
}
#endif



// ########################################################
//                       Destructor
// ########################################################
//

Parallel::~Parallel()
{

  // ########################################################
  // ########################################################




#ifdef USE_MPI
  int final;
  MPI_Finalized(&final);
  if (!final) {
    MPI_Finalize();
  }
#endif

}



