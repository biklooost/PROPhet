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





#ifndef __PARALLEL
#define __PARALLEL

#include <vector>
#include <string>

#ifdef USE_MPI
#include "mpi.h"
#else
#ifdef NOT_MD
#define MPI_Op int
enum {
  MPI_SUM,
  MPI_MIN,
  MPI_MAX,
};
#endif
#endif

#include "Common.h"


using namespace std;

class Parallel
{

public:

  Parallel();
  ~Parallel();



#ifdef USE_MPI
  int io_node();
  inline int rank()
  {
    return my_rank;
  }
  inline int Nprocs()
  {
    return my_Nprocs;
  }

  inline void Abort()
  {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
  }


  // scalar reduce via MPI
  template <typename type> inline type Reduce(type in, MPI_Op operation)
  {
    type receive;
    if (MPI_Reduce(&in, &receive, 1, MPI_Type(in), operation, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
      ERROR("MPI Reduce failed");
    }
    return receive;
  }

  // vector reduce via MPI
  template <typename type> inline vector<type> Reduce(vector<type> in, MPI_Op operation)
  {
    vector<type> receive(in.size());
    if (MPI_Reduce(&in.front(), &receive.front(), in.size(), MPI_Type(in.at(0)), operation, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
      ERROR("MPI Reduce failed");
    }
    return receive;
  }

  // vector send via MPI
  template <typename type> inline void Send(vector<type> in, int destination, int tag=0)
  {
    if (MPI_Send(&in.front(), in.size(), MPI_Type(in.at(0)), destination, tag, MPI_COMM_WORLD) != MPI_SUCCESS) {
      ERROR("MPI Send failed");
    }
  }
  // vector receive via MPI
  template <typename type> inline void Receive(vector<type> out, int size, int from, int tag=0)
  {
    out.resize(size);
    if (MPI_Recv(&out.front(), size, MPI_Type(out.at(0)), from, tag, MPI_COMM_WORLD) != MPI_SUCCESS) {
      ERROR("MPI Receive failed");
    }
  }

  // scalar Broadcast via MPI
  template <typename type> inline type Bcast(type in, int root=0)
  {
    if (MPI_Bcast(&in,1,MPI_Type(in),root,MPI_COMM_WORLD) != MPI_SUCCESS) {
      ERROR("MPI Bcast failed");
    }
    return in;
  }
  // Explicit specializtion for bool to work around openmpi convention changes
  template <bool> inline bool Bcast(bool in, int root=0)
  {
    int var = static_cast<int>(in);
    if (MPI_Bcast(&var,1,MPI_INT,root,MPI_COMM_WORLD) != MPI_SUCCESS) {
      ERROR("MPI Bcast failed");
    }
    return static_cast<bool>(var);
  }

  // vector Broadcast via MPI
  template <typename type> inline vector<type> Bcast(vector<type> in, int size, int root=0)
  {
    in.resize(size);
    if (MPI_Bcast(&in.front(), in.size(), MPI_Type(in.at(0)), root, MPI_COMM_WORLD) != MPI_SUCCESS) {
      ERROR("MPI Bcast failed");
    }
    return in;
  }

  //vector Gather (for when all procs are sending the same size object)
  template <typename type> vector<type> Gather(vector<type> &in, int size, int root=0)
  {
    type* receive_ptr;
    int N=this->Nprocs()*size;
    if (this->rank() == root) {
      receive_ptr = new type [N];
    }
    if (MPI_Gather(&in.front(), size, MPI_Type(in.at(0)), receive_ptr, N, MPI_Type(in.at(0)), root, MPI_COMM_WORLD) != MPI_SUCCESS) {
      ERROR("MPI Gather failed");
    }
    vector<type> output;
    if (this->rank()==root) {
      output.assign(receive_ptr, receive_ptr+N);
      delete [] receive_ptr;
    }
    return output;
  }

  // vector Gatherv (for different size objects from each proc)
  template <typename type> vector<type> Gatherv(vector<type> &in, int root=0)
  {
    type* receive_ptr;
    int*  sizes;
    int*  strides;

    if (this->rank()==root) {
      sizes = new int [this->Nprocs()];
      strides = new int [this->Nprocs()];
    }
    int my_size = in.size();
    MPI_Gather(&my_size, 1, MPI_INT, sizes, 1, MPI_INT, root, MPI_COMM_WORLD);

    int total=0;
    if (this->rank() == root) {
      total = sizes[0];
      strides[0] = 0;
      for (int i=1; i<this->Nprocs(); i++) {
        total += sizes[i];
        strides[i] = strides[i-1]+sizes[i-1];
      }
      receive_ptr = new type [total];
    }

    if (MPI_Gatherv(&in.front(), in.size(), MPI_Type(in.at(0)), receive_ptr, sizes, strides, MPI_Type(in.at(0)), root, MPI_COMM_WORLD) != MPI_SUCCESS) {
      ERROR("MPI_Gather failed");
    }

    vector<type> output;
    if (this->rank() == root) {
      output.assign(receive_ptr,receive_ptr+total);
      delete [] receive_ptr;
      delete [] sizes;
      delete [] strides;
    }

    return output;

  }


#else

  int io_node();
  inline int rank()
  {
    return 0;
  }
  inline int Nprocs()
  {
    return 1;
  }

  inline void Abort()
  {
    exit(1);
  }

  template <typename type > inline type Reduce(type in, int operation)
  {
    return in;
  }

  template <typename type > inline vector<type> Reduce(vector<type> in, int operation)
  {
    return in;
  }

  template <typename type > inline type Bcast(type in, int root=0)
  {
    return in;
  }

  template <typename type > inline vector<type> Bcast(vector<type> in, int root=0)
  {
    return in;
  }

  template <typename type > inline vector<type> Gather(vector<type> &in, int size, int root=0)
  {
    return in;
  }

  template <typename type > inline vector<type> Gatherv(vector<type> &in, int root=0)
  {
    return in;
  }

#endif


private:

  int my_rank;
  int my_Nprocs;

#ifdef USE_MPI
  // These are the only data types supported by this wrapper.
  // Others can be used with standard MPI_... calls, or by adding
  // a fuction here to return their type.
  inline MPI_Datatype MPI_Type(int in)
  {
    return MPI_INT;
  }
  inline MPI_Datatype MPI_Type(double in)
  {
    return MPI_DOUBLE;
  }
  inline MPI_Datatype MPI_Type(float in)
  {
    return MPI_FLOAT;
  }
  inline MPI_Datatype MPI_Type(long unsigned int in)
  {
    return MPI_UNSIGNED_LONG;
  }
  inline MPI_Datatype MPI_Type(long double in)
  {
    return MPI_LONG_DOUBLE;
  }
  //  inline MPI_Datatype MPI_Type(bool in) { return MPI_BOOL; }
#endif


};



#endif
