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
// This is the class that actually handles training a network. It
// contains a number of training algorithms including several flavors
// of steepest descent, L-BFGS, and Rprop. It handles interprocess
// communication internally by having only the root process determine
// the new trial set of parameters (via the selected algorithm) then
// broadcasting the set to other processes.
// ####################################################################





#ifndef __OPTIMIZER
#define __OPTIMIZER

#include <vector>
#include <deque>

#include "Parallel.h"
#include "Network_node.h"
#include "NF_node.h"
#include "Common.h"
#include "Functional_params.h"


using namespace std;

class Optimizer
{

public:

  Optimizer(Parallel *mpi_ptr);
  Optimizer(Parallel *mpi_ptr, Functional_params params, int Nsystems);
  ~Optimizer();

  void init(vector<vector<Network_node*> >* in_nodes);
  void internal_init();
  void init(vector<Network_node*> in_nodes);

  void set_params(Functional_params params);
  void update_network(REAL Error, vector<REAL> grad, REAL SSE = 0);

  inline bool is_converged()
  {
    return my_is_converged;
  }
  inline bool is_checkpoint()
  {
    return my_is_checkpoint;
  }
  inline bool is_backup()
  {
    if (this->my_Tbackup) {
      return my_is_backup;
    } else {
      return false;
    }
  }
  inline bool is_line_min()
  {
    return line_min;
  }
  inline int iteration()
  {
    return iteration_counter;
  }

  inline void set_Ncheckpoint(int N)
  {
    this->my_Ncheckpoint = N;
  }
  inline void set_Niterations(int N)
  {
    this->my_max_steps = N;
  }
  inline void set_Nsystems(int N)
  {
    this->my_Nsystems = N;
  }
  inline void set_Nbackup(int N)
  {
    this->my_Nbackup = N;
  }
  inline void set_Tbackup(bool T)
  {
    this->my_Tbackup = T;
  }

  void set_training_algorithm(string algorithm);
  inline void set_threshold(REAL threshold)
  {
    my_threshold = threshold;
  }
  inline void set_debug(REAL new_debug)
  {
    debug = new_debug;
  }

  inline REAL alpha()
  {
    return my_alpha;
  }
  inline REAL dphi_dt()
  {
    return my_dphi_dt;
  }

  inline REAL get_error()
  {
    return Errors.back();
  }
  //routines for early stopping/val SEE
  inline void set_early_stop(bool early_stop)
  {
    this->early_stop = early_stop;
  }
  inline void set_val_sse(REAL SEE, int nval)
  {
    this->SSE_val = SEE;
    this->nval = nval;
  }
private:



  // Member objects
  vector<vector<Network_node*> > *nodes;
  Parallel *mpi;
  Functional_params F_params;

  //functions for early stopping
  bool early_stop;
  REAL SSE_val;
  int nval;


  REAL line_min_threshold;

  bool do_print;
  bool my_Tbackup;
  int my_Nbackup;

  // Member variables
  int iteration_counter;
  bool my_is_converged, my_is_checkpoint, my_is_backup;
  REAL lambda0, lambda;
  bool is_initialized;
  bool allow_training_switch;
  REAL debug;
  int my_Nparameters;
  int my_Nsystems;

  REAL my_alpha;
  // User specified parameters
  REAL my_threshold;
  int my_Ncheckpoint;
  int my_max_steps;
  bool global_optimizer;
  int t0;
  REAL current_alpha_max;
  REAL E_min, dalpha_dt;
  REAL alpha_max;
  REAL phase_shift;
  REAL my_dphi_dt;
  REAL E_best,g_best;
  int cycles, sign_prop;

  vector<REAL> x_best;

  // general storage
  deque<vector<REAL> > x;
  deque<vector<REAL> > gradient;
  deque<REAL> Errors;
  vector<REAL> search_direction;

  // Storage for RPROP
  vector<REAL> rprop_weights;
  vector<REAL> rprop_step;
  REAL rprop_weights_min;
  REAL rprop_weights_max;


  // Storage for L-BFGS
  int BFGS_Nmax;
  REAL convergence;
  vector<REAL> BFGS_last_x;
  vector<REAL> BFGS_last_gradient;
  vector<vector<REAL> > s;
  vector<vector<REAL> > y;
  deque<REAL> rho;
  bool line_min;


  // Storage for steepest descent
  REAL last_lambda;

  // Storage for Brent's method of line minimization
  deque<REAL> line_g;
  deque<REAL> line_alpha;
  REAL line_last_E;
  REAL line_last_g, line_last_alpha;
  REAL brent_c, brent_g_c, brent_d, brent_g_d;
  bool brent_flag;


  // Member functions
  vector<REAL> get_parameters();
  void set_parameters(vector<REAL> new_parameters);

  inline  vector<REAL> switch_to_rprop()
  {
    search_direction.clear();
    cout << "----- Switching to rprop -----"<<endl;
    this->line_min = true;
    for (int i=0; i<rprop_weights.size(); i++) {
      rprop_weights.at(i) = this->lambda;
    }
    this->get_next_x = &Optimizer::rprop;
    return rprop();
  }

  inline vector<REAL> switch_to_sd()
  {
    cout << "----- Switching to SD -----" << endl;
    this->line_min = true;
    search_direction.clear();
    this->get_next_x = &Optimizer::steepest_descent;
    return steepest_descent();
  }

  inline vector<REAL> switch_to_BFGS()
  {
    cout << "----- Switching to BFGS -----" << endl;
    this->search_direction.clear();
    this->get_next_x = &Optimizer::LBFGS;
    return LBFGS();
  }



  REAL alpha1, alpha2;


  //  vector<REAL> delta_x();

  // Pointer to current optimization algorithm
  vector<REAL> (Optimizer::*get_next_x)(void);

  // Optimization algorithms (get_next_x points to one of these)
  vector<REAL> steepest_descent();
  vector<REAL> steepest_descent_with_line_opt();
  vector<REAL> steepest_descent_with_momentum();
  vector<REAL> rprop();
  vector<REAL> LBFGS();
  REAL brent_line_min(REAL tolerance);

  // debugging
  vector<REAL> plot();
  vector<REAL> endpoints;

  // Helper funtions
  inline REAL dot_product(vector<REAL> v1, vector<REAL> v2)
  {
    int N=v1.size();
    if (N != v2.size()) {
      ERROR("Cannot take dot product of unequal length vectors");
    }
    REAL result = 0.0;
    for (int i=0; i<N; i++) {
      result += v1[i]*v2[i];
    }
    return result;
  }

  inline vector<REAL> vector_diff(vector<REAL> &v1, vector<REAL> &v2)
  {
    vector<REAL> output(v1.size(),0.0);
    for (int i=0; i<v1.size(); i++) {
      output.at(i) = v1.at(i) - v2.at(i);
    }
    return output;
  }

  inline vector<REAL> normalized(vector<REAL> &v)
  {
    vector<REAL> output = v;
    REAL norm = sqrt(this->dot_product(v,v));
    for (int i=0; i<v.size(); i++) {
      output[i] /= norm;
    }
    return output;
  }

  template <typename type> inline int sign(type value)
  {
    return (type(0) < value)-(value < type(0));
  }


};


#endif
