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



#include "Optimizer.h"
#include <iomanip>
#include <limits>
#include <cstdio>





// ########################################################
//                       Constructor
// ########################################################
//

Optimizer::Optimizer(Parallel *mpi_ptr) : mpi(mpi_ptr) { 
  
  this->nodes = NULL;
  this->is_initialized = false;
  this->get_next_x = &Optimizer::rprop;
  this->my_Nsystems = 1;
  this->do_print = true;

}

// ########################################################
// ########################################################




// ########################################################
//                       Constructor
// ########################################################
//

Optimizer::Optimizer(Parallel *mpi_ptr, Functional_params params, int Nsystems) : mpi(mpi_ptr),F_params(params) {
  
  this->nodes = NULL;
  this->is_initialized = false;
  this->get_next_x = &Optimizer::rprop;
  this->my_Nsystems = Nsystems;
  this->do_print = true;

  this->set_Nsystems(Nsystems);
  this->set_Niterations(params.Niterations());
  this->set_Ncheckpoint(params.Ncheckpoint());
  this->set_training_algorithm(params.training_algorithm());
  this->set_threshold(params.threshold());
  this->set_debug(params.debug());

  
}

// ########################################################
// ########################################################




// ########################################################
//                       Desctructor
// ########################################################
//

Optimizer::~Optimizer() {  

}

// ########################################################
// ########################################################





// ########################################################
//                       SET_PARAMS
// ########################################################
// Sets internal parameter set.

void Optimizer::set_params(Functional_params params) { 
  F_params = params; 
}

// ########################################################
// ########################################################




// ########################################################
//                       SET_TRAINING_ALGORITHM
// ########################################################
// Selects training algorithm.

void Optimizer::set_training_algorithm(string algorithm) {

  if (algorithm == "default") {
    
    this->get_next_x = &Optimizer::rprop;
    this->allow_training_switch = true;
    
  } else if (algorithm == "rprop") {
    
    this->get_next_x = &Optimizer::rprop;
    this->allow_training_switch = false;
    
  } else if (algorithm == "steepest_descent") {
    
    this->get_next_x = &Optimizer::steepest_descent;
    this->allow_training_switch = false;
    
  } else if (algorithm == "steepest_descent_with_line_opt") {

    this->get_next_x = &Optimizer::steepest_descent_with_line_opt;
    this->allow_training_switch = false;

  } else if (algorithm == "bfgs") {
    
    this->get_next_x = &Optimizer::LBFGS;
    this->allow_training_switch = false;

  } else {

    if (mpi->io_node()) { ERROR("Training algorithm \""+algorithm+"\" is not implemented");}

  }

}

// ########################################################
// ########################################################




// ########################################################
//                       INIT
// ########################################################
// Initializes optimizer

void Optimizer::init(vector<Network_node*> in_nodes) {
  
  if (!is_initialized) {
    nodes = new vector<vector<Network_node*> >(1, in_nodes);
    is_initialized = true;
    return;
  }
  
  nodes->push_back(in_nodes);
  
}

// ########################################################
// ########################################################



// ########################################################
//                       INIT
// ########################################################
// Initializes optimizer

void Optimizer::init(vector<vector<Network_node*> > *in_nodes) {

  if (is_initialized) { 
    this->my_is_checkpoint = false;
    return; 
  }

  is_initialized = true;
  
  this->nodes = in_nodes;
  this->internal_init();

}

// ########################################################
// ########################################################



  
// ########################################################
//                       INTERNAL_INIT
// ########################################################
// Set up optimizer defaults

void Optimizer::internal_init() {
  
  if (mpi->io_node()) {
    cout << "Iteration    \u0190 = \u03A3(E\u00b2)       |\u2207\u0190|"<<endl;
    cout << "-------------------------------------" << endl;
  }
  
  this->iteration_counter = 1;
  this->my_is_converged = false;
  this->lambda0 = 1e-4;
  this->line_min = true;
  this->lambda = this->lambda0;
  this->my_is_checkpoint = false;
  this->my_Nparameters = 0;

  this->alpha1 = lambda0;
  this->my_alpha = 0.0;
  this->global_optimizer = true;
  this->alpha_max = F_params.lambda_max();
  this->current_alpha_max = alpha_max;
  this->phase_shift = 0.0;
  this->dalpha_dt = 0;
  this->cycles = 0;
  this->E_best = 9e9;

  x.push_back(this->get_parameters());  
  my_Nparameters = x.at(0).size();
  this->set_parameters(mpi->Bcast(x.at(0),my_Nparameters));


  
  BFGS_Nmax = 77;
  convergence = 100*std::numeric_limits<REAL>::epsilon();
  
  for (int i=0; i<my_Nparameters; i++) {
    rprop_weights.push_back(lambda0);
    rprop_step.push_back(0);
  }
  rprop_weights_max = 0.5;
  rprop_weights_min = 1000.0*std::numeric_limits<REAL>::epsilon();
  
}

// ########################################################
// ########################################################




// ########################################################
//                       SET_PARAMETERS
// ########################################################
// Loops over nodes and sets their parameters. i.e. updates
// network for next training step.

void Optimizer::set_parameters(vector<REAL> new_parameters) {

  vector<REAL>::iterator it = new_parameters.begin();
  for (int layer=0; layer<nodes->size(); layer++) {
    for (int node=0; node<nodes->at(layer).size(); node++) {
      nodes->at(layer).at(node)->set_parameters(vector<REAL>(it,it+nodes->at(layer).at(node)->Nparameters()));
      it+=nodes->at(layer).at(node)->Nparameters();
    }
  }
  if (it != new_parameters.end()) {ERROR("Node parameters not set properly");}
}

// ########################################################
// ########################################################




// ########################################################
//                       GET_PARAMETERS
// ########################################################
// Loops over nodes and gets their current parameters.

vector<REAL> Optimizer::get_parameters() {

  vector<REAL> params;
  if (my_Nparameters) {params.reserve(my_Nparameters);}
  
  for (int layer=0; layer<nodes->size(); layer++) {
    for (int node=0; node<nodes->at(layer).size(); node++) {
      
      vector<REAL> temp = nodes->at(layer).at(node)->get_parameters();
      params.insert(params.end(), temp.begin(), temp.end());

    }
  }

  return params;

}

// ########################################################
// ########################################################




//// ########################################################
////                       
//// ########################################################
////
//
//vector<REAL> Optimizer::plot() {
//
//  static REAL delta = 0.0;
//  static int start_it = 0;
//  line_min = false;
//  int N = 50;
//  
//  if (delta==0) {
//    start_it = iteration_counter;
//    delta = (endpoints.back()-endpoints.front())/(REAL)(N-1);
//    cout <<setprecision(10)<<endl<<endl<<"--------------"<<endl;
//    vector<REAL> x1 = BFGS_last_x;
//    for (int i=0; i<x1.size(); i++) {
//      x1.at(i) -= endpoints.front()*search_direction.at(i);
//    }
//    return x1;
//  }
//  
//  vector<REAL> x1 = x.back();
//  
//  cout << setprecision(20)<<++iteration_counter-start_it<<"   "<<Errors.back()<<"    "<< this->dot_product(search_direction,gradient.back())<<endl;
//  
//  for (int i=0; i<x1.size(); i++) {
//    x1.at(i) -= delta*search_direction.at(i);
//  }
//  
//  if (iteration_counter == start_it+N) {
//    cout <<endl<<"--------------------"<<endl<<endl;
//  }
//  return x1;
//
//}

      
// ########################################################
//                       STEEPEST_DESCENT
// ########################################################
// Implementation of a steepest descent algorithm.

vector<REAL> Optimizer::steepest_descent() {
  
  vector<REAL> x1 = x.back();
  
  search_direction = this->normalized(gradient.back());
  
  int N = Errors.size();
  if (N==1) { 
    alpha1 = 5e-5;
    for (int i=0; i<x1.size(); i++) {
      x1.at(i) -= alpha1*search_direction.at(i);
    }
    return x1;
  }  else if (Errors.at(N-1) <= Errors.at(N-2)) {
    alpha1 *= 1.4;
  } else {
    alpha1 *= 0.7;
  }
  
  if (alpha1 > rprop_weights_max) {
    alpha1 = rprop_weights_max;
  } else if (alpha1 < rprop_weights_min) {
    alpha1 = rprop_weights_min;
  }
  

  if (allow_training_switch && !(iteration_counter%50)) {
    REAL dg = 0.0;
    REAL ave = sqrt(dot_product(gradient.at(0),gradient.at(0)));
    for (int i=1; i<gradient.size(); i++) {
      ave += sqrt(dot_product(gradient.at(i),gradient.at(i)));
      dg += sqrt(dot_product(gradient.at(i),gradient.at(i))) -
	sqrt(dot_product(gradient.at(i-1),gradient.at(i-1)));
    }
    dg /= ave;
    if (dg > -this->debug && dg < 0) {
      return this->switch_to_BFGS();
    } else if (dg > this->debug) {
      return this->switch_to_rprop();
    }
  }
  
  for (int i=0; i<x1.size(); i++) {
    x1.at(i) = x1.at(i) - alpha1*(search_direction.at(i) + F_params.sd_momentum()*(x.at(N-1).at(i)-x.at(N-2).at(i)));
  }
  
  return x1;
  
}

// ########################################################
// ########################################################




// ########################################################
//               STEEPEST_DESENT_WITH_LINE_OPT
// ########################################################
// Implements a steepest descent algorithm with full
// line optimization.  Not recommended.

vector<REAL> Optimizer::steepest_descent_with_line_opt() {

  if (this->search_direction.size() == 0) {
    vector<REAL> x1 = x.back();
    search_direction = this->normalized(gradient.back());
    BFGS_last_x = x.back();
    BFGS_last_gradient = this->gradient.back();

    REAL delta = brent_line_min(F_params.line_min_epsilon());

    for (int i=0; i<gradient.back().size(); i++) {
      x1.at(i) -= delta*search_direction.at(i);
    }

    return x1;
  }

  if (!line_min) {
    
    vector<REAL> x1 = BFGS_last_x;
    REAL delta = brent_line_min(F_params.line_min_epsilon());
        
    for (int i=0; i<search_direction.size(); i++) {
      x1.at(i) -= delta*search_direction.at(i);
    }
    return x1;
    
  } else {
    
    vector<REAL> x1 = x.back();
    search_direction = this->normalized(gradient.back());
    BFGS_last_x = x.back();
    BFGS_last_gradient = this->gradient.back();

    REAL delta = brent_line_min(F_params.line_min_epsilon());

    for (int i=0; i<gradient.back().size(); i++) {
      x1.at(i) -= delta*search_direction.at(i);
    }

    return x1;
  }
  
}

// ########################################################
// ########################################################


    
  
   
// ########################################################
//                       RPROP
// ########################################################
// Implements the resilient back-propagation (RProp) algorithm.

vector<REAL> Optimizer::rprop() {
  
  this->line_min = true;

  int N = gradient.size();
  vector<REAL> x1 = x.back();
  
  if (allow_training_switch && !(iteration_counter % 67)) {
    if (gradient.size() > 1 && this->dot_product(gradient.back(), gradient.back())<1e11 && search_direction.size() >0) {
    for (int i=0; i<gradient.back().size(); i++) {
      if ( (gradient.back().at(i) - gradient.at(gradient.size()-2).at(i))/rprop_step.at(i) < 0) {
	break;
      }
      if (i == gradient.back().size() - 1) { 
	for (int i=0; i<gradient.back().size(); i++) {
	  if (rprop_weights.at(i) < this->lambda) { this->lambda = rprop_weights.at(i); }
	}
	this->alpha1 = lambda0;
	return this->switch_to_sd();
      }
    }
    
    }
  }
  
  if (search_direction.size() == 0) {
    search_direction.assign(gradient.at(0).size(),1);
    for (int i=0; i<gradient.at(0).size(); i++) {
      rprop_step.at(i) = -rprop_weights.at(i)*sign(search_direction.at(i));
      x1.at(i) = x1.at(i) + rprop_step.at(i); 
    }
  } else if (N > 1) {
    for (int i=0; i<gradient.at(0).size(); i++) {
      int sign_change = search_direction.at(i)*sign(gradient.back().at(i)*gradient.at(N-2).at(i));
      
      if (sign_change < 0) {
	search_direction.at(i) = 0;
      } else if (sign_change > 0) {
	rprop_weights.at(i) *= 1.2;
	if (rprop_weights.at(i) > rprop_weights_max && i != gradient[0].size()-1) { rprop_weights.at(i) = rprop_weights_max; }
	rprop_step.at(i) = -sign(gradient.back().at(i))*rprop_weights.at(i);
	x1.at(i) = x1.at(i) + rprop_step.at(i);
      } else if (sign_change == 0 && N >= 3) {
	REAL D2 = (gradient.at(N-1).at(i) - gradient.at(N-3).at(i))/rprop_step.at(i);
	if (D2 > 0) {
	  rprop_weights.at(i) = abs(gradient.back().at(i))/D2;
	} else {
	  
	}
	if (rprop_weights.at(i) < rprop_weights_min) { rprop_weights.at(i) = rprop_weights_min; }
	if (rprop_weights.at(i) > rprop_weights_max && i != gradient[0].size()-1) { rprop_weights.at(i) = rprop_weights_max; }
	rprop_step.at(i) = -sign(gradient.back().at(i))*rprop_weights.at(i);
	x1.at(i) = x1.at(i) + rprop_step.at(i);
	search_direction.at(i) = 1;
      }
    }
  } else {
    ERROR("RPROP called with no gradient");
  }
  
  
  return x1;

}

// ########################################################
// ########################################################





// ########################################################
//                       LBFGS
// ########################################################
// Implemnts the lesser memory BFGS algorithm.

vector<REAL> Optimizer::LBFGS() { 
  
  if (this->search_direction.size() == 0) {
    vector<REAL> x1 = x.back();
    search_direction = this->normalized(gradient.back());
    BFGS_last_x = x.back();
    BFGS_last_gradient = this->gradient.back();
    
    REAL delta = brent_line_min(F_params.line_min_epsilon());
    
    for (int i=0; i<gradient.back().size(); i++) {
      x1.at(i) -= delta*search_direction.at(i);
    }
    
    return x1;
  }
  
  if (!line_min) {
    
    vector<REAL> x1 = BFGS_last_x;
    REAL delta = brent_line_min(F_params.line_min_epsilon());
    
    for (int i=0; i<search_direction.size(); i++) {
      x1.at(i) -= delta*search_direction.at(i);
    }
    return x1;

  } else {
    
    vector<REAL> alpha(s.size()+1,0.0);
    REAL beta;
    
    search_direction = gradient.back();
    
    
    if (!(iteration_counter % (8*BFGS_Nmax))) {
      s.clear();
      y.clear();
      rho.clear();
    }
    
    s.push_back(vector<REAL>(search_direction.size(),0.0));
    y.push_back(vector<REAL>(search_direction.size(),0.0));
    for (int i=0; i<search_direction.size(); i++) {
      s.back().at(i) = x.back().at(i) - BFGS_last_x.at(i);
      y.back().at(i) = gradient.back().at(i) - BFGS_last_gradient.at(i);
    }

    rho.push_back(1/this->dot_product(y.back(),s.back()));
    
    if (s.size() > BFGS_Nmax) {
      s.erase(s.begin());
      y.erase(y.begin());
      rho.pop_front();
    }
    
    for (int i=s.size()-1; i>=0; i--) {
      alpha.at(i) = rho.at(i)*this->dot_product(s.at(i),search_direction);
      for (int j=0; j<search_direction.size(); j++) {
	search_direction[j] -= alpha[i]*y[i][j];
      }
    }
    
    REAL H = this->dot_product(y.back(),s.back())/this->dot_product(y.back(),y.back());
    for (int j=0; j<search_direction.size(); j++) {
      search_direction.at(j) *= H;
    }
    
    for (int i=0; i<s.size(); i++) {
      beta = rho.at(i)*this->dot_product(y.at(i),search_direction);
      for (int j=0; j<search_direction.size(); j++) {
	search_direction[j] += s[i][j]*(alpha[i] - beta);
      }
    }

    search_direction = this->normalized(search_direction);
    BFGS_last_x = x.back();
    BFGS_last_gradient = gradient.back();

    REAL delta = brent_line_min(F_params.line_min_epsilon());
    vector<REAL> x1 = BFGS_last_x;
    
    for (int i=0; i<search_direction.size(); i++) {
      x1.at(i) -= delta*search_direction.at(i);
    }
    return x1;
  }    

}

// ########################################################
// ########################################################




// ########################################################
//                       BRENT_LINE_MIN
// ########################################################
// Implements Brent's method for line minimization.

REAL Optimizer::brent_line_min(REAL tolerance) {

  REAL g = this->dot_product(normalized(search_direction),gradient.back());

  if (line_min == true || line_g.empty() ) {
    this->line_min = false;
    
    line_g.assign(1,g);
    line_alpha.assign(1,0.0);
    
    this->alpha1 = lambda0;
    this->line_last_E = Errors.back();
    line_last_alpha = alpha1;
    line_last_g = g;

    return alpha1;
    
  }
  
  if (abs(g) < tolerance || (line_g.size()==2 && abs(line_alpha.front()-line_alpha.back())<20*std::numeric_limits<REAL>::epsilon())) {

    this->line_min = true;
    
    if (abs(line_alpha.front()-line_alpha.back())<20*std::numeric_limits<REAL>::epsilon() && abs(g)>tolerance) {
      warning("In BFGS: Line optimization could not reach the requested precision");
    } 
    
    return alpha1;
    
  }
  
  if (line_g.size() == 1) {
    
    int g_sign = this->sign(g*line_g.back());
    
    if (g_sign <= 0) {
      if (abs(g) <= abs(line_g.back())) {
        line_g.push_back(g);
        line_alpha.push_back(alpha1);
      } else {
        line_g.push_front(g);
        line_alpha.push_front(alpha1);
      }
      
      this->brent_c = line_alpha.front();
      this->brent_g_c = line_g.front();
      this->brent_flag = true;
      
    } else {
      
      if (Errors.back() <= line_last_E) {

	line_alpha.push_back(alpha1);
	line_g.push_back(g);
	
	line_last_g = g;
	line_last_alpha = alpha1;
	line_last_E = Errors.back();
	
	if (abs(g) >= abs(line_g.front())) {
	  alpha1 *= 2.0;
	} else {
	  alpha1 = alpha1 - (g)*(line_alpha.back()-line_alpha.front())/(line_g.back()-line_g.front());
	}
	
	line_alpha.pop_front();
	line_g.pop_front();
	
      } else {

	if (alpha1 < 10*lambda0) {
	  line_min = true;
	  
	  s.clear();
	  y.clear();
	  rho.clear();
	  search_direction = this->normalized(gradient.back());
	  BFGS_last_x = x.back();
	  BFGS_last_gradient = this->gradient.back();
	  
	} else {
	  
	  alpha1 = line_last_alpha + sign(line_last_alpha)*lambda0;
	  
	  line_alpha.front() = line_last_alpha;
	  line_g.front() = line_last_g;
	} 
	
      }
      
      return alpha1;
      
    }  
    
  } else if (line_g.size() == 2) {
    int g_sign = this->sign(g*line_g.front());
    if (g_sign < 0) { 
      line_g.pop_back();
      line_g.push_back(g);
      line_alpha.pop_back();
      line_alpha.push_back(alpha1);
    } else {
      line_g.pop_front();
      line_g.push_front(g);
      line_alpha.pop_front();
      line_alpha.push_front(alpha1);
    }
    
    if (abs(line_g.front()) < abs(line_g.back())) {
      line_g.push_back(line_g.front());
      line_g.pop_front();
      line_alpha.push_back(line_alpha.front());
      line_alpha.pop_front();
    }
    
  } else { ERROR("Something went wrong");}

  if (line_g.size() >2 ||line_alpha.size()>2) {ERROR("Incorrect BFGS");}

  REAL new_alpha;
  if (brent_g_c != line_g.front() && brent_g_c != line_g.back()) {
    new_alpha = line_alpha.front()*line_g.back()*brent_g_c
      /((line_g.front()-line_g.back())*(line_g.front()-brent_g_c))
      + line_alpha.back()*line_g.front()*brent_g_c
      /((line_g.back()-line_g.front())*(line_g.back()-brent_g_c))
      + brent_c*line_g.front()*line_g.back()
      /((brent_g_c-line_g.front())*(brent_g_c-line_g.back()));
  } else {
    new_alpha = (line_alpha.front()*line_g.back() - line_alpha.back()*line_g.front())
      / (line_g.back() - line_g.front());
  }
  REAL a = abs(new_alpha - line_alpha.back());
  REAL b = abs(line_alpha.back() - brent_c);
  REAL c = abs(brent_c - brent_d);
  if ( (4*new_alpha > (3*line_alpha.front()+line_alpha.back())&&(new_alpha > line_alpha.back()))
	 || (4*new_alpha < (3*line_alpha.front()+line_alpha.back())&&(new_alpha < line_alpha.back()))
       || (brent_flag && (2*a >= b))
       || (!brent_flag && (2*a >= c))
       || (brent_flag && (b < convergence))
       || (!brent_flag && (c < convergence) )) {
    
    new_alpha = 0.5*(line_alpha.front() + line_alpha.back());
    brent_flag = true;

  } else {
    brent_flag = false;
  }
  
  brent_d = brent_c;
  brent_c = line_alpha.back();
  brent_g_c = line_g.back();

  this->alpha1 = new_alpha;
  return new_alpha;

}

// ########################################################
// ########################################################





// ########################################################
//                       UPDATE_NETWORK
// ########################################################
// Takes the current network state, figures out the next
// trial step and performs some bookkeeping.

void Optimizer::update_network(REAL Error, vector<REAL> grad, REAL SSE) {
  
  Errors.push_back(mpi->Reduce(Error, MPI_SUM)/my_Nsystems);
  gradient.push_back(mpi->Reduce(grad, MPI_SUM));
  
  if (SSE > 0) {
    SSE = mpi->Reduce(SSE,MPI_SUM)/my_Nsystems;
  } else {
    SSE = Errors.back();
  }
  
  
  if (F_params.regularization()) {
    deque<vector<REAL> >::iterator grad = --gradient.end();
    deque<vector<REAL> >::iterator xx = --x.end();
    for (int i=0; i<my_Nparameters-1; i++) {
      grad->at(i) += 2*F_params.regularization()*xx->at(i)/my_Nparameters;
      Errors.back() += F_params.regularization()*xx->at(i)*xx->at(i)/my_Nparameters/my_Nsystems;
    }
    
  }
  
  REAL norm = sqrt(this->dot_product(gradient.back(),gradient.back()))/my_Nparameters/my_Nsystems;
  
  if (mpi->io_node()) {
    this->my_is_converged = false;
    if (iteration_counter >= my_max_steps) {
      this->my_is_converged = true;
      line_min = true;
    } else if (norm < this->my_threshold) {
      if (alpha_max != 0) {
	if (my_alpha == 0) {
	  if (global_optimizer) {
	    if (SSE < E_best) {
	      E_best = SSE;
	      g_best = norm;
	      x_best = x.back();
	    }
	    cout << "----- converged; turning on perturbation -----" << endl;

	    E_min = SSE;
	    dalpha_dt = 3.14159265358979e-04;
	    t0 = iteration_counter;
	  } else {
	    this->my_is_converged = true;
	    line_min = true;
	    if (SSE > E_best) {
	      cout << "===== resetting to previous best =====" <<endl;
	      SSE = E_best;
	      norm = g_best;
	      set_parameters(x_best);
	    } 
	  }
	} 
      } else {
	this->my_is_converged = true;
	this->line_min = true;
      }
    }
  }
  
  my_is_converged = mpi->Bcast(my_is_converged);
  
  
  vector<REAL> temp_x;

  if (!my_is_converged) {
    
    if (mpi->rank() == 0) {
      temp_x = (this->*get_next_x)();
    }
    
    x.push_back(mpi->Bcast(temp_x,gradient.back().size()));
    set_parameters(x.back());
    
  }
  
  if (x.size() > 20) {
    x.pop_front();
    gradient.pop_front();
    Errors.pop_front();
  }
  
  
  if (do_print) {
  if (mpi->io_node()) {
    if (this->line_min) {

      if ( !(iteration_counter % F_params.Nprint()) || my_is_converged){
#ifdef USE_LONG_DOUBLE
	printf("%-9d %#13.6Le %#13.6Le\n",iteration_counter++,sqrt(SSE),norm);
#else
	printf("%-9d %#13.6e %#13.6e\n",iteration_counter++,sqrt(SSE),norm);
#endif
      } else {
	iteration_counter++;
      }
      my_is_checkpoint = ((iteration_counter-1) % my_Ncheckpoint ? false : true);
      if (iteration_counter == 1) { my_is_checkpoint = false; }
    }
    
  }
  
  my_is_checkpoint = mpi->Bcast(my_is_checkpoint);
  
  if (dalpha_dt) {
    
    if (SSE < 0.985*E_min && phase_shift == 0) {
      cout << "----- transition detected -----" << endl;
      current_alpha_max = my_alpha;
      t0 = iteration_counter;
      cycles = 0;
      phase_shift = 1.57079632679490;
      dalpha_dt *= 2;
    } 

    my_alpha = current_alpha_max*sin(dalpha_dt*(iteration_counter - t0) + phase_shift);
    my_dphi_dt = 0.5*RAND::Normal(1,0.15)*0.0314159/5;
    
    if (my_alpha <= 0) {
      if (current_alpha_max == alpha_max) {
	cout << "----- No transition found; completing optimization -----" << endl;
	cycles++;
	if (cycles == F_params.Ncycles()) {
	  global_optimizer = false;
	}
	my_alpha = 0;
	dalpha_dt = 0;
      } else {
	cout << "----- continuing optimization -----"<<endl;
	phase_shift = 0;
	my_alpha = 0;
	current_alpha_max = alpha_max;
	dalpha_dt = 0;
      }
    }
    
  }
  
  my_alpha = mpi->Bcast(my_alpha);
  my_dphi_dt = mpi->Bcast(my_dphi_dt);
  } else {
 
    if (this->line_min) {
      if (mpi->io_node()) {
	
	iteration_counter++;
	
	my_is_checkpoint = ((iteration_counter-1) % my_Ncheckpoint ? false : true);
	if (iteration_counter == 1) { my_is_checkpoint = false; }

      }
      
      my_is_checkpoint = mpi->Bcast(my_is_checkpoint);

    }

  }

}

// ########################################################
// ########################################################


