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
// This class holds information about the structure of a given
// system. It is required whenever the user requests mapping the
// structure to some other variable.
// ####################################################################




#ifndef __Structure
#define __Structure

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>

#include "Common.h"
#include "Grid_data.h"
#include "Functional_params.h"
#include "Atom.h"
#include "Tables.h"

using namespace std;

class Structure
{
public:
  Structure();
  ~Structure();
  void print();
  void Calc_G(int atom);

  void Allocate_G(int Gtype, vector <REAL>& row);
  void Set_Neighbor(REAL Rcut);
  void Normalize();
  void Set_Convert();
  REAL train_Local(Functional_params *F,REAL totalenergy);
  REAL unravel_Energy(REAL localenergy);
  int NG()
  {
    return my_NG;
  }

  void Get_Forces(const vector<vector<REAL> > &dE_dG, REAL **f);
  bool periodic;
  int Natom, inum;
  map<string,REAL> FE;
  inline void print_FE()
  {
    map<string,REAL>::iterator it;
    for (it = FE.begin(); it != FE.end(); it++) {
      cout << it->first << " : " << it->second << endl << endl;
    }
  }
  map<int,string> lammps_conv; //converts lammps type to FE energy symbol
  vector<vector<REAL>*> init_G(Functional_params *F);
  vector < vector<REAL> > pos, cpos;
  vector < vector <int> > firstneigh;
  vector <int> ilist, numneigh;
  vector <REAL> a, b, c;
  vector <string> specie;
  vector<Atom> types;
  vector < vector <REAL> > G1p, G2p, G3p, G4p;
  vector <vector<REAL> > G;
  vector <vector <REAL> > Forces;
  vector<vector<REAL>* > G_ptr;

  inline bool is_initialized()
  {
    return my_is_initialized;
  }
  void set_structure(double** x, int* in_types, int Natoms, int Nghosts);


  vector<REAL> mean(int atomic_number);
  int count(int atomic_number);
  vector<REAL> variance(int atomic_number, const vector<REAL> &means);



  bool CART;
  bool NORM;

  string train;

protected:
  int my_NG;
  int my_NG_types;

  map <char,REAL> norm;
  REAL nx, ny, nz, lx, ly, lz;
  vector < vector<int> > neighborlist;
  vector <REAL> na, nb, nc;
  bool my_is_initialized;
  vector < vector <REAL> > AtoF, FtoA;
  vector<vector<int> > translations;
  vector<vector<vector<int> > > trans_indices;
  bool is_MD;

  map<int,int> atom_index;
  vector<vector<int> > atom_matrix;
  int Npair_Gs;

  map<int,vector<REAL> > my_sums;
  map<int,int> my_counts;

  REAL Rcut;

  static map<REAL,REAL> prefactor_A;

  REAL old_theta;
  REAL prefactor_B;

  static int call_count;

  bool recalc_exp1, recalc_exp2, recalc_fc1, recalc_fc2;
  REAL G_prefactor, G_angular, G_exp1, G_exp2, G_fc1, G_fc2;
  REAL old_delta, old_term;
  REAL old_R, old_Ru, old_row0, old_row1, old_row2;
  REAL old_R3, old_Ru3, old_row03, old_row13, old_row23;

  static Tables my_cos;
  static Tables my_pow;

  inline REAL fc(REAL frac, int fcswitch)
  {
    if(frac >= 1.0) {
      return 0.0;
    } else {

      if(fcswitch == 0) {
        REAL x = M_PI*frac;
        if (x < 1.047197551196598) {
          return (0.043539571484819*x - 0.273567195834312)*x*x + 1.0;
        } else if (x < 2.094395102393195) {
          x -= 1.047197551196598;
          return ((0.087079142969639*x - 0.136783597917156)*x - 0.429718346348117)*x + 0.75;
        } else {
          x -= 2.094395102393195;
          return ((0.043539571484819*x + 0.136783597917156)*x - 0.429718346348117)*x + 0.25;
        }
      } else {
        return pow(tanh(1-frac),3);
      }
    }
  }

  inline vector < vector <REAL> > invert(vector < vector <REAL> > cell)
  {
    REAL det = cell[0][ 0] * (cell[1][ 1] * cell[2][ 2] - cell[2][ 1] * cell[1][ 2]) -
               cell[0][ 1] * (cell[1][ 0] * cell[2][ 2] - cell[1][ 2] * cell[2][ 0]) +
               cell[0][ 2] * (cell[1][ 0] * cell[2][ 1] - cell[1][ 1] * cell[2][ 0]);
    if (det == 0 ) {
      ERROR("Non-invertible unit cell");
    }
    REAL invdet = 1 / det;
    vector <REAL> tmp(3,0.0);
    vector < vector <REAL> > ident;
    for(int i =0; i < 3; i++) {
      ident.push_back(tmp);
      ident[i][i] = 1.0;
    }
    ident[0][ 0] = (cell[1][ 1] * cell[2][ 2] - cell[2][ 1] * cell[1][ 2]) * invdet;
    ident[0][ 1] = (cell[0][ 2] * cell[2][ 1] - cell[0][ 1] * cell[2][ 2]) * invdet;
    ident[0][ 2] = (cell[0][ 1] * cell[1][ 2] - cell[0][ 2] * cell[1][ 1]) * invdet;
    ident[1][ 0] = (cell[1][ 2] * cell[2][ 0] - cell[1][ 0] * cell[2][ 2]) * invdet;
    ident[1][ 1] = (cell[0][ 0] * cell[2][ 2] - cell[0][ 2] * cell[2][ 0]) * invdet;
    ident[1][ 2] = (cell[1][ 0] * cell[0][ 2] - cell[0][ 0] * cell[1][ 2]) * invdet;
    ident[2][ 0] = (cell[1][ 0] * cell[2][ 1] - cell[2][ 0] * cell[1][ 1]) * invdet;
    ident[2][ 1] = (cell[2][ 0] * cell[0][ 1] - cell[0][ 0] * cell[2][ 1]) * invdet;
    ident[2][ 2] = (cell[0][ 0] * cell[1][ 1] - cell[1][ 0] * cell[0][ 1]) * invdet;
    return ident;
    /*
    REAL prefactor;
    vector <REAL> tmp(3,0.0);
    vector < vector <REAL> > ident;
    for(int i =0; i < 3; i++) {
        ident.push_back(tmp);
        ident[i][i] = 1.0;
    }

    for (int i = 0; i < cell.size(); i++) {
        prefactor = cell[i][i];
        for (int j = 0; j < cell[i].size(); j++) {
            cell[i][j] /= prefactor;
            ident[i][j] /= prefactor;
        }
        for(int j =0; j<3; j++) {
            if (j != i) {
                prefactor = cell[j][i];
                for (int z = 0; z < 3; z++ ) {
                    cell[j][z] = cell[j][z]  - cell[i][z] * prefactor;
                    ident[j][z] = ident[j][z] - cell[i][z] * prefactor;
                }
            }
        }
    }
     */
    return ident;
  }

  inline vector < vector <REAL> > transpose(vector < vector <REAL> > cell)
  {
    //int i,j,k;
    vector <REAL> tmp(3,0.0);
    vector < vector <REAL> > ident;
    for(int i =0; i < 3; i++) {
      ident.push_back(tmp);
    }
    for(int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        ident[j][i] = cell[i][j];
      }
    }
    return ident;
  }


  // This functions calculates (1+lambda*cos(theta))^xi.
  // The algorithm first checks if xi is an integer and,
  // if so, calls int_power.
  // Otherwise, it performs the Taylor series approximation
  // of (1+lambda*cos(theta))^xi. It uses the C++ "pow"
  // function as a fallback if needed.

  inline REAL angular_term(REAL theta, vector<REAL> &row)
  {
    theta *= row[4];
    if (theta <= -1) {
      return 0;
    }
    if (row[3] == (int)(row[3])) {
      return int_power(1+theta,(int)(row[3]));
    }
    //return pow(1+theta,row[3]);

    //optional taylor series approximation to pow(1+theta)^xi. Comment out return
    // statment on previous line to use this (EXPERIMENTAL)

    REAL new_xi = row[3];
    REAL alpha = new_xi*theta;
    int term = 1;
    REAL result = 1.0;
    while(abs(alpha)>1e-8) {
      result += alpha;
      alpha *= taylor[term]*theta*(new_xi-term);
      ++term;
    }
    if (term >= 250) {
      return pow(1+theta,row[3]);
    } else {
      return (result+alpha);
    }
  }

  long int Nterms;
  long int my_call_count;

  // Effectively pow for integer exponent. It uses
  // squaring by bit shifts. Much faster than
  // native pow for integer exponents.
  inline REAL int_power(REAL base, int exponent)
  {

    base = (exponent >= 0 ? base : 1.0/base);
    unsigned int power = abs(exponent);
    REAL result = 1;

    while(power) {
      if (power & 1) {
        result *= base;
      }
      power >>= 1;
      base *= base;
    }
    return result;
  }

  vector<REAL> exp_integer;
  vector<REAL> taylor;

  inline REAL exp_term(REAL x)
  {
    return exp(x);
  }


  inline void To_Cart(vector <REAL>& proj )
  {
      /*
    proj[0] = this->FtoA[0][0]*proj[0] + this->FtoA[0][1]*proj[1] + this->FtoA[0][2]*proj[2];
    proj[1] = this->FtoA[1][1]*proj[1] + this->FtoA[1][2]*proj[2];
    proj[2] = this->FtoA[2][2]*proj[2];
       */
      REAL t1,t2,t3;
    t1 = this->FtoA[0][0]*proj[0] + this->FtoA[0][1]*proj[1] + this->FtoA[0][2]*proj[2];
    t2 = this->FtoA[1][0]*proj[0] + this->FtoA[1][1]*proj[1] + this->FtoA[1][2]*proj[2];
    t3 = this->FtoA[2][0]*proj[0] + this->FtoA[2][1]*proj[1] + this->FtoA[2][2]*proj[2];
    proj[0] = t1;
    proj[1] = t2;
    proj[2] = t3;
  }

  inline void To_Frac(vector <REAL>& proj)
  {
    /*
    proj[0] = this->AtoF[0][0]*proj[0] + this->AtoF[0][1]*proj[1] + this->AtoF[0][2]*proj[2];
    proj[1] = this->AtoF[1][1]*proj[1]+ this->AtoF[1][2]*proj[2];
    proj[2] = this->AtoF[2][2]*proj[2];
     */
      REAL t1, t2, t3;
    t1 = this->AtoF[0][0]*proj[0] + this->AtoF[0][1]*proj[1] + this->AtoF[0][2]*proj[2];
    t2 = this->AtoF[1][0]*proj[0] + this->AtoF[1][1]*proj[1] + this->AtoF[1][2]*proj[2];
    t3 = this->AtoF[2][0]*proj[0] + this->AtoF[2][1]*proj[1] + this->AtoF[2][2]*proj[2];
    proj[0] = t1;
    proj[1] = t2;
    proj[2] = t3;
  }
  inline REAL G1(vector <REAL>& row,REAL R)
  {
    return fc(R/row[0],row[1]);
  }
  inline REAL G2(vector <REAL>& row,REAL R)
  {
    return exp_term(-row[2] *pow(R - row[3],2))*fc(R/row[0],row[1]);
    //return 1.0;
  }
  inline REAL G3(vector <REAL>& row,REAL R,REAL Ru, REAL Rjk, REAL theta)
  {
      //cout << 1 + (row[4] * theta) << endl;
      REAL term = 1 + (row[4] * theta);
      if (term <= 0) { return 0.0; }
    return Structure::prefactor_A.at(row[3])*pow(1+(row[4]*theta),row[3])*exp(-row[2]*(R*R+Ru*Ru+Rjk*Rjk))*fc(R/row[0],row[1])*fc(Ru/row[0],row[1])*fc(Rjk/row[0],row[1]);
  }
  inline REAL G4(vector <REAL>& row,REAL R,REAL Ru, REAL Rjk,REAL theta)
  {

    G_prefactor = Structure::prefactor_A.at(row[3]);

    REAL term = (1.0+row[4]*theta);
    if (term <= 0) {
      return 0;
    }

    G_angular = angular_term(theta, row);

    recalc_exp1 = recalc_exp2 = recalc_fc1 = recalc_fc2 = false;

    if (R != old_R) {
      recalc_exp1 = true;
      recalc_fc1 = true;
      old_R = R;
    }
    if (Ru != old_Ru) {
      recalc_exp2 = true;
      recalc_fc2 = true;
      old_Ru = Ru;
    }
    if (row[0] != old_row0 || row[1] != old_row1) {
      recalc_fc1 = true;
      recalc_fc2 = true;
      old_row0 = row[0];
      old_row1 = row[1];
    }
    if (row[2] != old_row2) {
      recalc_exp1 = true;
      recalc_exp2 = true;
      old_row2 = row[2];
    }

    if (recalc_exp1) {
      G_exp1 = exp_term(-row[2]*(R*R));
    }
    if (recalc_exp2) {
      G_exp2 = exp_term(-row[2]*(Ru*Ru));
    }
    if (recalc_fc1) {
      G_fc1 = fc(R/row[0],row[1]);
    }
    if (recalc_fc2) {
      G_fc2 = fc(Ru/row[0],row[1]);
    }

    return G_prefactor*G_angular*G_exp1*G_exp2*G_fc1*G_fc2;

  }


  inline REAL d_fc(REAL Rc, REAL Rij, int fcswitch)
  {
    if(Rij/Rc >= 1.0) {
      return 0.0;
    } else {
      if(fcswitch == 0) {
        REAL x = M_PI*Rij/Rc;
        if (x < 1.047197551196598) {
          return x*(0.410350793751469*x - 1.718873385392471)/Rc;
        } else if (x < 2.094395102393195) {
          x -= 1.047197551196598;
          return (x*(0.820701587502935*x - 0.859436692696234) - 1.35)/Rc;
        } else {
          x -= 2.094395102393195;
          return (x*(0.410350793751469*x + 0.859436692696233) - 1.35)/Rc;
        }

      } else {
        return (-3*tanh(1-Rij/Rc)*tanh(1-Rij/Rc)*(1/cosh(1-Rij/Rc))*(1/cosh(1-Rij/Rc)))/(2*Rij*Rc);
      }
    }

  }
  inline void d_G1(vector <REAL>& force, vector <REAL>& row, REAL R, vector <REAL>& del)
  {
    for(int i = 0; i<3; i++) {
      force[i] = del[i]*d_fc(row[0],R,row[1]);
    }
  }
  inline void d_G2(vector <REAL>& force,vector <REAL>& row, REAL R, vector <REAL>& del)
  {
    for(int i =0; i<3; i++) {
      force[i] = exp(-row[2]*pow(R-row[3],2))*(d_fc(row[0],R,row[1])*del[i]-del[i]*fc(R/row[0],row[1])*((2*row[2]*(R-row[3]))));
    }
  }
  inline void d_G3(vector <REAL>& force,vector <REAL>& row,REAL R,REAL Ru, REAL Rjk, REAL theta,vector <REAL>& Dj, vector<REAL>& Dk)
  {
    REAL dot =  Dj[0]*Dk[0]+Dj[1]*Dk[1]+Dj[2]*Dk[2];
    REAL exponent = exp(-row[2]*(R*R+Ru*Ru+Rjk*Rjk));
    REAL lamcos = pow(1+row[4]*theta,row[3]);
    for(int i =0; i <3; i++) {
      REAL tmp;
      tmp = (row[3])*pow((1+row[4]*theta),row[3]-1)*row[4]*2*((Dj[i]+Dk[i])/(R*Ru))*exponent*fc(R/row[0],row[1])*fc(Ru/row[0],row[1]);
      tmp+= lamcos*exponent*(d_fc(row[0],R,row[1])*Dj[i]*fc(Ru/row[0],row[1])*fc(Rjk/row[0],row[1]));
      tmp+= lamcos*exponent*(fc(R/row[0],row[1])*d_fc(row[0],R,row[1])*Dk[i]);
      tmp -= lamcos*exponent*2*row[2]*(Dj[i]+Dk[i]);
      tmp *= pow(2,1-row[3]);
      force[i] = tmp;
    }
  }
  inline void d_G4(vector <REAL>& force,vector <REAL>& row,REAL R,REAL Ru, REAL Rjk, REAL theta,vector <REAL>& Dj, vector<REAL>& Dk)
  {
    REAL dot =  Dj[0]*Dk[0]+Dj[1]*Dk[1]+Dj[2]*Dk[2];
    REAL exponent = exp(-row[2]*(R*R+Ru*Ru));
    REAL lamcos = pow(1+row[4]*theta,row[3]);
    for(int i =0; i <3; i++) {
      REAL tmp;
      tmp = (row[3])*pow((1+row[4]*theta),row[3]-1)*row[4]*2*((Dj[i]+Dk[i])/(R*Ru))*exponent*fc(R/row[0],row[1])*fc(Ru/row[0],row[1]);
      tmp+= lamcos*exponent*(d_fc(row[0],R,row[1])*Dj[i]*fc(Ru/row[0],row[1]));
      tmp+= lamcos*exponent*fc(R/row[0],row[1])*d_fc(row[0],R,row[1])*Dk[i];
      tmp -= lamcos*exponent*2*row[2]*(Dj[i]+Dk[i]);
      tmp *= pow(2,1-row[3]);
      force[i] = tmp;
    }
  }


  void init_type_G_map(vector<int> atom_types);

  inline int G_index(int G_number, int N)
  {
      //cout << "atom_index " << N << " " << atom_index[N] << endl;
    return Npair_Gs*(atom_index[N])+G_number;

  }

  inline int G_index(int G_number, int N1, int N2)
  {
    return Npair_Gs*atom_index.size()+atom_matrix[atom_index[N1]][atom_index[N2]]+G_number;
  }
  
  //AMP Helper-functions
  /*  
  inline REAL kd(int i, int j) {
      if (i == j) { return 1.0; }
      else { return 0.0; }
  }
  
  inline vector <REAL> delRdelR_vec (int i, int j, int m, int l) {
      vector <REAL> dir_(3,1.0);
      REAL prefactor = kd(m,j) - kd(m,i);
      for (int i = 0; i < 3; i++) {
          dir_[i] = prefactor*kd(i,l);
      }
      return dir_;
  }
  
  inline REAL delRdelR (int i, int j, int m,REAL del_l, REAL Rij) {
      return (kd(m,j) - kd(m,i))*del_l/Rij;
  }
  
  template <typename T>
  inline T dot (vector <T> i, vector<T> j) {
      T r; 
      r = i[0]*j[0] + i[1]*i[1] + i[2]*i[2];
      return r;
  }
  //i,j,k,m are just indices
  // REAL Rij, Rik, R are the magnitudes of delij, delik, and R is the location of atom m
  // delij, delik are not unit vectors 
  inline REAL delCosdelR(int i, int j, int m, int l, int k, REAL Rij, REAL Rik, vector <REAL> delij, vector <REAL> delik) {
      vector <REAL> delRij_delRml_vec, delRik_delRml_vec;
      REAL delRij_delml, delRik_delml;
      REAL dot_1, dot_2, dot_3, dot_4;
      delRij_delRml_vec = delRdelR_vec(i,j,m,l);
      delRik_delRml_vec = delRdelR_vec(i,k,m,l);
      delRij_delml = delRdelR(i,j,m,delij[l],Rij);
      delRik_delml = delRdelR(i,k,m,delik[l],Rik);
      dot_1 = dot(delRij_delRml_vec,delik);
      dot_2 = dot(delij,delRik_delRml_vec);
      dot_3 = dot(delij,delik);
      dot_4 = dot(delij,delik);
      return (1/(Rij*Rik))*dot_1 + (1/(Rij*Rik))*dot_2 - (dot_3/(Rij*Rij*Rik))*delRij_delml - (dot_4/(Rij*Rik*Rik))*delRik_delml;
      
  }
  */

};




#endif
