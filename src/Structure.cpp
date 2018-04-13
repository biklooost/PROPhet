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



#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <math.h>

#include "Structure.h"
#include "Functional_params.h"


Tables Structure::my_cos = Tables("cos");

map<REAL,REAL> Structure::prefactor_A = std::map<REAL,REAL>();
int Structure::call_count = 0;



// ########################################################
//                       Constructor
// ########################################################
//

Structure::Structure()
{

  my_is_initialized = false;

  old_R = old_Ru = old_row0 = old_row1 = old_row2 = NOT_SET;
  old_R3 = old_Ru3 = old_row03 = old_row13 = old_row23 = NOT_SET;
  Nterms = 0;

  exp_integer.assign(100, 0.0);
  for (int i=0; i<100; i++) {
    exp_integer[i] = exp(-i);
  }

  taylor.assign(500,0.0);
  for(int i=0; i<300; i++) {
    taylor[i] = 1.0/(REAL)(i+1);
  }

  is_MD = false;
  this->NORM = false;
  this->train = "";

}

// ########################################################
// ########################################################




// ########################################################
//                       Destructor
// ########################################################
//

Structure::~Structure()
{
}

// ########################################################
// ########################################################

//#########################################################
//                      train_local
// Allows user to train to per-atom energies as opposed to
// total energies.
//
//#########################################################
REAL Structure::train_Local(Functional_params *F,REAL totalenergy)
{
  map<string,REAL> FE = F->FE();
  /*
  map<string,int> count;
  REAL numerator = 0., denominator = 1.,ENERGY = 0.0;
  for(int i = 0; i < pos.size(); i++) {
      count[types[i].atomic_symbol()] += 1;
  }

  map<string,int>::iterator it;
  for (it = count.begin(); it != count.end(); it++){
      string str = it->first;
      std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
      str.erase(end_pos, str.end());
      if(FE.count(str) == 0){
          ERROR("No Free Energy supplied for " + it->first);
      } else {
          numerator += it->second*FE[str];
          denominator *= it->second;
      }
  }
  ENERGY = (totalenergy - numerator);///denominator; //note: removed this for rational training.
   */
  REAL sum = 0.0;
  for(int i = 0; i < pos.size(); i++) {
    string str = types[i].atomic_symbol();
    std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
    str.erase(end_pos, str.end());
    if(FE.count(str) == 0) {
      ERROR("No Free Energy supplied for " + str);
    } else {
      sum += FE[str];
    }
  }

  this->NORM = true;
  //ENERGY = totalenergy - sum;
  return (totalenergy - sum);
}


REAL Structure::unravel_Energy(REAL localenergy)
{
  map<string,int> count;
  /*
  REAL numerator = 0., denominator = 1.,ENERGY = 0.0;
  for(int i = 0; i < pos.size(); i++) {
      if (!this->lammps_conv.empty()) {
          count[this->lammps_conv[types[i].atomic_number()]] += 1;
      }else {
          count[types[i].atomic_symbol()] += 1;
      }
  }
  map<string,int>::iterator it;
  for (it = count.begin(); it != count.end(); it++){
      string str = it->first;
      std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
      str.erase(end_pos, str.end());
      if(this->FE.count(str) == 0){
          ERROR("No Free Energy supplied for " + it->first);
      } else {
          numerator += it->second*this->FE[str];
          denominator *= it->second;
      }
  } */
  REAL sum = 0.0;
  for(int i = 0; i < pos.size(); i++) {
    string str;
    if(!this->lammps_conv.empty()) {
      str = this->lammps_conv[types[i].atomic_number()];
    } else {
      str = types[i].atomic_symbol();
    }
    std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
    str.erase(end_pos, str.end());
    if(FE.count(str) == 0) {
      ERROR("No Free Energy supplied for " + str);
    } else {
      sum += FE[str];
    }
  }
  //ENERGY = localenergy*denominator + numerator;
  //ENERGY = localenergy + numerator;
  return localenergy + sum;

}


// ########################################################
//                       GET_FORCES
// ########################################################
// Calculates the forces for the current structure.

void Structure::Get_Forces(const vector<vector<REAL> > &dE_dG, REAL **f)
{
  vector <REAL> ir(3,0), jr(3,0), ju(3,0);
  vector <REAL> del(3,0.0), del2(3,0.0), del3(3,0.0);
  vector <REAL> row(this->inum, 0.0);

  REAL R,Ru,Rjk, G1t, G2t, G3t, G4t;
  vector<REAL> force(3,0.0);
  REAL prefactor, term, fc_R, fc_Ru, fc_Rjk, Fij, Fij2, Fjk, dist_prod;
  REAL dE_dR, dG_dR, dG_dRu, dG_dRjk, dG_dcos;
  vector<REAL> G_vec;

  if(this->periodic == true) {
    for (int ii=0; ii<this->ilist.size(); ii++) {
      int i = this->ilist[ii];
      for(int z=0; z<3; z++) {
        ir[z] = this->pos[i][z];
      }
      vector <int> jlist = this->firstneigh[i];
      int jnum = this->numneigh[i];

      for(int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];

        vector<int> T = trans_indices[ii][jj];

        for (int tt=0; tt<T.size(); tt++) {
          for(int z=0; z<3; z++) {
            jr[z] = this->pos[j][z] + translations[T[tt]][z];
            del[z] = jr[z] - ir[z];
          }
          this->To_Cart(del);
          R = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);

          if (R < 0.01) {
            continue;
          }
          for (int z=0; z<3; z++) {
            del[z] /= R;
          }

          if(!G1p.empty()) {
            for (int r1=0; r1<G1p.size(); r1++) {
              if (R >= G1p[r1][0]) {
                continue;
              }
              REAL force = dE_dG[i][G_index(r1,types[j].atomic_number())]*d_fc(G1p[r1][0],R,G1p[r1][1]);
              for (int dir=0; dir<3; dir++) {
                f[i][dir] += force*del[dir];
                f[j][dir] -= force*del[dir];
              }
            }
          }
          if(!G2p.empty()) {
            for (int r1=0; r1<G2p.size(); r1++) {
              if (R >= G2p[r1][0]) {
                continue;
              }
              REAL delta_R = R-G2p[r1][3];
              dE_dR = dE_dG[i][G_index(r1+G1p.size(),types[j].atomic_number())]*exp_term(-G2p[r1][2]*delta_R*delta_R)
                      * (d_fc(G2p[r1][0],R,G2p[r1][1]) - 2*G2p[r1][2]*delta_R*fc(R/G2p[r1][0],G2p[r1][1]));
              for (int dir=0; dir<3; dir++) {
                f[i][dir] += dE_dR*del[dir];
                f[j][dir] -= dE_dR*del[dir];
              }
            }
          }
          if (!G3p.empty() || !G4p.empty()) {
            for(int jj2 = 0; jj2< jnum; jj2++) {
              int j2 = firstneigh[i][jj2];

              vector<int> T2 = trans_indices[ii][jj2];

              for (int tt2=0; tt2<T2.size(); tt2++) {
                for(int z=0; z<3; z++) {
                  ju[z] = this->pos[j2][z] + translations[T2[tt2]][z];
                  del2[z] = ju[z] - ir[z];
                  del3[z] = ju[z] - jr[z];
                }

                Ru = sqrt(del2[0]*del2[0] + del2[1]*del2[1] + del2[2]*del2[2]);
                Rjk = sqrt(del3[0]*del3[0]+del3[1]*del3[1]+del3[2]*del3[2]);

                if (Ru < 0.0 || Rjk < 0.01) {
                  continue;
                }
                if (Ru > Rcut) {
                  continue;
                }

                for (int z=0; z<3; z++) {
                  del2[z] /= Ru;
                  del3[z] /= Rjk;
                }

                REAL theta = 0;
                for (int temp = 0; temp < 3; temp++) {
                  theta += (del2[temp])*(del[temp]);
                }

                if(!G3p.empty()) {
                    cout << "In G3 Force Calculation\n";
                  if (Rjk > Rcut) {
                    continue;
                  }

                  for (int r1=0; r1<G3p.size(); r1++) {


                    if (R >= G3p[r1][0] || Ru >= G3p[r1][0] || Rjk >= G3p[r1][0]) {
                      continue;
                    }

                    term = (1+G3p[r1][4]*theta);
                    if (term <= 0) {
                      continue;
                    }
                    prefactor = Structure::prefactor_A.at(G3p[r1][3])
                              *dE_dG[i][G_index(r1,types[j].atomic_number(),types[j2].atomic_number())]*angular_term(theta,G3p[r1]);

                    G_exp1 = exp_term(-G3p[r1][2]*(R*R + Ru*Ru + Rjk*Rjk));

                    fc_R = fc(R/G3p[r1][0],G3p[r1][1]);
                    fc_Ru = fc(Ru/G3p[r1][0],G3p[r1][1]);
                    fc_Rjk = fc(Rjk/G3p[r1][0],G3p[r1][1]);

                    dG_dR = G_exp1*fc_Ru*fc_Rjk*(d_fc(G3p[r1][0],R,G3p[r1][1]) - 2*R*G3p[r1][2]*fc_R);
                    dG_dRu = G_exp1*fc_R*fc_Rjk*(d_fc(G3p[r1][0],Ru,G3p[r1][1]) - 2*Ru*G3p[r1][2]*fc_Ru);
                    dG_dRjk = G_exp1*fc_R*fc_Ru*(d_fc(G3p[r1][0],Rjk,G3p[r1][1]) - 2*Rjk*G3p[r1][2]*fc_Rjk);
                    dG_dcos = G_exp1*(G3p[r1][3]*G3p[r1][4]/term)*fc_R*fc_Ru*fc_Rjk/(R*Ru);
                    for (int dir=0; dir<3; dir++) {
                        Fjk = prefactor*dG_dRjk*del3[dir];
                        Fij  = prefactor*(dG_dR*del[dir] + Ru*dG_dcos*(del2[dir] - theta*del[dir]));
                        Fij2 = prefactor*(dG_dRu*del2[dir] + R*dG_dcos*(del[dir] - theta*del2[dir]));
                        f[j][dir]  -= Fij + Fjk;
                        f[j2][dir] -= Fij2 - Fjk;
                        f[i][dir]  += Fij + Fij2 ;
                    } 
                    /*
                    prefactor = prefactor_A.at(G3p[r1][3])*angular_term(theta,G3p[r1])*exp_term(-G3p[r1][2]*(R*R + Ru*Ru + Rjk*Rjk)); //pow(term,G3p[r1][3])*exp(-G3p[r1][2]*(R*R + Ru*Ru + Rjk*Rjk));
                    prefactor *= dE_dG[i][G_index(r1,types[j].atomic_number(),types[j2].atomic_number())];
                    fc_R = fc(R/G3p[r1][0],G3p[r1][1]);
                    fc_Ru = fc(Ru/G3p[r1][0],G3p[r1][1]);
                    fc_Rjk = fc(Rjk/G3p[r1][0],G3p[r1][1]);
                    dist_prod = R*Ru;

                    //dG_dR = prefactor*fc_Ru*fc_Rjk*(d_fc(G3p[r1][0],R,G3p[r1][1])/R - 2*G3p[r1][2]*fc_R);
                    //dG_dRu = prefactor*fc_R*fc_Rjk*(d_fc(G3p[r1][0],Ru,G3p[r1][1])/Ru - 2*G3p[r1][2]*fc_Ru);
                    //dG_dRjk = prefactor*fc_R*fc_Ru*(d_fc(G3p[r1][0],Rjk,G3p[r1][1])/Rjk - 2*G3p[r1][2]*fc_Rjk);
                    dG_dR = prefactor*fc_Ru*fc_Rjk*(d_fc(G3p[r1][0],R,G3p[r1][1]) - 2*G3p[r1][2]*fc_R);
                    dG_dRu = prefactor*fc_R*fc_Rjk*(d_fc(G3p[r1][0],Ru,G3p[r1][1]) - 2*G3p[r1][2]*fc_Ru);
                    dG_dRjk = prefactor*fc_R*fc_Ru*(d_fc(G3p[r1][0],Rjk,G3p[r1][1]) - 2*G3p[r1][2]*fc_Rjk);
                    dG_dcos = prefactor*(G3p[r1][3]*G3p[r1][4]/term)*fc_R*fc_Ru*fc_Rjk;

                    for (int dir=0; dir<3; dir++) {
                      Fij = dG_dR*del[dir] - dG_dRjk*del3[dir]
                            + dG_dcos*(del2[dir]/dist_prod - theta*del[dir]/(dist_prod*dist_prod));
                      Fij2 = dG_dRu*del2[dir] + dG_dRjk*del3[dir]
                             + dG_dcos*(del[dir]/dist_prod - theta*del2[dir]/(dist_prod*dist_prod));
                      //f[j][dir] -= Fij;
                      //f[j2][dir] -= Fij2;
                      //f[i][dir] += Fij + Fij2;
                    }
                     */
                  }
                }

                if (!G4p.empty()) {
                  for (int r1=0; r1<G4p.size(); r1++) {
                    if (R >= G4p[r1][0] || Ru >= G4p[r1][0]) {
                      continue;
                    }

                    term = (1+G4p[r1][4]*theta);
                    if (term <= 0) {
                      continue;
                    }
                    if (r1 == 0) {
                      prefactor = Structure::prefactor_A.at(G4p[r1][3])
                                  *dE_dG[i][G_index(r1+G3p.size(),types[j].atomic_number(),types[j2].atomic_number())]*angular_term(theta,G4p[r1]);

                      G_exp1 = exp_term(-G4p[r1][2]*(R*R + Ru*Ru));

                      fc_R = fc(R/G4p[r1][0],G4p[r1][1]);
                      fc_Ru = fc(Ru/G4p[r1][0],G4p[r1][1]);
                      dG_dR = G_exp1*fc_Ru*(d_fc(G4p[r1][0],R,G4p[r1][1]) - 2*R*G4p[r1][2]*fc_R);
                      dG_dRu = G_exp1*fc_R*(d_fc(G4p[r1][0],Ru,G4p[r1][1]) - 2*Ru*G4p[r1][2]*fc_Ru);
                      dG_dcos = G_exp1*(G4p[r1][3]*G4p[r1][4]/term)*fc_R*fc_Ru/(R*Ru);
                    } else {
                      prefactor = Structure::prefactor_A.at(G4p[r1][3])
                                  *dE_dG[i][G_index(r1+G3p.size(),types[j].atomic_number(),types[j2].atomic_number())]*angular_term(theta,G4p[r1]);
                      if (old_row2 != G4p[r1][2]) {
                        if (old_row0 != G4p[r1][0] || old_row1 != G4p[r1][1]) {
                          fc_R = fc(R/G4p[r1][0],G4p[r1][1]);
                          fc_Ru = fc(Ru/G4p[r1][0],G4p[r1][1]);
                          old_row0 = G4p[r1][0];
                          old_row1 = G4p[r1][1];
                        }
                        G_exp1 = exp_term(-G4p[r1][2]*(R*R + Ru*Ru));
                        dG_dR = G_exp1*fc_Ru*(d_fc(G4p[r1][0],R,G4p[r1][1]) - 2*R*G4p[r1][2]*fc_R);
                        dG_dRu = G_exp1*fc_R*(d_fc(G4p[r1][0],Ru,G4p[r1][1]) - 2*Ru*G4p[r1][2]*fc_Ru);
                        old_row2 = G4p[r1][2];
                      } else if (old_row0 != G4p[r1][0] || old_row1 != G4p[r1][1]) {
                        fc_R = fc(R/G4p[r1][0],G4p[r1][1]);
                        fc_Ru = fc(Ru/G4p[r1][0],G4p[r1][1]);
                        old_row0 = G4p[r1][0];
                        old_row1 = G4p[r1][1];
                      }
                      dG_dcos = G_exp1*(G4p[r1][3]*G4p[r1][4]/term)*fc_R*fc_Ru/(R*Ru);
                    }
                    for (int dir=0; dir<3; dir++) {
                      Fij  = prefactor*(dG_dR*del[dir]  + Ru*dG_dcos*(del2[dir] - theta*del[dir]));
                      Fij2 = prefactor*(dG_dRu*del2[dir] + R*dG_dcos*(del[dir] - theta*del2[dir]));
                      f[j][dir]  -= Fij;
                      f[j2][dir] -= Fij2;
                      f[i][dir]  += Fij + Fij2;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } else {
    int count = 0;
    for (int ii=0; ii<this->ilist.size(); ii++) {
      int i = this->ilist[ii];
      int jnum = this->numneigh[i];
      for(int jj = 0; jj < jnum; jj++) {
        int j = firstneigh[i][jj];
        if(j == i) {
          continue;
        }
        for(int z=0; z<3; z++) {
          del[z] = pos[j][z] - pos[i][z];
        }

        R = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);
        if (R > Rcut) {
          continue;
        }

        for (int z=0; z<3; z++) {
          del[z] /= R;
        }

        if(!G1p.empty()) {
          for (int r1=0; r1<G1p.size(); r1++) {
            if (R >= G1p[r1][0]) {
              continue;
            }
            REAL force = dE_dG[i][G_index(r1,types[j].atomic_number())]*d_fc(G1p[r1][0],R,G1p[r1][1]);

            for (int dir=0; dir<3; dir++) {
              f[i][dir] += force*del[dir];
              f[j][dir] -= force*del[dir];
            }

          }
        }
        if(!G2p.empty()) {
          for (int r1=0; r1<G2p.size(); r1++) {
            if (R >= G2p[r1][0]) {
              continue;
            }
            REAL delta_R = R-G2p[r1][3];
            dE_dR = dE_dG[i][G_index(r1+G1p.size(),types[j].atomic_number())]*exp_term(-G2p[r1][2]*delta_R*delta_R)
                    * (d_fc(G2p[r1][0],R,G2p[r1][1]) - 2*G2p[r1][2]*delta_R*fc(R/G2p[r1][0],G2p[r1][1]));

            for (int dir=0; dir<3; dir++) {
              f[i][dir] += dE_dR*del[dir];
              f[j][dir] -= dE_dR*del[dir];
            }
          }
        }

        if (!G3p.empty() || !G4p.empty()) {
          for(int jj2 = 0; jj2< jnum; jj2++) {
            int j2 = firstneigh[i][jj2];
            if((j2 == j) || (j2 == i)) {
              continue;
            }

            for(int z=0; z<3; z++) {
              del2[z] = pos[j2][z] - pos[i][z];
              del3[z] = pos[j][z] - pos[j2][z];
            }
            Ru = sqrt(del2[0]*del2[0] + del2[1]*del2[1] + del2[2]*del2[2]);
            Rjk = sqrt(del3[0]*del3[0]+del3[1]*del3[1]+del3[2]*del3[2]);

            if (Ru > Rcut) {
              continue;
            }

            for (int z=0; z<3; z++) {
              del2[z] /= Ru;
              del3[z] /= Rjk;
            }

            REAL theta = 0;
            for (int temp = 0; temp < 3; temp++) {
              theta += (del2[temp])*(del[temp]);
            }

            if(!G3p.empty()) {
                //cout << "In non-periodic G3 Force Calculation\n";
              if (Rjk > Rcut) {
                continue;
              }

              for (int r1=0; r1<G3p.size(); r1++) {

                if (R >= G3p[r1][0] || Ru >= G3p[r1][0] || Rjk >= G3p[r1][0]) {
                  continue;
                }
                term = (1+G3p[r1][4]*theta);
                if (term <= 0) {
                  continue;
                }
                /*
                prefactor = prefactor_A.at(G3p[r1][3])*pow(term,G3p[r1][3])*exp(-G3p[r1][2]*(R*R + Ru*Ru + Rjk*Rjk));
                prefactor *= dE_dG[i][G_index(r1,types[j].atomic_number(),types[j2].atomic_number())];
                fc_R = fc(R/G3p[r1][0],G3p[r1][1]);
                fc_Ru = fc(Ru/G3p[r1][0],G3p[r1][1]);
                fc_Rjk = fc(Rjk/G3p[r1][0],G3p[r1][1]);
                dist_prod = R*Ru;

                //dG_dR = prefactor*fc_Ru*fc_Rjk*(d_fc(G3p[r1][0],R,G3p[r1][1])/R - 2*G3p[r1][2]*fc_R);
                //dG_dRu = prefactor*fc_R*fc_Rjk*(d_fc(G3p[r1][0],Ru,G3p[r1][1])/Ru - 2*G3p[r1][2]*fc_Ru);
                //dG_dRjk = prefactor*fc_R*fc_Ru*(d_fc(G3p[r1][0],Rjk,G3p[r1][1])/Rjk - 2*G3p[r1][2]*fc_Rjk);
                dG_dR = prefactor*fc_Ru*fc_Rjk*(d_fc(G3p[r1][0],R,G3p[r1][1]) - 2*G3p[r1][2]*fc_R);
                dG_dRu = prefactor*fc_R*fc_Rjk*(d_fc(G3p[r1][0],Ru,G3p[r1][1]) - 2*G3p[r1][2]*fc_Ru);
                dG_dRjk = prefactor*fc_R*fc_Ru*(d_fc(G3p[r1][0],Rjk,G3p[r1][1]) - 2*G3p[r1][2]*fc_Rjk);
                dG_dcos = prefactor*(G3p[r1][3]*G3p[r1][4]/term)*fc_R*fc_Ru*fc_Rjk;
                
                for (int dir=0; dir<3; dir++) {
                  Fij = dG_dR*del[dir] 
                        + dG_dcos*(del2[dir]/dist_prod - theta*del[dir]/(dist_prod*dist_prod));
                  Fij2 = dG_dRu*del2[dir] 
                         + dG_dcos*(del[dir]/dist_prod - theta*del2[dir]/(dist_prod*dist_prod));
                  f[j][dir] -= Fij - dG_dRjk*del3[dir];
                  f[j2][dir] -= Fij2 + dG_dRjk*del3[dir];
                  f[i][dir] += Fij + Fij2;
                } */
                /* This is the old force fix:*/ 
                //if (r1 == 0) {
                  prefactor = Structure::prefactor_A.at(G3p[r1][3])
                              *dE_dG[i][G_index(r1,types[j].atomic_number(),types[j2].atomic_number())]*angular_term(theta,G3p[r1]);

                  G_exp1 = exp_term(-G3p[r1][2]*(R*R + Ru*Ru + Rjk*Rjk));

                  fc_R = fc(R/G3p[r1][0],G3p[r1][1]);
                  fc_Ru = fc(Ru/G3p[r1][0],G3p[r1][1]);
                  fc_Rjk = fc(Rjk/G3p[r1][0],G3p[r1][1]);

                  dG_dR = G_exp1*fc_Ru*fc_Rjk*(d_fc(G3p[r1][0],R,G3p[r1][1]) - 2*R*G3p[r1][2]*fc_R);
                  dG_dRu = G_exp1*fc_R*fc_Rjk*(d_fc(G3p[r1][0],Ru,G3p[r1][1]) - 2*Ru*G3p[r1][2]*fc_Ru);
                  dG_dRjk = G_exp1*fc_R*fc_Ru*(d_fc(G3p[r1][0],Rjk,G3p[r1][1]) - 2*Rjk*G3p[r1][2]*fc_Rjk);
                  dG_dcos = G_exp1*(G3p[r1][3]*G3p[r1][4]/term)*fc_R*fc_Ru*fc_Rjk/(R*Ru);

                /*} else {

                  prefactor = Structure::prefactor_A.at(G3p[r1][3])
                              *dE_dG[i][G_index(r1,types[j].atomic_number(),types[j2].atomic_number())]*angular_term(theta,G3p[r1]);
                  if (old_row23 != G3p[r1][2]) {

                    if (old_row03 != G3p[r1][0] || old_row13 != G3p[r1][1]) {
                      fc_R = fc(R/G3p[r1][0],G3p[r1][1]);
                      fc_Ru = fc(Ru/G3p[r1][0],G3p[r1][1]);
                      fc_Rjk = fc(Rjk/G3p[r1][0],G3p[r1][1]);
                      old_row03 = G3p[r1][0];
                      old_row13 = G3p[r1][1];
                    }

                    G_exp1 = exp_term(-G3p[r1][2]*(R*R + Ru*Ru + Rjk*Rjk));

                    dG_dR = G_exp1*fc_Ru*fc_Rjk*(d_fc(G3p[r1][0],R,G3p[r1][1]) - 2*R*G3p[r1][2]*fc_R);
                    dG_dRu = G_exp1*fc_R*fc_Rjk*(d_fc(G3p[r1][0],Ru,G3p[r1][1]) - 2*Ru*G3p[r1][2]*fc_Ru);
                    dG_dRjk = G_exp1*fc_R*fc_Ru*(d_fc(G3p[r1][0],Rjk,G3p[r1][1]) - 2*Rjk*G3p[r1][2]*fc_Rjk);
                    old_row23 = G3p[r1][2];

                  } else if (old_row03 != G3p[r1][0] || old_row13 != G3p[r1][1]) {
                    fc_R = fc(R/G3p[r1][0],G3p[r1][1]);
                    fc_Ru = fc(Ru/G3p[r1][0],G3p[r1][1]);
                    fc_Rjk = fc(Rjk/G3p[r1][0],G3p[r1][1]);
                    old_row03 = G3p[r1][0];
                    old_row13 = G3p[r1][1];
                  }

                  dG_dcos = G_exp1*(G3p[r1][3]*G3p[r1][4]/term)*fc_R*fc_Ru*fc_Rjk/(R*Ru);//*Rjk);

                }*/

                for (int dir=0; dir<3; dir++) {
                  Fjk = prefactor*dG_dRjk*del3[dir];
                  //Fjk = prefactor*(dG_dRjk*del3[dir] + Rjk*dG_dcos*(del[dir] - theta*del3[dir]));
                  Fij  = prefactor*(dG_dR*del[dir] + Ru*dG_dcos*(del2[dir] - theta*del[dir]));
                  Fij2 = prefactor*(dG_dRu*del2[dir] + R*dG_dcos*(del[dir] - theta*del2[dir]));
                  f[j][dir]  -= Fij + Fjk;
                  f[j2][dir] -= Fij2 - Fjk;
                  f[i][dir]  += Fij + Fij2 ;
                } 
                //Here we are implementing the AMP-notation
                /*
                vector <REAL> delij(3,0.0), delik(3,0.0), deljk(3,0.0);
                for (int zz = 0; zz < 3; zz++) {
                    delij[zz] = R*del[zz];
                    delik[zz] = Ru*del2[zz];
                    deljk[zz] = Rjk*del3[zz];
                }
                vector <REAL> G3p_(G3p[r1]);
                G3p_[3] = G3p_[3] - 1;
                prefactor = Structure::prefactor_A.at(G3p[r1][3])
                              *dE_dG[i][G_index(r1,types[j].atomic_number(),types[j2].atomic_number())]*angular_term(theta,G3p_);

                G_exp1 = exp_term(-G3p[r1][2]*(R*R + Ru*Ru + Rjk*Rjk));
                prefactor *= G_exp1; 
                fc_R = fc(R/G3p[r1][0],G3p[r1][1]);
                fc_Ru = fc(Ru/G3p[r1][0],G3p[r1][1]);
                fc_Rjk = fc(Rjk/G3p[r1][0],G3p[r1][1]);
                REAL dfc_ij, dfc_ik, dfc_jk; 
                REAL delRij, delRik, delRjk,dcos;
                dfc_ij = d_fc(G3p[r1][0],R,G3p[r1][1]);
                dfc_ik = d_fc(G3p[r1][0],Ru,G3p[r1][1]);
                dfc_jk = d_fc(G3p[r1][0],Ru,G3p[r1][1]);
                //G3p vector: Rcut cutoff_type eta zeta lambda
                //delRdelR (int i, int j, int m,REAL del_l, REAL Rij
                vector <REAL> tmp_f(3,0.0);
                for (int dir=0; dir<3; dir++) {
                    dcos = delCosdelR(i,j,i,dir,j2,R,Ru,deljk,delik);
                    delRij = delRdelR(i,j,i,delij[dir],R);
                    delRik = delRdelR(i,j2,i,delik[dir],Ru);
                    delRjk = delRdelR(j,j2,j,deljk[dir],Rjk);
                    tmp_f[dir] = prefactor*fc_R*fc_Ru*fc_Rjk*(G3p[r1][4]*G3p[r1][3]*dcos - 2*G3p[r1][2]*term*(R*delRij + Ru*delRik + Rjk*delRjk))
                            + prefactor*term*(dfc_R*delRij*fc_Ru*fc_Rjk+fc_R*)
                }*/
              }
            }

            if (!G4p.empty()) {
              for (int r1=0; r1<G4p.size(); r1++) {
                if (R >= G4p[r1][0] || Ru >= G4p[r1][0]) {
                  continue;
                }

                term = (1+G4p[r1][4]*theta);
                if (term <= 0) {
                  continue;
                }

                if (r1 == 0) {
                  prefactor = Structure::prefactor_A.at(G4p[r1][3])
                              *dE_dG[i][G_index(r1+G3p.size(),types[j].atomic_number(),types[j2].atomic_number())]*angular_term(theta,G4p[r1]);

                  G_exp1 = exp_term(-G4p[r1][2]*(R*R + Ru*Ru));

                  fc_R = fc(R/G4p[r1][0],G4p[r1][1]);
                  fc_Ru = fc(Ru/G4p[r1][0],G4p[r1][1]);

                  dG_dR = G_exp1*fc_Ru*(d_fc(G4p[r1][0],R,G4p[r1][1]) - 2*R*G4p[r1][2]*fc_R);
                  dG_dRu = G_exp1*fc_R*(d_fc(G4p[r1][0],Ru,G4p[r1][1]) - 2*Ru*G4p[r1][2]*fc_Ru);
                  dG_dcos = G_exp1*(G4p[r1][3]*G4p[r1][4]/term)*fc_R*fc_Ru/(R*Ru);

                } else {

                  prefactor = Structure::prefactor_A.at(G4p[r1][3])
                              *dE_dG[i][G_index(r1+G3p.size(),types[j].atomic_number(),types[j2].atomic_number())]*angular_term(theta,G4p[r1]);
                  if (old_row2 != G4p[r1][2]) {

                    if (old_row0 != G4p[r1][0] || old_row1 != G4p[r1][1]) {
                      fc_R = fc(R/G4p[r1][0],G4p[r1][1]);
                      fc_Ru = fc(Ru/G4p[r1][0],G4p[r1][1]);
                      old_row0 = G4p[r1][0];
                      old_row1 = G4p[r1][1];
                    }

                    G_exp1 = exp_term(-G4p[r1][2]*(R*R + Ru*Ru));

                    dG_dR = G_exp1*fc_Ru*(d_fc(G4p[r1][0],R,G4p[r1][1]) - 2*R*G4p[r1][2]*fc_R);
                    dG_dRu = G_exp1*fc_R*(d_fc(G4p[r1][0],Ru,G4p[r1][1]) - 2*Ru*G4p[r1][2]*fc_Ru);

                    old_row2 = G4p[r1][2];

                  } else if (old_row0 != G4p[r1][0] || old_row1 != G4p[r1][1]) {
                    fc_R = fc(R/G4p[r1][0],G4p[r1][1]);
                    fc_Ru = fc(Ru/G4p[r1][0],G4p[r1][1]);
                    old_row0 = G4p[r1][0];
                    old_row1 = G4p[r1][1];
                  }

                  dG_dcos = G_exp1*(G4p[r1][3]*G4p[r1][4]/term)*fc_R*fc_Ru/(R*Ru);

                }

                for (int dir=0; dir<3; dir++) {
                  Fij  = prefactor*(dG_dR*del[dir]  + Ru*dG_dcos*(del2[dir] - theta*del[dir]));
                  Fij2 = prefactor*(dG_dRu*del2[dir] + R*dG_dcos*(del[dir] - theta*del2[dir]));
                  f[j][dir]  -= Fij;
                  f[j2][dir] -= Fij2;
                  f[i][dir]  += Fij + Fij2;
                }

              }
            }
          }
        }
      }
    }
  }
  /*
  if(!this->FE.empty()) {
      for (int ii=0; ii<this->ilist.size(); ii++) {
          int i = this->ilist[ii];
          for (int dir=0; dir<3; dir++) {
              f[i][dir] *= this->pos.size();
          }
      }
  }
   */
}

// ########################################################
// ########################################################



// ########################################################
//                       PRINT
// ########################################################
// Prints the current structure to STDOUT.

void Structure::print()
{
  for (std::vector<REAL>::iterator it = this->a.begin(); it != a.end(); ++it) {
    cout << *it << "\t";
  }
  cout << endl;
  for (std::vector<REAL>::iterator it = this->b.begin(); it != b.end(); ++it) {
    cout << *it << "\t";
  }
  cout << endl;
  for (std::vector<REAL>::iterator it = this->c.begin(); it != c.end(); ++it) {
    cout << *it << "\t";
  }
  cout << endl << endl;


  if(pos.empty()) {
    cout<< "POS empty \n";
  } else {
    cout << "Specie\tx\ty\tz [A]\n";
    for(int i = 0; i < pos.size(); i++) {
      cout << types[i].atomic_symbol() << "\t";
      for(int j =0 ; j < 3; j++) {
        cout << pos[i][j] << "\t";
      }
      cout << "\n";
    }
  }
}

// ########################################################
// ########################################################




// ########################################################
//                       INIT_G
// ########################################################
// Initializes the mapping functions (G's) for the current
// potential.

vector<vector<REAL>*> Structure::init_G(Functional_params* F)
{
  this->is_MD = F->is_MD;
  this->FE = F->FE();
  this->Rcut = F->Rcut();
  vector<int> atom_types = F->atomic_numbers();

  /*
   * G1: row[0] = Rc, row[1] = Fctype
   * G2: row[0] = Rc, row[1] = Fctype, row[2] = eta, row[3] = Rs
   * G3: row[0] = Rc, row[1] = Fctype, row[2] = eta, row[3] = xi , row[4] = lambda
   * G4: row[0] = Rc, row[1] = Fctype, row[2] = eta, row[3] = xi , row[4] = lambda
   */

  G1p.clear();
  G2p.clear();
  G3p.clear();
  G4p.clear();

  for (int i=0; i<F->G1.size(); i++) {
    Allocate_G(1,F->G1.at(i));
  }
  for (int i=0; i<F->G2.size(); i++) {
    Allocate_G(2,F->G2.at(i));
  }
  for (int i=0; i<F->G3.size(); i++) {
    Allocate_G(3,F->G3.at(i));
    if (!Structure::prefactor_A.count(F->G3.at(i).at(3))) {
      Structure::prefactor_A.insert(pair<REAL,REAL>(F->G3.at(i).at(3),pow(2,1-F->G3.at(i).at(3))));
    }
  }
  for (int i=0; i<F->G4.size(); i++) {
    Allocate_G(4,F->G4.at(i));
    if (!Structure::prefactor_A.count(F->G4.at(i).at(3))) {
      Structure::prefactor_A.insert(pair<REAL,REAL>(F->G4.at(i).at(3),pow(2,1-F->G4.at(i).at(3))));
    }
  }


  if (F->G1.size()+F->G2.size()+F->G3.size()+F->G4.size() == 0) {
    if (F->Nradial()) {
      REAL stepsize = (F->Rcut())/(F->Nradial()-1);
      for (double rs=0.0; rs<=F->Rcut(); rs += stepsize) {
        vector<REAL> temp;
        temp.push_back(F->Rcut());
        temp.push_back(0);
        temp.push_back(11.5129254649702/(4*stepsize*stepsize));
        temp.push_back(rs);
        Allocate_G(2,temp);
        F->G2.push_back(temp);
      }
    }

    if (F->Nangular()) {
      double stepsize = 2*15.0/(F->Nangular()-2.0);
      for (double xi=1.0; xi<=16.001; xi+=stepsize) {
        vector<REAL> temp;
        temp.push_back(F->Rcut());
        temp.push_back(0);
        temp.push_back(4.60517018598809/(F->Rcut()*F->Rcut()));
        temp.push_back(xi);
        if (!Structure::prefactor_A.count(xi)) {
          Structure::prefactor_A.insert(pair<REAL,REAL>(xi,pow(2,1-xi)));
        }
        temp.push_back(1);
        Allocate_G(4,temp);
        F->G4.push_back(temp);
        temp.back() = -1;
        Allocate_G(4,temp);
        F->G4.push_back(temp);
      }
    }
  }

  if(this->periodic == true) {
    this->Set_Convert();
    this->Normalize();
  }

  this->my_NG_types = G1p.size()+G2p.size()+G3p.size()+G4p.size();
  this->my_NG = atom_types.size()*(G1p.size()+G2p.size())
                + 0.5*atom_types.size()*(atom_types.size()+1)*(G3p.size()+G4p.size());

  F->Ninput_nodes(my_NG);

  G.assign(my_NG,vector<REAL>(1,0.0));


  if (!F->is_MD) {

    if (!my_is_initialized) {
      this->Set_Neighbor(F->Rcut());
      init_type_G_map(atom_types);
    }

    vector<vector<REAL> > all_G(my_NG,vector<REAL>(this->Natom,0.0));
    my_sums.clear();
    my_counts.clear();
    for (int i=0; i<types.size(); i++) {
      my_sums.insert(pair<int,vector<REAL> >(types[i].atomic_number(),vector<REAL>(my_NG,0.0)));
      my_counts.insert(pair<int,int> (types[i].atomic_number(), 0));
    }



    for (int i=0; i<this->Natom; i++) {
      Calc_G(i);

      my_counts.at(types[i].atomic_number())++;
      for (int j=0; j<my_NG; j++) {
        all_G.at(j).at(i) = G.at(j).at(0);
        my_sums.at(types[i].atomic_number()).at(j) += G.at(j).at(0);
      }
    }
    G = all_G;

    G_ptr.clear();
    for (int i=0; i<my_NG; i++) {
      G_ptr.push_back(&G.at(i));
    }

  } else {

    this->Forces.assign(this->Natom, vector<REAL>(3,0.0));

    init_type_G_map(atom_types);

    Calc_G(0);
    G_ptr.clear();
    for (int i=0; i<my_NG; i++) {
      G_ptr.push_back(&G.at(i));
    }
  }

  my_is_initialized = true;

  return G_ptr;

}

// ########################################################
// ########################################################



// ########################################################
//                       ALLOCATE_G
// ########################################################
// Performs more set up needed for the mapping functions.

void Structure::Allocate_G(int Gtype, vector <REAL>& row)
{


  if(Gtype ==1) {
    if(row.size() == 2) {
      G1p.push_back(row);
    } else {
      ERROR("Improper number of parameters");
    }
  } else if(Gtype ==2) {
    if(row.size() == 4) {
      G2p.push_back(row);
    } else {
      ERROR("Improper number of parameters");
    }
  } else if(Gtype == 3) {
    if(row.size() == 5) {
      G3p.push_back(row);
    } else {
      ERROR("Improper number of parameters");
    }
  } else if(Gtype == 4) {
    if(row.size() == 5) {
      G4p.push_back(row);
    } else {
      ERROR("Improper number of parameters");
    }
  } else {
    ERROR("Improper functional type" + Gtype);
  }
}

// ########################################################
// ########################################################




// ########################################################
//                       CALC_G
// ########################################################
// Calculates the values of the mapping functions for a
// given atom in the current structure.

void Structure::Calc_G(int ii)
{

  vector <REAL> ir(3,0), jr(3,0), ju(3,0);
  vector <REAL> del(3,0.0), del2(3,0.0), del3(3,0.0);
  REAL R,Ru,Rjk, G1t, G2t, G3t, G4t;


  for (int i=0; i<G.size(); i++) {
    G[i][0] = 0.0;
  }
  if(this->periodic == true) {
    int i = this->ilist[ii];
    for(int z=0; z<3; z++) {
      ir[z] = this->pos[i][z];
    }
    vector <int> jlist = this->firstneigh[i];
    int jnum = this->numneigh[i];
    for(int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];

      vector<int> T;
      if (this->trans_indices.empty()) {
        T = vector<int>(3,0.0);
      } else {
        T = this->trans_indices[ii][jj];
      }


      for (int tt=0; tt<T.size(); tt++) {

        for(int z=0; z<3; z++) {
          jr[z] = this->pos[j][z] + translations[T[tt]][z];
          del[z] = jr[z] - ir[z];
        }

        this->To_Cart(del);
        R = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);
        if (R < 0.01) {
          continue;
        }
        if(!G1p.empty()) {
          for(int r1 = 0; r1<G1p.size(); r1++ ) {
            this->G[G_index(r1,types[j].atomic_number())][0] += G1(G1p[r1],R);
          }
        }
        if(!G2p.empty()) {
          for (int r1=0; r1<G2p.size(); r1++) {
            this->G[G_index(r1+G1p.size(),types[j].atomic_number())][0] += G2(G2p[r1],R);
          }
        }

        if (!G3p.empty() || !G4p.empty()) {
          for(int jj2 = 0; jj2< jnum; jj2++) {
            int j2 = jlist[jj2];
            vector<int> T2;
            if (this->trans_indices.empty()) {
              T2 = vector<int>(3,0.0);
            } else {
              T2 = this->trans_indices[ii][jj2];
            }

            for (int tt2=0; tt2<T2.size(); tt2++) {

              for(int z=0; z<3; z++) {
                del2[z] = (pos[j2][z]+translations[T2[tt2]][z]) - ir[z];
                del3[z] = (pos[j2][z]+translations[T2[tt2]][z]) - jr[z];
              }
              this->To_Cart(del2);
              this->To_Cart(del3);

              Ru =  sqrt(del2[0]*del2[0] + del2[1]*del2[1] + del2[2]*del2[2]);
              Rjk = sqrt(del3[0]*del3[0] + del3[1]*del3[1] + del3[2]*del3[2]);

              if (Ru < 0.01 || Rjk < 0.01) {
                continue;
              }

              REAL theta = 0;
              for (int temp = 0; temp < 3; temp++) {
                theta += (del2[temp])*(del[temp]);
              }
              theta = theta/(R*Ru);
              if (!G3p.empty()) {
                for (int r1=0; r1<G3p.size(); r1++) {
                  this->G[G_index(r1,types[j].atomic_number(),types[j2].atomic_number())][0] += G3(G3p[r1],R,Ru,Rjk,theta);
                }
              }
              if (!G4p.empty()) {
                for (int r1=0; r1<G4p.size(); r1++) {
                  this->G[G_index(r1+G3p.size(),types[j].atomic_number(),types[j2].atomic_number())][0] += G4(G4p[r1],R,Ru,Rjk,theta);
                }
              }
            }
          }
        }
      }
    }

  } else {

    int i = this->ilist[ii];
    for(int z=0; z<3; z++) {
      ir[z] = this->pos[i][z];
    }
    vector <int> jlist = this->firstneigh[i];
    int jnum = this->numneigh[i];

    for(int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];

      if(j == i) {
        continue;
      }

      for(int z=0; z<3; z++) {
        jr[z] = this->pos[j][z];
        del[z] = jr[z] - ir[z];
      }
      R = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);
      if (R > Rcut) {
        continue;
      }

      if(!G1p.empty()) {
        for(int r1 = 0; r1<G1p.size(); r1++ ) {
          this->G[G_index(r1,types[j].atomic_number())][0] += G1(G1p[r1],R);
        }
      }
      if(!G2p.empty()) {
        for (int r1=0; r1<G2p.size(); r1++) {
          this->G[G_index(r1+G1p.size(),types[j].atomic_number())][0] += G2(G2p[r1],R);
        }
      }
      if (!G3p.empty() || !G4p.empty() ) {
        for(int jj2 = 0; jj2< jnum; jj2++) {
          int j2 = jlist[jj2];
          if((j2 == j) || (j2 == i )) {
            continue;
          }
          for(int z=0; z<3; z++) {
            ju[z] = this->pos[j2][z];
            del2[z] = ju[z] - ir[z];
            del3[z] = jr[z] - ju[z];
          }
          Ru = sqrt(del2[0]*del2[0] + del2[1]*del2[1] + del2[2]*del2[2]);
          if (Ru > Rcut) {
            continue;
          }
          Rjk = sqrt(del3[0]*del3[0]+del3[1]*del3[1]+del3[2]*del3[2]);
          REAL theta = 0;
          for (int temp = 0; temp < 3; temp++) {
            theta += (del2[temp])*(del[temp]);
          }
          theta = theta/(R*Ru);
          if(!G3p.empty()) {
            if (Rjk > Rcut) {
              continue;
            }
            for (int r1=0; r1<G3p.size(); r1++) {
              this->G[G_index(r1,types[j].atomic_number(),types[j2].atomic_number())][0] += G3(G3p[r1],R,Ru,Rjk,theta);
            }
          }
          if(!G4p.empty()) {
            for (int r1=0; r1<G4p.size(); r1++) {
              this->G[G_index(r1+G3p.size(),types[j].atomic_number(), types[j2].atomic_number())][0] += G4(G4p[r1],R,Ru,Rjk,theta);
            }
          }
        }
      }
    }
  }

}

// ########################################################
// ########################################################



// ########################################################
//                       SET_CONVERT
// ########################################################
// Sets up matrices to move between cartesian/lattice
// cordinates.
void Structure::Set_Convert()
{
  if (!this->periodic) {
    return;
  }
  REAL Volume = this->a[0]*(this->b[1]*this->c[2]-this->b[2]*this->c[1])+this->a[1]*(this->b[2]*this->c[0]-this->b[0]*this->c[2])+this->a[2]*(this->b[0]*this->c[1]-this->b[1]*this->c[0]);
  vector <REAL> row(3,0.0);
  REAL norma = sqrt(pow(this->a[0],2) + pow(this->a[1],2) + pow(this->a[2],2));
  this->norm['a'] = norma;
  REAL normb = sqrt(pow(this->b[0],2) + pow(this->b[1],2) + pow(this->b[2],2));
  this->norm['b'] = normb;
  REAL normc = sqrt(pow(this->c[0],2) + pow(this->c[1],2) + pow(this->c[2],2));
  this->norm['c'] = normc;
  this->AtoF.clear();
  this->FtoA.clear();
  for(int i = 0; i<3; i++) {
    this->AtoF.push_back(row);
    this->FtoA.push_back(row);
  }
  vector < vector <REAL> > cell;
  cell.push_back(this->a);
  cell.push_back(this->b);
  cell.push_back(this->c);
  this->FtoA = transpose(cell);

  this->AtoF = invert(transpose(cell));
  /*
  cout << "CELL" << endl;
  for (int i =0; i<3; i++){
      for (int j = 0; j<3; j++) {
          cout << cell[i][j] << " ";
      }
      cout << endl;
  }
  cout << "AtoF" << endl;
  for (int i =0; i<3; i++){
      for (int j = 0; j<3; j++) {
          cout << this->AtoF[i][j] << " ";
      }
      cout << endl;
  }
  cout << "FtoA" << endl;
  for (int i =0; i<3; i++){
      for (int j = 0; j<3; j++) {
          cout << this->FtoA[i][j] << " ";
      }
      cout << endl;
  } */
  /*
  //The following are actually the cosine version of the angles
  REAL alpha = (this->b[0]*this->c[0]+this->b[1]*this->c[1]+this->b[2]*this->c[2])/(this->norm['b']*this->norm['c']);
  REAL beta = (this->a[0]*this->c[0]+this->a[1]*this->c[1]+this->a[2]*this->c[2])/(this->norm['a']*this->norm['c']);
  REAL gamma = (this->b[0]*this->a[0]+this->b[1]*this->a[1]+this->b[2]*this->a[2])/(this->norm['b']*this->norm['a']);
  //One sin version
  REAL sgamma = sin(acos(gamma));
  this->AtoF[0][0] = 1/this->norm['a'];
  this->AtoF[0][1] = -gamma/(this->norm['a']*sgamma);
  this->AtoF[0][2] = (((this->norm['b']*this->norm['c']*gamma*(alpha-beta*gamma))/sgamma) - this->norm['b']*this->norm['c']*beta*sgamma)/Volume;
  this->AtoF[1][1] = 1/(this->norm['b']*sgamma);
  this->AtoF[1][2] = - (this->norm['a']*this->norm['c']*(alpha-beta*gamma))/(Volume*sgamma);
  this->AtoF[2][2] = this->norm['a']*this->norm['b']*sgamma/Volume;
  this->FtoA[0][0] = this->norm['a'];
  this->FtoA[0][1] = this->norm['b']*gamma;
  this->FtoA[0][2] = this->norm['c']*beta;
  this->FtoA[1][1] = this->norm['b']*sgamma;
  this->FtoA[1][2] = this->norm['c']*(alpha-beta*gamma)/sgamma;
  this->FtoA[2][2] = Volume/(this->norm['a']*this->norm['b']*sgamma);
   */
}

// ########################################################
// ########################################################



// ########################################################
//                       NORMALIZE
// ########################################################
// Normalizes lattice vectors.

void Structure::Normalize()
{
  if ( this->periodic == true) {
    if( this->CART == true) {
      REAL norma = sqrt(pow(this->a[0],2) + pow(this->a[1],2) + pow(this->a[2],2));
      this->norm['a'] = norma;
      REAL normb = sqrt(pow(this->b[0],2) + pow(this->b[1],2) + pow(this->b[2],2));
      this->norm['b'] = normb;
      REAL normc = sqrt(pow(this->c[0],2) + pow(this->c[1],2) + pow(this->c[2],2));
      this->na.clear();
      this->nb.clear();
      this->nc.clear();

      this->norm['c'] = normc;
      for(int i = 0; i<3; i++) {
        this->na.push_back(this->a[i]/norma);
      }
      for(int i = 0; i<3; i++) {
        this->nb.push_back(this->b[i]/normb);
      }
      for(int i = 0; i<3; i++) {
        this->nc.push_back(this->c[i]/normc);
      }
      for (unsigned i=0; i< this->pos.size(); i++) {
        vector <REAL> proj;
        for(int j=0; j<3; j++) {
          proj.push_back(pos[i][j]);
        }
        this->To_Frac(proj);
        for(int j = 0; j<3; j++) {
          if (proj[j] < 0) {
            proj[j]+=1;
          } else if (proj[j] >1) {
            proj[j]-=1;
          }
          this->pos[i][j] = proj[j];
        }
      }
    }
  }
}
// ########################################################
// ########################################################
// ########################################################
//                       SET_NEIGHBOR
// ########################################################
// Determines the neighbor list for the required cutoff
// distance.
void Structure::Set_Neighbor(REAL Rcut)
{
  REAL delx, dely, delz;
  vector <REAL> ir(3,0.0), jr(3,0.0), del(3,0.0);
  vector <int> row;
  this->translations.clear();
  this->trans_indices.clear();
  this->ilist.clear();
  this->firstneigh.clear();
  this->numneigh.clear();

  if(this->periodic == true) {
    int max_lattice_X, max_lattice_Y, max_lattice_Z;
    max_lattice_X = ceil(1.1*Rcut/this->norm['a']);
    max_lattice_Y = ceil(1.1*Rcut/this->norm['b']);
    max_lattice_Z = ceil(1.1*Rcut/this->norm['c']);

    vector<int> R(3,0);
    for (int ix=-max_lattice_X; ix<=max_lattice_X; ix++) {
      for (int iy = -max_lattice_Y; iy<=max_lattice_Y; iy++) {
        for (int iz = -max_lattice_Z; iz<=max_lattice_Z; iz++) {
          R[0] = ix;
          R[1] = iy;
          R[2] = iz;
          this->translations.push_back(R);
        }
      }
    }

    for (int i=0; i<this->pos.size(); i++) {
      int tnumneigh = 0;

      this->ilist.push_back(i);
      row.clear();
      vector<vector<int> > i_trans_indices;
      for(int z=0; z<3; z++) {
        ir[z] = this->pos[i][z];
      }

      for(int j=0; j<this->pos.size(); j++) {

        vector<int> j_trans_indices;
        for (int t=0; t<translations.size(); t++) {

          for(int z=0; z<3; z++) {
            jr[z] = this->pos[j][z] + translations[t][z];
            del[z] = jr[z] - ir[z];
          }
          this->To_Cart(del);
          REAL distance = sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]);
          if (distance > 0.01 && distance <= Rcut) {
            if (j_trans_indices.empty()) {
              row.push_back(j);
              tnumneigh++;
            }
            j_trans_indices.push_back(t);
          }
        }
        if (!j_trans_indices.empty()) {
          i_trans_indices.push_back(j_trans_indices);
        }
      }
      this->firstneigh.push_back(row);
      this->numneigh.push_back(tnumneigh);
      this->trans_indices.push_back(i_trans_indices);
    }
    this->inum = this->pos.size();
  } else {
    for (int i=0; i<this->pos.size(); i++) {
      int tnumneigh = 0;
      this->ilist.push_back(i);
      row.clear();
      for(int z=0; z<3; z++) {
        ir[z] = this->pos[i][z];
      }
      for(int j=0; j<this->pos.size(); j++) {
        if(i==j) {
          continue;
        }
        for(int z=0; z<3; z++) {
          jr[z] = this->pos[j][z];
          del[z] = jr[z] - ir[z];
        }
        if(sqrt(del[0]*del[0]+del[1]*del[1]+del[2]*del[2]) < Rcut) {
          row.push_back(j);
          tnumneigh++;
        }
      }
      this->firstneigh.push_back(row);
      this->numneigh.push_back(tnumneigh);
    }
    this->inum = this->pos.size();
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       SET_STRUCTURE
// ########################################################
// Set the current structure.  Used by LAMMPS interface.

void Structure::set_structure(double** x, int* in_types, int Natoms, int Nghosts)
{

  this->Natom = Natoms;
  int Ntotal = Natoms + Nghosts;
  pos.assign(Ntotal,vector<REAL>(3,0.0));
  types.assign(Ntotal,::Atom());

  for (int i=0; i<Ntotal; i++) {
    pos[i][0] = x[i][0];
    pos[i][1] = x[i][1];
    pos[i][2] = x[i][2];
    types[i].set_type(in_types[i]);
  }

}

// ########################################################
// ########################################################




// ########################################################
//                       COUNT
// ########################################################
// Counts atoms of a specific type.

int Structure::count(int atomic_number)
{
  if (my_counts.count(atomic_number)) {
    return my_counts.at(atomic_number);
  } else {
    return 0;
  }

}

// ########################################################
// ########################################################



// ########################################################
//                       MEAN
// ########################################################
// Finds sum of each mapping function for the given atom
// type.

vector<REAL> Structure::mean(int atomic_number)
{
  if (my_sums.count(atomic_number)) {
    return my_sums.at(atomic_number);
  } else {
    return vector<REAL>(my_NG,0.0);
  }
}

// ########################################################
// ########################################################



// ########################################################
//                       VARIANCE
// ########################################################
// Finds sum of the squared difference between the mapping
// function values and their means.

vector<REAL> Structure::variance(int atomic_number, const vector<REAL> &means)
{
  vector<REAL> var(means.size(), 0.0);
  for (int atom=0; atom<Natom; atom++) {
    if (types[atom].atomic_number() == atomic_number) {
      for (int i=0; i<var.size(); i++) {
        var[i] += pow(G[i][atom] - means[i],2);
      }
    }
  }

  return var;

}

// ########################################################
// ########################################################




// ########################################################
//                       INIT_TYPE_G_MAP
// ########################################################
// Creates a map from atomic numbers of neighbor atoms
// to the correct mapping function associated with the
// interaction.

void Structure::init_type_G_map(vector<int>atom_types)
{

  this->atom_index.clear();
  this->atom_matrix.clear();
  for (int i=0; i<atom_types.size(); i++) {
    this->atom_index.insert(pair<int,int>(atom_types[i],i));
  }

  int index = 0;
  this->Npair_Gs = G1p.size()+G2p.size();
  this->atom_matrix.assign(atom_types.size(),vector<int>(atom_types.size(),0));

  for (int i=0; i<atom_types.size(); i++) {
    for (int j=i; j<atom_types.size(); j++) {
      this->atom_matrix[i][j] = index;
      this->atom_matrix[j][i] = index;
      index += G3p.size()+G4p.size();
    }
  }

}

// ########################################################
// ########################################################



