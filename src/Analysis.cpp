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
// This class performs some basic analysis to help determine the 
// quality of a fit. It implements histograms and an approximate
// Kolmogorov-Smirnov test for the similarity of two distributions.
// ####################################################################




#include "Analysis.h"
#include <algorithm>

// ########################################################
//                       Construtor
// ########################################################
//

Analysis::Analysis() { 

}

// ########################################################
// ########################################################



// ########################################################
//                       Constructor
// ########################################################
//

Analysis::~Analysis() {

}

// ########################################################
// ########################################################




// ########################################################
//                       BOUNDS
// ########################################################
// Finds the max and min of the given distribution

pair<REAL,REAL> Analysis::bounds(const vector<REAL> &data) {
  REAL min = data[0];
  REAL max = data[0];

  for (int i=1; i<data.size(); i++) {
    if (data[i] < min) { min = data[i]; }
    if (data[i] > max) { max = data[i]; }

  }

  return pair<REAL,REAL>(min,max);

}

// ########################################################
// ########################################################




// ########################################################
//                       HISTOGRAM
// ########################################################
// Generate a histogram of the distribution. (does not plot)

vector<int> Analysis::histogram(const vector<REAL> &data, int Nbins, REAL min, REAL max) {
  
  if (Nbins == 0) {
    Nbins = (int)((REAL)data.size()/10.0);
  }
  
  if (min == NOT_SET) {
    pair<REAL,REAL> B = this->bounds(data);
    min = B.first;
    if (max == NOT_SET) {
      max = B.second;
    }
  }
  
  REAL width = (max-min)/(Nbins-1e-6);
  vector<int> hist(Nbins,0);
  
  for (int i=0; i<data.size(); i++) {
    hist.at((int)((data[i]-min)/width))++;
  }
  
  return hist;
  
}

// ########################################################
// ########################################################





// ########################################################
//                       KS
// ########################################################
// Routine to calculate the APPROXIMATE Kolmogorov-Smirnov
// statistic, which gives indications about the simialrity
// of 2 distributions.  Here it is used to determine the 
// fidelity of a training set during a "validate" run.

REAL Analysis::KS(const vector<REAL> &data1, const vector<REAL> &data2) {
  
  pair<REAL,REAL> B1 = this->bounds(data1);
  pair<REAL,REAL> B2 = this->bounds(data2);
  
  REAL min = (B1.first < B2.first ? B1.first : B2.first);
  REAL max = (B1.second > B2.second ? B1.second : B2.second);

  int N = data1.size();

  vector<int> H1 = this->histogram(data1, N, min, max);
  vector<int> H2 = this->histogram(data2, N, min, max);
  
  REAL D = 0;
  REAL v1 = 0.0;
  REAL v2 = 0.0;
  
  for (int i=0; i<N; i++) {
    v1 += (REAL)H1[i]/data1.size();
    v2 += (REAL)H2[i]/data2.size();
    if (abs(v1-v2) > D) {
      D = abs(v1-v2);
    }
  }
  
  if (D == 0) {
    this->p = 1.0;
    return 0.0;
  }
    
  
  if (N > 30) {
    this->p = 0.0;
    for (int i=1; i<=501; i++) {
      p += pow(-1,i-1)*exp(-i*i*(REAL)N*D*D);
    }
    p *= 2;
  } else {
    this->p = 0.0;
  }

  
  return D;
  
}

// ########################################################
// ########################################################





