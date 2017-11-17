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




#ifndef __ANALYSIS__
#define __ANALYSIS__

#include <vector>


#include "Common.h"



using namespace std;

class Analysis
{

public:

  Analysis();
  ~Analysis();

  vector<int> histogram(const vector<REAL> &data, int Nbins = 0, REAL min=NOT_SET, REAL max=NOT_SET);
  REAL KS(const vector<REAL> &data1, const vector<REAL> &data2);
  pair<REAL,REAL> bounds(const vector<REAL> &data);

  REAL p;

private:



};


#endif

