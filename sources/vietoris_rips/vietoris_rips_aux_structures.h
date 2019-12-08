/*
    This file is part of the Stellar library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Stellar library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Stellar library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Stellar library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef VIETORIS_RIPS_AUX_STRUCTURES_H
#define VIETORIS_RIPS_AUX_STRUCTURES_H

#include <iostream>

#include "utilities/basic_wrappers.h"
#include "stellar_tree/stellar_tree.h"

using namespace std;

class local_generation_parameters
{
public:
    local_generation_parameters(double eps, Stellar_Tree* tree, bool debug)
    {
        this->eps = eps;
        this->tree = tree;
        this->debug = debug;
        maxV = maxOuterV = 0;
        maxR = maxP = maxX = 0;
        maxE = maxOuterE = 0;
    }

    inline void print_local_stats()
    {
        cerr<<"max-vertices: "<<maxV<<" -- max-entries:"<<maxE<<endl;
        cerr<<"max-outer-vertices: "<<maxOuterV<<" -- max-outer-entries:"<<maxOuterE<<endl;
        cerr<<"max-setR: "<<maxR<<" -- max-setP: "<<maxP<<" -- max-setX: "<<maxX<<endl;
    }

    double eps;
    Stellar_Tree* tree;
    bool debug;
    unsigned maxV, maxOuterV;
    long long maxE, maxOuterE;
    int maxR, maxP, maxX;
};

#endif // VIETORIS_RIPS_AUX_STRUCTURES_H
