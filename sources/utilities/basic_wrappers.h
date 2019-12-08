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

#ifndef TOPOLOGICAL_AUXILIARY_STRUCTURES_H
#define TOPOLOGICAL_AUXILIARY_STRUCTURES_H

#include <vector>
#include <set>
#include <queue>
#include <map>

typedef std::queue<int> iqueue;
typedef std::vector<int> ivect;
typedef ivect::iterator ivect_iter;
typedef ivect::const_iterator ivect_const_iter;

typedef std::vector<double> dvect;
typedef dvect::iterator dvect_iter;
typedef dvect::const_iterator dvect_const_iter;

typedef std::set<int> iset;
typedef iset::iterator iset_iter;
typedef iset::const_iterator iset_const_iter;

typedef std::vector<ivect> VT;
typedef std::vector<ivect> ET;

typedef std::vector<VT> leaf_VT;
typedef std::vector<iset> leaf_VV;
typedef std::map<ivect,ET> leaf_ET;

typedef std::set<ivect> leaf_p_faces;

#endif // TOPOLOGICAL_AUXILIARY_STRUCTURES_H
