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

#include "connectedness_validator_skeleton.h"

/// IA* ///
void connVal_skeleton::extract_top_skeleton(Top_CP_Cell &t)
{
    ivect edge; edge.assign(2,-1);
    for(int i=0; i<t.get_edges_num(); i++)
    {
        t.TE(edge,i);
        this->check_edge(edge);
    }
}

/// IA* ///
void connVal_skeleton::extract_top_skeleton(Top_Simplex &t)
{
    ivect edge; edge.assign(2,-1);
    // extraction of all the edges of a top simplex
    // not good-looking but it is the most efficient way to extract these iteratively
    for(int j=0; j<t.get_vertices_num(); j++)
    {
        for(int i=j+1; i<t.get_vertices_num(); i++)
        {
            t.TE(edge,j,i);
            this->check_edge(edge);
        }
    }
}

/// IA* ///
void connVal_skeleton::check_edge(ivect &edge)
{
    pair<map<ivect, int>::iterator, bool> ret = this->insert_edge(edge);
    if(ret.second) /// inserted
    {
        this->add_to_ve(edge[0]-1,ret.first->second);
            this->add_to_ve(edge[1]-1,ret.first->second);
    }
}
