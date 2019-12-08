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

#include "topological_query_extractor.h"

namespace topological_queries
{

void extract_VTop_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool,int> &p)
{
    extract_VTop(n,mesh,p);
}

void extract_VTop_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, pair<bool,int> &p)
{
    extract_VTop(n,mesh,p);
}

void extract_VV_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool,int> &p)
{
    extract_VV(n,mesh,p);
}

void extract_VV_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, pair<bool,int> &p)
{
    extract_VV(n,mesh,p);
}

void extract_links_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<bool,int,int,long> &p)
{
    extract_links(n,mesh,p);
}

void extract_links_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, tuple<bool,int,int,long> &p)
{
    extract_links(n,mesh,p);
}

void extract_v_degrees_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, ivect &p)
{
    extract_v_degree(n,mesh,p);
}

void extract_v_degrees_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, ivect &p)
{
    extract_v_degree(n,mesh,p);
}

void extract_VTop_stats_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, vector</*tuple<int,double,int>*/variance> &stats)
{
    extract_VTop_stats(n,mesh,stats);
}

void extract_VTop_stats_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, vector<variance> &stats)
{
    extract_VTop_stats(n,mesh,stats);
}

void extract_VTop_VE_sizes_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, vector<pair<int,int> > &v_stats)
{
    extract_VTop_VE_sizes(n,mesh,v_stats);
}

void extract_VTop_VE_sizes_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, vector<pair<int,int> > &v_stats)
{
    extract_VTop_VE_sizes(n,mesh,v_stats);
}

void extract_p_faces_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<int, bool, unsigned> &t)
{
    extract_p_faces(n,mesh,t);
}

void extract_p_faces_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, tuple<int, bool, unsigned> &t)
{
    extract_p_faces(n,mesh,t);
}

void extract_p_faces_Simplicial_wrapper_v2(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<int, bool, unsigned> &t)
{
    extract_p_faces_v2(n,mesh,t);
}

void extract_p_faces_CP_wrapper_v2(Node_Stellar &n, CP_Mesh &mesh, tuple<int, bool, unsigned> &t)
{
    extract_p_faces_v2(n,mesh,t);
}

}
