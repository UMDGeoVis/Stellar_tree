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

#ifndef TOPOLOGICALQUERIES_H
#define TOPOLOGICALQUERIES_H

#include <vector>
#include <list>
#include <map>
#include <tuple>
#include <omp.h>

#include <stellar_tree/box.h>
#include <stellar_tree/mesh.h>
#include <stellar_decomposition/run_iterator.h>
#include <stellar_decomposition/node_stellar.h>
#include <utilities/basic_wrappers.h>
#include <topological_ds/links_aux_structures.h>

using namespace std;

struct variance
{
    variance()
    {
        min = INT_MAX;
        avg = 0;
        max = 0;
        v_counter = 0;
    }

    int min;
    double avg;
    int max;
    int v_counter;
};

namespace topological_queries
{
/**
 * @brief A procedure defining a wrapper to extract the Vertex-Top relations on simplicial meshes in a given leaf block
 *
 * @param n the current leaf block
 * @param mesh the mesh indexed
 * @param p a pair containing debugging symbols, the first is a boolean enabling a debuggin session, while the other is a counter for the elements in all the VTop relations of the block
 */
void extract_VTop_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool,int> &p);
/**
 * @brief A procedure defining a wrapper to extract the Vertex-Top relations on CP meshes in a given leaf block
 *
 * @param n the current leaf block
 * @param mesh the mesh indexed
 * @param p a pair containing debugging symbols, the first is a boolean enabling a debuggin session, while the other is a counter for the elements in all the VTop relations of the block
 */
void extract_VTop_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, pair<bool,int> &p);
/**
 * @brief A procedure defining a wrapper to extract the Vertex-Vertex relations on simplicial meshes in a given leaf block
 *
 * @param n the current leaf block
 * @param mesh the mesh indexed
 * @param p a pair containing debugging symbols, the first is a boolean enabling a debuggin session, while the other is a counter for the elements in all the VTop relations of the block
 */
void extract_VV_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool,int> &p);
/**
 * @brief A procedure defining a wrapper to extract the Vertex-Vertex relations on CP meshes in a given leaf block
 *
 * @param n the current leaf block
 * @param mesh the mesh indexed
 * @param p a pair containing debugging symbols, the first is a boolean enabling a debuggin session, while the other is a counter for the elements in all the VTop relations of the block
 */
void extract_VV_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, pair<bool,int> &p);
/**
 * @brief A procedure defining a wrapper to extract the links relations on simplicial meshes in a given leaf block
 *
 * @param n the current leaf block
 * @param mesh the mesh indexed
 * @param p a pair containing debugging symbols, the first is a boolean enabling a debuggin session, while the others are counters for the elements encoded by the local links
 */
void extract_links_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<bool,int,int,long> &p);
/**
 * @brief A procedure defining a wrapper to extract the links relations on simplicial meshes in a given leaf block
 *
 * @param n the current leaf block
 * @param mesh the mesh indexed
 * @param p a pair containing debugging symbols, the first is a boolean enabling a debuggin session, while the others are counters for the elements encoded by the local links
 */
void extract_links_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, tuple<bool, int, int, long> &p);

void extract_v_degrees_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, ivect &p);
void extract_v_degrees_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, ivect &p);

void extract_VTop_stats_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, vector<variance> &stats);
void extract_VTop_stats_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, vector</*tuple<int, double, int>*/variance> &stats);

void extract_VTop_VE_sizes_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, vector<pair<int,int> > &v_stats);
void extract_VTop_VE_sizes_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, vector<pair<int,int> > &v_stats);

void extract_ETop_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool,int> &p);
void extract_FTop_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool,int> &p);

void extract_VE_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool,int> &p);
void extract_VF_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool,int> &p);

void extract_EF_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool,int> &p);


/**
 * @brief A procedure that extracts the local Vertex-Top relations for all the indexed vertices in a given leaf block
 *
 * @param n the current leaf block
 * @param mesh the mesh indexed
 * @param p a pair containing debugging symbols, the first is a boolean enabling a debuggin session, while the other is a counter for the elements in all the VTop relations of the block
 */
template<class C, class T> void extract_VTop(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p);

template<class C, class T> void extract_VTop(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p)
{
    leaf_VT allVTres;
    n.extract_local_VTop(mesh,allVTres);
    int isolated_vertices = 0;
    if(p.first)
    {
        int list_summation = 0;
        for(unsigned i=0; i<allVTres.size(); i++)
        {
            bool isolated = true;
            for(unsigned j=0; j<allVTres[i].size(); j++)
            {
                list_summation += allVTres[i][j].size();
                if(allVTres[i][j].size() > 0)
                    isolated = false;
            }
            if(isolated)
                isolated_vertices++;
        }

        if(p.second < list_summation)
            p.second = list_summation;
        if(isolated_vertices > 0)
            cout<<n<<"ISOLATED VERTICES: "<<isolated_vertices<<endl;
    }
}

template<class C, class T> void extract_VTop_stats(Node_Stellar &n, Mesh<C,T> &mesh, vector<variance/*tuple<int,double,int>*/ > &stats);

template<class C, class T> void extract_VTop_stats(Node_Stellar &n, Mesh<C,T> &mesh, vector<variance> &stats)
{
    leaf_VT allVTres;
    n.extract_local_VTop(mesh,allVTres);

    for(unsigned v=0; v<allVTres.size(); v++)
    {
        VT &vt = allVTres[v];
        int sum =0;
        for(unsigned d=0; d<vt.size(); d++)
        {
            int cob_size = vt[d].size();
            if(cob_size == 0) /// temporary
                continue;
            sum += cob_size;
            if(stats[d].min > cob_size)
                stats[d].min = cob_size;
            if(stats[d].max < cob_size)
                stats[d].max = cob_size;
            stats[d].avg += cob_size;
            if(cob_size > 0)
                stats[d].v_counter++;
        }

        if(stats[stats.size()-1].min > sum)
            stats[stats.size()-1].min = sum;
        if(stats[stats.size()-1].max < sum)
            stats[stats.size()-1].max = sum;
        stats[stats.size()-1].avg += sum;
    }
}

template<class C, class T> void extract_VTop_VE_sizes(Node_Stellar &n, Mesh<C,T> &mesh, vector<pair<int,int> > &v_stats);

template<class C, class T> void extract_VTop_VE_sizes(Node_Stellar &n, Mesh<C,T> &mesh, vector<pair<int,int> > &v_stats)
{
    leaf_VT allVTres;
    n.extract_local_VTop(mesh,allVTres);

    leaf_VV vvs;
    n.extract_local_VV(mesh,vvs);

    for(unsigned v=0; v<allVTres.size(); v++)
    {
        int global_v_pos = v + n.get_v_start() - 1;

        VT &vt = allVTres[v];
        int vtop_sum =0;
        for(unsigned d=0; d<vt.size(); d++)
        {
            int cob_size = vt[d].size();
            vtop_sum += cob_size;
        }

        int ve_sum = 0;
        ve_sum += vvs[v].size();

        pair<int,int> p = make_pair(vtop_sum,ve_sum);
        v_stats[global_v_pos] = p;
    }
}


/**
 * @brief A procedure that extracts the local Vertex-Vertex relations for all the indexed vertices in a given leaf block
 *
 * @param n the current leaf block
 * @param mesh the mesh indexed
 * @param p a pair containing debugging symbols, the first is a boolean enabling a debuggin session, while the other is a counter for the elements in all the VTop relations of the block
 */
template<class C, class T> void extract_VV(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p);

template<class C, class T> void extract_VV(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p)
{

    leaf_VV vvs;
    n.extract_local_VV(mesh,vvs);

    if(p.first)
    {
        int list_summation = 0;
        for(unsigned i=0; i<vvs.size(); i++)
        {
            list_summation += vvs[i].size();
        }

        if(p.second < list_summation)
            p.second = list_summation;
    }
}

/**
 * @brief A procedure that extract the links of the vertices and edges indexed by the leaf block
 *
 * @param n the current leaf block
 * @param mesh the mesh indexed
 * @param p a pair containing debugging symbols, the first is a boolean enabling a debuggin session, while the others are counters for the elements encoded by the local links
 */
template<class C, class T> void extract_links(Node_Stellar &n, Mesh<C,T> &mesh, tuple<bool,int,int,long> &p);

template<class C, class T> void extract_links(Node_Stellar &n, Mesh<C,T> &mesh, tuple<bool,int,int,long> &p)
{
    local_links links;
    n.extract_local_links(mesh,links);

    if(get<0>(p))
    {
        long reference_summation = 0;
        int vertices_num, edges_num = 0;

        vertices_num = n.get_real_v_array_size();
        for(vertices_link_list::iterator it=links.begin_vertices(); it!=links.end_vertices(); ++it)
        {
            reference_summation += it->get_references_number();
        }

        edges_num = links.edges_num();
        for(edge_link_map::iterator it=links.begin_edges(); it!=links.end_edges(); ++it)
        {
            reference_summation += it->second.get_references_number();
        }

        if(get<1>(p) < vertices_num)
            get<1>(p) = vertices_num;
        if(get<2>(p) < edges_num)
            get<2>(p) = edges_num;
        if(get<3>(p) < reference_summation)
            get<3>(p) = reference_summation;
    }
}

template<class C, class T> void extract_ETop(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p);

template<class C, class T> void extract_ETop(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p)
{
    leaf_ET allETres;
    n.extract_local_ETop(mesh,allETres);
}

template<class C, class T> void extract_FTop(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p);

template<class C, class T> void extract_FTop(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p)
{
    leaf_ET allETres;
    n.extract_local_FTop(mesh,allETres);
}

template<class C, class T> void extract_VE(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p);

template<class C, class T> void extract_VE(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p)
{
    leaf_Vi allVEres;
    n.extract_local_VE(mesh,allVEres);
}

template<class C, class T> void extract_VF(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p);

template<class C, class T> void extract_VF(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p)
{
    leaf_Vi allVFres;
    n.extract_local_VF(mesh,allVFres);
}

template<class C, class T> void extract_EF(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p);

template<class C, class T> void extract_EF(Node_Stellar &n, Mesh<C,T> &mesh, pair<bool,int> &p)
{
    leaf_ij allEFres;
    n.extract_local_EF(mesh,allEFres);
}

/**
 * @brief A procedure that extracts the degree for all the indexed vertices in a given leaf block
 *
 * @param n the current leaf block
 * @param mesh the mesh indexed
 * @param p an integer vector (with size equal to the vertices number) containing the vertices degree
 */
template<class C, class T> void extract_v_degree(Node_Stellar &n, Mesh<C,T> &mesh, ivect &v_degree);

template<class C, class T> void extract_v_degree(Node_Stellar &n, Mesh<C,T> &mesh, ivect &v_degree)
{
    leaf_VV vvs;
    n.extract_local_VV(mesh,vvs);

    int real_pos;
    for(unsigned i=0; i<vvs.size(); i++)
    {
        int v_id = i + n.get_v_start();
        if(!mesh.is_vertex_removed(v_id))
        {
            real_pos = v_id -1;
            v_degree[real_pos] = vvs[i].size();
        }
    }
}

void extract_p_faces_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<int, bool, unsigned> &t);
void extract_p_faces_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, tuple<int, bool, unsigned> &t);

/**
 *
 *
 */
template<class C, class T> void extract_p_faces(Node_Stellar &n, Mesh<C,T> &mesh, tuple<int, bool, unsigned> &t);

template<class C, class T> void extract_p_faces(Node_Stellar &n, Mesh<C,T> &mesh, tuple<int, bool, unsigned> &t)
{    
    pair<int,leaf_p_faces> p = make_pair(get<0>(t),leaf_p_faces());
    n.extract_local_p_faces(mesh,p);

    if(get<1>(t))
    {
        if(get<2>(t) < p.second.size())
            get<2>(t) = p.second.size();
    }
}


void extract_p_faces_Simplicial_wrapper_v2(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<int, bool, unsigned> &t);
void extract_p_faces_CP_wrapper_v2(Node_Stellar &n, CP_Mesh &mesh, tuple<int, bool, unsigned> &t);

/**
 *
 *
 */
template<class C, class T> void extract_p_faces_v2(Node_Stellar &n, Mesh<C,T> &mesh, tuple<int, bool, unsigned> &t);

template<class C, class T> void extract_p_faces_v2(Node_Stellar &n, Mesh<C,T> &mesh, tuple<int, bool, unsigned> &t)
{
    pair<int,leaf_p_faces> p = make_pair(get<0>(t),leaf_p_faces());
    n.extract_local_p_faces_v2(mesh,p);

    if(get<1>(t))
    {
        if(get<2>(t) < p.second.size())
            get<2>(t) = p.second.size();
    }
}

}



#endif // TOPOLOGICALQUERIES_H
