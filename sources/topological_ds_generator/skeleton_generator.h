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

#ifndef SKELETON_GENERATOR_H
#define SKELETON_GENERATOR_H

#include "stellar_tree/stellar_tree.h"
#include "topological_ds/skeleton.h"

namespace skeleton_generator
{

/// this method extracts a "directed" 1-skeleton, in which each edge is encoded only once
/// it is fine for visualization purpose
/// not for computing maximal cliques!!
void extract_global_directed_1skeleton_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, skeleton_graph &graph);
void extract_global_directed_1skeleton_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, skeleton_graph &graph);

template<class C, class T> void extract_global_directed_1skeleton(Node_Stellar &n, Mesh<C,T> &mesh, skeleton_graph &graph);

template<class C, class T> void extract_global_directed_1skeleton(Node_Stellar &n, Mesh<C,T> &mesh, skeleton_graph &graph)
{
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                T& t = mesh.get_top_cell(d,*t_id);
                ivect f;
                for(int v=0; v<t.get_vertices_num(); v++)
                {
                    int v_id = abs(t.TV(v));

                    if(n.indexes_vertex(v_id))
                    {
                        t.TF(f,v);
                        graph.conditional_insert(v_id,f,1);
                    }
                }
            }
        }
    }
}

template<class C, class T> void extract_local_directed_1skeleton(Node_Stellar &n, Mesh<C,T> &mesh, skeleton_graph &graph);

template<class C, class T> void extract_local_directed_1skeleton(Node_Stellar &n, Mesh<C,T> &mesh, skeleton_graph &graph)
{
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                T& t = mesh.get_top_cell(d,*t_id);
                ivect f;
                for(int v=0; v<t.get_vertices_num(); v++)
                {
                    int v_id = abs(t.TV(v));

                    if(n.indexes_vertex(v_id))
                    {
                        t.TF(f,v);
                        graph.conditional_insert(v_id,f,n.get_v_start());
                    }
                }
            }
        }
    }
}

void extract_global_undirected_1skeleton_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, skeleton_graph &graph);
void extract_global_undirected_1skeleton_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, skeleton_graph &graph);

template<class C, class T> void extract_global_undirected_1skeleton(Node_Stellar &n, Mesh<C,T> &mesh, skeleton_graph &graph);

template<class C, class T> void extract_global_undirected_1skeleton(Node_Stellar &n, Mesh<C,T> &mesh, skeleton_graph &graph)
{
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                T& t = mesh.get_top_cell(d,*t_id);

                for(int v=0; v<t.get_vertices_num(); v++)
                {
                    int v_id = abs(t.TV(v));

                    for(int v2=v+1; v2<t.get_vertices_num(); v2++)
                    {
                        int v2_id = abs(t.TV(v2));

                        if(n.indexes_vertex(v_id))
                            graph.insert(v_id-1,v2_id);

                        if(n.indexes_vertex(v2_id))
                            graph.insert(v2_id-1,v_id);
                    }
                }
            }
        }
    }
}

}

#endif // SKELETON_GENERATOR_H
