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

#ifndef LINKS_VALIDATOR_H
#define LINKS_VALIDATOR_H

/// ONLY FOR TETRAHEDRAL, TRIANGLE AND QUAD MESHES ///
#include <boost/algorithm/minmax.hpp>

#include "stellar_tree/mesh.h"
#include "stellar_decomposition/node_stellar.h"
#include "topological_ds/links_aux_structures.h"
#include "topological_ds/adjacency_aux_structure.h"

class Links_Validator
{
public:
    Links_Validator() {  }

    inline static void validate_links_CP_wrapper(Node_Stellar& n, CP_Mesh& mesh, edge_link_map &faulty_links)
    {
        Links_Validator::validate_link_conditions(n,mesh,faulty_links);
    }
    inline static void validate_links_Simplicial_wrapper(Node_Stellar& n, Simplicial_Mesh& mesh, edge_link_map &faulty_links)
    {
        Links_Validator::validate_link_conditions(n,mesh,faulty_links);
    }
    static void print_validation_result(edge_link_map &faulty_links, bool debug);

private:

    template<class C, class T> static void validate_link_conditions(Node_Stellar& n, Mesh<C,T>& mesh, edge_link_map &faulty_links);

    template<class C, class T> static void get_links_from_tet(Node_Stellar &n, int dim, int t_id, local_links& links, Mesh<C,T>& mesh);
    template<class C, class T> static void get_links_from_tri(Node_Stellar &n, int dim, int t_id, local_links& links, Mesh<C,T>& mesh);
    template<class C, class T> static void get_links_from_quad(Node_Stellar &n, int dim, int t_id, local_links& links, Mesh<C,T>& mesh);

    static void check_link_conditions(Node_Stellar &n, local_links& links, edge_link_map &faulty_links);
    static bool valid_entities_num(s_link &link);
    static bool valid_connected_components_num(s_link &link);
};


template<class C, class T> void Links_Validator::validate_link_conditions(Node_Stellar& n, Mesh<C,T>& mesh, edge_link_map &faulty_links)
{
    local_links links;

    links.init(n,mesh.get_top_cell(0,1).get_sub_types_num());

    for(int d=0; d<n.get_num_top_cells_encoded(); d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            //I consider the current tetra only if it is not removed
            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                if(mesh.get_type(d)==TRIANGLE)
                    Links_Validator::get_links_from_tri(n,d,*t_id,links,mesh);
                if(mesh.get_type(d)==TETRA)
                    Links_Validator::get_links_from_tet(n,d,*t_id,links,mesh);
                if(mesh.get_type(d)==QUAD)
                    Links_Validator::get_links_from_quad(n,d,*t_id,links,mesh);
            }
        }
    }

#pragma omp critical
    {
        //controllo le relazione di connettività fra le entità della foglia
        Links_Validator::check_link_conditions(n,links,faulty_links);
    }
}

template<class C, class T>
void Links_Validator::get_links_from_tet(Node_Stellar &n, int dim, int t_id, local_links& links, Mesh<C,T> &mesh)
{
    T& t = mesh.get_top_cell(dim,t_id);

    ivect f;
    ivect e;
    ivect edge_key;

    for(int j=0; j<t.get_vertices_num(); j++)
    {
        int real_v_index = abs(t.TV(j));
        int v_pos = real_v_index - n.get_v_start();

        if(n.indexes_vertex(real_v_index))
        {
            //popolo l'insieme delle facce che sono nel link di un vertice
            t.TF(f,j);
            links.add_simplex_to_vertex_link(v_pos,f);


            for(int v=1; v<t.get_vertices_num(); v++)
            {
                int other_v = abs(t.TV((j+v)%t.get_vertices_num()));
                //popolo l'insieme dei vertici che sono nel link di un vertice
                links.add_v_to_vertex_link(v_pos,other_v);

                //popolo l'insieme degli edge che sono nel link di un vertice
                for(int v1=v+1; v1<t.get_vertices_num(); v1++)
                {
                    e.clear();
                    e.push_back(other_v);
                    e.push_back(abs(t.TV((j+v1)%t.get_vertices_num())));
                    sort(e.begin(),e.end());
                    links.add_simplex_to_vertex_link(v_pos,e);
                }

                //recupero le informazioni necessarie al controllo sugli edge
                //considero solo gli edge che hanno il vertice interno alla foglia con indice minore
                //in questo modo controllo ogni edge una sola volta

                if(real_v_index < other_v)
                {
                    edge_key.clear();
                    edge_key.push_back(real_v_index);
                    edge_key.push_back(other_v);

                    int other_two_v[2];
                    int count=0;

                    for(int v1=1; v1<4; v1++)
                    {
                        int current = abs(t.TV((j+v1)%4));
                        if(current != other_v)
                        {
                            other_two_v[count] = current;
                            count++;
                        }
                    }

                    //inserisco i quattro vertici del tetraedro (da fare meglio)
                    edge_link_map::iterator it = links.get_edge_iterator(edge_key,t.get_sub_types_num()-1);
                    links.add_v_to_edge_link(it,real_v_index);
//                    edge_link_map::iterator it = links.add_v_to_edge_link(edge_key,real_v_index,t.get_sub_types_num()-1);
                    links.add_v_to_edge_link(it,other_v);
                    links.add_v_to_edge_link(it,other_two_v[0]);
                    links.add_v_to_edge_link(it,other_two_v[1]);

                    //inserisco i quattro edge incidenti in uno dei due estremi
                    for(int i=0; i<2; i++)
                    {
                        e.clear();
                        e.push_back(real_v_index);
                        e.push_back(other_two_v[i]);
                        sort(e.begin(),e.end());
                        links.add_simplex_to_edge_link(it,e);

                        e.clear();
                        e.push_back(other_v);
                        e.push_back(other_two_v[i]);
                        sort(e.begin(),e.end());
                        links.add_simplex_to_edge_link(it,e);
                    }

                    //inserisco l'edge opposto a quello di riferimento
                    e.clear();
                    e.push_back(other_two_v[0]);
                    e.push_back(other_two_v[1]);
                    sort(e.begin(),e.end());
                    links.add_simplex_to_edge_link(it,e);

                    //inserisco le due facce del tetra incidenti nei due vertici che compongono l'edge ma non nell'edge stesso
                    e.clear();
                    e.push_back(real_v_index);
                    e.push_back(other_two_v[0]);
                    e.push_back(other_two_v[1]);
                    sort(e.begin(),e.end());
                    links.add_simplex_to_edge_link(it,e);

                    e.clear();
                    e.push_back(other_v);
                    e.push_back(other_two_v[0]);
                    e.push_back(other_two_v[1]);
                    sort(e.begin(),e.end());
                    links.add_simplex_to_edge_link(it,e);
                }
            }
        }
    }
}

template<class C, class T>
void Links_Validator::get_links_from_tri(Node_Stellar &n, int dim, int t_id, local_links& links, Mesh<C,T> &mesh)
{
    T& t = mesh.get_top_cell(dim,t_id);

    ivect f;

    for(int j=0; j<t.get_vertices_num(); j++)
    {
        int real_v_index = abs(t.TV(j));
        int v_pos = real_v_index - n.get_v_start();

        if(n.indexes_vertex(real_v_index))
        {
            //popolo l'insieme degli edge che sono nel link di un vertice
            t.TF(f,j);
            links.add_simplex_to_vertex_link(v_pos,f);

            for(int v=1; v<t.get_vertices_num(); v++)
            {
                int other_v = abs(t.TV((j+v)%t.get_vertices_num()));
                //popolo l'insieme dei vertici che sono nel link di un vertice
                links.add_v_to_vertex_link(v_pos,other_v);
            }
        }
    }
}

template<class C, class T>
void Links_Validator::get_links_from_quad(Node_Stellar &n, int dim, int t_id, local_links& links, Mesh<C,T> &mesh)
{
    T& t = mesh.get_top_cell(dim,t_id);

    ivect e;

    for(int j=0; j<t.get_vertices_num(); j++)
    {
        int real_v_index = abs(t.TV(j));
        int v_pos = real_v_index - n.get_v_start();

        if(n.indexes_vertex(real_v_index))
        {
            for(int v=1; v<t.get_vertices_num(); v++)
            {
                int other_v = abs(t.TV((j+v)%t.get_vertices_num()));
                //popolo l'insieme dei vertici che sono nel link di un vertice
                links.add_v_to_vertex_link(v_pos,other_v);
            }

            //popolo l'insieme degli edge che sono nel link di un vertice
            t.TF(e,j+1);
            links.add_simplex_to_vertex_link(v_pos,e);
//            e.clear();
            t.TF(e,j+2);
            links.add_simplex_to_vertex_link(v_pos,e);
        }
    }
}

#endif // LINKS_VALIDATOR_H
