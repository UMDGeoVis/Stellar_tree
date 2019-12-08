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

#ifndef HALFEDGE_GENERATOR_H
#define HALFEDGE_GENERATOR_H

#include <vector>
#include <map>
#include <boost/dynamic_bitset.hpp>
#include <queue>

#include "stellar_tree/mesh.h"
#include "stellar_decomposition/node_stellar.h"
#include "utilities/basic_wrappers.h"
#include "topological_ds/adjacency_aux_structure.h"
#include "topological_ds/halfedge.h"

using namespace std;

class HalfEdge_Generator
{
public:
    HalfEdge_Generator() { max_half_edges = 0; max_half_vertices = 0; max_half_faces = 0;} /// for LOCAl approach
    template<class C, class T> HalfEdge_Generator(Mesh<C,T>& mesh) /// for GLOBAL approach
    {
        this->half_vertices.assign(mesh.get_vertices_num(),-1);
        this->half_faces.assign(mesh.get_top_cells_num(0),-1);
        max_half_edges = 0;
        max_half_vertices = 0;
        max_half_faces = 0;
        max_edges_list_size = 0;
    }

    template<class M> void build_HalfEdge(Node_Stellar& n, M &mesh);
    template<class M> void build_global_HalfEdge(Node_Stellar& n, M &mesh, bool debug);

    void print_half_edges();
    void print_stats();

private:
    /// these variables are used by the GLOBAL approach only
    vector<HalfEdge> half_edges;
    ivect half_vertices;
    ivect half_faces;

    /// for stats gathering
    int max_half_edges, max_half_vertices, max_half_faces, max_edges_list_size;

    template<class C, class T> void build_HalfEdge_leaf(Node_Stellar& n, Mesh<C,T> &mesh);
    template<class M> void build_global_HalfEdge_visit(Node_Stellar& n, M &mesh);
    template<class C, class T> void build_global_HalfEdge_leaf(Node_Stellar& n, Mesh<C,T> &mesh);

    inline void set_half_edge(HalfEdge &he, int v, int t_id, int pos, int &current_id, int num_faces)
    {
        he.set_pointed_vertex(v); /// set pointed vertex
        he.set_reference_face(t_id); /// set reference face

        current_id = half_edges.size();

        /// this if...else works only on the global approach.. and in the local if we consider all the edges of the top simplices..
        /// but if we discretize the external edges.. then the indices obtained can point to an uninitialized entry..
        if(pos+1 == num_faces)
            he.set_next_half_edge(current_id-pos);
        else
            he.set_next_half_edge(current_id+1);

        half_edges.push_back(he);
    }

    inline void set_local_half_vertex(Node_Stellar &n, int v, int current_id)
    {
        /// set an internal vertex from which the half edge starts..
        if(n.indexes_vertex(v) && half_vertices[v-n.get_v_start()] == -1)
            half_vertices[v-n.get_v_start()] = current_id;
    }
    inline void set_global_half_vertex(Node_Stellar &n, int v, int current_id)
    {
        /// set an internal vertex from which the half edge starts..
        if(n.indexes_vertex(v) && half_vertices[v-1] == -1)
            half_vertices[v-1] = current_id;
    }

    inline bool is_unset_half_face(int t_id) { return (half_faces[t_id-1] == -1); }
    inline void set_half_face(int t_id, int current_id) { half_faces[t_id-1] = current_id; }
    inline int get_half_face(int t_id) { return half_faces[t_id-1]; }

    inline void add_half_top(face_top_vector &all_edges, ivect &edge, /*int t_id, */int current_id)
    {
        face_top ht;
        ht.set_and_sort_face(edge);
        ht.set_t_id(current_id);
        all_edges.push_back(ht);
    }

    inline void set_opposite_edges(face_top_vector &all_edges)
    {
        sort(all_edges.begin(),all_edges.end());

        for(unsigned i=0; i<all_edges.size(); i++)
        {
            if(i+1<all_edges.size() && all_edges[i] == all_edges[i+1])
            {
                half_edges[all_edges[i].get_t_id()].set_opposite_half_edge(all_edges[i+1].get_t_id());
                half_edges[all_edges[i+1].get_t_id()].set_opposite_half_edge(all_edges[i].get_t_id());
                i++;
            }

            if(i>=all_edges.size())
                break;
        }
    }

    void fix_local_half_edges_ordering();
    void fix_global_half_edges_ordering();
    void fix_half_edges_ordering(iqueue &he_queue, boost::dynamic_bitset<> &visited, unsigned &count);

    void test_local_half_edge_structure();
    void test_global_half_edge_structure();
    void test_half_edge_structure(iqueue &he_queue, boost::dynamic_bitset<> &visited);
    void print_test_result(int conn_components, boost::dynamic_bitset<> &visited);

    inline void reset_variables()
    {
        half_edges.clear();
        half_vertices.clear();
        half_faces.clear();
    }
};

/// LOCAL VERSION ///
template<class M>  void HalfEdge_Generator::build_HalfEdge(Node_Stellar &n, M &mesh)
{
    if (n.is_leaf())
    {
        this->build_HalfEdge_leaf(n,mesh);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->build_HalfEdge(**it,mesh);
        }
    }
}

template<class C, class T>  void HalfEdge_Generator::build_HalfEdge_leaf(Node_Stellar &n, Mesh<C, T> &mesh)
{
    /// we use some global variable for local purpose
    this->half_vertices.assign(n.get_v_end()-n.get_v_start(),-1);

    /// for top faces we cannot use the global variables because the range is not consecutive
    map<int,int> local_half_faces;

    ivect edge_vector;
    HalfEdge current_he;
    int current_id;

    face_top_vector all_edges;

    int num_faces = -1;

    /// extract
    for(int d=0; d<n.get_num_top_cells_encoded(); d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                T &top = mesh.get_top_cell(d,*t_id);

                if(num_faces == -1)
                    num_faces = top.get_dfaces_num();

                for(int i=0; i<top.get_dfaces_num(); i++)
                {
                    top.TF_unsorted(edge_vector,i);

                    set_half_edge(current_he,edge_vector.back(),*t_id,i,current_id,num_faces);
                    set_local_half_vertex(n,edge_vector.front(),current_id);

                    /// if the current simplex is not set then add the current half edge to it
                    local_half_faces.insert(make_pair(*t_id,current_id));

                    add_half_top(all_edges,edge_vector,current_id);
                }
            }
        }
    }

    set_opposite_edges(all_edges);

    fix_local_half_edges_ordering();

    /// gathering execution statistics
    if((int)half_edges.size() > max_half_edges)
        max_half_edges = half_edges.size();
    if((int)half_vertices.size() > max_half_vertices)
        max_half_vertices = half_vertices.size();
    if((int)local_half_faces.size() > max_half_faces)
        max_half_faces = local_half_faces.size();
    if((int)all_edges.size() > max_edges_list_size)
        max_edges_list_size = all_edges.size();

    /// resetting the global variables used locally
    reset_variables();
}

/// GLOBAL VERSION ///
template<class M>  void HalfEdge_Generator::build_global_HalfEdge(Node_Stellar &n, M &mesh, bool debug)
{
    build_global_HalfEdge_visit(n,mesh);
    fix_global_half_edges_ordering();
    
    /// gathering execution statistics
    max_half_edges = half_edges.size();
    max_half_vertices = half_vertices.size();
    max_half_faces = half_faces.size();

    ///for debug
    if(debug)
    {
        test_global_half_edge_structure();
    }

    reset_variables();
}
template<class M>  void HalfEdge_Generator::build_global_HalfEdge_visit(Node_Stellar &n, M &mesh)
{
    if (n.is_leaf())
    {
        this->build_global_HalfEdge_leaf(n,mesh);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->build_global_HalfEdge_visit(**it,mesh);
        }
    }
}

template<class C, class T>  void HalfEdge_Generator::build_global_HalfEdge_leaf(Node_Stellar &n, Mesh<C, T> &mesh)
{
    ivect edge_vector;
    HalfEdge current_he;
    int current_id;

    face_top_vector all_edges;

    int num_faces = -1;

    /// extract
    for(int d=0; d<n.get_num_top_cells_encoded(); d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            if(!mesh.is_top_cell_removed(d,*t_id))
            {                
                T &top = mesh.get_top_cell(d,*t_id);

                if(num_faces == -1)
                    num_faces = top.get_dfaces_num();

                if(this->is_unset_half_face(*t_id))
                {
                    for(int i=0; i<top.get_dfaces_num(); i++)
                    {
                        top.TF_unsorted(edge_vector,i);

                        set_half_edge(current_he,edge_vector.back(),*t_id,i,current_id,num_faces);
                        set_global_half_vertex(n,edge_vector.front(),current_id);

                        /// set it one time
                        if(i==0)
                            set_half_face(*t_id,current_id);

                        ///I put into the local queue
                        if(n.indexes_cell(edge_vector))
                            add_half_top(all_edges,edge_vector,current_id);
                    }
                }
                else
                {
                    int next = this->get_half_face(*t_id);

                    do
                    {
                        HalfEdge &he = half_edges[next];

                        if(he.get_opposite_half_edge() == -1 && n.indexes_vertex(he.get_pointed_vertex()))
                        {
                            /// we need to extract the edge encoded by the current half edge (the one which has the second entry = to he.pointed_vertex
                            for(int i=0; i<top.get_dfaces_num(); i++)
                            {
                                top.TF_unsorted(edge_vector,i);
                                if(edge_vector[1] == he.get_pointed_vertex())
                                    break;
                            }

                            add_half_top(all_edges,edge_vector,next);
                        }

                        next = he.get_next_half_edge();
                    }
                    while(this->get_half_face(*t_id) != next);
                }
            }
        }
    }

    set_opposite_edges(all_edges);

    if((int)all_edges.size() > max_edges_list_size)
        max_edges_list_size = all_edges.size();
}

#endif // HALFEDGE_GENERATOR_H
