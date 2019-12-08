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

#ifndef IASTAR_GENERATOR_H
#define IASTAR_GENERATOR_H

#include <fstream>
#include <bm/bm.h>
#include <queue>
#include <map>
#include <set>
#include <boost/tuple/tuple.hpp>

#include "stellar_tree/mesh.h"
#include "utilities/timer.h"
#include "stellar_decomposition/run_iterator.h"
#include "stellar_decomposition/node_stellar.h"
#include "utilities/container_utilities.h"
#include "topological_ds/iastar.h"
#include "iastar_statistics.h"

/// this class builds the generalized indexed mesh with adjacencies (IA*)
class IAstar_Generator
{
public:
    IAstar_Generator() {} /// dummy constructor

    /// WRAPPER FUNCTIONS FOR SIMPLICIAL AND CP-MESHES ///
    static void global_generation_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, tuple<IAstar,bool,IAstar_stats> &p)
    { IAstar_Generator::global_generation(n,mesh,get<0>(p),get<1>(p),get<2>(p)); }
    static void local_generation_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, pair<bool,IAstar_stats> &p)
    { IAstar_Generator::local_generation(n,mesh,p.first,p.second); }
    static void global_generation_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<IAstar,bool,IAstar_stats> &p)
    { IAstar_Generator::global_generation(n,mesh,get<0>(p),get<1>(p),get<2>(p)); }
    static void local_generation_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool,IAstar_stats> &p)
    { IAstar_Generator::local_generation(n,mesh,p.first,p.second); }
    /// ****** ///

private:
    template<class N, class C, class T> static void local_generation(N& n, Mesh<C,T> &mesh, bool debug, IAstar_stats &stats);
    template<class N, class C, class T> static void global_generation(N& n, Mesh<C,T> &mesh, IAstar &gia, bool debug, IAstar_stats &stats);

    /// extract topological relations
    template<class N, class C, class T> static void extract_topological_relations_G(N &n, vector<Vertex_Coboundary> &total_vt, vector< face_top_v_vector > &all_faces,
                                                                                    Mesh<C,T> &mesh, IAstar &gia);
    template<class N, class C, class T> static void extract_topological_relations_L(N &n, vector<Vertex_Coboundary> &total_vt, vector< face_top_v_vector > &all_faces,
                                                                                    Mesh<C,T> &mesh, IAstar_local &gia);
    template<class N, class C, class T> static void get_top_topological_relations(N& n, int dim, int t_id, int local_t_id, vector<Vertex_Coboundary> &total_vt,
                                                                           vector< face_top_vector > &all_faces, Mesh<C,T> &mesh);
    template<class N, class C, class T> static void get_top_topological_relations(N& n, int dim, int t_id, vector<Vertex_Coboundary> &total_vt,
                                                                           vector< face_top_vector > &all_faces, Mesh<C,T> &mesh, IAstar &gia);
    ///set the adjacency relations
    template<class I> static void set_adjacencies(vector< face_top_v_vector > &all_faces, I &gia, int &face_count);

    ///set coboundary relations
    // GLOBAL
    template<class N, class C, class T> static void set_vertices_coboundaries_G(N &n, vector<Vertex_Coboundary> &total_vt, Mesh<C,T> &mesh, IAstar &gia);
    static void visit_adj_queue_G(int top_id, int d, int real_v_index, bm::bvector<> &bv, Simplicial_Mesh &mesh, IAstar &gia);
    static void visit_adj_queue_G(int top_id, int d, int real_v_index, bm::bvector<> &bv, CP_Mesh &mesh, IAstar &gia);
    // LOCAL
    template<class N, class C, class T> static void set_vertices_coboundaries_L(N &n, vector<Vertex_Coboundary> &total_vt, Mesh<C,T> &mesh, IAstar_local &gia);
    static void visit_adj_queue_L(int top_id, int d, int real_v_index, bm::bvector<> &bv, Simplicial_Mesh &mesh, IAstar_local &gia);
    static void visit_adj_queue_L(int top_id, int d, int real_v_index, bm::bvector<> &bv, CP_Mesh &mesh, IAstar_local &gia);
};

/// ---- LOCAL VERSION ----
template<class N, class C, class T> void IAstar_Generator::local_generation(N &n, Mesh<C,T> &mesh, bool debug, IAstar_stats &stats)
{
    Timer time;
    //DEBUG
    int face_count = 0;
    //DEBUG

    if(debug)
        time.start();

    IAstar_local local_gia = IAstar_local(n,mesh);

    vector<Vertex_Coboundary> total_vt; total_vt.assign(n.get_v_end()-n.get_v_start(),Vertex_Coboundary(n.get_num_top_cells_encoded()));

    vector< face_top_v_vector > all_faces;
    vector< face_top_vector > tmp_vect; tmp_vect.assign(n.get_v_end()-n.get_v_start(),face_top_vector());
    all_faces.assign(n.get_num_top_cells_encoded(),tmp_vect);

    if(debug)
    {
        time.stop();
#pragma omp atomic
        stats.time_init += time.get_elapsed_time();
        time.start();
    }

    /// (1) first we extract the VT relation for the internal vertices
    ///     and the association between the d-1 cells and the top cells
    IAstar_Generator::extract_topological_relations_L(n,total_vt,all_faces,mesh,local_gia);

    if(debug)
    {
        time.stop();
#pragma omp atomic
        stats.time_gath += time.get_elapsed_time();
        time.start();
    }

    /// two possible adjacency cases
    /// (1) if I have a manifold adjacency, I set the adjacency between the two k-simplices
    /// (2) otherwise, if I have a non-manifold adjacency, we set the co-boundary relation in the k-1 face that is not top,
    ///     then, we save each the position index of that face in each top k-simplex incident in it
    IAstar_Generator::set_adjacencies(all_faces,local_gia,face_count);

    if(debug)
    {
        time.stop();
#pragma omp atomic
        stats.time_tt += time.get_elapsed_time();
        time.start();
    }

    /// now we have set the adjacencies for all top simplices with the minimum vertex indexed by the leaf block
    /// we set the coboundary relation for all the internal vertices, checking the presence of the connected components incident in each vertex
    IAstar_Generator::set_vertices_coboundaries_L(n,total_vt,mesh,local_gia);

    if(debug)
    {
        time.stop();
#pragma omp atomic
        stats.time_vt += time.get_elapsed_time();
    }

    //DEBUG
    if(debug)
    {
#pragma omp critical
        {
            if(stats.max_queue_size < face_count)
                stats.max_queue_size = face_count;
            local_gia.set_stats(stats.max_local_vt,stats.max_local_tt,stats.max_local_nm);
        }
    }
    //DEBUG
}

/// ---- GLOBAL VERSION ----
template<class N, class C, class T> void IAstar_Generator::global_generation(N &n, Mesh<C,T> &mesh, IAstar &gia, bool debug, IAstar_stats &stats)
{
    Timer time;

    if(debug)
    {
        time.start();
    }

    vector<Vertex_Coboundary> total_vt; total_vt.assign(n.get_v_end()-n.get_v_start(),Vertex_Coboundary(n.get_num_top_cells_encoded()));

    vector< face_top_v_vector > all_faces;
    vector< face_top_vector > tmp_vect; tmp_vect.assign(n.get_v_end()-n.get_v_start(),face_top_vector());
    all_faces.assign(n.get_num_top_cells_encoded(),tmp_vect);

    if(debug)
    {
        time.stop();
#pragma omp atomic
        stats.time_init += time.get_elapsed_time();
        time.start();
    }

    /// (1) first we extract the VT relation for the internal vertices
    ///     and the association between the d-1 cells and the top cells
    IAstar_Generator::extract_topological_relations_G(n,total_vt,all_faces,mesh,gia);

    //DEBUG
    int face_count = 0;
    //DEBUG

    if(debug)
    {
        time.stop();
#pragma omp atomic
        stats.time_gath += time.get_elapsed_time();
        time.start();
    }

    /// two possible adjacency cases
    /// (1) if I have a manifold adjacency, I set the adjacency between the two k-simplices
    /// (2) otherwise, if I have a non-manifold adjacency, we set the co-boundary relation in the k-1 face that is not top,
    ///     then, we save each the position index of that face in each top k-simplex incident in it
    IAstar_Generator::set_adjacencies(all_faces,gia,face_count);

    if(debug)
    {
        time.stop();
#pragma omp atomic
        stats.time_tt += time.get_elapsed_time();
        time.start();
    }

    /// now we have set the adjacencies for all top simplices with the minimum vertex indexed by the leaf block
    /// we set the coboundary relation for all the internal vertices, checking the presence of the connected components incident in each vertex
    IAstar_Generator::set_vertices_coboundaries_G(n,total_vt,mesh,gia);

    if(debug)
    {
        time.stop();
#pragma omp atomic
        stats.time_vt += time.get_elapsed_time();
    }

#pragma omp critical
    {
        //DEBUG
        if(stats.max_queue_size < face_count)
            stats.max_queue_size = face_count;
        //DEBUG
    }
}

template<class N, class C, class T> void IAstar_Generator::extract_topological_relations_G(N &n, vector<Vertex_Coboundary> &total_vt,
                                                                                         vector< face_top_v_vector > &all_faces, Mesh<C,T> &mesh, IAstar &gia)
{
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            IAstar_Generator::get_top_topological_relations(n,d,*t_id,total_vt,all_faces[d],mesh,gia);
        }
    }
}

template<class N, class C, class T> void IAstar_Generator::extract_topological_relations_L(N &n, vector<Vertex_Coboundary> &total_vt,
                                                                                         vector< face_top_v_vector > &all_faces, Mesh<C,T> &mesh, IAstar_local &gia)
{
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        int local_t_counter = 1;
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            IAstar_Generator::get_top_topological_relations(n,d,*t_id,local_t_counter,total_vt,all_faces[d],mesh);
            gia.set_index(d,local_t_counter,*t_id);
            local_t_counter++;
        }
    }
}

template<class N, class C, class T> void IAstar_Generator::get_top_topological_relations(N &n, int dim, int t_id, int local_t_id, vector<Vertex_Coboundary> &total_vt,
                                                                                         face_top_v_vector &all_faces, Mesh<C,T> &mesh/*, IAstar_local &local_gia*/)
{
    T& t = mesh.get_top_cell(dim,t_id);
    face_top ft;
    ft.set_t_id(local_t_id);
    int v_id=0;

    for(int j=0; j<t.get_vertices_num(); j++) /// these type of meshes have the same number of vertices and faces
    {
        int real_index = abs(t.TV(j));

        if(n.indexes_vertex(real_index))
        {
            total_vt[real_index-n.get_v_start()].add_incident_top(dim,local_t_id);
        }

        if(j<t.get_dfaces_num()) /// needed for the hexahedra as #F != #V
        {
            t.TF(ft.get_face(),j);
            ft.set_f_pos(j);
            if(n.indexes_cell(ft.get_face(),v_id))
            {
                all_faces[v_id-n.get_v_start()].push_back(ft);
            }
        }
    }
}

template<class N, class C, class T> void IAstar_Generator::get_top_topological_relations(N &n, int dim, int t_id, vector<Vertex_Coboundary> &total_vt,
                                                                                         face_top_v_vector &all_faces, Mesh<C,T> &mesh, IAstar &gia)
{
    T& t = mesh.get_top_cell(dim,t_id);
    face_top ft;
    ft.set_t_id(t_id);
    int v_id=0;

    for(int j=0; j<t.get_vertices_num(); j++) /// these type of meshes have the same number of vertices and faces
    {
        int real_index = abs(t.TV(j));

        if(n.indexes_vertex(real_index))
        {
            total_vt[real_index-n.get_v_start()].add_incident_top(dim,t_id);
        }

        if(j<t.get_dfaces_num()) /// needed for the hexahedra as #F != #V
        {
            /// probably it is an unset adjacency..
            if(gia.get_TT(dim,t_id,j) == UNSETADJ)
            {
                t.TF(ft.get_face(),j);
                ft.set_f_pos(j);

                if(n.indexes_cell(ft.get_face(),v_id))
                {
                    all_faces[v_id-n.get_v_start()].push_back(ft);
                }
            }
        }
    }
}

template<class I> void IAstar_Generator::set_adjacencies(vector< face_top_v_vector > &all_faces, I &gia, int &face_count)
{
    for(unsigned d=0; d<all_faces.size(); d++)
    {
        face_top_v_vector &d_faces = all_faces[d];

        for(unsigned i=0; i<d_faces.size(); i++)
        {
            if(d_faces[i].size() > 0)
            {
                //DEBUG
                face_count += d_faces[i].size();
                //DEBUG

                sort_container(d_faces[i]);

                for(unsigned j=0; j<d_faces[i].size(); j++)
                {
                    if(j+1<d_faces[i].size() && d_faces[i][j] == d_faces[i][j+1])
                    {
                        /// if I have a manifold adjacency, then adj_counter = 1, otherwise, it contains the number of top simplices
                        int adj_counter = 1;
                        while((j+1+adj_counter)<d_faces[i].size() && d_faces[i][j] == d_faces[i][j+1+adj_counter] )
                            adj_counter++;

                        if(adj_counter==1) /// manifold adjacency
                        {
                            #pragma omp critical
                            {
                                gia.set_TT(d,d_faces[i][j].get_t_id(),d_faces[i][j].get_f_pos(),d_faces[i][j+1].get_t_id());
                                gia.set_TT(d,d_faces[i][j+1].get_t_id(),d_faces[i][j+1].get_f_pos(),d_faces[i][j].get_t_id());
                            }
                        }
                        else /// non manifold adjacency
                        {
                            #pragma omp critical
                            {
                                /// the k-1 face must be not top, thus, a simplex not explicitly stored
                                gia.set_nmTT(d_faces[i],j,adj_counter,d);
                            }
                        }

                        j+=adj_counter;
                    }
                    if(j>=d_faces[i].size())
                        break;
                }
            }
        }
    }
}

template<class N, class C, class T> void IAstar_Generator::set_vertices_coboundaries_G(N &n, vector<Vertex_Coboundary> &total_vt, Mesh<C,T> &mesh, IAstar &gia)
{
    for(unsigned i=0; i<total_vt.size(); i++)
    {
        int real_v_index = i+n.get_v_start();
        Vertex_Coboundary &v = total_vt[i];

        for(int d=0; d<v.get_incident_top_types_num(); d++)
        {
            bm::bvector<> bv;
            for(int s=0; s<v.get_incident_top_number(d); s++)
            {
                int top_id = v.get_incident_top(d,s);

                if(mesh.get_top_cell(d,top_id).get_vertices_num() == 2) /// it is an edge
                {
                    #pragma omp critical
                    {
                        /// then simply push in the partial vt of the vertex
                        gia.set_VTstar(n,real_v_index,d,top_id);
                    }
                }
                else
                {
                    if(!bv[top_id])
                    {

                        #pragma omp critical
                        {
                            ///push into the partial vt of the current vertex
                            gia.set_VTstar(n,real_v_index,d,top_id);
                        }

                        /// I must use a queue to check if, with the adjacencies relation, I can reach all the top k-simplices
                        bv.set(top_id);
                        visit_adj_queue_G(top_id,d,real_v_index,bv,mesh,gia);

                        /// if we visit all the top simplices, then we have all the entries in the bit-vector set to 1
                        /// -> then, we have a single connected component
                        /// otherwise, we have some top simplices that belong to a different component
                    }
                }
            }
        }
    }
}



template<class N, class C, class T> void IAstar_Generator::set_vertices_coboundaries_L(N &n, vector<Vertex_Coboundary> &total_vt, Mesh<C,T> &mesh, IAstar_local &gia)
{
    for(unsigned i=0; i<total_vt.size(); i++)
    {
        int real_v_index = i+n.get_v_start();
        Vertex_Coboundary &v = total_vt[i];

        for(int d=0; d<v.get_incident_top_types_num(); d++)
        {
            bm::bvector<> bv;
            for(int s=0; s<v.get_incident_top_number(d); s++)
            {
                int top_id = v.get_incident_top(d,s);

                if(mesh.get_top_cell(d,gia.get_index(d,top_id)).get_vertices_num() == 2) /// it is an edge
                {
                    #pragma omp critical
                    {
                        /// then simply push in the partial vt of the vertex
                        gia.set_VTstar(n,real_v_index,d,top_id);
                    }
                }
                else
                {
                    if(!bv[top_id])
                    {
                        #pragma omp critical
                        {
                            ///push into the partial vt of the current vertex
                            gia.set_VTstar(n,real_v_index,d,top_id);
                        }

                        /// I must use a queue to check if, with the adjacencies relation, I can reach all the top k-simplices
                        bv.set(top_id);
                        visit_adj_queue_L(top_id,d,real_v_index,bv,mesh,gia);

                        /// if we visit all the top simplices, then we have all the entries in the bit-vector set to 1
                        /// -> if so, we have a single connected component
                        /// otherwise, we have some top simplices that belong to a different component
                    }
                }
            }
        }
    }
}



#endif // IASTAR_GENERATOR_H
