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

#ifndef CONNECTEDNESS_VALIDATOR_H
#define CONNECTEDNESS_VALIDATOR_H

#include <boost/algorithm/minmax.hpp>
#include <algorithm>

#include "stellar_tree/mesh.h"
#include "utilities/timer.h"
#include "io/writer.h"
#include "stellar_decomposition/node_stellar.h"
#include "topological_ds/iastar.h"
#include "connectedness_validator_skeleton.h"
#include "connectedness_validator_adjs.h"
#include "connectedness_validator_stats.h"

class Connectedness_Validator
{
public:
    Connectedness_Validator() {}

    /// WRAPPER FUNCTIONS FOR SIMPLICIAL AND CP-MESHES ////
    static inline void validate_0_connectedness_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, pair<Connectedness_Validator_stats,bool> &params)
    { Connectedness_Validator::check_0_connectness(n,mesh,params.first,params.second); }
    static inline void validate_0_connectedness_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<Connectedness_Validator_stats,bool> &params)
    { Connectedness_Validator::check_0_connectness(n,mesh,params.first,params.second); }
    static inline void validate_0_connectednessVVs_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, pair<Connectedness_Validator_stats,bool> &params)
    { Connectedness_Validator::check_0_connectness_with_VVs(n,mesh,params.first,params.second); }
    static inline void validate_0_connectednessVVs_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<Connectedness_Validator_stats,bool> &params)
    { Connectedness_Validator::check_0_connectness_with_VVs(n,mesh,params.first,params.second); }
    static inline void validate_d_connectedness_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, pair<Connectedness_Validator_stats,bool> &params)
    { Connectedness_Validator::check_d_connectness(n,mesh,params.first,params.second); }
    static inline void validate_d_connectedness_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, pair<Connectedness_Validator_stats,bool> &params)
    { Connectedness_Validator::check_d_connectness(n,mesh,params.first,params.second); }
    /// ****** ///

    template<class C, class T> void validate_connectedness_IAstar(IAstar &gia, Mesh<C,T>& mesh, Connectedness_Validator_stats &stats); /// IA*

private:
    /// FUNCTIONS USED to CHECK CONNECTNESS
    template<class C, class T> static void check_0_connectness(Node_Stellar &n, Mesh<C,T> &mesh, Connectedness_Validator_stats &stats, bool debug);
    template<class C, class T> static void check_0_connectness_with_VVs(Node_Stellar &n, Mesh<C,T> &mesh, Connectedness_Validator_stats &stats, bool debug);
    /// skeleton-based       
    static void visit_skeleton(Node_Stellar& n, connVal_skeleton& skeleton, Connectedness_Validator_stats &stats);
    static void visit_edge_queue(Node_Stellar& n, connVal_skeleton &skeleton, iqueue &q, int &flag, Connectedness_Validator_stats &stats);
    void visit_skeleton(connVal_skeleton& skeleton, Connectedness_Validator_stats &stats); /// IA*
    void visit_edge_queue(connVal_skeleton &skeleton, iqueue &q, int &flag, Connectedness_Validator_stats &stats); /// IA*
    /// VV-based
    static void visit_VVs(Node_Stellar &n, leaf_VV &vvs, Connectedness_Validator_stats &stats);
    template<class C, class T> void visit_VVs(IAstar& gia, Mesh<C,T> &mesh, Connectedness_Validator_stats &stats); /// IA*

    /// functions used to check d-1 connectness
    template<class C, class T> static void check_d_connectness(Node_Stellar &n, Mesh<C,T> &mesh, Connectedness_Validator_stats &stats, bool debug);
    static void visit_local_adjacencies(connVal_d_1_adj& local, Connectedness_Validator_stats &stats);
    static void visit_top_queue(connVal_d_1_adj& local, iqueue &q, int &flag, Connectedness_Validator_stats &stats);
    ///
    template<class C, class T> void visit_adjacencies(IAstar &gia, Mesh<C,T>& mesh, Connectedness_Validator_stats &stats); /// IA*
    template<class C, class T> void visit_top_queue(iqueue &q, int &flag, IAstar &gia, Mesh<C,T>& mesh, boost::dynamic_bitset<> &is_visited, Connectedness_Validator_stats &stats); /// IA*

    static void set_flag(int id, int &flag, Connectedness_Validator_stats &stats);
};

template<class C, class T> void Connectedness_Validator::check_0_connectness(Node_Stellar &n, Mesh<C,T>& mesh, Connectedness_Validator_stats &stats, bool debug)
{
    Timer t;

    if(debug)
        t.start();

    connVal_skeleton skeleton(n.get_v_end() - n.get_v_start());
    skeleton.extract_local_skeleton(n,mesh);

    if(debug)
    {
        t.stop();
#pragma omp critical
        stats.time_extraction += t.get_elapsed_time();
        t.start();
    }

    skeleton.init_edge_vector();

    if(debug)
    {
        t.stop();
#pragma omp critical
        stats.time_init += t.get_elapsed_time();
    }

    if(debug)
        t.start();

#pragma omp critical
    visit_skeleton(n,skeleton,stats);

    if(debug)
    {
        t.stop();
#pragma omp critical
        stats.time_visit += t.get_elapsed_time();
        skeleton.get_stats(stats.edges_num,stats.ve_size);
    }
}



template<class C, class T> void Connectedness_Validator::check_0_connectness_with_VVs(Node_Stellar &n, Mesh<C,T> &mesh, Connectedness_Validator_stats &stats, bool debug)
{
//    cout<<n<<endl;
    Timer t;

    if(debug)
        t.start();

    leaf_VV vvs;
    n.extract_local_VV(mesh,vvs);

    if(debug)
    {
        t.stop();
#pragma omp critical
        stats.time_extraction += t.get_elapsed_time();
        t.start();
    }

#pragma omp critical
    visit_VVs(n,vvs,stats);

    if(debug)
    {
        t.stop();
#pragma omp critical
        stats.time_visit += t.get_elapsed_time();
    }
}

template<class C, class T> void Connectedness_Validator::check_d_connectness(Node_Stellar &n, Mesh<C,T>& mesh, Connectedness_Validator_stats &stats, bool debug)
{
    Timer t;
    if(debug)
        t.start();

    ///
    connVal_d_1_adj local(n,n.get_num_top_cells_encoded() - 1,mesh);

    if(debug)
    {
        t.stop();
#pragma omp critical
        stats.time_init += t.get_elapsed_time();
        t.start();
    }

    local.extract_d_1_faces(n,mesh,n.get_num_top_cells_encoded()-1);

    if(debug)
    {
        t.stop();
#pragma omp critical
        stats.time_extraction += t.get_elapsed_time();
        t.start();
    }

    local.set_adjacencies();

    if(debug)
    {
        t.stop();
#pragma omp critical
        stats.time_simpl_gath += t.get_elapsed_time();
        t.start();
    }

#pragma omp critical
    Connectedness_Validator::visit_local_adjacencies(local,stats);

    if(debug)
    {
        t.stop();
#pragma omp critical
        stats.time_visit += t.get_elapsed_time();
    }
}

/// //////////////////////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////////////
/// IA* code /////////////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////////////

template<class C, class T> void Connectedness_Validator::validate_connectedness_IAstar(IAstar &gia, Mesh<C,T>& mesh, Connectedness_Validator_stats &stats)
{
    Timer t;

    /// (1) check 0-connectedness
    stats.cc.assign(mesh.get_vertices_num(),0);
    t.start();
    this->visit_VVs(gia,mesh,stats);
    t.stop();
    t.print_elapsed_time("[validator] check_connectness (through VVs): ");

    stats.print_validation_result();
    stats.reset();

    /// (2) check d-connectedness and pseudo-manifoldness
    stats.cc.assign(mesh.get_top_cells_num(mesh.get_top_cells_types()-1),0); //init an array for the highest dimensional top simplexes
    t.start();
    visit_adjacencies(gia,mesh,stats);
    t.stop();
    t.print_elapsed_time("[validator] check_d_connecteness: ");

    stats.print_validation_result();
    stats.print_pseudo_manifold_result();
    stats.reset();

    stats.cc.assign(mesh.get_vertices_num(),0);
    t.start();
    connVal_skeleton skeleton = connVal_skeleton(mesh.get_vertices_num());
    skeleton.extract_global_skeleton(mesh);
    skeleton.init_edge_vector();
    visit_skeleton(skeleton,stats);
    t.stop();
    t.print_elapsed_time("[validator] check_connectness: ");

    stats.print_validation_result();
    skeleton.get_stats();
    skeleton.reset();
    stats.reset();
}



template<class C, class T> void Connectedness_Validator::visit_VVs(IAstar &gia, Mesh<C,T> &mesh, Connectedness_Validator_stats &stats)
{
    iqueue q;
    int cc_size = 0;

    for(int v=1; v<=mesh.get_vertices_num(); v++)
    {
        if(stats.cc[v-1]==0)
        {
            stats.cc[v-1]=v;
            stats.connected_components.insert(make_pair(v,deque<int>()));
            cc_size++;

            iset vv = gia.VV(v,mesh);
            for(iset_iter it=vv.begin(); it!=vv.end(); ++it)
            {
                if(stats.cc[*it-1]==0)
                {
                    stats.cc[*it-1]=v;
                    q.push(*it);
                    cc_size++;
                }
            }

            while(!q.empty())
            {
                int v_id = q.front();
                q.pop();

                iset vv2 = gia.VV(v_id,mesh);
                for(iset_iter it2=vv2.begin(); it2!=vv2.end(); ++it2)
                {
                    if(stats.cc[*it2-1]==0)
                    {
                        stats.cc[*it2-1]=v;
                        q.push(*it2);
                        cc_size++;
                    }
                }
            }

            if(cc_size == mesh.get_vertices_num())
                return;
        }
    }
}

template<class C, class T> void Connectedness_Validator::visit_adjacencies(IAstar &gia, Mesh<C,T>& mesh, Connectedness_Validator_stats &stats)
{
    int flag;
    iqueue top_queue;

    int max_dim = mesh.get_top_cells_types() -1;
    boost::dynamic_bitset<> is_visited = boost::dynamic_bitset<>(mesh.get_top_cells_num(max_dim));

    for(int t=1; t<=mesh.get_top_cells_num(max_dim); t++)
    {
        if(!is_visited[t-1])
        {
            if(stats.cc[t-1] == 0)
            {
                flag = t;

                /// every time we start a graph visit we have a new a connected components,
                /// we flag each connected components with an integer, i.e. the index of the first vertex
                stats.connected_components.insert(make_pair(flag,deque<int>()));
            }
            else
            {
                flag = stats.cc[t-1];
            }


            is_visited.set(t-1);
            top_queue.push(t);

            /// visit the queue
            visit_top_queue(top_queue,flag,gia,mesh,is_visited,stats);

            /// if we visit all the edges then we have done within this leaf
            if(is_visited.count() == is_visited.size())
            {
                return;
            }
        }
    }
}

template<class C, class T> void Connectedness_Validator::visit_top_queue(iqueue &q, int &flag, IAstar &gia, Mesh<C,T>& mesh, boost::dynamic_bitset<> &is_visited, Connectedness_Validator_stats &stats)
{
    int d = mesh.get_top_cells_types() - 1;

    while(!q.empty())
    {
        int t = q.front();
        set_flag(t,flag,stats);
        q.pop();

        for(int adj=0; adj<gia.get_TT_Num(d,t); adj++)
        {
            int adj_id = gia.get_TT(d,t,adj);

            if(gia.is_manifold_adj(d,t,adj))
            {
                if(adj_id > UNSETADJ && !is_visited[adj_id-1]) /// exists an adjacent top cell
                {
                    is_visited.set(adj_id-1);
                    q.push(adj_id);
                }
            }
            else
            {
                /// if I have an adjacency with more than two top simplices incident then we do not have a pseudo manifold
                stats.is_pseudo_manifold = false;

                for(int ad=0; ad<gia.get_nmTT_Num(d,adj_id); ad++)
                {
                    int nm_adj_id = gia.get_nmTT(d,adj_id,ad);

                    if(!is_visited[nm_adj_id-1])
                    {
                        is_visited.set(nm_adj_id-1);
                        q.push(nm_adj_id);
                    }
                }
            }
        }
    }
}

#endif // CONNECTEDNESS_VALIDATOR_H
