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

#ifndef CONNECTEDNESS_VALIDATOR_SKELETON_H
#define CONNECTEDNESS_VALIDATOR_SKELETON_H

#include <vector>
#include <map>
#include <boost/dynamic_bitset.hpp>

#include "stellar_tree/mesh.h"
#include "stellar_decomposition/run_iterator.h"

using namespace std;

/// ***** ///
/// \brief The edge_skeleton class
///
class connVal_skeleton
{
public:
    connVal_skeleton(int num_v) { local_ve.assign(num_v,ivect()); }

    inline void init_edge_vector()
    {
        indexed_edges.assign(edge_to_id_map.size(),ivect());
        is_visited_edge = boost::dynamic_bitset<>(edge_to_id_map.size(),0);

        for(map<ivect,int>::iterator it=edge_to_id_map.begin(); it!=edge_to_id_map.end(); ++it)
        {
            indexed_edges[it->second] = it->first;
        }

        edge_to_id_map.clear();
    }

    inline pair<map<ivect, int>::iterator, bool> insert_edge(ivect &e)
    {
        return edge_to_id_map.insert(make_pair(e,edge_to_id_map.size()));
    }

    inline long long get_edges_num() { return indexed_edges.size(); }
    inline ivect& get_edge(int i) { return indexed_edges[i]; }

    inline bool is_visited(int i) { return is_visited_edge[i]; }
    inline void set_visited(int i) { is_visited_edge[i] = 1; }
    inline bool all_visited() { return is_visited_edge.count() == is_visited_edge.size(); }

    inline void add_to_ve(int v_pos, int e_id) { local_ve[v_pos].push_back(e_id); }
    inline int get_ve_num() { return local_ve.size(); }
    inline ivect& get_ve(int i) { return local_ve[i]; }

    /// **** ///
    template<class N, class C, class T> void extract_local_skeleton(N& n, Mesh<C,T>& mesh);
    template<class C, class T> void extract_global_skeleton(Mesh<C,T>& mesh); /// IA*
    /// **** ///

    inline void reset()
    {
        this->edge_to_id_map.clear();
        this->indexed_edges.clear();
        this->is_visited_edge.reset(0);
        this->local_ve.clear();
    }

    inline void get_stats()
    {
        cerr<<"edges number: "<<this->get_edges_num()<<" ("<<((4*2*this->get_edges_num()) / (1024.0*1024.0))<<" MBs)"<<endl;
        long long counter = 0;
        for(int i=0; i<this->get_ve_num(); i++)
            counter += this->get_ve(i).size();
        cerr<<"VE entries: "<<counter<<" ("<<((4*counter) / (1024.0*1024.0))<<" MBs)"<<endl;
        cerr<<"==> TOT skeleton storage: "<<((4*2*this->get_edges_num()+4*counter) / (1024.0*1024.0))<<" MBs"<<endl;
    }

    inline void get_stats(int &edges_num, int &ve_size)
    {
        if(edges_num < this->get_edges_num())
            edges_num = this->get_edges_num();
        int counter = 0;
        for(int i=0; i<this->get_ve_num(); i++)
            counter += this->get_ve(i).size();
        if(ve_size < counter)
            ve_size = counter;
    }

    inline void print_skeleton(int start)
    {
        cout<<"local skeleton:"<<endl;
        for(unsigned i=0; i<indexed_edges.size(); i++)
        {
            cout<<i<<"] ";
            for(unsigned j=0; j<indexed_edges[i].size(); j++)
                cout<<indexed_edges[i][j]<<" ";
            cout<<endl;
        }
        cout<<"local VE:"<<endl;
        for(unsigned i=0; i<local_ve.size(); i++)
        {
            cout<<i+start<<"] ";
            for(ivect_iter it=local_ve[i].begin(); it!=local_ve[i].end(); ++it)
                cout<<*it<<" ";
            cout<<endl;
        }
        int a; cin>>a;
    }

private:
    map<ivect,int> edge_to_id_map;
    vector<ivect > indexed_edges;
    boost::dynamic_bitset<> is_visited_edge;
    vector<ivect > local_ve; /// from vector to set..

    /// **** ///
    template<class N> void extract_top_skeleton(N &n, Top_CP_Cell &t);
    template<class N> void extract_top_skeleton(N &n, Top_Simplex &t);
    template<class N> void check_edge(N &n, ivect &edge);
    void extract_top_skeleton(Top_CP_Cell &t); /// IA*
    void extract_top_skeleton(Top_Simplex &t); /// IA*
    void check_edge(ivect &edge); /// IA*
    /// **** ///
};

template<class N, class C, class T> void connVal_skeleton::extract_local_skeleton(N &n, Mesh<C,T>& mesh)
{
    for(int d=0; d<n.get_num_top_cells_encoded(); d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            extract_top_skeleton(n,mesh.get_top_cell(d,*t_id));
        }
    }
}

template<class N> void connVal_skeleton::extract_top_skeleton(N &n, Top_CP_Cell &t)
{
    ivect edge; edge.assign(2,-1);

    for(int i=0; i<t.get_edges_num(); i++)
    {
        t.TE(edge,i);
        check_edge(n,edge);
    }
}

template<class N> void connVal_skeleton::extract_top_skeleton(N &n, Top_Simplex &t)
{
    ivect edge; edge.assign(2,-1);
    // extraction of all the edges of a top simplex
    // not good-looking but it is the most efficient way to extract these iteratively
    for(int j=0; j<t.get_vertices_num(); j++)
    {
        for(int i=j+1; i<t.get_vertices_num(); i++)
        {
            t.TE(edge,j,i);
            check_edge(n,edge);
        }
    }
}

template<class N> void connVal_skeleton::check_edge(N &n, ivect &edge)
{
    if(n.indexes_vertex(edge[0]))
    {
        pair<map<ivect, int>::iterator, bool> ret = this->insert_edge(edge);
        if(ret.second) /// inserted
        {
            this->add_to_ve(edge[0]-n.get_v_start(),ret.first->second);
            if(n.indexes_vertex(edge[1]))
                this->add_to_ve(edge[1]-n.get_v_start(),ret.first->second);
        }
    }
}

/// IA* ///
template<class C, class T> void connVal_skeleton::extract_global_skeleton(Mesh<C,T>& mesh)
{
    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        for(int i=1; i<=mesh.get_top_cells_num(d); i++)
        {
            T &t = mesh.get_top_cell(d,i);
            extract_top_skeleton(t);
        }
    }
}

#endif // CONNECTEDNESS_VALIDATOR_SKELETON_H
