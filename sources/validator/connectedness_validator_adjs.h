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

#ifndef CONNECTEDNESS_VALIDATOR_ADJS_H
#define CONNECTEDNESS_VALIDATOR_ADJS_H

#include <queue>
#include <map>
#include <boost/dynamic_bitset.hpp>

#include "stellar_tree/mesh.h"
#include "stellar_decomposition/run_iterator.h"
#include "utilities/container_utilities.h"
#include "topological_ds/adjacency_aux_structure.h"

/// **** ///
typedef queue<int> iqueue;
typedef vector<Top_nmAdj>::iterator iterator_TT;
///

/// class used to check d-1 connecteness ///
/// \brief The local_adj class
class connVal_d_1_adj
{
public:
    template<class N, class C, class T> connVal_d_1_adj(N &n, int d, Mesh<C,T> &m)
    {
        T &t = m.get_top_cell(d,1);

        all_faces.assign(n.get_v_end() - n.get_v_start(),/*tmp_v*/face_top_vector());

        Top_nmAdj tmp = Top_nmAdj(t.get_dfaces_num());
        local_TT.assign(n.get_real_t_array_size(d),tmp);

        local_to_global_indexes.assign(n.get_real_t_array_size(d),0);

        is_visited_top = boost::dynamic_bitset<>(n.get_real_t_array_size(d));
    }

    inline Top_nmAdj& get_adj(int key) { return local_TT[key-1]; }
    inline int get_adj_size() { return this->local_TT.size(); }

    inline bool is_visited(int i) { return is_visited_top[i-1]; } // -1 needed only by boost::dynamic_bitset
    inline void set_visited(int i) { is_visited_top.set(i-1); } // -1 needed only by boost::dynamic_bitset
    inline bool all_visited() { return is_visited_top.count() == is_visited_top.size(); }

    inline void set_index(int pos, int id) { this->local_to_global_indexes[pos-1] = id; }
    inline int get_index(int pos) { return this->local_to_global_indexes[pos-1]; }
    inline int get_index_array_size() { return this->local_to_global_indexes.size(); }

    template<class N, class C, class T> void extract_d_1_faces(N& n, Mesh<C,T>& mesh, int pos);
    void set_adjacencies();

    inline Top_nmAdj convert_to_global_adjacencies(int pos)
    {
        Top_nmAdj adj = this->get_adj(pos);
        for(int f=0; f<adj.get_adj_num(); f++)
        {
            for(int i=0; i<adj.get_adj_num(f); i++)
                adj.update_adj(f,i,this->get_index(adj.get_adj(f,i)));
        }
        return adj;
    }

    inline void print_adjacencies()
    {
        for(unsigned i=0; i<local_TT.size(); i++)
        {
            cout<<i+1<<") ["<<get_index(i+1)<<"] "<<local_TT[i]<<endl;
        }
    }

private:
    boost::dynamic_bitset<> is_visited_top;
    face_top_v_vector all_faces;

    ivect local_to_global_indexes;
    vector<Top_nmAdj> local_TT;

    template<class N, class C, class T> void extract_top_d_1_faces(N &n, int dim, int t_id, int local_t_id, Mesh<C,T>& mesh);
    inline void add_face(int pos, face_top &ft) { all_faces[pos].push_back(ft); }
};

template<class N, class C, class T> void connVal_d_1_adj::extract_d_1_faces(N &n, Mesh<C,T>& mesh, int pos)
{
    int local_t_counter = 1;
    for(RunIteratorPair itPair = n.make_t_array_iterator_pair(pos); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        // we consider a top cell only if it is not flagged for removal
        if(!mesh.is_top_cell_removed(pos,*t_id))
        {
            this->extract_top_d_1_faces(n,pos,*t_id,local_t_counter,mesh);
            this->set_index(local_t_counter,*t_id);
            local_t_counter++;
        }
    }
}

template<class N, class C, class T> void connVal_d_1_adj::extract_top_d_1_faces(N &n, int dim, int t_id, int local_t_id, Mesh<C,T>& mesh)
{
    T& t = mesh.get_top_cell(dim,t_id);

    face_top ft;
    ft.set_t_id(local_t_id);
    int v_id=0;

    for(int j=0; j<t.get_dfaces_num(); j++)
    {
        t.TF(ft.get_face(),j);
        ft.set_f_pos(j);

        if(n.indexes_cell(ft.get_face(),v_id))
        {
            this->add_face(v_id-n.get_v_start(),ft);
        }
    }
}

#endif // CONNECTEDNESS_VALIDATOR_ADJS_H
