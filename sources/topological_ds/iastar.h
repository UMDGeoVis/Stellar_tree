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

#ifndef IASTAR_H
#define IASTAR_H

#include "adjacency_aux_structure.h"

#include <stellar_tree/mesh.h>
#include <stellar_decomposition/node_stellar.h>

#include <map>
#include <bm/bm.h>
#include <queue>
#include <set>
#include <fstream>

using namespace std;

/// in the case of a non-manifold adjacency we have to create an entry in face_to_index and in face_coboundary
/// ---> the index of a k-1 face is computed starting from the number of top k-simplex, |Tk|, the number of entries in face_to_index, plus 1
///      |Tk| + |face_to_index| + 1
/// finally it is make negative
///
/// NOTA: this face position_index is encoded as a positive index within each top k-simplex in its co-boundary
///

class IAstar
{
public:

    IAstar() {}
    template<class C, class T> IAstar(Mesh<C,T> &m) // GLOBAL CONSTRUCTOR
    {
        for(int d=0; d<m.get_top_cells_types(); d++)
        {
            vector<Top_Adj> tmp_vect;

            if(m.get_top_cells_num(d) > 0)
            {
                T &t = m.get_top_cell(d,1);
                Top_Adj tmp = Top_Adj(t.get_dfaces_num());
                tmp_vect.assign(m.get_top_cells_num(d),tmp);
            }

            TT_rel.push_back(tmp_vect);
        }
        Vertex_Coboundary v_tmp(m.get_top_cells_types());
        this->VTstar_rel.assign(m.get_vertices_num(),v_tmp);

        this->face_coboundary.assign(m.get_top_cells_types(),map<int,ivect >());
    }


    template<class N> inline void set_VTstar(N &, int v, int dim, int t) { this->VTstar_rel[v-1].add_incident_top(dim,t); }
    inline Vertex_Coboundary& get_VTstar(int pos) { return this->VTstar_rel[pos-1]; }
    inline int get_VTstar_Num() { return this->VTstar_rel.size(); }

    inline void set_nmTT(face_top_vector &tuples, int start, int end, int d)
    {
        int f_id = face_coboundary[d].size() + 1;
        ivect cob; cob.assign(end+1,0);

        for(unsigned i=0; i<cob.size(); i++)
        {
            /// and we set the TT with the negative face index
            this->set_TT(d,tuples[i+start].get_t_id(),tuples[i+start].get_f_pos(),-f_id);
            ///set up co-boundary!
            cob[i] = tuples[i+start].get_t_id();
        }

        face_coboundary[d].insert(make_pair(f_id,cob));
    }
    inline int get_nmTT_Num() { return face_coboundary.size(); }
    inline int get_nmTT_Num(int dim) { return face_coboundary[dim].size(); }
    inline int get_nmTT_Num(int dim, int f_id) { return face_coboundary[dim][abs(f_id)].size(); }
    inline int get_nmTT(int dim, int f_id, int pos) { return face_coboundary[dim][abs(f_id)][pos];}
    inline ivect& get_nmTT(int dim, int f_id) { return face_coboundary[dim][abs(f_id)];}

    inline void set_TT(int dim, int t, int pos, int tt) { this->TT_rel[dim][t-1].set_adj(pos,tt); }
    inline int get_TT(int dim, int t, int pos) { return this->TT_rel[dim][t-1].get_adj(pos); }
    inline Top_Adj& get_TT(int dim, int t) { return this->TT_rel[dim][t-1]; }

    inline int get_TT_Num(int d, int t) { return this->TT_rel[d][t-1].size(); }
    inline int get_TT_Num(int d) { return this->TT_rel[d].size(); }
    inline int get_TT_Num() { return this->TT_rel.size(); }

    inline int is_manifold_adj(int d, int t, int pos) { return (this->TT_rel[d][t-1].get_adj(pos) >= UNSETADJ); }

    inline void check_adj(int adj, int d, iqueue &q, bm::bvector<> &bv)
    {
        if(adj > UNSETADJ && !bv[adj]) /// manifold adjacency
        {
            bv.set(adj);
            q.push(adj);
        }
        else if(adj < UNSETADJ)/// non-manifold adjacency
        {
            for(int ad=0; ad<this->get_nmTT_Num(d,adj); ad++)
            {
                if(!bv[this->get_nmTT(d,adj,ad)])
                {
                    bv.set(this->get_nmTT(d,adj,ad));
                    q.push(this->get_nmTT(d,adj,ad));
                }
            }
        }
    }

    /// this function outputs only the co-boundary and adjacency relations of the IA* data structure
    void save_IAstar(string filename, ivect &original_v_pos, vector<ivect > &original_t_pos);
    /// this function outputs the complete IA* data structure
    /// === BUGGY IMPLEMENTATION TO CHECK ===
//    template<class C, class T> void save_IAstar_full(string filename, Mesh<C,T> &m/*, ivect &original_v_pos, vector<ivect > &original_t_pos*/);
    void print_non_manifold_adjacencies();
    void print_stats();

    template<class C, class T> inline void extract_all_VTop(Mesh<C,T> &mesh)
    {
        for(int v=1; v<=mesh.get_vertices_num(); v++)
        {
            this->VTop(v,mesh);
        }
    }

     /// this function extract the relation R_{0,k}
    template<class C, class T> vector<ivect > VTop(int v_id, Mesh<C,T> &mesh);
    void get_res_and_next_VTop(int v_id, int current, Top_Simplex &top, int d, iqueue &q, bm::bvector<> &bv, VT &res);
    void get_res_and_next_VTop(int v_id, int current, Top_CP_Cell &top, int d, iqueue &q, bm::bvector<> &bv, VT &res);

    template<class C, class T> inline void extract_all_VV(Mesh<C,T> &mesh)
    {
        for(int v=1; v<=mesh.get_vertices_num(); v++)
        {
            this->VV(v,mesh);
        }
    }

    /// this function extract the relation R_{0,0}
    template<class C, class T> iset VV(int v_id, Mesh<C,T> &mesh);
    void get_res_and_next_VV(int pos, int current, Top_Simplex &top, int d, iqueue &q, bm::bvector<> &bv, iset &res);
    void get_res_and_next_VV(int pos, int current, Top_CP_Cell &top, int d, iqueue &q, bm::bvector<> &bv, iset &res);

    ///
    template<class C, class T> void compute_storage(Mesh<C,T> &mesh);

protected:
    /// R0,k
    /// Rk,k
    /// Rk-1,k --> for each k-1 non-manifold simplex (NON top)
    /// Mesh * ???
    ///
    vector<vector<Top_Adj> > TT_rel; /// to init consistently

    vector<Vertex_Coboundary> VTstar_rel;
    vector<map<int,ivect > > face_coboundary;

    void check_adj(int adj, int d, iqueue &q, bm::bvector<> &bv, vector<ivect > &res);
};

template<class C, class T> vector<ivect > IAstar::VTop(int v_id, Mesh<C,T> &mesh)
{
    VT res;
    res.assign(mesh.get_top_cells_types(),ivect());

    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        bm::bvector<> bv;

        for(int k=0; k<this->get_VTstar(v_id).get_incident_top_number(d); k++)
        {
            int vtstar = this->get_VTstar(v_id).get_incident_top(d,k);

            iqueue q;
            q.push(vtstar);

            if(!bv[vtstar]) /// needed only for the edges
            {
                bv.set(vtstar);
                res[d].push_back(vtstar);
            }

            while(!q.empty())
            {
                int current = q.front();

                T& top = mesh.get_top_cell(d,current);

                get_res_and_next_VTop(v_id,current,top,d,q,bv,res);

                q.pop();
            }
        }
    }
    return res;
}

template<class C, class T> iset IAstar::VV(int v_id, Mesh<C,T> &mesh)
{
    iset res;

    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        bm::bvector<> bv;

        for(int k=0; k<this->get_VTstar(v_id).get_incident_top_number(d); k++)
        {
            int vtstar = this->get_VTstar(v_id).get_incident_top(d,k);

            iqueue q;
            q.push(vtstar);

            if(!bv[vtstar]) /// needed only for the edges
            {
                bv.set(vtstar);
            }

            int pos = -1;

            while(!q.empty())
            {
                pos = -1;
                int current = q.front();

                T& top = mesh.get_top_cell(d,current);
                pos = top.vertex_index(v_id);

                get_res_and_next_VV(pos,current,top,d,q,bv,res);

                q.pop();
            }
        }
    }

    return res;
}

template<class C, class T> void IAstar::compute_storage(Mesh<C,T> &mesh)
{
    long sum_v_coboundary = 0, tot_sum_t_adj = 0, tot_sum_nm_adj = 0, tot_sum_t_bound = 0;
    vector<long> sum_t_bound; sum_t_bound.assign(this->get_TT_Num(),0);
    vector<long> sum_t_adj; sum_t_adj.assign(this->get_TT_Num(),0);
    vector<long> sum_nm_adj; sum_nm_adj.assign(this->get_nmTT_Num(),0);

    for(int i=1; i<=this->get_VTstar_Num(); i++)
        sum_v_coboundary += 4 * this->get_VTstar(i).count_entries(); /// 4 * because this are integer references
    cerr<<"TOT storage for vertices coboundary: "<< sum_v_coboundary / (1024.0*1024.0) <<" MBs"<<endl;

    for(int i=0; i<mesh.get_top_cells_types(); i++)
    {
        /// 4 * because this are integer references
        sum_t_bound[i] += 4 * mesh.get_top_cell(i,1).get_vertices_num() * mesh.get_top_cells_num(i); /// TV
        cerr<<"   storage for "<<mesh.print_type(i)<<" boundaries: "<< sum_t_bound[i] / (1024.0*1024.0) <<" MBs"<<endl;
        tot_sum_t_bound += sum_t_bound[i];

        sum_t_adj[i] += 4 * mesh.get_top_cell(i,1).get_dfaces_num() * mesh.get_top_cells_num(i); /// TT
        cerr<<"   storage for "<<mesh.print_type(i)<<" adjacencies: "<< sum_t_adj[i] / (1024.0*1024.0) <<" MBs"<<endl;
        tot_sum_t_adj += sum_t_adj[i];
    }
    cerr<<"TOT storage for boundary relations: "<< tot_sum_t_bound / (1024.0*1024.0) <<" MBs"<<endl;
    cerr<<"TOT storage for adjacency relations: "<< tot_sum_t_adj / (1024.0*1024.0) <<" MBs"<<endl;

    for(int i=0; i<this->get_nmTT_Num(); i++)
    {
        map<int,ivect > &m = this->face_coboundary[i];

        for(map<int,ivect >::iterator it=m.begin(); it!=m.end(); ++it)
        {
            sum_nm_adj[i] += 4 * it->second.size(); /// 4 * because this are integer references
        }
        cerr<<"   storage for i-th nm adjacencies: "<< sum_nm_adj[i] / (1024.0*1024.0) <<" MBs"<<endl;
        tot_sum_nm_adj += sum_nm_adj[i];
    }
    cerr<<"TOT storage for nm adjacencies: "<< tot_sum_nm_adj / (1024.0*1024.0) <<" MBs"<<endl;
    cerr<<"==> TOT storage for IA*: "<< (sum_v_coboundary + tot_sum_t_adj + tot_sum_nm_adj + tot_sum_t_bound) / (1024.0*1024.0) <<" MBs"<<endl;
}

//template<class C, class T> void IAstar::save_IAstar_full(string filename, Mesh<C,T> &m/*, ivect &original_v_pos, vector<ivect > &original_t_pos*/)
//{
////    stringstream ss; ss<<filename<<".ia";
//    ofstream output (filename);
//    output.unsetf( std::ios::floatfield ); // floatfield not set
//    output.precision(15);

//    //number of vertices and number of different types of top cells
//    output << m.get_vertices_num() << " " << m.get_top_cells_types() << endl;

//    map<int,int> pos_to_dim_map;

//    //for each type of top simplex, its real index in the vector of top cells
//    for(int tdim=0; tdim < m.get_top_cells_types(); tdim++)
//    {
//        output << m.get_top_cell(tdim,1).get_vertices_num()-1 << " " << tdim << endl;
//        pos_to_dim_map[tdim] = m.get_top_cell(tdim,1).get_vertices_num()-1;
//    }
//    output << endl;

//    cout<<"vertices"<<endl;

//    //for each vertex
//    for(int vid=1; vid<=m.get_vertices_num(); vid++)
//    {
////        cout<<vid<<endl;

////        auto real_pos = vid+1;
//        Vertex<COORDBASETYPE> v = m.get_vertex(vid);
//        //number of coordinates and values
//        output << v.get_dimension() << " ";
//        for(int j=0; j<v.get_dimension(); j++)
//            output << v.getC(j) << " ";
//        output << endl;
////        cout << v << endl;

//        // partial vertex co-boundary
//        Vertex_Coboundary &vcob = this->get_VTstar(vid);

////        cout<<vcob<<endl;
//        // co-boundary size
//        output << vcob.get_incident_top_types_num() << endl;
//        for(int d=0; d<vcob.get_incident_top_types_num(); d++)
//        {
//            // dimension of top + co-boundary size
//            output << pos_to_dim_map[d] << " " << vcob.get_incident_top_number(d) << " ";
//            // followed by the tops indices
//            for(int i=0; i<vcob.get_incident_top_number(d); i++)
//                output << vcob.get_incident_top(d,i)-1 << " ";
//            output << endl;
//        }
//    }

//    cout<<"workaround"<<endl;

//    // temporary workaroud.. to output something compatible with the IA* available on github
//    // we have to create a mapping for the non-manifold faces
//    int f_counter = 0;
//    vector<map<int,int> > nm_f_map;
//    nm_f_map.assign(m.get_top_cells_types(),map<int,int>());
//    for(int tdim=0; tdim < m.get_top_cells_types(); tdim++)
//    {
////        int key_counter = 1;
//        for(auto tops : this->face_coboundary[tdim])
//        {
//            nm_f_map[tdim].insert(make_pair(tops.first,f_counter));
//            f_counter++;
////            key_counter++;
//        }
////        nm_f_map[tdim].insert(make_pair(this->face_coboundary[tdim]->first,f_counter));
////        f_counter++;
//    }

//    cout<<"tops"<<endl;

//    // for each top cells type
//    for(int tdim=0; tdim < m.get_top_cells_types(); tdim++)
//    {
//        output << m.get_top_cells_num(tdim) << endl;
////        cout << m.get_top_cells_num(tdim) << endl;
//        //for each top cell of that type
//        for(int tid = 1; tid<=m.get_top_cells_num(tdim); tid++)
//        {
//            T &top = m.get_top_cell(tdim,tid);
//            Top_Adj t_adj = this->get_TT(tdim,tid);

//            // workaround: the github implementation sorts the TV list - buggy...
////            for(int i=0; i<top.get_vertices_num(); i++)
////            {
////                for(int j=i+1; j<top.get_vertices_num(); j++)
////                {
////                    if(top.TV(i) > top.TV(j))
////                    {
////                        int tmp = top.TV(i);
////                        top.setTV(i,top.TV(j));
////                        top.setTV(j,tmp);

////                        tmp = t_adj.get_adj(i);
////                        t_adj.set_adj(i,t_adj.get_adj(j));
////                        t_adj.set_adj(j,tmp);
////                    }
////                }
////            }

//            // boundary vertices
//            output << top.get_vertices_num() << endl;
//            for(int i=0; i<top.get_vertices_num(); i++)
//                output << top.TV(i)-1 << " ";
//            output << endl;
//            //adjacencies
//            // temporary workaroud.. to output something compatible with the IA* available on github
//            for(int i=0; i<t_adj.size(); i++)
//            {
//                if(t_adj.get_adj(i) == UNSETADJ) // null adjacency
//                    output << INT_MAX -1 << " ";
//                else if(t_adj.get_adj(i) < 0) // non-manifold adjacency
////                    output << t_adj.get_adj(i) - m.get_top_cells_num(tdim) -1 << " "; // -1 as the the first f_id must be 0
//                    output << nm_f_map[tdim][abs(t_adj.get_adj(i))] << " ";
//                else // manifold adjacency
//                    output << -t_adj.get_adj(i) << " ";
//            }
//            output << endl;
//        }
//    }

//    cout<<"nm"<<endl;

//    int nm_f_counter = 0;
//    for(int tdim=0; tdim < m.get_top_cells_types(); tdim++)
//        for(auto tops : this->face_coboundary[tdim])
//            nm_f_counter++;

//    if(nm_f_counter == 0)
//        output << 0 << endl;
//    else
//    {
//        output << this->face_coboundary.size() << endl;
//        for(int tdim=0; tdim < m.get_top_cells_types(); tdim++)
//        {
//            for(auto tops : this->face_coboundary[tdim])
//            {
//                ivect tvect = tops.second;
//                output << tvect.size() << " ";
//                for(auto t : tvect)
//                    output << t << " ";
//                output << endl;
//            }
//    //        ivect &tops = this->face_coboundary[tdim][fid];
//        }
//    }

//    output.close();
//}

// here for efficiency we use a local indexing for the adjacencies, thus,
// we have to keep the mapping between the local indexing of the leaf block and the global position indexes of the top cells in the mesh
class IAstar_local : public IAstar
{
public:
    IAstar_local() {}
    template<class C, class T> IAstar_local(Node_Stellar &n, Mesh<C,T> &m) // LOCAL CONSTRUCTOR
    {
        VTstar_rel.assign(n.get_v_end()-n.get_v_start(),Vertex_Coboundary(n.get_num_top_cells_encoded()));
        for(int d=0; d<n.get_num_top_cells_encoded(); d++)
        {
            T &t = m.get_top_cell(d,1);
            Top_Adj tmp = Top_Adj(t.get_dfaces_num());
            vector<Top_Adj> tmp_vect; tmp_vect.assign(n.get_real_t_array_size(d),tmp);
            TT_rel.push_back(tmp_vect);

            ivect tmp_ivect; tmp_ivect.assign(n.get_real_t_array_size(d),0);
            this->local_to_global_indexes.push_back(tmp_ivect);
        }
        this->face_coboundary.assign(n.get_num_top_cells_encoded(),map<int,ivect >());
    }

    template<class N> inline void set_VTstar(N &n, int v, int dim, int t) { this->VTstar_rel[v-n.get_v_start()].add_incident_top(dim,t); }
    inline Vertex_Coboundary& get_VTstar(int pos) { return this->VTstar_rel[pos]; }

    inline void set_index(int d, int pos, int id) { this->local_to_global_indexes[d][pos-1] = id; }
    inline int get_index(int d, int pos) { return this->local_to_global_indexes[d][pos-1]; }

    inline void set_stats(int &max_local_vt, ivect &max_local_tt, ivect &max_local_nm)
    {
        if(max_local_vt < this->get_VTstar_Num())
            max_local_vt = this->get_VTstar_Num();

        for(int i=0; i<this->get_TT_Num(); i++)
        {
            if(max_local_tt[i] < this->get_TT_Num(i))
                max_local_tt[i] = this->get_TT_Num(i);
        }

        for(int i=0; i<this->get_nmTT_Num(); i++)
        {
            if(max_local_nm[i] < this->get_nmTT_Num(i))
                max_local_nm[i] = this->get_nmTT_Num(i);
        }
    }

private:
    vector< ivect > local_to_global_indexes; // we must keep the linking between the local indexing and the global

};


#endif // IASTAR_H
