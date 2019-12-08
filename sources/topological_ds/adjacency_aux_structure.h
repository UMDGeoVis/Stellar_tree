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

#ifndef ADJACENCY_AUX_STRUCTURE_H
#define ADJACENCY_AUX_STRUCTURE_H

#include <vector>
#include <algorithm>

#include <stellar_decomposition/top_simplex.h>
#include <utilities/basic_wrappers.h>

using namespace std;

#define UNSETADJ 0

/// a container used to extract the adjacency relation for a generic top simplex
class face_top
{
protected:
    /// the generic face
    ivect f; /// nota the face index must be sorted using the ascending ordering
    /// the top index
    int t;
    /// the face position in t
    ushort f_pos;

public:
    face_top() { t = 0; f_pos = 0; }

    inline ivect& get_face() { return f; }
    inline int get_first_face_vertex() { return f[0]; }

    inline int get_t_id() { return t; }
    inline ushort get_f_pos() { return f_pos; }

    inline void set_and_sort_face(ivect &s)
    {
        f = s;
        sort(f.begin(),f.end());
    }
    inline void set_t_id(int t) { this->t = t; }
    inline void set_f_pos(ushort p) { this->f_pos = p; }

    inline bool has_not(int v_ind)
    {
        for(unsigned i=0; i<f.size(); i++)
            if(v_ind == f[i])
                return false;
        return true;
    }

    inline bool operator<(const face_top& p) const { return f < p.f; }
    inline bool operator==(const face_top& s) const { return f==s.f; }
    inline bool operator!=(const face_top& s) const { return !(*this==s); }

    inline friend ostream& operator<<(ostream& out, const face_top& a)
    {
        out <<"F[ ";
        for(ivect_const_iter i=a.f.begin(); i!=a.f.end(); ++i)
            out << *i << " ";
        out << "] T[ "<< a.t <<" ] is in position: "<<a.f_pos;
        return out;
    }
};

typedef vector<face_top> face_top_vector;
typedef vector<face_top_vector> face_top_v_vector; // the vector associate each tuple to the list associated to the vertex with the minimum position index

/// a class that stores the adjacency relations of a generic top simplex/cell
/// it can be used in both manifold a non-manifold models
/// ================== NOTA =============================
/// this class is sufficient for both IA and IA* data structures
/// =====================================================
class Top_Adj
{
protected:
    ivect adj;

public:
    Top_Adj() {}
    Top_Adj(const Top_Adj &orig) { this->adj = orig.adj; }
    Top_Adj(int adj_num) { adj.assign(adj_num,UNSETADJ); }/// zero and not -1 because the negative indexes refer to a face of the top simplex
    inline void set_adj(int pos, int t_id) { adj[pos] = t_id; }
    inline int get_adj(int pos) const { return adj[pos]; }
    inline int size() { return this->adj.size(); }

    inline friend ostream& operator<<(ostream& out, const Top_Adj& a)
    {
        for(unsigned i=0; i<a.adj.size(); i++)
            out<<a.get_adj(i)<<" ";
        return out;
    }
};

/// VERBOSE encoding of the adjacencies (o non-manifold structures, i.e. on manifold one has the same requirements of Top_adj)
/// use it when the storage requirements are not an issue
class Top_nmAdj
{
protected:
    vector<ivect> adj;

public:
    Top_nmAdj() {}
    Top_nmAdj(int adj_num) { adj.assign(adj_num,ivect()); } // zero and not -1 because the negative indexes refer to a face of the top simplex

    inline void set_adj(int pos, int t_id) { adj[pos].push_back(t_id); }
    inline ivect& get_adj(int fpos) { return adj[fpos]; }
    inline int get_adj(int fpos, int pos) { return adj[fpos][pos]; }
    inline int get_adj_num() { return adj.size(); }
    inline int get_adj_num(int fpos) { return adj[fpos].size(); }

    inline void update_adj(int fpos, int pos, int t_id) { adj[fpos][pos] = t_id; }
    inline void update_unset_adj(int fpos, ivect_iter it) { adj[fpos].erase(it); }
    inline bool unset_adj(int fpos) { return (adj[fpos].size()==0); }
    inline bool has_unset_adj()
    {
        for(unsigned d=0; d<adj.size();d++)
        {
            if(unset_adj((d)))
                return true;
        }
        return false;
    }

    inline void merge(Top_nmAdj &v)
    {
        for(unsigned d=0; d<adj.size();d++)
        {
            if(unset_adj((d)) && !v.unset_adj(d))
            {
                for(int i=0; i<v.get_adj_num(d); i++)
                    this->set_adj(d,v.get_adj(d,i));
            }
        }
    }

    inline friend ostream& operator<<(ostream& out, const Top_nmAdj& a)
    {
        out<<"ADJ: ";
        for(unsigned d=0; d<a.adj.size();d++)
        {
            out<<"F"<<d<<"[";
            for(unsigned i=0; i<a.adj[d].size(); i++)
            {
                out<<a.adj[d][i];
                if(i+1 < a.adj[d].size())
                    out<<" ";
            }
            out<<"] ";
        }
        return out;
    }
};

/// class used during validation (during link link condition)
class Top_nmAdj_and_boundary : public Top_Simplex, public Top_nmAdj
{
public:
    Top_nmAdj_and_boundary() {}
    Top_nmAdj_and_boundary(const ivect &vertices) { this->vertices = vertices; adj.assign(this->vertices.size(),ivect()); }

    inline void init(const ivect &vertices) { this->vertices = vertices; adj.assign(this->vertices.size(),ivect()); }

    inline friend ostream& operator<<(ostream& out, const Top_nmAdj_and_boundary& p)
    {
        out << "s_dim["<<p.get_vertices_num()<<"] v[";
        for(int i=0; i<p.get_vertices_num(); i++)
        {
            out << p.TV(i) << " ";
        }
        out << "] ";

        out<<"ADJ: ";
        for(unsigned d=0; d<p.adj.size();d++)
        {
            out<<"F"<<d<<"[";
            for(unsigned i=0; i<p.adj[d].size(); i++)
            {
                out<<p.adj[d][i];
                if(i+1 < p.adj[d].size())
                    out<<" ";
            }
            out<<"] ";
        }
        return out;
    }

private:
};

class Vertex_Coboundary
{
private:
    VT partial_VT;

public:
    Vertex_Coboundary(int num_top) { partial_VT.assign(num_top,ivect()); }
    Vertex_Coboundary(const Vertex_Coboundary &orig) { this->partial_VT = orig.partial_VT; }

    inline void add_incident_top(int dim, int t_id) { partial_VT[dim].push_back(t_id); }
    inline int get_incident_top(int dim, int pos) { return partial_VT[dim][pos]; }
    inline int get_incident_top_types_num() { return partial_VT.size(); }
    inline int get_incident_top_number(int dim) { return partial_VT[dim].size(); }
    inline vector<ivect >& get_partial_VT() { return this->partial_VT; }

    inline bool non_manifold_vertex()
    {
        int counter = 0;
        for(unsigned i=0; i<partial_VT.size(); i++)
        {
            if(partial_VT[i].size() > 1) /// in dimension i is non-manifold
                return true;
            else if(partial_VT[i].size() == 1)
                counter++;
        }

        if(counter > 1) /// I have a vertex with more than one top simplex type incident in it
            return true;

        return false;
    }

    inline int count_entries()
    {
        int counter = 0;
        for(unsigned i=0; i<partial_VT.size(); i++)
            counter += partial_VT[i].size();
        return counter;
    }

    inline friend ostream& operator<<(ostream& out, const Vertex_Coboundary& a)
    {
        for(unsigned d=0; d<a.partial_VT.size();d++)
        {
            for(unsigned i=0; i<a.partial_VT[d].size(); i++)
            {
                out<<a.partial_VT[d][i];
                if(i+1 < a.partial_VT[d].size())
                    out<<" ";
            }
            out<<" # ";
        }

        return out;
    }
};

#endif // ADJACENCY_AUX_STRUCTURE_H
