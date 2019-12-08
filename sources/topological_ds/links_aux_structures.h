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

#ifndef LINKS_AUX_STRUCTURES_H
#define LINKS_AUX_STRUCTURES_H

#include <set>
#include <queue>
#include <map>

#include <utilities/basic_wrappers.h>
#include "stellar_decomposition/top_simplex.h"
#include "stellar_decomposition/top_cp_cell.h"

using namespace std;

/**
 * @brief A class representing the complete link of a cell/simplex
 * The class can represent both the cell and simplices, and in that way we have defined specific insertion procedures
 */
class s_link
{
private:
    iset v_in_link;
    vector< set<ivect > > s_in_link;

public:
    s_link() {}
    s_link(int n) { s_in_link.assign(n,set<ivect >()); }

    inline pair<iset_iter,bool> insert_vertex(int v) { return v_in_link.insert(v); }
    inline void insert_vertices(ivect &v_ids) { v_in_link.insert(v_ids.begin(),v_ids.end()); }

    // the offset = 2 works only on simplices (if we encode all the sub-simplices)
    inline pair<set<ivect>::iterator,bool> insert_simplex(ivect &s, int offset = 2)
    {
        return s_in_link[s.size()-offset].insert(s);
    }
    inline pair<set<ivect>::iterator,bool> insert_cell(ivect &s) // if we are managing a CP-complex
    {
        if(s.size() == 2) // we are inserting an edge
            return s_in_link[0].insert(s);
        else if(s.size() == 4) // we are inserting a quad (can happen only on hexahedral meshes)
            return s_in_link[1].insert(s);
        else
        {
            cerr<<"[ERROR] s_link -> insert_cell undefined behavior"<<endl;
            int a; cin>>a;
            return s_in_link[s.size()-2].insert(s);
        }
    }

    inline iset& get_v_in_link() { return this->v_in_link; }
    inline vector< set<ivect> >& get_s_in_link() { return this->s_in_link; }
    inline set<ivect>& get_s_in_link(int d) { return this->s_in_link[d]; }

    inline int get_v_in_link_size() { return this->v_in_link.size(); }
    inline int get_s_in_link_size() { return this->s_in_link.size(); }
    inline int get_s_in_link_size(int d) { return this->s_in_link[d].size(); }

    inline set<ivect >::iterator get_s_in_link_begin_iterator(int d) { return this->s_in_link[d].begin(); }
    inline set<ivect >::iterator get_s_in_link_end_iterator(int d) { return this->s_in_link[d].end(); }

    inline void clear_s_link() { this->s_in_link.clear(); }
    inline void clear_s_link(int d) { this->s_in_link[d].clear(); }

    // NOTA: v_in_v_link must be sorted
    template<class T> inline void populate_link(ivect &v_in_v_link, T &t)
    {
        //(1) we add the vertices to the link of v_id
        insert_vertices(v_in_v_link);

        if(v_in_v_link.size() > 1)
        {
            //(2) we add the simplex composed by those vertices to the link of v_id
            insert_simplex(v_in_v_link);

            if(v_in_v_link.size() > 2)
            {
                this->extract_subs(v_in_v_link,t);
            }
        }
    }

    void print_link();

    // for statistical purpose
    // the procedure counts the number of integer references
    inline int get_references_number()
    {
        int counter = 0;
        counter += this->v_in_link.size();

        for(vector< set<ivect > >::iterator it=this->s_in_link.begin(); it!=s_in_link.end(); ++it)
            for(set<ivect >::iterator it2=it->begin(); it2!=it->end(); ++it2)
                counter += it2->size();

        return counter;
    }

private:
    void extract_subs(ivect &v_in_v_link, Top_Simplex &);
    void extract_subs(ivect &v_in_v_link, Top_CP_Cell &);
};

typedef map<ivect, s_link> edge_link_map;
typedef vector<s_link> vertices_link_list;

/**
 * @brief A class representing the local links of the vertices and edges encoded in a leaf block
 */
class local_links
{
private:
    vertices_link_list vertices_links;
    edge_link_map edges_links; //used only from simplicial 3-simplices and above

public:

    local_links() {}

    template<class N> inline void init(N &n, int s_type_num)
    {
        int internal_vertices = n.get_v_end()-n.get_v_start();
        vertices_links.assign(internal_vertices,s_link(s_type_num-1)); // -1 to have the correct number of init entries
    }

    inline vertices_link_list::iterator begin_vertices() { return vertices_links.begin(); }
    inline vertices_link_list::iterator end_vertices() { return vertices_links.end(); }
    inline int vertices_num() { return vertices_links.size(); }

    inline s_link& get_vertex_link(int i) { return vertices_links[i]; }
    inline s_link& get_edge_link(ivect &e) { return edges_links[e]; }

    inline void add_v_to_vertex_link(int v_pos, int v_id) { vertices_links[v_pos].insert_vertex(v_id); }
    inline void add_v_ids_to_vertex_link(int v_pos, ivect &v_ids) { vertices_links[v_pos].insert_vertices(v_ids); }
    inline void add_simplex_to_vertex_link(int v_pos, ivect &s) { vertices_links[v_pos].insert_simplex(s); }

    //the below functions are used only from simplicial 3-simplices and above
    inline edge_link_map::iterator begin_edges() { return edges_links.begin(); }
    inline edge_link_map::iterator end_edges() { return edges_links.end(); }
    inline int edges_num() { return edges_links.size(); }

    inline edge_link_map::iterator get_edge_iterator(ivect &e, int s_type_num)
    {
        edge_link_map::iterator it = edges_links.find(e);
        if(it == edges_links.end())
        {
            s_link new_one = s_link(s_type_num);
            pair<edge_link_map::iterator,bool> ret = edges_links.insert(make_pair(e,new_one));
            return ret.first;
        }
        else
            return it;
    }

    inline void add_v_to_edge_link(edge_link_map::iterator &it, int v_id) { it->second.insert_vertex(v_id); }
    inline void add_simplex_to_edge_link(edge_link_map::iterator &it, ivect &s) { it->second.insert_simplex(s); }

    template<class T> inline void populate_vertex_link(int local_pos, ivect &v_in_v_link, T &t){ this->vertices_links[local_pos].populate_link(v_in_v_link,t); }
    template<class T> inline void populate_edge_link(ivect &e, ivect &v_in_e_link, T &t)
    {
        edge_link_map::iterator it = get_edge_iterator(e,t.get_sub_types_num()-1);
        it->second.populate_link(v_in_e_link,t);
    }
};

#endif // LINKS_AUX_STRUCTURES_H
