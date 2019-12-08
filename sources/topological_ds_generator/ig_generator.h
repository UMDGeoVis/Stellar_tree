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

#ifndef IG_GENERATOR_H
#define IG_GENERATOR_H

#include "topological_ds/ig.h"
#include "stellar_tree/mesh.h"
#include "stellar_decomposition/run_iterator.h"

// IMPLEMENTATION NOTES ///
// in the LOCAL approach we have to give a coherent view (i.e., we have to extract the IG every time we visit a top simplex
// conversely, in the GLOBAL approach we can visit each top simplex once, i.e., in the leaf block that index the vertex with the lower position index

class IG_Generator{
public:
    IG_Generator() {}
    IG_Generator(int num_simpl_types)
    {
        stats = IG_stats(num_simpl_types);
    }

    template<class N, class C, class T> void build_IG(N& n, Mesh<C,T> &mesh, bool debug);
    template<class N, class C, class T> void build_IG(N& n, Mesh<C,T> &mesh, IG &ig);

    inline void print_global_stats(IG &ig) { ig.print_stats(stats); }
    inline void print_local_stats() { stats.print_stats(); }

private:
    IG_stats stats;

    template<class N, class C, class T> void build_IG_Leaf(N& n, Mesh<C,T> &mesh, bool debug);
    template<class N, class C, class T> void build_IG_Leaf(N& n, Mesh<C,T> &mesh, IG &ig);

    void init_IG(int dim, int t_id, Simplicial_Mesh &mesh, IG &ig);
    void init_IG(int dim, int t_id, CP_Mesh &mesh, IG &ig);

    template<class T> void visit_sub_simplices(T &t, IG_node *father_node, IG &ig);
    template<class T> void visit_hexahedron(T &t, IG_node *top_node, IG &ig);
    void set_coboundaries(IG_node *father_node);
};

template<class N, class C, class T> void IG_Generator::build_IG(N &n, Mesh<C, T> &mesh, bool debug)
{
    if (n.is_leaf())
    {
        this->build_IG_Leaf(n,mesh,debug);
    }
    else
    {
        for(typename N::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->build_IG(**it,mesh,debug);
        }
    }
}

template<class N, class C, class T> void IG_Generator::build_IG(N &n, Mesh<C, T> &mesh, IG &ig)
{
    if (n.is_leaf())
    {
        this->build_IG_Leaf(n,mesh,ig);
    }
    else
    {
        for(typename N::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->build_IG(**it,mesh,ig);
        }
    }
}

template<class N, class C, class T> void IG_Generator::build_IG_Leaf(N& n, Mesh<C,T> &mesh, bool debug)
{
    IG local_ig = IG(mesh.get_implicitly_and_explicitly_encoded_cells_num());

    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            this->init_IG(d,*t_id,mesh,local_ig);
        }
    }

    if(debug)
    {
        local_ig.update_stats(stats);
    }

    local_ig.clear_IG();
}

template<class N, class C, class T> void IG_Generator::build_IG_Leaf(N& n, Mesh<C,T> &mesh, IG &ig)
{
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            /// in this way we visit each top simplex exactly once
            if(n.indexes_vertex(mesh.get_top_cell(d,*t_id).min_v_index()))
            {
                this->init_IG(d,*t_id,mesh,ig);
            }
        }
    }
}

template<class T> void IG_Generator::visit_sub_simplices(T &t, IG_node *father_node, IG &ig)
{
    IG_node* son_node = new IG_node();
    if(t.get_vertices_num() > 1) /// recursively call the function on the d-1 simplices
    {
        for(int f=0; f<t.get_dfaces_num(); f++)
        {
            ivect verts;
            t.TF(verts,f);
            T son(verts);
            this->visit_sub_simplices(son,son_node,ig);
        }
    }
    else /// t is a vertex
    {
        son_node->set_vertex(t.TV(0));
    }
    ig.add_node(t.get_vertices_num()-1,son_node);

    father_node->set_boundary(son_node);
}

template<class T> void IG_Generator::visit_hexahedron(T &t, IG_node *top_node, IG &ig)
{
    /// quad-faces
    for(int f=0; f<t.get_dfaces_num(); f++)
    {
        IG_node* q_node = new IG_node();
        ivect eids;
        t.edges_in_face(f,eids);
        for(unsigned e=0; e<eids.size(); e++)
        {
            IG_node* e_node = new IG_node();

            ivect edge;
            t.TE(edge,eids[e]);

            for(unsigned v=0; v<edge.size(); v++)
            {
                IG_node* v_node = new IG_node();
                v_node->set_vertex(edge[v]);

                ig.add_node(0,v_node); /// NOTA: we force the positioning
                e_node->set_boundary(v_node);
            }

            ig.add_node(1,e_node); /// NOTA: we force the positioning
            q_node->set_boundary(e_node);
        }

        ig.add_node(2,q_node); /// NOTA: we force the positioning
        top_node->set_boundary(q_node);
    }

    ig.add_node(3,top_node);
}


#endif // IG_GENERATOR_H
