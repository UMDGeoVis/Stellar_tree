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

#include "ig_generator.h"

void IG_Generator::init_IG(int dim, int t_id, Simplicial_Mesh &mesh, IG &ig)
{
    Top_Simplex& top = mesh.get_top_cell(dim,t_id);
    IG_node* top_node = new IG_node();

    /// (1) first init the boundary relations
    for(int f=0; f<top.get_dfaces_num(); f++)
    {
        ivect verts;
        top.TF(verts,f);
        Top_Simplex son(verts);
        this->visit_sub_simplices(son,top_node,ig);
    }
    ig.add_node(top.get_vertices_num()-1,top_node);  /// this works only on simplices

    /// (2) then the co-boundary
    this->set_coboundaries(top_node);
}

void IG_Generator::init_IG(int dim, int t_id, CP_Mesh &mesh, IG &ig)
{
    Top_CP_Cell& top = mesh.get_top_cell(dim,t_id);
    IG_node* top_node = new IG_node();

    if(mesh.get_type(dim) == HEX)
    {
        this->visit_hexahedron(top,top_node,ig);
    }
    else if(mesh.get_type(dim) == QUAD)
    {
        for(int f=0; f<top.get_dfaces_num(); f++)
        {
            ivect verts;
            top.TF(verts,f);
            Top_CP_Cell son(verts);
            this->visit_sub_simplices(son,top_node,ig);
        }
        ig.add_node(2,top_node);
    }

    /// (2) then the co-boundary
    this->set_coboundaries(top_node);
}


void IG_Generator::set_coboundaries(IG_node *father_node)
{
    IG_node::IG_nodes_set &b = father_node->get_boundary();
    for(IG_node::IG_nodes_iterator it=b.begin(); it!=b.end(); ++it)
    {
        (*it)->set_coboundary(father_node);
        set_coboundaries(*it);
    }
}
