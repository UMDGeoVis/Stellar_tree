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

#include "stellar_tree.h"

void Stellar_Tree::init_leaves_list(Node_Stellar& n)
{
    if(!n.is_leaf())
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
            {
                if((*it)->is_leaf())
                {
                    this->leaves.push_back(*it);
                }
                else
                {
                    this->init_leaves_list(**it);
                }
            }
        }
    }
}

void Stellar_Tree::init_tops_lists(Node_Stellar& n, int num_top_cells_types)
{
    if(!n.is_leaf())
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
            {
                if((*it)->is_leaf())
                {
                    (*it)->init_sub_cells_vectors(num_top_cells_types);
                }
                else
                {
                    this->init_tops_lists(**it,num_top_cells_types);
                }
            }
        }
    }
}

void Stellar_Tree::get_leaf_indexing_vertex(Node_Stellar &n, int v_id, Node_Stellar *&res)
{
    if (n.is_leaf())
    {
        res = &n;
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL && (*it)->indexes_vertex(v_id)) // NOTA
                this->get_leaf_indexing_vertex(**it,v_id,res);
        }
    }
}

void Stellar_Tree::get_leaves_indexing_cell(Node_Stellar &n, ivect &v_ids, vector<Node_Stellar*> &res)
{
    if (n.is_leaf())
    {
        res.push_back(&n);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL && (*it)->indexes_cell(v_ids)) // NOTA
                this->get_leaves_indexing_cell(**it,v_ids,res);
        }
    }
}

//void Stellar_Tree::get_leaves_indexing_vertices(Node_Stellar &n, ivect &v_ids, unordered_map<int,Node_Stellar*> &res)
//{
//    if (n.is_leaf())
//    {
//        for(ivect_iter it=v_ids.begin(); it!=v_ids.end(); ++it)
//        {
//            if(n.indexes_vertex(*it))
//                res.insert(make_pair(*it,&n));
//        }
//    }
//    else
//    {
//        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
//        {
//            if(*it != NULL && (*it)->indexes_cell(v_ids)) // NOTA
//                this->get_leaves_indexing_vertices(**it,v_ids,res);
//        }
//    }
//}

void Stellar_Tree::compact_and_update_top_cells_lists(Node_Stellar &n, vector<ivect > &new_top_positions, boost::dynamic_bitset<> &all_deleted)
{
    if (n.is_leaf())
    {
        n.update_and_compress_top_cells_arrays(new_top_positions,all_deleted);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->compact_and_update_top_cells_lists(**it,new_top_positions,all_deleted);
        }
    }
}

void Stellar_Tree::compact_and_update_top_d_cells_lists(Node_Stellar &n, pair<int,ivect > &p)
{
    if (n.is_leaf())
    {
        n.update_and_compress_top_d_cells_array(p.first,p.second);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->compact_and_update_top_d_cells_lists(**it,p);
        }
    }
}

void Stellar_Tree::update_tree(Node_Stellar &n, ivect &new_v_positions, vector<ivect > &new_top_positions, boost::dynamic_bitset<> &all_deleted)
{
    if (n.is_leaf())
    {
        n.update_vertex_indices(new_v_positions);
        if(new_top_positions.size()!=0) // if not all the top simplices have been removed
            n.update_and_compress_top_cells_arrays(new_top_positions,all_deleted);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->update_tree(**it,new_v_positions,new_top_positions,all_deleted);
        }
    }
}

//void Stellar_Tree::clear_top_d_cells_lists(Node_Stellar &n, int d)
//{
//    if (n.is_leaf())
//    {
//        n.clear_top_d_cells_array(d);
//    }
//    else
//    {
//        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
//        {
//            if(*it != NULL)
//                this->clear_top_d_cells_lists(**it,d);
//        }
//    }
//}

//void Stellar_Tree::clear_top_cells_lists(Node_Stellar &n, boost::dynamic_bitset<> &all_deleted)
//{
//    if (n.is_leaf())
//    {
//        n.clear_top_cells_lists(all_deleted);
//    }
//    else
//    {
//        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
//        {
//            if(*it != NULL)
//                this->clear_top_cells_lists(**it,all_deleted);
//        }
//    }
//}
