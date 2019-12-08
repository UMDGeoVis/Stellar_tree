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

#include "reindexer.h"

void Reindexer::compress_tree_representation_vertices(Node_Stellar &n, bool save_original_positions)
{
    if (n.is_leaf())
    {
        int start = vertCounter;
        for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& v_id = itPair.first;
            new_v_positions[*v_id-1] = vertCounter;

            if(save_original_positions)
                original_v_positions[vertCounter - 1] = *v_id;

            vertCounter++;
        }
        int end = vertCounter;
        n.clear_v_list();
        n.set_v_range(start,end);
    }
    else
    {        
        int start = vertCounter;
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                compress_tree_representation_vertices(**it,save_original_positions);
        }
        int end = vertCounter;
        n.set_v_range(start,end);
    }
}

void Reindexer::compress_tree_representation_top(Node_Stellar &n)
{
    if (n.is_leaf())
    {
        compress_tree_representation_top_leaf(n);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                compress_tree_representation_top(**it);
        }
    }
}

void Reindexer::compress_tree_representation_top_leaf(Node_Stellar& n)
{
    for(int i=0; i<n.get_num_top_cells_encoded(); i++)
    {
         ivect new_t_list;
         for(RunIteratorPair itPair = n.make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
         {
             RunIterator const& t_id = itPair.first;
             new_t_list.push_back(new_tops_positions[i][*t_id-1]);
         }
         n.clear_top_d_cells_array(i);

         if(new_t_list.size()>0)
         {
             n.compress_top_d_cell_array(i,new_t_list);
         }
    }
}

void Reindexer::get_top_tuple_key(Node_Stellar &n, ivect &v_ids, ivect &key)
{
    if (n.is_leaf())
    {
        // when I'm entering a leaf the v of a top are for sure indexed..
        key.push_back(n.get_v_start());
        // clean-up phase for v_ids array
        for(ivect_iter it=v_ids.begin(); it!=v_ids.end();)
        {
            if(n.indexes_vertex(*it))
            {
                v_ids.erase(it);
            }
            else
                ++it;
        }
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
            {
                for(ivect_iter itv=v_ids.begin(); itv!=v_ids.end(); ++itv)
                {
                    if((*it)->indexes_vertex(*itv))
                    {
                        get_top_tuple_key(**it,v_ids,key);
                        break;
                    }
                }
            }
        }
    }
}

void Reindexer::get_tops_reordered_indexes(bool save_original_positions)
{
    for(unsigned i=0; i<new_tops_positions.size(); i++)
    {
        if(new_tops_positions[i].size() > 0)
        {
            ivect I;
            I.assign(leaf_tuples_array[i].size(),0);
            int counter = 1;


            //get the prefix sum of the leaf_tuple_array counts by iterating it
            //this will give us, for each grouped set of top cells, the initial index of this group in the reindexed array
            for(map<ivect,pair<int,int> >::iterator it=leaf_tuples_array[i].begin(); it!=leaf_tuples_array[i].end(); ++it)
            {
                I[it->second.first] = counter;
                counter += it->second.second;
            }
            leaf_tuples_array[i].clear();

            //we updated the values in reordered_tops_ids_array setting the new index value for each top simplex
            for(unsigned j=0; j<new_tops_positions[i].size(); j++)
            {
                int leaf_key = new_tops_positions[i][j];
                new_tops_positions[i][j] = I[leaf_key];

                if(save_original_positions)
                    original_tops_positions[i][I[leaf_key]-1] = j+1; // WARNING KEEP TRACK OF THE INDICES!!

                I[leaf_key]++;
            }
        }
    }
}
