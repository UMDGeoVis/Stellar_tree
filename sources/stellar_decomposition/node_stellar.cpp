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

#include "node_stellar.h"

void Node_Stellar::compress_top_d_cell_array(int pos, ivect &new_t_list)
{
    sort(new_t_list.begin(),new_t_list.end());

    int count=0;
    int start_t_id = new_t_list[0];

    if(new_t_list.size()==1)
    {
       this->add_top_cell(pos,start_t_id);
       return;
    }
    //otherwise visit new_t_list

    //now obtain the new encoding
    for(unsigned i=0; i<new_t_list.size(); i++)
    {
        if((i+1<new_t_list.size()) && (new_t_list[i]+1) == new_t_list[i+1]) //I have a consecutive range in the t_list
        {
            count++;
        }
        else //found a possible new range
        {
            if(count > 1) //more than two consecutive tetrahedra
            {
                //the range should be from start_t_id to start_t_id + count
                //example: 12 - 13 - 14 - 15
                //is encoded: -12 - 3
                //the first index is implicit
                this->add_top_cell(pos,-start_t_id);
                this->add_top_cell(pos,count);
            }
            else //less or equal to two
            {
                this->add_top_cell(pos,start_t_id);
                if(count==1)
                    this->add_top_cell(pos,start_t_id+count);
            }
            //re-init
            count = 0;
            start_t_id = new_t_list[i+1];
        }
    }
}

void Node_Stellar::extract_top_link(Top_Simplex &t, local_links &links)
{
    ivect &v_arr = t.get_vertices_array();
    ivect e; e.assign(2,0);

    for(int i=0; i<t.get_vertices_num(); i++)
    {
        int v_id = t.TV(i);
        int real_index = abs(v_id);
        if(this->indexes_vertex(real_index))
        {
            int local_pos = real_index - this->get_v_start();

            //(1) we extract the vertices of t that are in the link of v_id
            ivect v_in_v_link(t.get_vertices_num()-1);
            copy_if(v_arr.begin(),v_arr.end(),v_in_v_link.begin(), [v_id](int v){return (v!=v_id);} );
            sort(v_in_v_link.begin(),v_in_v_link.end());
            //(2) we populate the link of v_id
            links.populate_vertex_link(local_pos,v_in_v_link,t);
        }        

        // (3) then we process the edges of t
        for(int j=i+1; j<t.get_vertices_num(); j++)
        {
            t.TE(e,i,j);
            if(this->indexes_vertex(e[1])) // we process an edge only if it has alle the extrema already processed
            {
                //(4) we extract the vertices of t that are in the link of e
                int v_id1 = e[0], v_id2 = e[1];
                ivect v_in_e_link(t.get_vertices_num()-2);
                copy_if(v_arr.begin(),v_arr.end(),v_in_e_link.begin(), [v_id1,v_id2](int v){return (v!=v_id1 && v!=v_id2);} );
                sort(v_in_e_link.begin(),v_in_e_link.end());
                //(5) we populate the link of edge e
                links.populate_edge_link(e,v_in_e_link,t);
            }
        }
    }
}

void Node_Stellar::extract_top_link(Top_CP_Cell &t, local_links &links)
{    
    for(int j=0; j<t.get_vertices_num(); j++)
    {
        int v_id = t.TV(j);
        int real_index = abs(v_id);
        if(this->indexes_vertex(real_index))
        {
            int local_pos = real_index - this->get_v_start();

            //(1) we extract the d-1 faces id that are in the link of v_id
            ivect not_in_v;
            t.faces_NOT_in_vertex(j,not_in_v);

            for(ivect_iter it=not_in_v.begin(); it!=not_in_v.end(); ++it)
            {
                ivect face;
                t.get_d_cell(face,t.get_sub_types_num()-1,*it); // here we have to extrat quads (in 3D) or edges (in 2D)
                //(2) we add the face to the link of v_id
                links.populate_vertex_link(local_pos,face,t);
            }
        }
    }

    // (3) then we process the edges of t
    ivect e;

    for(int j=0; j<t.get_edges_num(); j++)
    {
        t.TE(e,j);

        if(this->indexes_vertex(e[1])) // we process an edge only if it has alle the extrema already processed
        {
            //(4) we extract the vertices of t that are in the link of e
            int v_pos1 = t.vertex_index(e[0]);
            int v_pos2 = t.vertex_index(e[1]);

            ivect not_in_v1;
            t.faces_NOT_in_vertex(v_pos1,not_in_v1);
            sort(not_in_v1.begin(),not_in_v1.end());

            ivect not_in_v2;
            t.faces_NOT_in_vertex(v_pos2,not_in_v2);
            sort(not_in_v2.begin(),not_in_v2.end());

            ivect f_in_link;
            set_intersection(not_in_v1.begin(),not_in_v1.end(),not_in_v2.begin(),not_in_v2.end(),std::back_inserter(f_in_link));

            for(ivect_iter it=f_in_link.begin(); it!=f_in_link.end(); ++it)
            {
                ivect face;
                t.get_d_cell(face,t.get_sub_types_num()-1,*it); // here we have to extrat quads (in 3D) or edges (in 2D)
                //(5) we add the face to the link of edge e
                links.populate_edge_link(e,face,t);
            }
        }
    }
}

void Node_Stellar::extract_top_ETop(int d, int t_id, Top_Simplex &t, leaf_ET &ets, ET &empty_et)
{
    ivect e;

    for(int i=0; i<t.get_vertices_num(); i++)
    {
        for(int j=i+1; j<t.get_vertices_num(); j++)
        {
            t.TE(e,i,j);
            if(this->indexes_vertex(e[1])) // we process an edge only if it has alle the extrema already processed
            {
                leaf_ET::iterator it = ets.find(e);
                if(it!=ets.end())
                {
                    (it->second)[d].push_back(t_id);
                }
                else
                {
                    pair<leaf_ET::iterator,bool> p = ets.insert(make_pair(e,empty_et));
                    ((p.first)->second)[d].push_back(t_id);
                }
            }
        }
    }
}

void Node_Stellar::extract_top_ETop(int d, int t_id, Top_CP_Cell &t, leaf_ET &ets, ET &empty_et)
{
    ivect e;

    // this loop is compatible on both simplicial and CP-complexes
    // but it is inefficient on simplicial ones
    // as we can loop directly on the vertices
    for(int j=0; j<t.get_edges_num(); j++)
    {
        t.TE(e,j);

        if(this->indexes_vertex(e[1])) // we process an edge only if it has alle the extrema already processed
        {
            leaf_ET::iterator it = ets.find(e);
            if(it!=ets.end())
            {
                (it->second)[d].push_back(t_id);
            }
            else
            {
                pair<leaf_ET::iterator,bool> p = ets.insert(make_pair(e,empty_et));
                ((p.first)->second)[d].push_back(t_id);
            }
        }
    }
}

void Node_Stellar::update_vertex_indices(ivect &new_v_indices)
{
    ivect new_v_list;
    for(RunIteratorPair itPair = make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& v_id = itPair.first;
        if(new_v_indices[*v_id-1] != -1)
            new_v_list.push_back(new_v_indices[*v_id-1]);
    }
    this->clear_v_list();
    this->set_v_array(new_v_list);
}

void Node_Stellar::update_and_compress_top_cells_arrays(vector<ivect > &new_top_positions, boost::dynamic_bitset<> &all_deleted)
{
    for(int d=0; d<this->get_num_top_cells_encoded(); d++)
    {
        if(all_deleted[d])
            this->clear_top_d_cells_array(d);
        else
            this->update_and_compress_top_d_cells_array(d,new_top_positions[d]);
    }
}

void Node_Stellar::update_and_compress_top_d_cells_array(int d, ivect &new_indices)
{
    ivect t_list;

    for(RunIteratorPair itPair = this->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        if(new_indices[*t_id-1] != -1) // the top d-cell still exists
        {
            t_list.push_back(new_indices[*t_id-1]);
        }
    }

    this->clear_top_d_cells_array(d);
    if(t_list.size() > 0)
        this->compress_top_d_cell_array(d,t_list);
}

void Node_Stellar::compute_run_histogram(IndexStatistics &is)
{
    for(unsigned d=0; d<top_cells.size(); d++)
    {
        for(unsigned t=0; t<top_cells[d].size(); t++)
        {
            if(top_cells[d][t] <0) //we encounter a run
            {
                int run_size = top_cells[d][t+1] + 1;
                t++;
                is.set_run_histogram(run_size);
            }
            else
            {
                is.set_run_histogram(1);
            }
        }
    }
}

void Node_Stellar::print_top_cells_arrays()
{
    for(int d=0; d<this->get_num_top_cells_encoded();d++)
    {
        this->print_top_cells_array(d);
    }
}

void Node_Stellar::print_top_cells_array(int d)
{
    if(this->get_t_array_size(d)>0)
        cout<<d<<"-position simplex"<<endl;
    for(RunIteratorPair itPair = this->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        cout<<*t_id<<" ";
    }
    if(this->get_t_array_size(d)>0)
        cout<<endl;
}

bool Node_Stellar::check_duplicate_top_cells_array(int d, int t_id)
{
    for(RunIteratorPair itPair = this->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t = itPair.first;
        if(*t==t_id)
            return true;
    }
    return false;
}

//void Node_Stellar::extract_local_p_faces_v2(Simplicial_Mesh &mesh, pair<int,leaf_p_faces> &p)
//{
//    ivect pcell;

//    /// we cannot place d instead of 0 because we do not know here if we are encoding the tops using a verbose encoding or a compact one
//    for(int i=0; i<this->get_num_top_cells_encoded();i++)
//    {
//        for(RunIteratorPair itPair = this->make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
//        {
//            RunIterator const& runIt = itPair.first;
//            Top_Simplex& top = mesh.get_top_cell(i,*runIt);

//            if(top.get_vertices_num() <= p.first)
//                break; /// we stop this for as we are analyzing a top simplex that is lower dimensional than the one that we want to extract

//            ivect top_verts;
//            top.get_positive_sorted_vertices(top_verts);
////            if(this->completely_indexes_top(top))
//            if(this->indexes_vertex(top_verts.front()) && this->indexes_vertex(top_verts.back()))
//            {
//                int num_s = top.get_sub_types_num(p.first);
////                if(num_s == -1)
////                {
////                    break; /// we stop this for as we are analyzing a top simplex that is lower dimensional than the one that we want to extract
////                }
//                /// we extract all the p-faces without checking
//                for(int s=0; s<num_s; s++)
//                {
//                    top.get_d_cell(pcell,p.first,s);
//                    p.second.insert(pcell);
//                    pcell.clear();
//                }
//            }
//            else
//            {
////                ivect top_verts = top.get_vertices_array();
//                Top_Simplex t_aux;
//                for(auto it=top_verts.begin(); it != top_verts.end();)
//                {
//                    int current_vertex = *it;
//                    if(this->indexes_vertex(current_vertex))
//                    {
//                        it = top_verts.erase(it); /// I remove the current vertex
//                        t_aux = Top_Simplex(top_verts);
//                        /// we extract all the p-1 faces of t_aux
//                        /// and then we add the indexed vertex of b
//                        int num_s_1 = t_aux.get_sub_types_num(p.first-1);
////                        if(num_s_1 == -1)
////                        {
////                            break; /// we stop this for as we are analyzing a top simplex that is lower dimensional than the one that we want to extract
////                        }
//                        for(int s=0; s<num_s_1; s++)
//                        {
//                            t_aux.get_d_cell(pcell,p.first-1,s);
//                            pcell.push_back(current_vertex);
//                            sort(pcell.begin(),pcell.end());
//                            p.second.insert(pcell);
//                            pcell.clear();
//                        }
//                    }
//                    /// the array is sorted.. once I found the first vertex on the "right" of the current block I can avoid to iterate
//                    else if(this->is_before_vertex(current_vertex))
//                        break;
//                    else
//                    {
//                        ++it;
//                    }
//                }
//            }

////            for(int s=0; s<num_s; s++)
////            {
////                top.get_d_cell(pcell,p.first,s);
////                if(this->indexes_cell(pcell))
////                {
////                    p.second.insert(pcell);
////                }
////                pcell.clear();
////            }
//        }
//    }
//}

//void Node_Stellar::extract_local_p_faces_v2(CP_Mesh &mesh, pair<int,leaf_p_faces> &p)
//{
//    ivect pcell;

//    /// we cannot place d instead of 0 because we do not know here if we are encoding the tops using a verbose encoding or a compact one
//    for(int i=0; i<this->get_num_top_cells_encoded();i++)
//    {
//        for(RunIteratorPair itPair = this->make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
//        {
//            RunIterator const& runIt = itPair.first;
//            Top_CP_Cell& top = mesh.get_top_cell(i,*runIt);

//            int num_s = top.get_sub_types_num(p.first);
//            if(num_s == -1)
//            {
//                break; /// we stop this for as we are analyzing a top simplex that is lower dimensional than the one that we want to extract
//            }

//            if(this->completely_indexes_top(top))
//            {
//                /// we extract all the p-faces without checking
//                for(int s=0; s<num_s; s++)
//                {
//                    top.get_d_cell(pcell,p.first,s);
//                    p.second.insert(pcell);
//                    pcell.clear();
//                }
//            }
//            else
//            {
//                for(int s=0; s<num_s; s++)
//                {
//                    top.get_d_cell(pcell,p.first,s);
//                    if(this->indexes_cell(pcell))
//                    {
//                        p.second.insert(pcell);
//                    }
//                    pcell.clear();
//                }
//            }
//        }
//    }
//}
