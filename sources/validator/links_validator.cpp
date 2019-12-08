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

#include <utilities/container_utilities.h>
#include "links_validator.h"

void Links_Validator::check_link_conditions(Node_Stellar &n, local_links& links, edge_link_map &faulty_links)
{
    /// visit all the i-simplices
    for(int i=0; i<links.vertices_num(); i++)
    {
        s_link &link = links.get_vertex_link(i);

        if(!Links_Validator::valid_entities_num(link) || !Links_Validator::valid_connected_components_num(link))
        {
            ivect v;
            v.push_back(i+n.get_v_start());
            faulty_links.insert(make_pair(v,link));
        }
    }

    for(edge_link_map::iterator it = links.begin_edges(); it!=links.end_edges(); ++it)
    {
        s_link &link = it->second;

        if(!Links_Validator::valid_entities_num(link) || !Links_Validator::valid_connected_components_num(link))
        {
            faulty_links.insert(make_pair(it->first,link));
        }
    }
}

bool Links_Validator::valid_entities_num(s_link &link)
{
    if(link.get_s_in_link_size() == 1) /// triangle or quad mesh case
        return (link.get_v_in_link_size() == link.get_s_in_link_size(0));
    else if(link.get_s_in_link_size() == 2) /// tetrahedral mesh
    {
        int sum = link.get_v_in_link_size() - link.get_s_in_link_size(0) + link.get_s_in_link_size(1);
        return (sum == 1 || sum == 2);
    }
    else
    {
        cout<<"s_in_link_size: "<<link.get_s_in_link_size()<<endl;
        cout<<"SOMETHING UNKNOWN"<<endl;
        int a; cin>>a;
        return false;
    }
}

bool Links_Validator::valid_connected_components_num(s_link &link)
{
    /// firstly we get the indexed representation of the higher dimensional faces
    vector<Top_nmAdj_and_boundary> indexed_s;
    indexed_s.reserve(link.get_s_in_link_size());

    int max_d = link.get_s_in_link_size() - 1; /// we have to check the connected components on the d-1 simplices

    for(set<ivect>::iterator it=link.get_s_in_link_begin_iterator(max_d); it!=link.get_s_in_link_end_iterator(max_d); ++it)
    {
        Top_nmAdj_and_boundary gt = Top_nmAdj_and_boundary(*it);
        indexed_s.push_back(gt);
    }

    /// then we extract the faces
    face_top_vector all_faces;
    face_top ft;

    int dim = -1, num_sub_simpl = -1;

    for(unsigned i=0; i<indexed_s.size(); i++)
    {
        Top_nmAdj_and_boundary &gt = indexed_s[i];

        ft.set_t_id(i);

        /// here we are guarantee that all the simplex have the same dimension
        if(dim==-1)
            dim = gt.get_sub_types_num() - 1;
        if(num_sub_simpl==-1)
            num_sub_simpl = gt.get_sub_types_num(dim);

        for(int f=0; f<num_sub_simpl; f++)
        {
            gt.get_d_cell(ft.get_face(),dim,f);
            ft.set_f_pos(f);
            all_faces.push_back(ft);
        }
    }

    /// and we pair adjacent simplices
    sort(all_faces.begin(),all_faces.end());
    for(unsigned j=0; j<all_faces.size(); j++)
    {
        if(j+1<all_faces.size() && all_faces[j] == all_faces[j+1])
        {
            /// se ho un'adiacenza manifold adj_counter e' 1, altrimenti contiene il numero di top
            int adj_counter = 1;
            while((j+1+adj_counter)<all_faces.size() && all_faces[j] == all_faces[j+1+adj_counter] )
                adj_counter++;

            if(adj_counter==1)
            {
                indexed_s[all_faces[j].get_t_id()].set_adj(all_faces[j].get_f_pos(),all_faces[j+1].get_t_id());
                indexed_s[all_faces[j+1].get_t_id()].set_adj(all_faces[j+1].get_f_pos(),all_faces[j].get_t_id());
            }
            else
            {
                for(unsigned ad=j; ad<=(j+adj_counter); ad++)
                {
                    for(unsigned bd=ad+1; bd<=(j+adj_counter); bd++)
                    {
                        indexed_s[all_faces[ad].get_t_id()].set_adj(all_faces[ad].get_f_pos(),all_faces[bd].get_t_id());
                        indexed_s[all_faces[bd].get_t_id()].set_adj(all_faces[bd].get_f_pos(),all_faces[ad].get_t_id());
                    }
                }
            }

            j+=adj_counter;
        }
        if(j>=all_faces.size())
            break;
    }

    all_faces.clear();

    /// finally we count the connected components
    boost::dynamic_bitset<> visited_face(indexed_s.size());

    for(unsigned i=0; i<indexed_s.size(); i++)
    {
        if(!visited_face[i])
        {
            iqueue q;
            visited_face.set(i);
            q.push(i);

            while(!q.empty())
            {
                int current = q.front();
                q.pop();

                for(int a=0; a<indexed_s[current].get_adj_num(); a++)
                {
                    for(int j=0; j<indexed_s[current].get_adj_num(a); j++)
                    {
                        if(!visited_face[indexed_s[current].get_adj(a,j)])
                        {
                            visited_face.set(indexed_s[current].get_adj(a,j));
                            q.push(indexed_s[current].get_adj(a,j));
                        }
                    }
                }
            }

            if(visited_face.count() == visited_face.size())
                return true;
            else
            {
                return false;
            }
        }
    }


    cerr<<"[valid_connected_components_num] should never reach this point.."<<endl;
    link.print_link();
    int a; cin>>a;
    return false;
}

void Links_Validator::print_validation_result(edge_link_map &faulty_links, bool debug)
{
    if(faulty_links.empty())
        cerr<<"The model has VALID link conditions"<<endl;
    else
    {
        cerr<<"The model has the following simplices that violate the link conditions: "<<faulty_links.size()<<endl;

        /// for debug only ---- STANDARD OUTPUT
        if(debug)
        {
            for(edge_link_map::iterator it = faulty_links.begin(); it!=faulty_links.end(); ++it)
            {
                const ivect &simplex = it->first;
                s_link &link = it->second;

                cout<<"simplex that violate the link conditions: ";
                for(unsigned i=0; i<simplex.size(); i++)
                    cout<<simplex[i]<<" ";
                cout<<endl;
                link.print_link();
            }
            cout<<endl;
        }
    }
}
