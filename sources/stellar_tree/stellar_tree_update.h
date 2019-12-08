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

template<class C, class T> void Stellar_Tree::compact_vertices_lists(Node_Stellar &n, Mesh<C,T> &mesh, ivect &surviving_vertices)
{
    if (n.is_leaf())
    {
        n.compact_vertices_array(mesh,surviving_vertices);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->compact_vertices_lists(**it,mesh,surviving_vertices);
        }
    }
}

template<class C, class T> int Stellar_Tree::visit_and_unify(Node_Stellar &n, Mesh<C,T> &mesh)
{    
//    cout<<n<<endl;
    int vertex_counter = 0;
    if(n.is_leaf())
    {
        vertex_counter = this->count_indexed_vertices(n,mesh);
//        if(vertex_counter==0)
//        {
//            cout<<"nessun vertice contenuto in "<<n<<endl;
//        }
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
            {
                int local_num_v = this->visit_and_unify(**it,mesh);
                if(local_num_v == 0)
                {
//                    cout<<"nessun vertice contenuto cancello "<<**it<<endl;
                    delete *it;
                    *it = NULL;
                }
                else
                    vertex_counter += local_num_v;
            }
        }

        if(vertex_counter <= this->kv)
        {
            vector<iset> internal_top_cells;
            internal_top_cells.assign(mesh.get_top_cells_types(),iset());
            pair<iset_iter,bool> coppia;

            n.clear_v_list(); // we have to reset the vertices lists as it contains the internal range of vertices
            n.init_sub_cells_vectors(mesh.get_top_cells_types()); // and we have to init the top cells lists

            for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
            {
                if(*it != NULL)
                {
                    Node_Stellar &son = **it;
                    // we reinsert the vertices
                    for(RunIteratorPair itPair = son.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
                    {
                        RunIterator const& v_id = itPair.first;
                        if(!mesh.is_vertex_removed(*v_id))
                            n.add_vertex(*v_id);
                    }

                    // and we reinsert the top cells
                    for(int d=0; d<son.get_num_top_cells_encoded(); d++)
                    {
                        for(RunIteratorPair itPair = son.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
                        {
                            RunIterator const& t_id = itPair.first;
                            if(!mesh.is_top_cell_removed(d,*t_id))
                            {
                                coppia = internal_top_cells[d].insert(*t_id);
                                if(coppia.second == true)
                                {
                                    n.add_top_cell(d,*t_id);
                                }
                            }
                        }
                    }
                }
            }

            n.delete_sons();
        }
    }

    return vertex_counter;
}

template<class C,class T> int Stellar_Tree::count_indexed_vertices(Node_Stellar &n, Mesh<C, T> &mesh)
{
    int internal_vertices = 0;

    for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& v_id = itPair.first;
        if(!mesh.is_vertex_removed(*v_id))
            internal_vertices++;
    }

    return internal_vertices;
}
