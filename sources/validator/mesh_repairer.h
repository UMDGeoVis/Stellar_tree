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

#ifndef MESH_REPAIRER_H
#define MESH_REPAIRER_H

#include "io/writer.h"
#include "stellar_tree/mesh_updater.h"

template<class C, class T> class Mesh_Repairer
{
public:
    Mesh_Repairer() { duplicated_vertices = 0; degenerate_top = 0; duplicated_top = 0; isolated_vertices=0; }

    template<class N> void repair_mesh(N& n, Mesh<C,T>& mesh);
    inline void print_repair_output(string mesh_name, Mesh<C,T>& mesh)
    {
        if(duplicated_vertices == 0 && duplicated_top == 0 && degenerate_top == 0 && isolated_vertices == 0)
            cerr<<"The model was fine.."<<endl;
        else
        {
            cerr<<"The model has :"<<endl;
            cerr<<" "<<duplicated_vertices<<" pairs of duplicated vertices"<<endl;
            cerr<<" "<<isolated_vertices<<" isolated vertices"<<endl;
            cerr<<" "<<degenerate_top<<" degenerate top simplexes"<<endl;
            cerr<<" "<<duplicated_top<<" duplicated top simplexes"<<endl;
            Mesh_Updater mu;
            cerr<<"Fixing mesh..."<<endl;
            mu.clean_mesh(mesh);
            mesh.print_mesh_stats(cerr);
            Writer::write_mesh(mesh_name,"repair",mesh); /// OFF format
        }
    }

private:
    int duplicated_vertices, degenerate_top, duplicated_top, isolated_vertices;

    template<class N> int check_duplicated_vertices(N& n, int v_id, Vertex<C> &v, Mesh<C,T>& mesh);
    int update_top_cells(int v_new, int v_old, int dim, ivect &vt_new, ivect &vt_old, Mesh<C, T> &mesh);
};

template<class C, class T> template<class N> void Mesh_Repairer<C,T>::repair_mesh(N& n, Mesh<C,T>& mesh)
{
    if (n.is_leaf())
    {
        /// (1) get the internal VT relations
//        cout<<"get VTs"<<endl;
        leaf_VT allVTres;
        n.extract_local_VTop(mesh,allVTres);

        /// NOTICE: temporarily disabled steps from 2 to 4
//        /// (2) check for duplicated vertices
////        cout<<"check duplicate vertices"<<endl;
//        for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
//        {
//            if(!mesh.is_vertex_removed(*itPair.first))
//            {
//                Vertex<C> &v = mesh.get_vertex(*itPair.first);
//                int ind = this->check_duplicated_vertices(n,*itPair.first,v,mesh);
//                if(ind != -1)
//                {
//                    mesh.remove_vertex(ind);
//                    this->degenerate_top += update_top_cells(*itPair.first,ind,0,allVTres[*itPair.first-n.get_v_start()][0],allVTres[ind-n.get_v_start()][0],mesh);
//                    this->duplicated_vertices++;
//                }
//            }
//        }
//        /// ------
//        ///
//        /// (3) check for duplicated top simplexes
////        cout<<"check duplicate top cells"<<endl;
//        for(int i=0; i<n.get_num_top_cells_encoded();i++)
//        {
//            for(RunIteratorPair itPair = n.make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
//            {
//                RunIterator const& t_id = itPair.first;
//                T& t = mesh.get_top_cell(i,*t_id);

//                for(RunIteratorPair itPair2 = n.make_t_array_iterator_pair(i); itPair2.first != itPair2.second; ++itPair2.first)
//                {
//                    RunIterator const& t_id2 = itPair2.first;
//                    T& t2 = mesh.get_top_cell(i,*t_id2);

//                    if(*t_id != *t_id2)
//                    {
//                        if(t == t2)
//                        {
//                            mesh.remove_top_cell(i,*t_id2);
//                            this->duplicated_top++;
//                        }
//                    }
//                }
//            }
//        }

//        /// (4) check for degenerate top simplexes
////        cout<<"check degenerate top cells"<<endl;
//        for(int i=0; i<n.get_num_top_cells_encoded();i++)
//        {
//            for(RunIteratorPair itPair = n.make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
//            {
//                RunIterator const& t_id = itPair.first;
//                T& t = mesh.get_top_cell(i,*t_id);

//                if(t.is_degenerate())
//                {
//                    mesh.remove_top_cell(i,*t_id);
//                    this->degenerate_top++;
//                }
//            }
//        }

        /// (5) check for isolated vertices
//        cout<<"check isolated vertices"<<endl;
        for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
        {
            if(!mesh.is_vertex_removed(*itPair.first))
            {
                int local_id = *itPair.first - n.get_v_start();
                int local_counter = 0;
                for (int id = 0; id < allVTres[local_id].size(); id++)
                {
                    local_counter += allVTres[local_id][id].size();
                }
//                cout << local_counter << " [" << (local_counter == 0) << "] " << flush;
                if(local_counter == 0)
                {
                    mesh.remove_vertex(*itPair.first);
                    this->isolated_vertices++;
                }
            }
        }
//        cout<<endl;
    }
    else
    {
        for(typename N::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
            {
                this->repair_mesh(**it,mesh);
            }
        }
    }
}

template<class C, class T> template<class N> int Mesh_Repairer<C,T>::check_duplicated_vertices(N& n, int v_id, Vertex<C> &v, Mesh<C, T> &mesh)
{
    int ind = -1;
    for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        if(!mesh.is_vertex_removed(*itPair.first))
        {
            Vertex<C> &vert = mesh.get_vertex(*itPair.first);
            if(v_id != *itPair.first && v == vert)
            {
                ind = *itPair.first;
                break;
            }
        }
    }
    return ind;
}

template<class C, class T> int Mesh_Repairer<C,T>::update_top_cells(int v_new, int v_old, int dim, ivect &vt_new, ivect &vt_old, Mesh<C, T> &mesh)
{
    int deg_top = 0;

    for(unsigned i=0; i<vt_old.size(); i++)
    {
        /// we update only the tops that are not deleted
        if(!mesh.is_top_cell_removed(dim,vt_old[i]))
        {
            T& t = mesh.get_top_cell(dim,vt_old[i]);

            int old_id = t.vertex_index(v_old);
            int new_id = t.vertex_index(v_new);

            /// if the top simplex is degenerate (i.e. it contains already v_new)
            /// then remove it.. (EXPERIMENTAL)
            if(new_id != -1)
            {
                mesh.remove_top_cell(dim,vt_old[i]);
                deg_top++;
            }
            else
            {
                /// otherwise we set v_new at old_id position
                /// keeping the sign of the index
                if(t.TV(old_id)<0)
                    t.setTV(old_id,-v_new);
                else
                    t.setTV(old_id,v_new);

                vt_new.push_back(vt_old[i]);
            }
        }
    }

    return deg_top;
}

#endif // MESH_REPAIRER_H
