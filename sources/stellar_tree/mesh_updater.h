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

#ifndef MESH_UPDATER_H
#define MESH_UPDATER_H

#include "mesh.h"

/**
 * @brief The Mesh_Updater class represents an interface for updating a simplicial (or CP) mesh.
 *
 */
class Mesh_Updater
{
public:
    Mesh_Updater() {}

    /**
     * @brief A public method that compresses the indexed mesh representation
     * The procedure removes the vertices and top cells flagged as deleted during a simplification procedure
     *
     * NOTA: during the procedure are also updated the boundary relations of the top cells
     *
     * @param mesh a Mesh& representing the indexed mesh
     */
    template<class C, class T> void clean_mesh(Mesh<C,T> &mesh);
    /**
     * @brief A public method that compresses the vertices array of the mesh
     * The procedure removes the vertices flagged as deleted during a simplification procedure
     * and returns the new position indexes of the vertices not deleted during the simplification.
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers that at the end of the procedure contains the new position indexes of the vertices
     *
     * NOTA: the size of the new_v_positions array is equal to the size of the original array of the mesh, and the
     * entries corresponding to deleted elements are flagged with -1 value
     *
     */
    template<class C, class T> void clean_vertices_array(Mesh<C,T> &mesh, ivect &new_v_positions);
    /**
     * @brief A public method that compresses the vertices array of the mesh
     * The procedure removes the vertices flagged as deleted during a simplification procedure
     * and returns the new position indexes of the vertices not deleted during the simplification.
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers that at the end of the procedure contains the new position indexes of the vertices
     * @param surviving_vertices a vector of integers containing the non-deleted vertices
     *
     * NOTA: the size of the new_v_positions array is equal to the size of the original array of the mesh, and the
     * entries corresponding to deleted elements are flagged with -1 value.
     */
    template<class C, class T> void clean_vertices_array(Mesh<C,T> &mesh, ivect &new_v_positions, ivect &surviving_vertices);
    /**
     * @brief A public method that compresses the top cells arrays of the mesh
     * The procedure removes the top cells flagged as deleted during a simplification procedure
     * and returns, for each top k-cell array, if it has been completely erased.
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     * @param new_top_positions a vector of vectors that, at the end of the procedure, contains the new position indexes of the top cells
     * @param simplification_counters a vector referring to the number of top cells removed during a simplification procedure
     * @return a bit-vector sized as the number of top cells encoded, returning the top cells type that have been completely erased by a simplification procedure
     *
     * NOTA: the size of the new_v_positions and new_top_positions arrays is equal to the size of the original arrays of the mesh, and the
     * entries corresponding to deleted elements are flagged with -1 value.
     *
     * NOTA2: during the procedure are also updated the boundary relations of the top cells
     */
    template<class C, class T> boost::dynamic_bitset<> update_and_clean_top_cells_arrays(Mesh<C,T> &mesh, ivect &new_v_positions,
                                                                                         vector<ivect> &new_top_positions, ivect &simplification_counters);

    /**
     * @brief A public method that updates the boundary relations of the top cells, by assigning the new position indexes of the corresponding vertices,
     * and resorting the vertices array in the mesh in order to be coherent to the new vertices ordering (calls reorder_vertices_array)
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     */
    template<class C, class T> void update_tops_boundary(Mesh<C,T> &mesh, ivect &new_v_positions);

    /**
     * @brief A private method that update the vertices array of the mesh by assigning the new position indexes.
     * The procedure execute the swaps inline an thus it does not require extra storage.
     * The number of swap operation is exactly |V|, where |V| is the number of vertices of the mesh
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers that at the end of the procedure contains the new position indexes of the vertices
     */
    template<class C, class T> void reorder_vertices_array(Mesh<C,T> &mesh, ivect &new_v_positions);
    /**
     * @brief A public method that update the top cells arrays of the mesh by assigning the new position indexes.
     * The procedure execute the swaps inline an thus it does not require extra storage.
     * The number of swap operation is exactly |T|, where |T| is the number of top cells of the mesh
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_top_positions a vector of vectors that, at the end of the procedure, contains the new position indexes of the top cells
     */
    template<class C, class T> void reorder_tops_array(Mesh<C,T> &mesh, vector<ivect> &new_tops_positions);

private:
    /**
     * @brief A public method that compresses the top d-cells arrays of the mesh.
     * The procedure removes the top d-cells flagged as deleted during a simplification procedure
     * and returns if the corresponding array has been completely erased.
     * @param mesh a Mesh& representing the indexed mesh
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     * @param new_top_positions a vector that, at the end of the procedure, contains the new position indexes of the top d-cells
     * @param d an integer representing the dimension of the top d-cells in the mesh arrays
     * @param counter an integer representing the number of the top d-cells deleted during the simplification procedure
     * @return a boolean, true if all the top d-cells have been deleted, false otherwise
     */
    template<class C, class T> bool update_and_clean_top_d_cells_array(Mesh<C,T> &mesh, ivect &new_v_positions, ivect &new_top_positions, int d, int counter);
    /**
     * @brief A public method that compresses the top d-cells array of the mesh.
     * The procedure removes the top d-cell deleted during a simplification procedure
     * and returns the new position indexes of the top d-cells not deleted during the simplification.
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param d an integer representing the dimension of the top d-cells in the mesh arrays
     * @param new_d_top_positions a vector of integers that at the end of the procedure contains the new position indexes of the top d-cells
     *
     * NOTA: the size of the new_d_top_positions array is equal to the size of the original array of the mesh, and the
     * entries corresponding to deleted elements are flagged with -1 value.
     */
    template<class C, class T> void clean_top_d_cells_array(Mesh<C,T> &mesh, int d, ivect &new_d_top_positions);
    /**
     * @brief A public procedure that updates the boudary relation for all the top d-cells of the mesh
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param d an integer representing the dimension of the top d-cells in the mesh arrays
     * @param new_v_positions a vector of integers containing the updated position indexes of the vertices
     */
    template<class C, class T> void update_d_tops_boundary(Mesh<C,T> &mesh, int d, ivect &new_v_positions);
    /**
     * @brief A public procedure that coherently reorder all the top d-cells of the mesh
     *
     * @param mesh a Mesh& representing the indexed mesh
     * @param d an integer representing the dimension of the top d-cells in the mesh arrays
     * @param new_top_positions a vector that, at the end of the procedure, contains the new position indexes of the top d-cells
     */
    template<class C, class T> void reorder_d_tops_array(Mesh<C,T> &mesh, int d, ivect new_d_tops_positions);
};

template<class C, class T> void Mesh_Updater::clean_mesh(Mesh<C,T> &mesh)
{
    int v_counter = 1, v = 0; // v_counter keeps the new vertex indexing while v is needed to set the new index
    ivect new_v_positions;
    new_v_positions.assign(mesh.get_vertices_num(),-1);

    for(typename Mesh<C,T>::v_array_iter it=mesh.v_array_begin(); it!=mesh.v_array_end(); )
    {
        if(mesh.is_vertex_removed(*it))
        {
            it = mesh.erase_vertex(it);
        }
        else
        {
            new_v_positions[v] = v_counter;
            v_counter++;
            ++it;
        }
        v++;
    }

    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        for(typename Mesh<C,T>::top_d_array_iter it=mesh.top_d_array_begin(d); it!=mesh.top_d_array_end(d); )
        {
            if(mesh.is_top_cell_removed(*it))
            {
                it = mesh.erase_top_cell(d,it);
            }
            else
            {
                for(int v=0; v<it->get_vertices_num(); v++)
                    it->setTV(v, new_v_positions[abs(it->TV(v))-1]);
                ++it;
            }
        }
    }
}

template<class C, class T> void Mesh_Updater::clean_vertices_array(Mesh<C,T> &mesh, ivect &new_v_positions)
{
    int v_counter = 1, v = 0; /// vertexCounter keeps the new vertex indexing while v is needed to set the new index
    new_v_positions.assign(mesh.get_vertices_num(),-1);

    for(typename Mesh<C,T>::v_array_iter it=mesh.v_array_begin(); it!=mesh.v_array_end(); )
    {
        if(mesh.is_vertex_removed(*it))
        {
            it = mesh.erase_vertex(it);
        }
        else
        {
            new_v_positions[v] = v_counter;
            v_counter++;
            ++it;
        }
        v++;
    }
}

template<class C, class T> void Mesh_Updater::clean_vertices_array(Mesh<C,T> &mesh, ivect &new_v_positions, ivect &surviving_vertices)
{
    int v_counter = 1; // v_counter keeps the new vertex indexing while v is needed to set the new index
    new_v_positions.assign(mesh.get_vertices_num(),-1);

    vector<Vertex<C> > old_list = mesh.get_vertices_array();
    mesh.reset_vertices_vector();
    mesh.reserve_vertices_space(surviving_vertices.size());

    for(ivect_iter it=surviving_vertices.begin(); it!=surviving_vertices.end(); ++it)
    {
        mesh.add_vertex(old_list[*it-1]);
        new_v_positions[*it-1] = v_counter;
        v_counter++;
    }

    old_list.clear();
}

template<class C, class T> boost::dynamic_bitset<> Mesh_Updater::update_and_clean_top_cells_arrays(Mesh<C,T> &mesh, ivect &new_v_positions, vector<ivect > &new_top_positions,
                                                                                                   ivect &simplification_counters)
{
    boost::dynamic_bitset<> all_deleted = boost::dynamic_bitset<>(mesh.get_top_cells_types());
    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        ivect d_top_pos;
        all_deleted[d] = update_and_clean_top_d_cells_array(mesh,new_v_positions,d_top_pos,d,simplification_counters[d]);
        new_top_positions.push_back(d_top_pos);
    }
    return all_deleted;
}

template<class C, class T> bool Mesh_Updater::update_and_clean_top_d_cells_array(Mesh<C,T> &mesh, ivect &new_v_positions, ivect &new_top_positions, int d, int counter)
{
    // step 0: if we have removed all the top d-cells, then, simply clear the corresponding array
    if(counter==mesh.get_top_cells_num(d))
    {
        mesh.clear_top_d_cell_list(d);
        return true;
    }
    else
    {
        clean_top_d_cells_array(mesh,d,new_top_positions);
        update_d_tops_boundary(mesh,d,new_v_positions);
    }
    return false;
}

template<class C, class T> void Mesh_Updater::clean_top_d_cells_array(Mesh<C,T> &mesh, int d, ivect &new_d_top_positions)
{
    int t_counter = 1, t = 0; // t_counter keeps the new top indexing while t is needed to set the new index

    // (1) init the tops position indices array
    new_d_top_positions.assign(mesh.get_top_cells_num(d),-1);
    // (2) get the new position indices for the top cells that are not deleted
    for(typename Mesh<C,T>::top_d_array_iter it=mesh.top_d_array_begin(d); it!=mesh.top_d_array_end(d); ++it)
    {
        if(!mesh.is_top_cell_removed(*it))
        {
            new_d_top_positions[t] = t_counter;
            t_counter++;
        }
        t++;
    }
    // (3) reorder (inline) the top cells arrays
    reorder_d_tops_array(mesh,d,new_d_top_positions);
    // (4) delete the positions from t_counter to the end of the array
    mesh.resize_top_d_cell_array(d,t_counter-1);
}

template<class C, class T> void Mesh_Updater::update_tops_boundary(Mesh<C,T> &mesh, ivect &new_v_positions)
{
    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        update_d_tops_boundary(mesh,d,new_v_positions);
    }
}

template<class C, class T> void Mesh_Updater::update_d_tops_boundary(Mesh<C,T> &mesh, int d, ivect &new_v_positions)
{
    int t_id = 1;
    for(typename Mesh<C,T>::top_d_array_iter it=mesh.top_d_array_begin(d); it!=mesh.top_d_array_end(d); ++it)
    {
        if(!mesh.is_top_cell_removed(*it))
        {
            for(int v=0; v<it->get_vertices_num(); v++)
            {
                if(new_v_positions[abs(it->TV(v))-1]==-1)
                {
                    cout<<"UPDATING A TOP-SIMPLEX WITH A DELETED VERTEX: "<<abs(it->TV(v))-1<<endl;
                    cout<<t_id<<" "<<*it<<endl;
                    int a; cin>>a;
                }
                it->setTV(v, new_v_positions[abs(it->TV(v))-1]);
            }
        }
        else
            cerr<<"[WARNING] update_top_d_cells_array does not removes deleted top d-cells"<<endl;

        t_id++;
    }
}

template<class C, class T> void Mesh_Updater::reorder_vertices_array(Mesh<C,T> &mesh, ivect &new_v_positions)
{
    // (2) update global vertices array
    //     -1 identifies an already set vertex
    // ** the reordering is done inline the global vertex array ** ///
    for(int i=1;i<=mesh.get_vertices_num();i++)
    {
        int j = i -1; // -1 is needed to avoid array lookup error

        if(new_v_positions[j] == i || new_v_positions[j] < 0)
        {
            // mark the last entry visited...
            new_v_positions[j] = -1;
            continue;
        }

        while(new_v_positions[j] != i)
        {
            mesh.vertices_swap(i,new_v_positions[j]);
            int j_prime = new_v_positions[j] -1; // -1 is needed to avoid array lookup error
            new_v_positions[j] = -1;
            j = j_prime;
        }
        // mark the last entry visited...
        new_v_positions[j] = -1;
    }
}

template<class C, class T> void Mesh_Updater::reorder_tops_array(Mesh<C,T> &mesh, vector<ivect> &new_tops_positions)
{
    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        reorder_d_tops_array(mesh,d,new_tops_positions[d]);
    }
}

template<class C, class T> void Mesh_Updater::reorder_d_tops_array(Mesh<C,T> &mesh, int d, ivect new_d_tops_positions)
{
    for(int i=1; i<=mesh.get_top_cells_num(d); i++)
    {
        int j = i -1; // -1 is needed to avoid array lookup error

        if(new_d_tops_positions[j] == i || new_d_tops_positions[j] < 0)
        {
            // mark the last entry visited...
            new_d_tops_positions[j] = -1;
            continue;
        }

        while(new_d_tops_positions[j] != i && new_d_tops_positions[j] > 0)
        {
            mesh.tops_swap(d,i,new_d_tops_positions[j]);
            int j_prime = new_d_tops_positions[j] -1; // -1 is needed to avoid array lookup error
            new_d_tops_positions[j] = -1;
            j = j_prime;
        }

        // mark the last entry visited...
        new_d_tops_positions[j] = -1;
    }
}

#endif // MESH_UPDATER_H
