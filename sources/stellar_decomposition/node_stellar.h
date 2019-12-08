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

#ifndef NODE_STELLAR_H
#define NODE_STELLAR_H

#include "node_v.h"
#include "stellar_tree/mesh.h"
#include "utilities/basic_wrappers.h"
#include "topological_ds/links_aux_structures.h"
#include "statistics/index_statistics.h"

/**
 * @brief A class, extending Node_V, representing the node used by a Stellar tree.
 *
 * This nodes encodes also an array (or arrays) of top cells, separating each top cells type in a separate array.
 *
 */
class Node_Stellar : public Node_V<Node_Stellar>
{
public:
    /**
     * @brief A constructor method
     *
     */
    Node_Stellar() : Node_V<Node_Stellar>() { }
    /**
     * @brief A constructor method
     *
     * @param orig
     */
    Node_Stellar(const Node_Stellar& orig) : Node_V<Node_Stellar>(orig) { this->top_cells = orig.top_cells; }
    /**
     * @brief A destructor method
     *
     */
    virtual ~Node_Stellar() {}
    /**
     * @brief A public method that clears the two arrays encoded by the node
     *
     */
    inline void clear_lists()
    {
        this->vertices.clear();
        this->top_cells.clear();
    }
    /**
     * @brief A public method that clears the arrays of top cells
     *
     */
    inline void clear_top_d_cells_array() { this->top_cells.clear(); }
    /**
     * @brief A public method that clears the array of top cells at position pos
     *
     * @param pos is the position of the array to clear
     */
    inline void clear_top_d_cells_array(int pos) { this->top_cells[pos].clear(); }
    /**
     * @brief A public method that initializes the arrays of top cells in the node
     *
     * @param dim is the number of top cells types that cam be indexed by the node
     */
    inline void init_sub_cells_vectors(int dim) { top_cells.assign(dim,ivect()); }
    /**
     * @brief A public method that adds a top cell to a given
     *
     * @param dim is the position of the top cells array
     * @param ind is the position index of the top cell to add
     */
    inline void add_top_cell(int dim, int ind) { this->top_cells[dim].push_back(ind); }
    /**
     * @brief A public method that returns a pair of iterator to a top cells array at a given dim
     *
     * @param dim is the position of the top cells array
     * @return RunIteratorPair
     */
    inline RunIteratorPair make_t_array_iterator_pair(int dim) { return run_iterator<int>::make_run_iterator_pair(top_cells[dim]); }
    /**
     * @brief A public method that returns the iterator to the beginning of a top cells array at a given dimension
     *
     * @param dim is the position of the top cells array
     * @return RunIterator
     */
    inline RunIterator t_array_begin_iterator(int dim) { return run_iterator<int>(top_cells[dim].begin(),top_cells[dim].end()); }
    /**
     * @brief A public method that returns the iterator to the end of a top cells array at a given dimension
     *
     * @param dim is the position of the top cells array
     * @return RunIterator
     */
    inline RunIterator t_array_end_iterator(int dim) { return run_iterator<int>(top_cells[dim].end()); }
    /**
     * @brief A public method that returns the number of indexed top cells of a given dimension
     *
     * @param dim is the position of the top cells array
     * @return int
     *
     * NOTA: this method expand the runs and returns the real number of top d-cells indexed in the node
     */
    inline int get_real_t_array_size(int dim) { return run_iterator<int>(top_cells[dim].begin(),top_cells[dim].end()).elementCountFast(top_cells[dim]); }
    /**
     * @brief A public method that returns the size of the top cells array of a given dimension
     *
     * @param dim is the position of the top cells array
     * @return int
     */
    inline int get_t_array_size(int dim) { return this->top_cells[dim].size(); }
    /**
     * @brief A public method that returns the number of top d-cells types indexed in the node
     *
     * @return int
     */
    inline int get_num_top_cells_encoded() { return this->top_cells.size(); }
    /**
     * @brief A public method that update the list of top d-cells encoded by the leaf block
     *
     * @param d is the position of the top cells array to update
     * @return tl contains the updated top d-cells list
     *
     * NOTA: this update can invalidate the run compression for the top cells
     */
    inline void set_t_array(int d, ivect &tl) { this->top_cells[d] = tl; }
    /**
     * @brief A public method that set some statistics related to the run size
     *
     * @param dim is the position of the top cells array
     * @param run_counter is a global counter of the runs in the tree
     * @param min_run refers to the minimum run size that appears in the tree
     * @param max_run refers to the maximum run size that appears in the tree
     * @return an integer containing the summation of the top cells at a given dimension
     */
    inline int get_t_in_run_size(int dim, int &run_counter, int &min_run, int &max_run)
    {
        return run_iterator<int>(top_cells[dim].begin(),top_cells[dim].end()).runElementCount(top_cells[dim],run_counter,min_run,max_run);
    }

    // EXTRACT TOPOLOGICAL RELATIONS //

    /**
     * @brief A public method that extracts the local Vertex-Top-coboundary relation for the indexed vertices in the current node
     *
     * @param mesh is the mesh indexed by the tree
     * @param vts saves the local vertex-top-coboundary relations for the indexed vertices
     */
    template<class C, class T> void extract_local_VTop(Mesh<C,T> &mesh, leaf_VT &vts);
    /**
     * @brief A public method that extracts the local Vertex-Vertex relation for the indexed vertices in the current node
     *
     * @param mesh is the mesh indexed by the tree
     * @param vvs saves the local vertex-vertices relations for the indexed vertices
     */
    template<class C, class T> void extract_local_VV(Mesh<C,T> &mesh, leaf_VV &vvs);
    /**
     * @brief A public method that extracts the local links for the indexed vertices and edges in the current node
     *
     * @param mesh is the mesh indexed by the tree
     * @param links saves the local links for the indexed vertices and edges
     */
    template<class C, class T> void extract_local_links(Mesh<C,T> &mesh, local_links &links);
    /**
     * @brief A public method that extracts the local Edge-Top-coboundary relation for the edges in the current node
     *
     * @param mesh is the mesh indexed by the tree
     * @param ets saves the local edge-top-coboundary relations for the indexed vertices
     *
     * NOTA: the procedure extract the ETop relation for those edges that have maximum entreme indexed by the current leaf block.
     */
    template<class C, class T> void extract_local_ETop(Mesh<C,T> &mesh, leaf_ET &ets);

    /**
     *
     */
    template<class C, class T> void extract_local_p_faces(Mesh<C,T> &mesh, pair<int,leaf_p_faces> &p);
    /**
     *
     */
    template<class C, class T> void extract_local_p_faces_v2(Mesh<C,T> &mesh, pair<int,leaf_p_faces> &p);

    // UPDATE/COMPRESS THE VERTICES/TOPS ARRAYS //

    /**
     * @brief A public method that removes the vertices that have been deleted during a simplification process
     * The procedure saves in an auxiliary vector the position indexes of those not deleted
     *
     * @param mesh is the mesh indexed by the tree
     * @param surviving_vertices, a vector of integer containing the position indexes of the vertices not deleted
     */
    template<class C,class T> void compact_vertices_array(Mesh<C, T> &mesh, ivect &surviving_vertices);
    /**
     * @brief A public method that updates the position indexes of the vertices
     *
     * @param new_v_indices, a vector of integers containing the updated position indexes
     *
     * NOTA: the SRE encoding is not valid any more. A reindexing operation is needed after calling this procedure on all the leaf blocks
     */
    void update_vertex_indices(ivect &new_v_indices);
    /**
     * @brief A public method that updates the top cells list of the leaf block after a simplification process.

     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     */
    template<class C,class T> inline void compact_top_cell_arrays(Mesh<C, T> &mesh)
    {
        for(int d=0; d<this->get_num_top_cells_encoded(); d++)
            this->compact_top_d_cells_array(d,mesh);
    }
    /**
     * @brief A public method that updates the top cells encoded in the leaf block by updating the position indexes
     *
     * @param new_top_positions, a vector of vectors of integers that contains the new position indexes of the top cells
     * @param all_deleted, a bit-vector stating if a specific top d-cell type has been completely erased by a simplification procedure
     *
     * NOTA: this update removes the run compression!
     */
    void update_and_compress_top_cells_arrays(vector<ivect > &new_top_positions, boost::dynamic_bitset<> &all_deleted);
    /**
     * @brief A public method that updates the top d-cells encoded in the leaf block by updating the position indexes
     *
     * @param d an integer representing the list to update
     * @param new_indices, a vector of integers that contains the new position indexes of the top d-cells
     *
     * NOTA: this update removes the run compression!
     */
    void update_and_compress_top_d_cells_array(int d, ivect &new_indices);
    /**
     * @brief A public procedure that compress a top d-cells array of the leaf block
     *
     * The procedure compresses the top d-cells array using the Sequential-Range Encoding (SRE).
     *
     * @param pos an integer identifying the position index of the array to update/compress
     * @param new_t_list an integer vector containing the updated array
     */
    void compress_top_d_cell_array(int pos, ivect &new_t_list);

    // STATISTICS AND DEBUG PRINTS //

    /**
     * @brief A public method that computes the histogram of the runs on the current node
     *
     * @param is, an IndexStatistics variable, saves the run statistics
     */
    void compute_run_histogram(IndexStatistics &is);
    /**
     * @brief A public method that prints on the standard output the histogram of the runs on the current node
     *
     * @param is, an IndexStatistics variable, saves the run statistics
     */
    template<class C, class T> void print_top_cells_arrays(Mesh<C,T> &mesh);

    template<class C, class T> void print_top_cells_arrays(int d, Mesh<C,T> &mesh);
    /**
     * @brief print_top_cells_lists
     */
    void print_top_cells_arrays();
    /**
     * @brief print_top_cells_list
     * @param d
     */
    void print_top_cells_array(int d);

    bool check_duplicate_top_cells_array(int d, int t_id);

    void print_node_stats()
    {
        for(int d=0; d<this->get_num_top_cells_encoded(); d++)
        {
            int top_size = this->get_real_t_array_size(d);
            if(top_size > 0)
                cerr << "   top "<<d+1<<"-simplices num: "<<top_size<<endl;
        }
    }

protected:
    ///A private variable representing the list containing the to cells (of different types) indexed in the node
    vector<ivect > top_cells;

private:
    /**
     * @brief A public method that compact the top d-cells list of the leaf block after a simplification process.
     *
     * @param d an integer representing the list to update
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     */
    template<class C,class T> void compact_top_d_cells_array(int d, Mesh<C, T> &mesh);
    /**
     * @brief A public method that extracts the local links for the indexed vertices and edges from the current top simplex
     * This procedure defines a wrapper for simplicial meshes
     *
     * @param t, the current top simplex
     * @param links saves the local links for the indexed vertices and edges
     */
    void extract_top_link(Top_Simplex &t, local_links &links);
    /**
     * @brief A public method that extracts the local links for the indexed vertices and edges from the current top cp-cell
     * This procedure defines a wrapper for CP meshes
     *
     * @param t, the current top cp-cell
     * @param links saves the local links for the indexed vertices and edges
     */
    void extract_top_link(Top_CP_Cell &t, local_links &links);
    /**
     * @brief A public method that extracts the local Edge-Top-coboundary relation for the edges of the current top simplex
     * NOTA: the procedure extract the ETop relation for those edges that have maximum entreme indexed by the current leaf block.
     * NOTA2: the method defines a wrapper for simplicial meshes
     *
     * @param d, an integer representing the dimension of the current top simplex
     * @param t_id, the position index of the current top simplex
     * @param t, the current top simplex
     * @param ets saves the local Edge-Top-coboundary relations for the indexed vertices
     * @param empty_et, an empty ETop relation used to initialize the ETop relation of a 'firstly touched' edge
     */
    void extract_top_ETop(int d, int t_id, Top_Simplex &t, leaf_ET &ets, ET &empty_et);
    /**
     * @brief A public method that extracts the local Edge-Top-coboundary relation for the edges of the current top cp-cell
     * NOTA: the procedure extract the ETop relation for those edges that have maximum entreme indexed by the current leaf block.
     * NOTA2: the method defines a wrapper for CP meshes
     *
     * @param d, an integer representing the dimension of the current top cp-cell
     * @param t_id, the position index of the current top cp-cell
     * @param t, the current top cp-cell
     * @param ets saves the local Edge-Top-coboundary relations for the indexed vertices
     * @param empty_et, an empty ETop relation used to initialize the ETop relation of a 'firstly touched' edge
     */
    void extract_top_ETop(int d, int t_id, Top_CP_Cell &t, leaf_ET &ets, ET &empty_et);
};

template<class C, class T> void Node_Stellar::extract_local_VTop(Mesh<C,T> &mesh, leaf_VT &vts)
{
    VT tmp;
    tmp.assign(mesh.get_top_cells_types(),ivect());

    // we reserve the maximum size a VT can have (in order to avoid reallocation)
    for(int d=0; d<this->get_num_top_cells_encoded();d++)
        tmp[d].reserve(this->get_real_t_array_size(d));

    vts.assign(this->get_v_end()-this->get_v_start(),tmp);

    for(int d=0; d<this->get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = this->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                T& t = mesh.get_top_cell(d,*t_id);
                for(int j=0; j<t.get_vertices_num(); j++)
                {
                    int real_index = abs(t.TV(j));
                    if(this->indexes_vertex(real_index))
                    {
                        vts[real_index-this->get_v_start()][d].push_back(*t_id);
                    }
                }
            }
        }
    }
}

template<class C, class T> void Node_Stellar::extract_local_VV(Mesh<C,T> &mesh, leaf_VV &vvs)
{
    vvs.assign(this->get_v_end()-this->get_v_start(),iset());

    for(int d=0; d<this->get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = this->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                T& t = mesh.get_top_cell(d,*t_id);
                for(int j=0; j<t.get_vertices_num(); j++)
                {
                    int real_index = abs(t.TV(j));
                    if(this->indexes_vertex(real_index))
                    {
                        if(mesh.is_simplex(d)) // simplicial complexes
                        {
                            for(int i=1; i<t.get_vertices_num(); i++)
                            {
                                vvs[real_index-this->get_v_start()].insert(t.TV((int)(j+i)%t.get_vertices_num()));
                            }
                        }
                        else // CP-complexes
                        {
                            for(int i=0; i<t.get_vertices_num(); i++)
                            {
                                if(i != j)
                                    vvs[real_index-this->get_v_start()].insert(t.TV(i));
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class C, class T> void Node_Stellar::extract_local_links(Mesh<C,T> &mesh, local_links &links)
{
    links.init(*this,mesh.get_implicitly_encoded_cells_num());
    for(int d=0; d<this->get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = this->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                T& t = mesh.get_top_cell(d,*t_id);
                this->extract_top_link(t,links);
            }
        }
    }
}

template<class C, class T> void Node_Stellar::extract_local_ETop(Mesh<C,T> &mesh, leaf_ET &ets)
{
    ET empty_et;
    empty_et.assign(mesh.get_top_cells_types(),ivect());

    for(int d=0; d<this->get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = this->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                T& t = mesh.get_top_cell(d,*t_id);
                this->extract_top_ETop(d,*t_id,t,ets,empty_et);
            }
        }
    }
}

template<class C, class T> void Node_Stellar::extract_local_p_faces(Mesh<C,T> &mesh, pair<int,leaf_p_faces> &p)
{
    ivect pcell;

    /// we cannot place d instead of 0 because we do not know here if we are encoding the tops using a verbose encoding or a compact one
    for(int i=0; i<this->get_num_top_cells_encoded();i++)
    {
        for(RunIteratorPair itPair = this->make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& runIt = itPair.first;
            T& top = mesh.get_top_cell(i,*runIt);

            int num_s = top.get_sub_types_num(p.first);

            if(num_s == -1)
            {
                break; /// we stop this for as we are analyzing a top simplex that is lower dimensional than the one that we want to extract
            }

            if(this->completely_indexes_top(top))
            {
                /// we extract all the p-faces without checking
                for(int s=0; s<num_s; s++)
                {
                    top.get_d_cell(pcell,p.first,s);
                    p.second.insert(pcell);
                    pcell.clear();
                }
            }
            else
            {
                for(int s=0; s<num_s; s++)
                {
                    top.get_d_cell(pcell,p.first,s);
                    if(this->indexes_cell(pcell))
                    {
                        p.second.insert(pcell);
                    }
                    pcell.clear();
                }
            }
        }
    }
}

template<class C, class T> void Node_Stellar::extract_local_p_faces_v2(Mesh<C,T> &mesh, pair<int,leaf_p_faces> &p)
{
    ivect pcell;

    /// we cannot place d instead of 0 because we do not know here if we are encoding the tops using a verbose encoding or a compact one
    for(int i=0; i<this->get_num_top_cells_encoded();i++)
    {
        for(RunIteratorPair itPair = this->make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& runIt = itPair.first;
            T& top = mesh.get_top_cell(i,*runIt);

            if(!this->indexes_vertex(top.min_v_index())) /// we consider a top only in the first leaf block indexing it
                continue;

            int num_s = top.get_sub_types_num(p.first);

            if(num_s == -1)
            {
                break; /// we stop this for as we are analyzing a top simplex that is lower dimensional than the one that we want to extract
            }

            if(this->completely_indexes_top(top))
            {
                /// we extract all the p-faces without checking
                for(int s=0; s<num_s; s++)
                {
                    top.get_d_cell(pcell,p.first,s);
                    p.second.insert(pcell);
                    pcell.clear();
                }
            }
            else
            {
                for(int s=0; s<num_s; s++)
                {
                    top.get_d_cell(pcell,p.first,s);
                    if(this->indexes_cell(pcell))
                    {
                        p.second.insert(pcell);
                    }
                    pcell.clear();
                }
            }
        }
    }
}

template<class C, class T> void Node_Stellar::print_top_cells_arrays(Mesh<C,T> &mesh)
{
    for(int d=0; d<this->get_num_top_cells_encoded();d++)
    {
        if(this->get_t_array_size(d)>0)
            cout<<d<<"-position simplex"<<endl;
        for(RunIteratorPair itPair = this->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                cout<<*t_id<<"] "<<mesh.get_top_cell(d,*t_id);
                if(mesh.get_top_cell(d,*t_id).is_degenerate())
                    cout<<" --> is degenerate"<<endl;
                else
                    cout<<endl;
            }
        }
        if(this->get_t_array_size(d)>0)
            cout<<endl;
    }
}

template<class C, class T> void Node_Stellar::print_top_cells_arrays(int d, Mesh<C,T> &mesh)
{
//    for(int d=0; d<this->get_num_top_cells_encoded();d++)
//    {
    if(this->get_t_array_size(d)>0)
        cout<<d<<"-position simplex"<<endl;
    for(RunIteratorPair itPair = this->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        if(!mesh.is_top_cell_removed(d,*t_id))
        {
            cout<<*t_id<<"] "<<mesh.get_top_cell(d,*t_id);
            if(mesh.get_top_cell(d,*t_id).is_degenerate())
                cout<<" --> is degenerate "/*<<endl*/;
            else
                cout<</*endl*/" ";
        }
    }
    if(this->get_t_array_size(d)>0)
        cout<<endl;
//    }
}

template<class C,class T> void Node_Stellar::compact_vertices_array(Mesh<C, T> &mesh, ivect &surviving_vertices)
{
    ivect new_v_list;
    for(RunIteratorPair itPair = make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& v_id = itPair.first;
        if(!mesh.is_vertex_removed(*v_id))
        {
            new_v_list.push_back(*v_id);
            surviving_vertices.push_back(*v_id);
        }
    }
    this->clear_v_list();
    this->set_v_array(new_v_list);
}

template<class C,class T> void Node_Stellar::compact_top_d_cells_array(int d, Mesh<C, T> &mesh)
{
    ivect t_list;

    for(RunIteratorPair itPair = this->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& t_id = itPair.first;
        // NEW the top cell must be into the leaf block as well
        // (for reusing position indices)
        // NOTA -> this completely removes the run-lenght encoding
        if(!mesh.is_top_cell_removed(d,*t_id) && this->indexes_top(mesh.get_top_cell(d,*t_id)))
        {
            t_list.push_back(*t_id);
        }
    }

    this->clear_top_d_cells_array(d);
    if(t_list.size() > 0)
        this->compress_top_d_cell_array(d,t_list); // this sets again the SRE encoding
}

#endif // NODE_STELLAR_H
