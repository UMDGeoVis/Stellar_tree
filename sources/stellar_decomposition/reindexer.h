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

#ifndef REINDEXER_H
#define REINDEXER_H

#include <vector>
#include <set>
#include <iostream>
#include <map>

#include "stellar_tree/vertex.h"
#include "stellar_tree/mesh.h"
#include "stellar_tree/stellar_tree.h"
#include "stellar_tree/mesh_updater.h"

/**
 * @brief A class that exploits the spatial coherence of the indexed enties of a Stellar tree and resort accordingly these in the indexed mesh and Stellar tree representation
 * This class also compresses the tree representation using the Compressed encoding.
 */
class Reindexer
{
public:
    /**
     * @brief A constructor method
     *
//     * @param type a ReindexType that tells how much reorganize and compress the tree and the indexed mesh.
     * ReindexType can be VERTICES or ALL_SAVE if we want to reoder only the vertices or every entities in the mesh/tree.
     */
    Reindexer(/*ReindexType type*/)
    {
        this->vertCounter=1;
//        this->reindexing = type;
    }

    /**
     * @brief A public method that exploits the spatial coherence of vertices and top cells and compresses their representation within the tree
     *
     * @param tree a Stellar_tree& argument representing the tree to reindex and compress
     * @param mesh a Mesh& argument representing the indexed mesh to reindex
     * @param save_original_positions a boolean that enables the backup of the vertices and top cells original position indexes in two auxiliary arrays
     */
    template<class C, class T> void reorganize_index_and_mesh(Stellar_Tree& tree, Mesh<C,T> &mesh, bool save_original_positions);
    /**
     * @brief A public method that exploits the spatial coherence of vertices and compresses their representation within the tree
     *
     * @param tree a Stellar_tree& argument representing the tree to reindex and compress
     * @param mesh a Mesh& argument representing the indexed mesh to reindex
     * @param save_original_positions a boolean that enables the backup of the vertices original position indexes in a auxiliary array
     */
    template<class C, class T> void reorganize_mesh_vertices(Stellar_Tree &tree, Mesh<C,T> &mesh, bool save_original_positions);
    /**
     * @brief A public method that exploits the spatial coherence of top cells and compresses their representation within the tree
     *
     * @param tree a Stellar_tree& argument representing the tree to reindex and compress
     * @param mesh a Mesh& argument representing the indexed mesh to reindex
     * @param save_original_positions a boolean that enables the backup of the top cells original position indexes in a auxiliary array
     */
    template<class C, class T> void reorganize_mesh_top_simplices(Stellar_Tree &tree, Mesh<C,T> &mesh, bool save_original_positions);

    /// the below operation is for debug purpose
    /**
     * @brief A public method that returns the array of original position indexes of the vertices
     *
     * @return ivect
     */
    inline ivect& get_original_vertices_ordering() { return original_v_positions; }
    /**
     * @brief A public method that returns the array of original position indexes of the top cells
     *
     * @return vector<ivect >
     */
    inline vector<ivect >& get_original_tops_ordering() { return original_tops_positions; }

    /**
     * @brief A public method that resets the auxiliary class variables associated to the vertices
     *
     */
    inline void reset_vertices_variables()
    {
        this->new_v_positions.clear();
        this->vertCounter=1;
    }

    /**
     * @brief A public method that clear the arrays of original position indexes of vertices and top cells
     *
     */
    inline void reset_original_ordering_variables()
    {
        this->original_v_positions.clear();
        this->original_tops_positions.clear();
    }
    /**
     * @brief A public method that resets the auxiliary variables of the class
     *
     */
    inline void reset()
    {
        this->new_tops_positions.clear();
        this->new_v_positions.clear();
        this->vertCounter=1;
        this->leaf_tuples_array.clear();
    }

private:
//    ///A private variable representing the type or reindexing and compression adopted for the tree and the mesh. ReindexType can be VERTICES or ALL_SAVE if we want to reoder only the vertices or every entity in the mesh/tree.
//    ReindexType reindexing;
    ///A private variable encoding the new position indexes of the vertices
    ivect new_v_positions;
    ///A private variable that keeps track of the next position index available
    int vertCounter;
    ///A private variable encoding the new position indexes of the vertices
    vector<ivect> new_tops_positions;
    ///A private variable that encodes for each top cells type an associative array, in which the key represents a leaves tuple and the value is the number of top cells indexed by the leaves tuple and the starting position index for the top cells group.
    vector<map<ivect,pair<int,int> > > leaf_tuples_array;

    ///A private variable encoding the original position indexes of the vertices
    ivect original_v_positions;
    ///A private variable encoding the original position indexes of the top cells
    vector<ivect> original_tops_positions;

    /**
     * @brief A private method that compress the vertices representation in the tree
     * The tree is visited recursively and in the end each vertices list is represented with just two integers.
     * Moreover, each internal node encodes the vertices range of the vertices encoded in the subtree rooted in the internal node
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param save_original_positions a boolean that enables the backup of the vertices original position indexes in a auxiliary array
     */
    void compress_tree_representation_vertices(Node_Stellar& n, bool save_original_positions);
    /**
     * @brief A private method that inits the auxiliary variables needed for the exploiting of the spatial coherence for the top cells
     *
     * @param mesh a Mesh& argument representing the indexed mesh to reindex
     * @param save_original_positions a boolean that enables the backup of the top cells original position indexes in a auxiliary array
     */
    template<class C, class T> void init_tops_variables(Mesh<C,T> &mesh, bool save_original_positions);
    /**
     * @brief A private method that extract for each top cell indexed in the tree the leaves tuple indexing it (visiting procedure).
     * The procedure visits recursively the nodes of the tree, and for each leaf block calls extract_leaf_top_association_leaf
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param mesh a Mesh& argument representing the indexed mesh to reindex
     * @param root a Node_Stellar& argument representing the root of the tree
     */
    template<class C, class T> void extract_leaf_top_association(Node_Stellar& n, Mesh<C, T> &mesh, Node_Stellar &root);
    /**
     * @brief A private method that extract for each top cell indexed in the tree the leaves tuple indexing it.
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param mesh a Mesh& argument representing the indexed mesh to reindex
     * @param root a Node_Stellar& argument representing the root of the tree
     */
    template<class C, class T> void extract_leaf_top_association_leaf(Node_Stellar& n, Mesh<C,T>& mesh, Node_Stellar& root);
    /**
     * @brief A private method that visit the tree in order to extract the leaves tuple for a top cell indexed in more than a leaf block
     * The procedure visits recursively the nodes of the tree.
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param v_ids an integer vector containing the vertices of the top cell indexed outside the target node
     * @param key an integer vector encoding the leaves tuple
     */
    void get_top_tuple_key(Node_Stellar& n, ivect &v_ids, ivect &key);
    /**
     * @brief A private method that extracts the spatially coherent position indexes for the top cells of the mesh
     *
     * @param save_original_positions a boolean that enables the backup of the top cells original position indexes in a auxiliary array
     */
    void get_tops_reordered_indexes(bool save_original_positions);
    /**
     * @brief A private method that compress the top cell arrays in the tree (visiting procedure).
     *
     * The procedure visits recursively the nodes of the tree and call for each leaf block the compress_tree_representation_top_leaf procedure.
     *
     * @param n a Node_Stellar& argument representing the current node
     */
    void compress_tree_representation_top(Node_Stellar& n);
    /**
     * @brief A private method that compress the top cell arrays in the tree.
     *
     * The procedure extracts the new list with the updated top cells position indexes and then compress it using the Sequential-range encoding (SRE) (internally calls compress_t_list).
     *
     * @param n a Node_Stellar& argument representing the current node
     */
    void compress_tree_representation_top_leaf(Node_Stellar& n);
};

template<class C, class T> void Reindexer::reorganize_index_and_mesh(Stellar_Tree& tree, Mesh<C, T> &mesh, bool save_original_positions)
{
//    cerr<<"reorganize_mesh_vertices"<<endl;
    reorganize_mesh_vertices(tree,mesh,save_original_positions);
//    cerr<<"reorganize_mesh_top_simplices"<<endl;
    reorganize_mesh_top_simplices(tree,mesh,save_original_positions);
//    cerr<<"reset"<<endl;
    reset();
    return;
}

template<class C, class T> void Reindexer::reorganize_mesh_vertices(Stellar_Tree &tree, Mesh<C, T> &mesh, bool save_original_positions)
{
    new_v_positions.assign(mesh.get_vertices_num(),-1);
    if(save_original_positions)
        original_v_positions.assign(mesh.get_vertices_num(),-1);
    compress_tree_representation_vertices(tree.get_root(),save_original_positions);

    Mesh_Updater mu;
    mu.update_tops_boundary(mesh,this->new_v_positions);
    mu.reorder_vertices_array(mesh,this->new_v_positions);
}

template<class C, class T> void Reindexer::init_tops_variables(Mesh<C,T> &mesh, bool save_original_positions)
{
    // INIT PHASE
    //for the sub-simplexes (if any)
    for(int i=0; i<mesh.get_top_cells_types(); i++)
    {
        new_tops_positions.push_back(ivect());
        new_tops_positions[i].assign(mesh.get_top_cells_num(i),-1);

        if(save_original_positions)
        {
            original_tops_positions.push_back(ivect());
            original_tops_positions[i].assign(mesh.get_top_cells_num(i),-1);
        }

        leaf_tuples_array.push_back(map<ivect,pair<int,int> >());
    }
}

template<class C, class T> void Reindexer::reorganize_mesh_top_simplices(Stellar_Tree &tree, Mesh<C, T> &mesh, bool save_original_positions)
{
    if(/*this->reindexing != VERTICES &&*/ mesh.count_top_cells_num()>0)
    {
//        cerr<<"init_tops_variables"<<endl;
        init_tops_variables(mesh,save_original_positions);
//        cerr<<"extract_leaf_top_association"<<endl;
        extract_leaf_top_association(tree.get_root(),mesh,tree.get_root());
//        cerr<<"get_tops_reordered_indexes"<<endl;
        get_tops_reordered_indexes(save_original_positions);
//        cerr<<"compress_tree_representation_top"<<endl;
        compress_tree_representation_top(tree.get_root());
        Mesh_Updater mu;
//        cerr<<"reorder_tops_array"<<endl;
        mu.reorder_tops_array(mesh,this->new_tops_positions);
    }
    return;
}

///////////////////////***************************/////////////////////////////
template<class C, class T> void Reindexer::extract_leaf_top_association(Node_Stellar& n, Mesh<C, T> &mesh, Node_Stellar &root)
{
    if (n.is_leaf())
    {
        extract_leaf_top_association_leaf(n,mesh,root);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                extract_leaf_top_association(**it,mesh,root);
        }
    }
}

template<class C, class T>  void Reindexer::extract_leaf_top_association_leaf(Node_Stellar &n, Mesh<C, T> &mesh, Node_Stellar &root)
{
    ivect internal_t_key; internal_t_key.push_back(n.get_v_start());

    for(int i=0; i<n.get_num_top_cells_encoded(); i++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& top_id = itPair.first;
            if(n.has_lower_index_inside(mesh.get_top_cell(i,*top_id)))
            {
                ivect key = internal_t_key;
                pair<int,int> value = make_pair(this->leaf_tuples_array[i].size(),1);

                if(!n.completely_indexes_top(mesh.get_top_cell(i,*top_id)))
                {
                    ivect outer_v_ids;
                    n.get_outer_v_ids(mesh.get_top_cell(i,*top_id),outer_v_ids);
                    get_top_tuple_key(root,outer_v_ids,key);
                }

                pair<map<ivect,pair<int,int> >::iterator,bool> ret = this->leaf_tuples_array[i].insert(make_pair(key,value));
                if(ret.second) // inserted
                {
                    this->new_tops_positions[i][*top_id-1] = ret.first->second.first;
                }
                else
                {
                    this->new_tops_positions[i][*top_id-1] = ret.first->second.first;
                    ret.first->second.second++;
                }
            }
        }
    }
}

#endif // REINDEXER_H
