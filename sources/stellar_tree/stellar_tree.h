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

#ifndef STELLAR_TREE_H
#define STELLAR_TREE_H

#include "box.h"
#include "io/writer.h"
#include "mesh.h"
#include "cells_soup.h"
#include "utilities/timer.h"
#include "spatial_subdivision.h"
#include "stellar_decomposition/node_stellar.h"

#include <unordered_map>

/**
 * @brief A class representing the spatio-topological index, the Stellar tree,
 *
 * This class represents the true core of the Stellar tree of this framework.
 */
class Stellar_Tree
{
public:
    /**
     * @brief A constructor method
     *
     * @param vertices_per_leaf an integer argument, represents the maximum number of vertices admitted for a tree node
     */
    Stellar_Tree(int vertices_per_leaf)
    {
        this->kv = vertices_per_leaf;
        this->root = Node_Stellar();
    }
    /**
     * @brief A copy-constructor method
     *
     * @param orig
     */
    Stellar_Tree(const Stellar_Tree& orig)
    {
        this->kv = orig.kv;
        this->subdivision = orig.subdivision;
        this->root = orig.root;
    }
    /**
     * @brief A destructor method
     *
     */
    virtual ~Stellar_Tree() {}
    /**
     * @brief A public method that returns the root node of the tree
     *
     * @return Node_Stellar variable, representing the root node
     */
    inline Node_Stellar& get_root() { return this->root; }
    ///
    /**
     * @brief A public method that initializes the root node of the tree
     *
     */
    inline void init_root() { this->root = Node_Stellar(); }
    /**
     * @brief A public method that returns the division type associated to the tree
     *
     * @return Spatial_Subdivision& variable, representing the division type of the tree
     */
    inline Spatial_Subdivision& get_subdivision() { return this->subdivision; }
    /**
     * @brief A public method that sets the maximum number of vertices admitted for a node
     *
     * @param maxV an integer argument, represents the maximum number of vertices
     */
    inline void set_max_vertex(int maxV) { this->kv = maxV; }
    /**
     * @brief A public method that returns the maximum number of vertices admitted for a node
     *
     * @return int value representing the node vertices threshold
     */
    inline int get_max_vertex() { return this->kv; }
    /**
     * @brief A public method that initializes the spatial subdivision of the tree
     *
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     */
    template<class C, class T> inline void set_subdivision(Mesh<C,T> &mesh) { subdivision.set(mesh.get_domain().getMinPoint().get_dimension()); }
    /**
     * @brief A public method that initializes the spatial subdivision of the tree
     *
     * @param dim an integer argument representing the dimension in which the mesh is embedded
     */
    inline void set_subdivision(int dim) { subdivision.set(dim); }

    /**
     * @brief A public methods that builds the tree starting from a indexed mesh representation

     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     * @param filename a string argument representing the filename given to the file containing the Stellar tree
     * @param reindexer a R& argument representing the object that extract the spatial coherence of vertices and top cells and compress the tree representation
     * @param save_original_positions a bool argument that enables the saving of the original position indexes for vertices and top cells
     */
    template<class C, class T, class R> void build(Mesh<C,T> &mesh, string filename, R &reindexer, bool save_original_positions);
    /**
     * @brief A public method that builds the tree starting from a soup of cells and extract the indexed representation of the soup
     *
     * @param soup a Cells_Soup& argument representing the soup of cells indexed by the Stellar tree
     * @param mesh a Mesh& argument representing the indexed mesh representation of the soup
     * @param reindexer a R& argument representing the object that extract the spatial coherence of vertices and top cells and compress the tree representation
     * @param save_original_positions a bool argument that enables the saving of the original position indexes for vertices and top cells
     */
    template<class C, class T, class R> void build(Cells_Soup &soup, Mesh<C,T> &mesh, R &reindexer, bool save_original_positions);
    /**
     * @brief A public method that generate the starting spatial subdivision by inserting the vertices of the mesh and exploits the spatial coherence of them
     *
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     * @param reindexer a R& argument representing the object that extract the spatial coherence of vertices and top cells and compress the tree representation
     */
    template<class C, class T, class R> void generate_spatial_subdivision(Mesh<C,T> &mesh, R &reindexer);
    /**
     * @brief A public method that adds the top cells to a Stellar tree and exploit the spatial coherence of them
     *
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     * @param reindexer a R& argument representing the object that extract the spatial coherence of vertices and top cells and compress the tree representation
     */
    template<class C, class T, class R> void add_and_resort_top_cells(Mesh<C,T> &mesh, R &reindexer);

    /**
     * @brief A public method that initialize the array of leaf blocks references in order to enable the parallel visit of leaf blocks
     *
     * @param n a Node_Stellar& argument representing the current node
     *
     * NOTA: this procedure must be called after loading/building a Stellar tree
     */
    void init_leaves_list(Node_Stellar& n);
    /**
     * @brief A public method that returns the leaf block at i-th position in the leaves array
     *
     * @param i an integer argument representing the position index of the leaf
     * @return Node_Stellar* a pointer reference to the leaf block
     */
    inline Node_Stellar* get_leaf(unsigned i) { return this->leaves[i]; }
    /**
     * @brief A public method that returns the number of leaf blocks in the Stellar tree
     *
     * @return unsigned
     */
    inline unsigned get_leaves_number() { return this->leaves.size(); }

    // UPDATE FUNCTIONS //

    /**
     * @brief A public method that updates the Stellar tree after a simplification process.
     *
     * The procedure removes deleted vertices and top cells and updates the position indexes of vertices and top-cells
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param new_v_positions, a vector of integers that contains the new position indexes of the vertices
     * @param new_top_positions, a vector of vectors of integers that contains the new position indexes of the top cells
     */
    void update_tree(Node_Stellar &n, ivect &new_v_positions, vector<ivect > &new_top_positions, boost::dynamic_bitset<> &all_deleted);
    /**
     * @brief A public method that updates the top cell indices in the Stellar tree after a simplification process.
     *
     * The procedure removes deleted top cells and updates the position indexes of top-cells
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param new_top_positions, a vector of vectors of integers that contains the new position indexes of the top cells
     * @param
     */
    void compact_and_update_top_cells_lists(Node_Stellar &n, vector<ivect > &new_top_positions, boost::dynamic_bitset<> &all_deleted);
    /**
     * @brief A public method that updates the top d-cell indices in the Stellar tree after a simplification process.
     *
     * The procedure removes deleted top d-cells and updates the position indexes of top-cells
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param p, a pair containing in the first position the position of the top d-cells in the Stellar tree, while in the second position a vector of integers that contains the new position indexes of the top d-cells
     */
    void compact_and_update_top_d_cells_lists(Node_Stellar &n, pair<int,ivect > &p);
    /**
     * @brief A public method that compacts the vertices lists in the Stellar tree after a simplification process.
     *
     * The procedure removes the deleted vertices and returns a list of vertices that are still existing
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     * @param surviving_vertices, a vector that in the end contains the position indexes of the valid vertices
     */
    template<class C, class T> void compact_vertices_lists(Node_Stellar &n, Mesh<C,T> &mesh, ivect &surviving_vertices);
    /**
     * @brief A public method that recursively visits the Stellar tree and unify the tree branches that index less than kv vertices.
     *
     * The procedure returns the number of vertices indexed in a tree branch.
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     * @return int an integer containing the number of vertices indexed in a tree branch
     */
    template<class C, class T> int visit_and_unify(Node_Stellar& n, Mesh<C, T> &mesh);

    // VISIT FUNCTIONS //

    /**
     * @brief A public function that recursively visits the nodes of a Stellar tree
     *     
     * @param function_name a F argument representing the function name to call within each leaf block
     * @param n a Node_Stellar& argument representing the current node
     * @param mesh a M& template argument representing the mesh indexed by the Stellar tree
     * @param parameters a P& argument representing the input/output parameters needed by the function function_name
     *
     * The function F can be a simple RAW function or a static function defined into a class.
     * The name of the function must be unique, otherwise the linker cannot understand correctly which function have to use
     * and, finally, the arguments of the functions must be all defined (i.e. templated parameter in the function definition are not allowed).
     * The function arguments are order dependent and each function F must have as first argument a Stellar tree node, then the indexed mesh and finally the parameters needed by F
     */
    template<typename F, class M, class P> void visit(F function_name, Node_Stellar &n, M &mesh, P &parameters);

    /**
     * @brief A public function that recursively visits the nodes of a Stellar tree
     *     
     * @param function_name a F argument representing the function name to call within each leaf block
     * @param n a Node_Stellar& argument representing the current node
     * @param parameters a P& argument representing the input/output parameters needed by the function function_name
     *
     * The function F can be a simple RAW function or a static function defined into a class.
     * The name of the function must be unique, otherwise the linker cannot understand correctly which function have to use
     * and, finally, the arguments of the functions must be all defined (i.e. templated parameter in the function definition are not allowed).
     * The function arguments are order dependent and each function F must have as first argument a Stellar tree node, then the parameters needed by F
     */
    template<typename F, class P> void visit(F function_name, Node_Stellar &n, P &parameters);

    /**
     * @brief A public function that recursively visits the nodes of a Stellar tree
     *
     * @param function_name a F argument representing the function name to call within each leaf block
     * @param n a Node_Stellar& argument representing the current node
     * @param dom, a Box& argument representing the domain of n
     * @param level, an integer representing the level of node n
     * @param mesh a M& template argument representing the mesh indexed by the Stellar tree
     * @param parameters a P& argument representing the input/output parameters needed by the function function_name
     *
     * The function F can be a simple RAW function or a static function defined into a class.
     * The name of the function must be unique, otherwise the linker cannot understand correctly which function have to use
     * and, finally, the arguments of the functions must be all defined (i.e. templated parameter in the function definition are not allowed).
     * The function arguments are order dependent and each function F must have as first argument a Stellar tree node, then the parameters needed by F
     */
    template<typename F, class M, class P> void visit(F function_name, Node_Stellar &n, Box<COORDBASETYPE> &dom, int level, M &mesh, P &parameters);

    /**
     * @brief A public function that parallel visits the leaf blocks of a Stellar tree
     *
     * @param function_name a F argument representing the function name to call within each leaf block
     * @param mesh a M& template argument representing the mesh indexed by the Stellar tree
     * @param parameters a P& argument representing the input/output parameters needed by the function function_name
     *
     * The parallel access of the leaf blocks is gathered with openMP framework
     *
     * The function F can be a simple RAW function or a static function defined into a class.
     * The name of the function must be unique, otherwise the linker cannot understand correctly which function have to use
     * and, finally, the arguments of the functions must be all defined (i.e. templated parameter in the function definition are not allowed)
     *
     * The function arguments are order dependent and each function F must have as first argument a Stellar tree leaf node, then the indexed mesh and finally the parameters needed by F
     */
    template<typename F, class M, class P> void parallel_visit(F function_name, M &mesh, P &parameters);

    /**
     * @brief A public function that parallel visits the leaf blocks of a Stellar tree
     *
     * @param function_name a F argument representing the function name to call within each leaf block
     * @param parameters a P& argument representing the input/output parameters needed by the function function_name
     *
     * The parallel access of the leaf blocks is gathered with openMP framework
     *
     * The function F can be a simple RAW function or a static function defined into a class.
     * The name of the function must be unique, otherwise the linker cannot understand correctly which function have to use
     * and, finally, the arguments of the functions must be all defined (i.e. templated parameter in the function definition are not allowed)
     *
     * The function arguments are order dependent and each function F must have as first argument a Stellar tree leaf node, then the parameters needed by F
     */
    template<typename F, class P> void parallel_visit(F function_name, P &parameters);

    /**
     * @brief A public function that recursively visits the nodes of a Stellar tree and saves in a LRU-cache the most recent used leaf blocks
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param function_name a F argument representing the function name to call within each leaf block
     * @param mesh a M& template argument representing the mesh indexed by the Stellar tree
     * @param cache an LRU-cache encoding the most used expanded leaf block representations
     * @param parameters a P& argument representing the input/output parameters needed by the function function_name
     *
     * The function F can be a simple RAW function or a static function defined into a class.
     * The name of the function must be unique, otherwise the linker cannot understand correctly which function have to use
     * and, finally, the arguments of the functions must be all defined (i.e. templated parameter in the function definition are not allowed)
     *
     * The function arguments are order dependent and each function F must have as first argument a Stellar tree node, the indexed mesh, the LRU-cache and finally the parameters needed by F
     */
    template<typename F, class M, class C, class P> void visit_with_cache(F function_name, Node_Stellar &n, M &mesh, C &cache, P &parameters);

    /**
     * @brief A public function that recursively visits the nodes of a Stellar tree and saves in a LRU-cache the most recent used leaf blocks
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param function_name a F argument representing the function name to call within each leaf block
     * @param cache an LRU-cache encoding the most used expanded leaf block representations
     * @param parameters a P& argument representing the input/output parameters needed by the function function_name
     *
     * The function F can be a simple RAW function or a static function defined into a class.
     * The name of the function must be unique, otherwise the linker cannot understand correctly which function have to use
     * and, finally, the arguments of the functions must be all defined (i.e. templated parameter in the function definition are not allowed)
     *
     * The function arguments are order dependent and each function F must have as first argument a Stellar tree node, the LRU-cache and finally the parameters needed by F
     */
    template<typename F, class C, class P> void visit_with_cache(F function_name, Node_Stellar &n, C &cache, P &parameters);

    /**
     * @brief A public function that parallel visits the leaf blocks of a Stellar tree and saves in a LRU-cache the most recent used leaf blocks
     *
     * @param function_name a F argument representing the function name to call within each leaf block
     * @param mesh a M& template argument representing the mesh indexed by the Stellar tree
     * @param cache an LRU-cache encoding the most used expanded leaf block representations
     * @param parameters a P& argument representing the input/output parameters needed by the function function_name
     *
     * The parallel access of the leaf blocks is gathered with openMP framework
     *
     * The function F can be a simple RAW function or a static function defined into a class.
     * The name of the function must be unique, otherwise the linker cannot understand correctly which function have to use
     * and, finally, the arguments of the functions must be all defined (i.e. templated parameter in the function definition are not allowed)
     *
     * The function arguments are order dependent and each function F must have as first argument a Stellar tree leaf node, the indexed mesh, the LRU-cache and finally the parameters needed by F
     */
    template<typename F, class M, class C, class P> void parallel_visit_with_cache(F function_name, M &mesh, C &cache, P &parameters);

    /**
     * @brief A public function that parallel visits the leaf blocks of a Stellar tree and saves in a LRU-cache the most recent used leaf blocks
     *
     * @param function_name a F argument representing the function name to call within each leaf block
     * @param cache an LRU-cache encoding the most used expanded leaf block representations
     * @param parameters a P& argument representing the input/output parameters needed by the function function_name
     *
     * The parallel access of the leaf blocks is gathered with openMP framework
     *
     * The function F can be a simple RAW function or a static function defined into a class.
     * The name of the function must be unique, otherwise the linker cannot understand correctly which function have to use
     * and, finally, the arguments of the functions must be all defined (i.e. templated parameter in the function definition are not allowed)
     *
     * The function arguments are order dependent and each function F must have as first argument a Stellar tree leaf node, the LRU-cache and finally the parameters needed by F
     */
    template<typename F, class C, class P> void parallel_visit_with_cache(F function_name, C &cache, P &parameters);

    /**
     * @brief A public function that identify the leaf block indexing the target vertex and returns the pointer to the leaf block
     *
     * @param n a Node_Stellar variable, representing the current node block
     * @param v_id an integer representing the position index of the vertex
     * @param res a Node_Stellar* variable that stores the pointer to the leaf block indexing v_id
     */
    void get_leaf_indexing_vertex(Node_Stellar &n, int v_id, Node_Stellar *&res);
    /**
     * @brief get_leaves_indexing_cell
     *
     * @param n
     * @param v_ids
     * @param res
     */
    void get_leaves_indexing_cell(Node_Stellar &n, ivect &v_ids, vector<Node_Stellar*> &res);

private:
    ///A private variable representing the maximum number of vertices admitted in a node
    int kv;
    ///A private variable representing the root node of the tree
    Node_Stellar root;
    ///A private variable representing the division type of the tree
    Spatial_Subdivision subdivision;
    ///A private variable that contains the list of pointers to the leaf block of the tree. This enables the PARALLELIZATION of the index lookup
    vector<Node_Stellar*> leaves;

    /**
     * @brief A private function that insert the vertices of the mesh to the Stellar tree
     *
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     */
    template<class C, class T> void add_vertices(Mesh<C,T> &mesh);
    /**
     * @brief A private function that insert the vertices of the cells soup in the Stellar tree and extract the indexed mesh representation of the soup
     *
     * @param soup
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     */
    template<class C, class T> void add_vertices(Cells_Soup &soup, Mesh<C,T> &mesh);
    /**
     * @brief A private function that insert the top cells of the mesh in the Stellar tree
     *
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     */
    template<class C, class T> void add_top_cells(Mesh<C,T> &mesh);    

    /**
     * @brief A private method that initializes the lists of tops cell accordingly with the number of top cells types encoded by the mesh
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     */
    void init_tops_lists(Node_Stellar& n, int num_top_cells_types);
    /**
     * @brief A private method that insert a vertex in the leaf block that contains it
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * @param v a Vertex& argument representing the vertex to insert
     * @param vertex_index
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     *
     * This method checks if a vertex is contained by a node, and then insert the vertex into
     * the vertices array if the node is a leaf, otherwise if the node is internal a recursive call is activated.
     * If the maximum number of vertices for a node is exceeded, a splitting operation is activated.
     */
    template<class C, class T> void insert_vertex(Node_Stellar& n, Box<C>& domain, int level, Vertex<C>& v, int vertex_index, Mesh<C,T> &mesh);
    /**
     * @brief A private method that insert a vertex from a cells soup in the leaf block that contains it
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * @param v a Vertex& argument representing the vertex to insert
     * @param vertex_index an integer argument that after the insertion of v saves the position index of v in the indexed mesh representation
     * @param indexed_simplex a vector of integers representing the indexed representation of the cell containing the vertex v
     * @param first_time a boolean that is true if we insert the vertex for the first time, false otherwise
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     *
     * This method checks if a vertex is contained by a node, and then insert the vertex into
     * the vertices array if the node is a leaf, otherwise if the node is internal a recursive call is activated.
     * If the maximum number of vertices for a node is exceeded, a splitting operation is activated.
     */
    template<class C, class T> void insert_vertex_from_soup(Node_Stellar& n, Box<C>& domain, int level, Vertex<C>& v, int vertex_index, ivect &indexed_cell, bool &first_time, Mesh<C,T> &mesh);
    /**
     * @brief A private method that adds a top cell to the Stellar tree
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param t a T& argument representing the top cell to insert
     * @param topsimplex_index an integer argument representing the position index of the top cell in the mesh array
     * @param vector_position an integer argument representing the position of the top cells array in which we have to insert the top cell
     *
     * This method checks if a top cell has one of the vertices indexed by a tree branch, and then insert the top cell into the node top cells array if the node is a leaf,
     * otherwise if the node is internal a recursive call is activated. The procedure must be used after exploiting the spatial coherence of the vertices.
     * In this way to understand if a top cell is into a tree branch is just executing an integer range comparison.
     */
    template<class T> void insert_top_cell(Node_Stellar& n, T& t, int topsimplex_index, int vector_position);
    /**
     * @brief A private method that split a node, creating the sons node, following the current division type
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param domain a Box& argument, represents the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     *
     * This method also reinsert all the vertices in the splitted node into the sons node.
     */
    template<class C, class T> void split(Node_Stellar& n, Box<C>& domain, int level, Mesh<C, T> &mesh);
    /**
     * @brief A private method that counts the number of vertices in the leaf block n
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     * @return int an integer containing the number of vertices indexed in leaf n
     */
    template<class C, class T> int count_indexed_vertices(Node_Stellar& n, Mesh<C, T> &mesh);
    /**
     * @brief A private method that checkss if a vertex is already into the vertices array of leaf n.
     *
     * The method returns the position index of the vertex if a vertex is already into the vertices array, -1 otherwise
     *
     * @param n a Node_Stellar& argument representing the current node
     * @param v a Vertex& argument representing the target vertex
     * @param mesh a Mesh& argument representing the mesh indexed by the Stellar tree
     * @return int the position index of the vertex if v is already into the array, -1 otherwise
     */
    template<class C, class T> int is_already_inserted(Node_Stellar& n, Vertex<C> &v, Mesh<C, T> &mesh);
};

#include "stellar_tree_generation.h"
#include "stellar_tree_update.h"
#include "stellar_tree_visit.h"

#endif // STELLAR_TREE_H
