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

#ifndef _STATISTICS_H
#define	_STATISTICS_H

#include <vector>
#include <string>
#include <set>
#include <sstream>
#include <cmath>

#include "io/writer.h"
#include "io/reader.h"
#include "stellar_tree/mesh.h"
#include "index_statistics.h"
#include "utilities/string_management.h"
#include "stellar_tree/stellar_tree.h"

using namespace std;

/**
 * @brief A class representing a statistical object used to compute statistics over a spatial index
 *
 */
class Statistics
{
public:
    /**
     * @brief A constructor method
     *
     */
    Statistics() { this->indexStats = IndexStatistics(); }    
    /**
     * @brief A constructor method
     *
     * @param orig
     */
    Statistics(const Statistics& orig) { this->indexStats = orig.indexStats; }
    /**
     * @brief A destructor method
     *
     */
    virtual ~Statistics() {}    
    /**
     * @brief A public method that computes index statistics over a given tree
     *
     * @param tree is the tree on which we want to compute the statistics
     * @param mesh is the mesh indexed by the tree
     * @param
     */
    template<class C, class T> void get_index_statistics(Stellar_Tree& tree, Mesh<C, T> &mesh, bool compute_run_histogram);
    ///
    /**
     * @brief A public method that computes the number of simplices (encoded or not in the tree and by the mesh) without duplicates
     *
     * @param tree is the tree on which we want to compute the statistics
     * @param mesh is the mesh indexed by the tree
     */
    template<class C,class T> void get_unique_cp_cells_counters(Stellar_Tree& tree, Mesh<C,T>& mesh);

protected:
    ///A protected variable in which the index statistics can be saved
    IndexStatistics indexStats;
    /**
     * @brief A protected method that simulates a tree visit and calculates the spatial index statistics
     *
     * @param n is the current node of the tree to visit
     * @param n_level is the level of the node in the tree
     * @param mesh is the mesh indexed by the tree
     */
    template<class C,class T> void visit_tree(Node_Stellar& n, int n_level, Mesh<C,T>& mesh, bool compute_run_histogram);
    /**
     * @brief A protected method that initializes the vectors used for the index statistics extraction
     *
     * @param mesh is the mesh indexed by the tree
     */
    template<class C,class T> void init_vector(Mesh<C,T>& mesh);
    /**
     * @brief A protected method that calculates index statistics which is not possible to compute during the tree visit
     *
     */
    void calc_remaining_index_stats();
    /**
     * @brief A protected method that checks common inconsistencies of index statistics operation
     *
     */
    void check_inconsistencies();
    ///
    /**
     * @brief A protected method that computes the statistics of a leaf block
     *
     * @param n represents the current leaf block
     * @param mesh is the mesh indexed by the tree
     */
    template<class C, class T> void get_leaf_stats(Node_Stellar& n, Mesh<C,T>& mesh, bool compute_run_histogram);
    /**
     * @brief A protected method that computes the unique cells number on a Constant-Polytope mesh
     *
     * @param tree is the tree on which we want to compute the statistics
     * @param mesh is the mesh indexed by the tree
     * @param p a pair containing the parameters needed by the inner wrapper procedures. first contains the dimension of the cells to extract, while second the unique number of these cells
     */
    void get_unique_cells_counter(Stellar_Tree &tree, CP_Mesh& mesh, pair<int,long long> &p);
    /**
     * @brief A protected method that computes the unique cells number on a Simplicial mesh
     *
     * @param tree is the tree on which we want to compute the statistics
     * @param mesh is the mesh indexed by the tree
     * @param p a pair containing the parameters needed by the inner wrapper procedures
     */
    void get_unique_cells_counter(Stellar_Tree& tree, Simplicial_Mesh& mesh, pair<int, long long> &p);
    /**
     * @brief A protected method that computes the unique cells number on a Simplicial mesh (wrapper function)
     *
     * @param n represents the current leaf block
     * @param mesh is the mesh indexed by the tree
     * @param p a pair containing the parameters needed by the inner wrapper procedures. first contains the dimension of the cells to extract, while second the unique number of these cells
     */
    static void get_unique_dcells_counter_CP_wrapper(Node_Stellar& n, CP_Mesh& mesh, pair<int,long long> &p)
    { Statistics::get_unique_dcells_counter(n,mesh,p.first,p.second); }
    /**
     * @brief A protected method that computes the unique cells number on a Simplicial mesh (wrapper function)
     *
     * @param n represents the current leaf block
     * @param mesh is the mesh indexed by the tree
     * @param p a pair containing the parameters needed by the inner wrapper procedures. first contains the dimension of the cells to extract, while second the unique number of these cells
     */
    static void get_unique_dcells_counter_Simplicial_wrapper(Node_Stellar& n, Simplicial_Mesh& mesh, pair<int,long long> &p)
    { Statistics::get_unique_dcells_counter(n,mesh,p.first,p.second); }
    /**
     * @brief A protected method that computes the unique cells number on a mesh (main function)
     *
     * @param n represents the current leaf block
     * @param mesh is the mesh indexed by the tree
     * @param d contains the dimension of the cells to extract
     * @param num contains the unique number of these cells
     */
    template<class C, class T> static void get_unique_dcells_counter(Node_Stellar& n, Mesh<C,T>& mesh, int d, long long &num);
};

#include "statistics_index.h"

#endif	/* _STATISTICS_H */


