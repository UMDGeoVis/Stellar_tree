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

#ifndef CONTRACTION_SIMPLIFIER_H
#define CONTRACTION_SIMPLIFIER_H

#include "stellar_tree/mesh.h"
#include "utilities/lru_cache.h"
#include "stellar_tree/stellar_tree.h"
#include "topological_ds/links_aux_structures.h"
#include "statistics/statistics.h"
#include "utilities/container_utilities.h"
#include "stellar_tree/mesh_updater.h"
#include "utilities/usage.h"
#include "utilities/cli_parameters.h"
#include "utilities/string_management.h"
#include "topological_queries/topological_query_extractor.h"

#include "simplification_aux_structures.h"

/// Contraction_Simplifier class ///
class Contraction_Simplifier_Weak
{
public:
    Contraction_Simplifier_Weak() {}
    /// this procedure simplify the simplicial complex without considering any weight for the edges
    void simplify(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli);
    void simplify_length(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli); /// JUST FOR TRIANGLE MESHES

    static void update_mesh_and_tree(Stellar_Tree &tree, Simplicial_Mesh &mesh,contraction_parameters &params);

protected:
    static void simplify_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params);
    static void simplify_tri_length_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params); /// JUST FOR TRIANGLE MESHES

    static void extract_edges_lengths(Node_Stellar &n, Simplicial_Mesh &mesh, contraction_parameters &params); /// JUST FOR TRIANGLE MESHES
    static void extract_edge_length_stats(Node_Stellar &n, Simplicial_Mesh &mesh, contraction_parameters &params); /// JUST FOR TRIANGLE MESHES

    /// initialize the queue
    static void extract_target_edges(Node_Stellar &n, Simplicial_Mesh &mesh, edge_queue &edges);
    static void extract_target_edges(Node_Stellar &n, Simplicial_Mesh &mesh, edge_pqueue &edges, w_contraction_parameters &params);
    static void extract_target_edges_length(Node_Stellar &n, Simplicial_Mesh &mesh, edge_pqueue &edges, contraction_parameters &params); /// FOR TRIANGLE AND TETRA MESHES

    /// skips an edge that has one of the two extreme deleted
    inline static bool skip_edge(const ivect &e, Simplicial_Mesh &mesh) /// the procedure checks if one of the two extreme is already deleted
    {
        return (mesh.is_vertex_removed(e[0]) || mesh.is_vertex_removed(e[1]));
    }
    /// this procedure checks the weak link condition as defined by [Attali et al., 2011]
    /// reference paper:
    ///
    /// for each p-dimensional simplex we extract the link of e and of its two extremes
    /// then we check the p-link condition
    /// if all the p-link conditions are valid then we can contract the edge
    /// (1) the top-down version first checkes the p-condition on the vertices then moves to
    ///     d-simplices and checks iteratively going down up to the 1-simplices
    static bool is_weak_link_condition_valid_top_down_version(ivect &e, ET &et, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, contraction_parameters &params);
//    /// (2) the bottom-up version checks the p-condition on the vertices, then on 1-simplices ... up to d-simplices (on average SLOWER)
//    static bool is_weak_link_condition_valid_bottom_up_version(ivect &e, ET &et, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, contraction_parameters &params);

    ///
    static void contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, leaf_VT &vtops, Node_Stellar &outer_v_block, edge_queue &edges,
                              Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params);
    ///
    static void contract_weighted_edge(ivect &e, double e_weight, ET &et, VT &vt0, VT &vt1, leaf_VT &vtops, Node_Stellar &outer_v_block, edge_pqueue &edges,
                              Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, w_contraction_parameters &params);
    static void contract_edge_length(ivect &e, ET &et, VT &vt0, VT &vt1, leaf_VT &vtops, Node_Stellar &outer_v_block, edge_pqueue &edges,
                                     Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params); /// JUST FOR TRIANGLE MESHES

    /// initializes the VTop and ETop of the edge and its two vertices
    /// plus saves the second leaf block if we are processing a cross-edge
    static void get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_Stellar *& outer_v_block,
                                   Node_Stellar &n, Simplicial_Mesh &mesh, leaf_VT &vtops, LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params);
    /// the VTop is always without removed top-simplices
    static VT* get_VTop(int v_id, Node_Stellar &n, Simplicial_Mesh &mesh, leaf_VT &vtops, LRU_Cache<int,leaf_VT> &cache,
                        Stellar_Tree &tree, Node_Stellar *& v_block, contraction_parameters &params);
    static void get_ETop(ivect &e, ET &et, Node_Stellar &n, Simplicial_Mesh &mesh, leaf_VT &vtops);
    static void get_link_lite(int v_id, VT &vt, s_link &link, Simplicial_Mesh &mesh); /// extracts the vertices and the d-1 simplices (in any dim) in a link of a vertex
    static void get_link_lite(const ivect &e, ET &et, s_link &link, Simplicial_Mesh &mesh); /// extracts the vertices and the d-1 simplices (in any dim) in a link of an edge
    ///
    static void extract_sub_d_simplices(int d, s_link &link);
    static void extract_d_simplices(int d, VT& cob, s_link &link, Simplicial_Mesh &mesh);

    /// the procedure removes from the VTop relation the deleted top d-simplices
    static void clean_coboundary(VT &cob, Simplicial_Mesh &mesh);

    static void check_and_remove_top_edge(const ivect &e, int &e_top_id, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, contraction_parameters &params);
    /// the procedure checks if some boundary d-1 faces have not an adjacent d-top simplex. if not, we add that face to the list of top simplices
    static void check_for_new_tops(int v, VT &vt, ET &et, Node_Stellar &n, Simplicial_Mesh &mesh,
                                   leaf_VT &vtops, LRU_Cache<int,leaf_VT> &cache, Stellar_Tree &tree, contraction_parameters &params);
    static bool exists_simplex_in_coboundary(int start_d, const ivect &s, VT &vt, Simplicial_Mesh &mesh);
    static bool exists_simplex_in_coboundary_v2(int start_d, ivect &s, vector<vector<ivect> > &explicit_vt);
    static void add_new_top_in_VTop(int v_id, int d, int t_id, Node_Stellar &n, vector<Node_Stellar *> &leaves, leaf_VT &vtops,
                                    LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params);

    /// the procedure updates
    /// (1) the VTop relation of the surviving vertex
    /// (2) the top simplices for which we change the v2 with v1
    /// (3) checks and updates the top list in the target leaf block n (if a top simplex is now incident in it)
    /// (4) and add the new edges to the edge queue
    static void update(const ivect &e, VT& vt, VT& difference, Node_Stellar &n, Node_Stellar &v_block, edge_queue &edges,
                                          Simplicial_Mesh &mesh, contraction_parameters &params);
    static bool update(int v1, int v2, double e_weight, VT& vt1, VT& vt2, Node_Stellar &target, Node_Stellar &origin, edge_pqueue &edges,
                                          Simplicial_Mesh &mesh, w_contraction_parameters &params);
    static void update_length(const ivect &e, VT& vt, VT& difference, Node_Stellar &n, Node_Stellar &v_block, edge_pqueue &edges,
                                          Simplicial_Mesh &mesh, contraction_parameters &params); /// JUST FOR TRIANGLE MESHES
    static void remove_from_mesh(int to_delete_v, int e_top_id, ET &et, Simplicial_Mesh &mesh, contraction_parameters &params);

//    /// this procedure removes duplicated top simplices obtained after a collapse
//    ///    this check is needed only by the procedure that do not preserve the homology
//    ///    as we can have "weird" situations that are avoided by using top-link condition
//    static void remove_duplicated_tops(VT& vt, Simplicial_Mesh &mesh);
//    static void remove_degenerate_tops(VT& vt, Simplicial_Mesh &mesh, bool print);

    /// just for debug
//    static void explicit_top_coboundary(VT &vt, vector<vector<ivect> > &explicit_top_cob, Simplicial_Mesh &mesh);

    /// FOR DEBUG ONLY
    static int get_top_edge_id(const ivect &e, VT &vt0, Simplicial_Mesh &mesh);
    static bool exists_degenerate_tops(Node_Stellar &n, Simplicial_Mesh &mesh);
    static bool exists_degenerate_tops(VT& vt, Simplicial_Mesh &mesh);
    static bool exists_duplicated_tops(VT& vt, Simplicial_Mesh &mesh);
    static bool everything_is_synced(const ivect &e, VT &vt0, Node_Stellar &n, Node_Stellar *&v_block, Stellar_Tree &tree, LRU_Cache<int,leaf_VT> &cache);
    static bool cache_is_synced(Node_Stellar &n, Simplicial_Mesh &mesh, Stellar_Tree &tree, LRU_Cache<int,leaf_VT> &cache);
    static void compute_cache_stats(LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params);
    static void print_verbose_coboundary_relation(VT &vt, Simplicial_Mesh &mesh);
    static bool is_degenerate_edge(const ivect &e,Node_Stellar &n, Simplicial_Mesh &mesh);
};

#endif // CONTRACTION_SIMPLIFIER_H
