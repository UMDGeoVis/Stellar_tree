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

#ifndef CONTRACTION_SIMPLIFIER_CLIQUES_H
#define CONTRACTION_SIMPLIFIER_CLIQUES_H

#include "contraction_simplifier_weak.h"

class Contraction_Simplifier_Top : public Contraction_Simplifier_Weak
{
public:
    Contraction_Simplifier_Top() {}
    /// this procedure simplify the simplicial complex without considering any weight for the edges
    void simplify(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli, ivect &v_orig, string base_file_name);
    /// this procedure simplify the simplicial complex considering a weight for the edges
    void simplify_weighted(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli, ivect &v_orig, string base_file_name);

//    void simplify_tri_qem(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli);
    void simplify_length(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli); /// JUST FOR TRIANGLE MESHES

protected:
    static void simplify_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params);
//    static void simplify_weighted_with_rejected_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, w_contraction_parameters &params);
    static void simplify_weighted_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, w_contraction_parameters &params);
//    static void simplifyW_no_homology_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, w_contraction_parameters &params);
//    static void simplify_tri_qem_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params);
    static void simplify_tri_length_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params); /// FOR TRIANGLE AND TETRA MESHES

//    static bool check_edge();
    /// the functions check an edge for applying a contraction on it
    /// returns 1 if the contraction happens
    /// -1 if the link condition is not verified
    /// 0 if the edge has been skipped (as one of the two extrema is deleted)
    static int check_weighted_edge(edge_weight &ew, edge_pqueue &edges, leaf_VT &vtops, Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, w_contraction_parameters &params);

    static bool is_tops_link_condition_valid(ivect &e, ET &et, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, contraction_parameters &params);
    static bool is_tops_link_condition_valid_v2(ivect &e, ET &et, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, contraction_parameters &params);

    static void get_vertices_in_link(int v_id, VT &vt, iset &v_in_link, Simplicial_Mesh &mesh);
    static void get_vertices_in_link(const ivect &e, ET &et, iset &v_in_link, Simplicial_Mesh &mesh);

    static void get_shared_simplices(VT &vt, iset &v_intersection, int d, Simplicial_Mesh &mesh, set<ivect> &shared_s);

};

#endif // CONTRACTION_SIMPLIFIER_CLIQUES_H
