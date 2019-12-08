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

#ifndef VIETORIS_RIPS_GENERATOR_H
#define VIETORIS_RIPS_GENERATOR_H

#include <map>
#include <unordered_map>
#include <boost/function.hpp>
#include <tuple>

#include "utilities/container_utilities.h"
#include "vietoris_rips_aux_structures.h"
#include "topological_ds/skeleton.h"

typedef map<int, list<ivect> > top_simplices_map;

class VietorisRips_Generator
{
public:
    VietorisRips_Generator() { }

    /// extract the maximal cliques from a neighborhood graph generated from a point cloud
    void extract_VietorisRips_complex_GLOBAL(Stellar_Tree &tree, Simplicial_Mesh &mesh, double eps, bool debug); /// (hierarchycal) global extraction of n-graph + global extraction of maximal cliques
    void extract_VietorisRips_complex_LOCAL(Stellar_Tree &tree, Simplicial_Mesh &mesh, double eps, bool debug); /// local extraction graph + of maximal cliques (leaf block scope)
    void extract_VietorisRips_complex_HYBRID(Stellar_Tree &tree, Simplicial_Mesh &mesh, double eps, bool debug, string mesh_name); /// global graph + local hierarchy based extraction of maximal cliques
    void extract_VietorisRips_complex_PARALLEL(Stellar_Tree &tree, Simplicial_Mesh &mesh, double eps); /// global graph + local parallel extraction of maximal cliques

    /// extract the maximal cliques from a 1-skeleton
    void extract_VietorisRips_complex_GLOBAL(Simplicial_Mesh &mesh, skeleton_graph &graph, bool debug); /// global extraction of maximal cliques
    void extract_VietorisRips_complex(Stellar_Tree &tree, Simplicial_Mesh &mesh, skeleton_graph &graph, bool debug); /// local hierarchy based extraction of maximal cliques
    void extract_VietorisRips_complex_PARALLEL(Stellar_Tree &tree, Simplicial_Mesh &mesh, skeleton_graph &graph); /// local parallel extraction of maximal cliques

    /// FOR DEBUG ONLY
    void validity_check(Stellar_Tree &tree, Simplicial_Mesh &mesh);

protected:
    static void extract_neighborhood_graph(Node_Stellar &n, Box<COORDBASETYPE> &n_domain, Simplicial_Mesh &mesh, pair<skeleton_graph, Stellar_Tree *> &p);    
    static void extract_neighborhood_graph_PARALLEL(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<skeleton_graph, Stellar_Tree *, map<int,Box<COORDBASETYPE> > > &p);

    static void get_leaf_block_domain(Node_Stellar &n, Box<COORDBASETYPE> &n_domain, Simplicial_Mesh &mesh, map<int,Box<COORDBASETYPE> > &domains);

    static void get_inner_nearest_vertices(int v_id, Vertex<COORDBASETYPE> &v, skeleton_graph &graph, Node_Stellar &n, Simplicial_Mesh &mesh, int offset);
    static void get_outer_nearest_vertices(Node_Stellar &pivot_n, Box<COORDBASETYPE> &pivot_domain, skeleton_graph &graph,
                                    Node_Stellar &n, Box<COORDBASETYPE> &dom, int level, Simplicial_Mesh &mesh, Spatial_Subdivision& subdiv, int offset);
    static void get_all_outer_nearest_vertices(Node_Stellar &pivot_n, Box<COORDBASETYPE> &pivot_domain, skeleton_graph &graph,
                                               Node_Stellar &n, Box<COORDBASETYPE> &dom, int level, Simplicial_Mesh &mesh, Spatial_Subdivision& subdiv, int offset,
                                               unordered_map<int, iset> &outer_vertices);

    void get_maximal_cliques(skeleton_graph &graph, Simplicial_Mesh &mesh, bool debug);
    static void get_maximal_cliques(skeleton_graph &graph, int offset, iset &setR, iset setP, iset &setX,
                             map<int, list<ivect > > &top_simplexes_local, bool debug, unsigned &maxX);

    static void extract_VR_complex(Node_Stellar &n, Box<COORDBASETYPE> &n_domain, Simplicial_Mesh &mesh, local_generation_parameters &p/*pair<double, Stellar_Tree *> &p*/);

    /// for local strategy
    static void get_local_maximal_cliques(skeleton_graph &graph, unordered_map<int,iset> &outer_vertices, Simplicial_Mesh &mesh, Node_Stellar &n, local_generation_parameters &p);
    static void get_local_maximal_cliques_rec(skeleton_graph &graph, unordered_map<int,iset> &outer_vertices, iset &setR, iset setP, iset &setX,
                                              map<int, list<ivect > > &top_simplexes_local, Node_Stellar &n, int offset, local_generation_parameters &p);
    /// for parallel and hybrid strategies
    static void get_local_maximal_cliques_v2(Node_Stellar &n, Simplicial_Mesh &mesh, skeleton_graph &graph);
    static void get_local_maximal_cliques_v2_rec(skeleton_graph &graph, iset &setR, iset setP, iset &setX,
                                                 map<int, list<ivect > > &top_simplexes_local, Node_Stellar &n, int offset);

    /// FOR DEBUG ONLY    
    static void extract_skeleton(Node_Stellar &n, Simplicial_Mesh &mesh, skeleton_graph &graph);
    void check_maximal_cliques(skeleton_graph &graph, Simplicial_Mesh &mesh);
};

#endif // VIETORIS_RIPS_GENERATOR_H

