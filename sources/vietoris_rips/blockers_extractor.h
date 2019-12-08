#ifndef BLOCKERS_EXTRACTOR_H
#define BLOCKERS_EXTRACTOR_H

#include "vietoris_rips_generator.h"
#include "topological_ds_generator/skeleton_generator.h"
#include "simplification/contraction_simplifier_top.h"

//typedef Top_Simplex Blocker; /// convinience alias
typedef set<Top_Simplex> blocker_set;
typedef vector<Top_Simplex> blocker_array;
typedef vector<set<ivect> > failing_array;

struct fs_stats
{

};

class Blockers_Extractor : public VietorisRips_Generator, public Contraction_Simplifier_Top
{
public:
    Blockers_Extractor() {}

    blocker_set extract_blockers(Stellar_Tree &tree, Simplicial_Mesh &mesh, bool extract_from_maximal_cliques, ivect &orig_v_pos, string blockers_file_name);
    void add_blockers(Stellar_Tree &tree, Simplicial_Mesh &mesh, vector<blocker_array> &b_arr);
    void extract_aux_structures_statistics(Stellar_Tree &tree, Simplicial_Mesh &mesh);

private:    
    static void extract_blockers_from_maximal_cliques(top_simplices_map &tops, leaf_VT &all_vt, Simplicial_Mesh &mesh, blocker_set &p); ///GLOBAL - [SLOW] - CLIQUES BASED
    static void extract_blockers(failing_array &failing_simplices, leaf_VT &all_vt, Simplicial_Mesh &mesh, blocker_set &p); ///GLOBAL - FAILING SIMPLICES BASED
    static void extract_leaf_blockers(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<blocker_set, leaf_VT, pair<int, int> > &tupla); /// BUGGY - hierarchy-based
    static void add_blockers_to_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, vector<ivect> &blockers_id); /// hierarchy-based
    static void get_blockers_stats(blocker_set &blockers);
    static void write_blockers_file(blocker_set &blockers, string file_name, ivect &orig_v_pos);

    static void extract_global_VTops(leaf_VT &all_vt, Simplicial_Mesh &mesh); ///GLOBAL

    /// these two procedures are needed just for cliques computation
    static void extract_undirected_1skeleton(Node_Stellar &n, Simplicial_Mesh &mesh, skeleton_graph &graph, unordered_map<int,iset> &ov);
    static void extract_undirected_1skeleton(Simplicial_Mesh &mesh, skeleton_graph &graph); ///GLOBAL

    static void get_local_maximal_cliques(skeleton_graph &graph, unordered_map<int, iset> &outer_vertices, top_simplices_map &tops, Node_Stellar &n);
    static void get_maximal_cliques(skeleton_graph &graph, top_simplices_map &tops, Simplicial_Mesh &mesh); ///GLOBAL
    static void get_maximal_cliques_stats(top_simplices_map &tops);

    static bool is_blocker(Top_Simplex &blocker_candidate, int current_d, Node_Stellar &n, Simplicial_Mesh &mesh, blocker_set &bs);
    static bool is_blocker(Top_Simplex &blocker_candidate, int current_d, leaf_VT &all_vt, Simplicial_Mesh &mesh, blocker_set &bs); ///GLOBAL - CLIQUES BASED
    static bool is_blocker(Top_Simplex &blocker_candidate, int current_d, leaf_VT &all_vt, Simplicial_Mesh &mesh, failing_array &failing_simplices, blocker_set &bs); ///GLOBAL - FAILING SIMPLICES BASED

    /// this procedure checks if f is in a boundary of a top d+1 simplex or if f is a top simplex itself
    static bool is_in_Sigma(ivect &f, int current_d, Node_Stellar &n, Simplicial_Mesh &mesh);
//    static bool is_in_Sigma(int_vect &f, int current_d, Simplicial_Mesh &mesh, bool debug); ///GLOBAL - CLIQUES BASED
    static bool is_in_Sigma(ivect &f, int current_d, leaf_VT &all_vt, Simplicial_Mesh &mesh, bool debug); ///GLOBAL

    /// for debug only
    static void check_blockers(blocker_set &bs, Simplicial_Mesh &mesh);

    //// this procedure extracts all the top-simplices
    static void extract_failing_simplices(skeleton_graph &graph, Simplicial_Mesh &mesh, leaf_VT &all_vt, vector<set<ivect> > &failing_simplices);
//    static void extract_failing_simplices(Node_Stellar &n, skeleton_graph &graph, Simplicial_Mesh &mesh, leaf_VT &all_vt, vector<set<int_vect> > &failing_simplices);
    static void extract_failing_simplices(ivect &e, ET &et, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, vector<set<ivect> > &failing_simplices);
    static void get_failing_simplices_stats(vector<set<ivect> > &failing_simplices);
    static void get_failing_simplices_stats(vector<set<ivect> > &failing_simplices, pair<int,int> &s);

};

#endif // BLOCKERS_EXTRACTOR_H
