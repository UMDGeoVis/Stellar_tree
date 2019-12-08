#ifndef CLI_PARAMETERS_H
#define CLI_PARAMETERS_H

#include <string>
using namespace std;

#define DEFAULT_X_PER_LEAF -1
#define DEFAULT "null"
#define DEFAULT_SIZE 100

struct cli_parameters
{
    string mesh_path, exe_name;
    bool get_index_stats, count_sub_simplices;
    int vertices_per_leaf;

    string kind_of_operation;
//    ReindexType reindex_type;

    bool debug_mode;
    bool debug_prints_mode;
    bool is_cp_mesh;

    bool verbose_encoding;

    bool save_original_indices;

    double eps;

    int cache_size;

    string weights_file;
//    bool no_homology_preservation;
    string extract_blockers;
    bool /*qem,*/ len, force3D;

    cli_parameters()
    {
        vertices_per_leaf = DEFAULT_X_PER_LEAF;
        get_index_stats = false;
        count_sub_simplices = false;
        kind_of_operation = DEFAULT;
//        reindex_type = ALL_SAVE;

        debug_mode = false;
        debug_prints_mode = false;
        is_cp_mesh = false;

        verbose_encoding = false;
        save_original_indices = false;
        cache_size = DEFAULT_SIZE;

        weights_file = DEFAULT;
//        no_homology_preservation = false;
        extract_blockers = DEFAULT;
//        qem = false;
        len = false;
        force3D = false;
    }

    inline bool exec_weighted_simplification() { return weights_file != DEFAULT; }
};

#endif // CLI_PARAMETERS_H
