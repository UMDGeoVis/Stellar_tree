#include "main_utility_functions.h"
#include "main_simplicial_functions.h"
#include "main_cp_functions.h"

template<class C, class T> void load_index(Stellar_Tree& tree, Mesh<C,T> &mesh, Statistics& stats, Reindexer &reindexes, cli_parameters &cli);
template<class M> void exec_topo_queries(Stellar_Tree& tree, M &mesh, cli_parameters &cli);

int main(int argc, char** argv)
{
    if(argc == 1)
    {
        print_help();
        return 0;
    }

    cli_parameters cli = cli_parameters();

    if (read_arguments(argc, argv,cli) == -1)
    {
        print_usage();
        return (EXIT_FAILURE);
    }

    /// if kv=INT_MAX we do not want to have a hierarchical decomposition
    int kv = (cli.vertices_per_leaf == 0) ? INT_MAX : cli.vertices_per_leaf;

    Stellar_Tree tree = Stellar_Tree(/*cli.vertices_per_leaf*/kv);
    Reindexer reindexer = Reindexer(/*variables.reindex_type*/);
    Statistics stats = Statistics();

    // SIMPLICIAL MESHES -> TRI - TETRA - SIERP - GMV - nDim
    if(is_simplicial_mesh(cli))
    {
        cerr << "==Simplicial_Mesh==" <<endl;

        Simplicial_Mesh mesh = Simplicial_Mesh();

        load_mesh(mesh,cli);
        cerr << "RAM peak for Indexing the mesh: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//        mesh.print_mesh_stats(cerr);
//        cerr << mesh.get_domain() << endl;
        load_index(tree,mesh,stats,reindexer,cli);
        cerr << "RAM peak for encoding the Stellar tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

        // TOPOLOGICAL RELATION EXTRACTION
        exec_topo_queries(tree,mesh,cli);
    }

    return (EXIT_SUCCESS);
}

template<class C, class T> void load_index(Stellar_Tree &tree, Mesh<C,T> &mesh, Statistics &stats, Reindexer &reindexer, cli_parameters &cli)
{
    // init the division variable in tree
    tree.set_subdivision(mesh);

    stringstream out;
    string kv = (cli.vertices_per_leaf == 0) ? "inf" : to_string(cli.vertices_per_leaf);
    // NOTA: the verbose encoding is needed/used by the contraction and reduction simplifiers
    out << get_path_without_file_extension(cli.mesh_path) << "_v_" << kv << "_" << get_file_extension(cli.mesh_path);
    if(cli.verbose_encoding)
        out << "_verbose";
    if(cli.force3D)
        out << "_force3D";


    out << "_.stellar";
    cout << out.str() << endl;

    // first we check if exists a basic Stellar tree already built using the given kv
    if (!Reader::read_tree(tree.get_root(), tree.get_subdivision(), out.str()))
    {
        // if not we have to build the index with the given threshold
        tree.build(mesh,out.str(),reindexer,cli.save_original_indices);
    }
    else
    {
        // we need to explicitly exploit the spatial coherence
        stringstream base_info;
        base_info << kv << " Stellar-tree [Reordering Index] ";
        Timer t;
        t.start();
        // if we want to build the global IA/IA* data structure
        // we need the original position indexes of vertices and of top cells in order to save them consistently in the output file
        reindexer.reorganize_index_and_mesh(tree,mesh,cli.save_original_indices);
        t.stop();
        t.print_elapsed_time(base_info.str());
    }

    tree.init_leaves_list(tree.get_root()); // NOTA: enables the parallel OpenMP-based processing of the leaf block

    if (cli.get_index_stats)
        stats.get_index_statistics(tree,mesh,cli.debug_mode);
    if (cli.count_sub_simplices)
        stats.get_unique_cp_cells_counters(tree,mesh);
}

template<class M> void exec_topo_queries(Stellar_Tree &tree, M &mesh, cli_parameters &cli)
{
    extract_VV_relations(tree,mesh);
    cerr << "RAM peak for extracting VV relations: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
    extract_VE_relations(tree,mesh);
    cerr << "RAM peak for extracting VE relations: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
    extract_VF_relations(tree,mesh);
    cerr << "RAM peak for extracting VF relations: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
    extract_VTop_relations(tree,mesh);
    cerr << "RAM peak for extracting VT relations: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    extract_EF_relations(tree,mesh);
    cerr << "RAM peak for extracting EF relations: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
    extract_ETop_relations(tree,mesh);
    cerr << "RAM peak for extracting ET relations: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    extract_FTop_relations(tree,mesh);
    cerr << "RAM peak for extracting FT relations: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
}
