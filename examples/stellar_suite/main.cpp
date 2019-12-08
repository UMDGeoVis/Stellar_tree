#include "main_utility_functions.h"
#include "main_simplicial_functions.h"
#include "main_cp_functions.h"

template<class C, class T> void load_index(Stellar_Tree& tree, Mesh<C,T> &mesh, Statistics& stats, Reindexer &reindexes, cli_parameters &cli);
template<class M> void exec_topo_queries(Stellar_Tree& tree, M &mesh, cli_parameters &cli);

// NOTA: old approach, the visit function is implemented into the Generator class
template<class C, class T> int extract_half_edges(Stellar_Tree& tree, Mesh<C,T>& mesh, cli_parameters &cli);
// NOTA: old approach, the visit function is implemented into the Generator class
template<class C, class T> int extract_incidence_graph(Stellar_Tree& tree, Mesh<C,T>& mesh, cli_parameters &cli);
// NOTA: old approach, the visit function is implemented into the Generator class
int extract_vietoris_rips_complex(Stellar_Tree& tree, Statistics& stats, Reindexer &reindexer, cli_parameters &cli);

//EXPERIMENTAL FEATURE..
// NOTA: old approach, the visit function is implemented into the Repair class
template<class T,class C, class Top> int repair_mesh(T& tree, Mesh<C,Top>& mesh, cli_parameters &cli);

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

    if(cli.kind_of_operation == "vrG" || cli.kind_of_operation == "vrL"
            || cli.kind_of_operation == "vrH" || cli.kind_of_operation == "vrP")
        // if we want to build a Vietoris-Rips complex, then we have a dedicated main
    {
        extract_vietoris_rips_complex(tree,stats,reindexer,cli);
    }
    // SIMPLICIAL MESHES -> TRI - TETRA - SIERP - GMV - nDim
    else if(is_simplicial_mesh(cli))
    {
        cerr << "==Simplicial_Mesh==" <<endl;

        Simplicial_Mesh mesh = Simplicial_Mesh();

        load_mesh(mesh,cli);
        cerr << "RAM peak for Indexing the mesh: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//        mesh.print_mesh_stats(cerr);
//        cerr << mesh.get_domain() << endl;
        load_index(tree,mesh,stats,reindexer,cli);
        cerr << "RAM peak for encoding the Stellar tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

        /// get the adjacency graph -- DISABLE AFTER DEBUGGING ///////////////////////////////
//        if
//        {
//            stringstream ss;
//            ss << get_file_name(cli.mesh_path) << "_kv_";
//            if(cli.vertices_per_leaf == 0)
//                ss << "inf";
//            else
//                ss << cli.vertices_per_leaf;
//    //            mesh.print_mesh_adj_graph(ss2.str());
////            /// prints the 1-skeleton indexed by each block
////            string f = ss.str();
////            tree.visit(Adjacency_Graph_Extractor::extract_local_1skeleton_graph_Simplicial_wrapper,tree.get_root(),mesh,f);

//            /// prints the 1-skeleton of the complex
//            cerr<<"[OUTPUT] compute the GLOBAL 1-skeleton within each leaf block of the INITIAL mesh"<<endl;
//            skeleton_graph n_graph = skeleton_graph(mesh.get_vertices_num());
//            tree.visit(skeleton_generator::extract_skeleton_Simplicial_wrapper,tree.get_root(),mesh,n_graph);
//            n_graph.write_gv_file(get_file_name(cli.mesh_path));

////            /// prints an adjacency graph for each leaf block
////            Local_Adjacency_Graph g(ss.str(),mesh);
////            tree.visit(Adjacency_Graph_Extractor::extract_local_adjacency_graph_Simplicial_wrapper,tree.get_root(),mesh,g);

//            /// prints an adjacency graph for the complete complex
//            cerr<<"[OUTPUT] computing the GLOBAL adjacency matrix of the INITIAL mesh"<<endl;
//            Adjacency_Graph gg(mesh);
//            tree.visit(Adjacency_Graph_Extractor::extract_lite_adjacency_graph_Simplicial_wrapper,tree.get_root(),mesh,gg);
//            gg.write_graph_file(get_file_name(cli.mesh_path),mesh);
//            gg.write_fixed_graph_file(get_file_name(cli.mesh_path),mesh);

//            /// TODO compute graphs statistics
//            cerr<<"===[TODO] COMPUTE GRAPHS STATISTICS==="<<endl;
//        }
        ///////////////////////////////////////////////////////////////////////////////////

        // TOPOLOGICAL RELATION EXTRACTION
        exec_topo_queries(tree,mesh,cli);
        // TOPOLOGICAL DATA STRUCTURE GENERATION
        extract_adjacencies(tree,mesh,reindexer,cli);
        if(get_file_extension(cli.mesh_path) == "tri")
            extract_half_edges(tree,mesh,cli);
        extract_incidence_graph(tree,mesh,cli);
        // VALIDATOR
        validate_connectivity(tree,mesh,cli);
        if(get_file_extension(cli.mesh_path) == "ts" || get_file_extension(cli.mesh_path) == "tri")
        {
            validate_link_conditions(tree,mesh,cli); // the "classical" link conditions make sense only on low dimensional simplices
            /// to-do: add gmv support
        }
        //SIMPLIFICATION
        if(get_file_extension(cli.mesh_path) == "ts" ||
                get_file_extension(cli.mesh_path) == "tri" ||
                get_file_extension(cli.mesh_path) == "off")
        {
            simplify(tree,reindexer,mesh,cli,stats);
        }        
        // FOR DEBUGGING THE VRIPS GENERATOR
        if(cli.kind_of_operation == "vrDebug")
        {
            VietorisRips_Generator vr;
            vr.validity_check(tree,mesh);
        }
        //EXPERIMENTAL FEATURE
        repair_mesh(tree,mesh,cli);
        //
        extract_global_skeleton(tree,mesh,cli);        
    }
    // CP MESHES -> QUAD - HEXA
    else if(is_cp_mesh(cli))
    {
        cerr << "==CP_Mesh==" <<endl;

        CP_Mesh mesh = CP_Mesh();

        load_mesh(mesh,cli);
        cerr << "RAM peak for Indexing the mesh: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        load_index(tree,mesh,stats,reindexer,cli);
        cerr << "RAM peak for encoding the Stellar tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        // TOPOLOGICAL RELATION EXTRACTION
        exec_topo_queries(tree,mesh,cli);
        // TOPOLOGICAL DATA STRUCTURE GENERATION
        extract_adjacencies(tree,mesh,reindexer,cli);
        if(get_file_extension(cli.mesh_path) == "off")
            extract_half_edges(tree,mesh,cli);
        extract_incidence_graph(tree,mesh,cli);
        // VALIDATOR
        validate_connectivity(tree,mesh,cli);
        if(get_file_extension(cli.mesh_path) == "off") // JUST FOR QUAD-MESHES
        {
            validate_link_conditions(tree,mesh,cli);
            //EXPERIMENTAL FEATURE
            repair_mesh(tree,mesh,cli);
        }
    }
    // GET INDEXED REPRESENTATION FROM SOUP + APPS
    else if(get_file_extension(cli.mesh_path) == "soup")
    {
        Cells_Soup soup;
        Simplicial_Mesh mesh;

        if (!Reader::read_soup(soup, cli.mesh_path))
        {
            cerr << "Error Loading soup file. Execution Stopped." << endl;
            return (EXIT_FAILURE);
        }

        tree.set_subdivision(soup.get_cell(0,1).get_vertex(0).get_dimension());
        tree.build(soup,mesh,reindexer,(cli.kind_of_operation == "adjG"));

        if (cli.get_index_stats)
            stats.get_index_statistics(tree,mesh,true);
        if (cli.count_sub_simplices)
            stats.get_unique_cp_cells_counters(tree,mesh);

        exec_topo_queries(tree,mesh,cli);

        extract_adjacencies(tree,mesh,reindexer,cli);
        extract_incidence_graph(tree,mesh,cli);

        validate_connectivity(tree,mesh,cli);
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
    if(cli.kind_of_operation == "vt")
        extract_VTop_relations(tree,mesh);
    if(cli.kind_of_operation == "vv")
        extract_VV_relations(tree,mesh);
    if(cli.kind_of_operation == "link")
        extract_links_relations(tree,mesh);
    /// pick one of these two strategies for computing p-faces ///
    if(cli.kind_of_operation == "pfaces")
    {
        extract_p_faces(tree,mesh);
    }
    if(cli.kind_of_operation == "pfaces2")
        extract_p_faces_v2(tree,mesh);
}

// NOTA: old approach, the visit function is implemented into the HalfEdge_Generator class
template<class C, class T> int extract_half_edges(Stellar_Tree& tree, Mesh<C,T>& mesh, cli_parameters &cli)
{
    HalfEdge_Generator he_gen;
    Timer time;

    if(cli.kind_of_operation == "heG")
    {
        he_gen = HalfEdge_Generator(mesh);
        time.start();
        he_gen.build_global_HalfEdge(tree.get_root(),mesh,cli.debug_mode);//,debug);
        time.stop();
        time.print_elapsed_time("Build global HalfEdge ");
        he_gen.print_stats();

    }
    else if(cli.kind_of_operation == "heL")
    {
        he_gen = HalfEdge_Generator();
        time.start();
        he_gen.build_HalfEdge(tree.get_root(),mesh);
        time.stop();
        time.print_elapsed_time("Build local HalfEdge ");
        he_gen.print_stats();

    }

    return (EXIT_SUCCESS);
}

// NOTA: old approach, the visit function is implemented into the IG_Generator class
template<class C, class T> int extract_incidence_graph(Stellar_Tree& tree, Mesh<C,T>& mesh, cli_parameters &variables)
{
    IG_Generator ig_generator(mesh.get_implicitly_and_explicitly_encoded_cells_num());
    Timer time;

    if(variables.kind_of_operation == "igG")
    {
        IG ig = IG(mesh.get_implicitly_and_explicitly_encoded_cells_num());
        time.start();
        ig_generator.build_IG(tree.get_root(),mesh,ig);
        time.stop();
        time.print_elapsed_time("Build global IG ");
        ig_generator.print_global_stats(ig);
        ig.clear_IG();
    }
    else if(variables.kind_of_operation == "igL")
    {
        time.start();
        ig_generator.build_IG(tree.get_root(),mesh,variables.debug_mode);
        time.stop();
        time.print_elapsed_time("Build local IG ");
        if(variables.debug_mode)
            ig_generator.print_local_stats();
    }

    return (EXIT_SUCCESS);
}

int extract_vietoris_rips_complex(Stellar_Tree& tree, Statistics& stats, Reindexer &reindexer, cli_parameters &cli)
{
    Simplicial_Mesh mesh = Simplicial_Mesh();
    VietorisRips_Generator vr_generator;

    Reader::read_vertices(mesh,cli.mesh_path);
    /// init the division variable in tree
    tree.set_subdivision(mesh);

    double diagonal = mesh.get_domain().get_diagonal();
//    for(int i=0; i< mesh.get_domain().getMinPoint().get_dimension(); i++)
//        cerr << mesh.get_domain().getMinPoint().getC(i) << " " << mesh.get_domain().getMaxPoint().getC(i) << " "
//             << i << "_side: " << mesh.get_domain().get_side(i) << endl;
    cerr << "[STATS] espilon: " << cli.eps << " diagonal: " << diagonal << " -> ratio: " << cli.eps / diagonal << endl;

    cerr << "RAM peak for loading the vertices: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    tree.generate_spatial_subdivision(mesh,reindexer);
    tree.init_leaves_list(tree.get_root());

    cerr << "RAM peak for generating the spatial subdivision + leaves list: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    ///for debug only
    if(cli.debug_mode)
    {
        mesh.get_vertices_distance_stats();
    }

    if(cli.eps > 0) /// we generate the Vietoris-Rips complex from a point cloud
    {
        if(cli.kind_of_operation == "vrG")
            // GLOBAL graph + GLOBAL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_GLOBAL(tree,mesh,cli.eps,cli.debug_mode);
        else if(cli.kind_of_operation == "vrL")
            // LOCAL graph + LOCAL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_LOCAL(tree,mesh,cli.eps,cli.debug_mode);
        else if(cli.kind_of_operation == "vrH")
            // GLOBAL graph + LOCAL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_HYBRID(tree,mesh,cli.eps,cli.debug_mode,get_path_without_file_extension(cli.mesh_path));
        else if(cli.kind_of_operation == "vrP")
            // GLOBAL graph + PARALLEL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_PARALLEL(tree,mesh,cli.eps);
    }
    else /// otherwise we extract the maximal cliques from a 1-skeleton
    {
        skeleton_graph graph;
        stringstream g_path;
        g_path << get_path_without_file_extension(cli.mesh_path) << ".edges";
//        cout<<"A"<<endl;
        Reader::read_skeleton(graph,mesh,g_path.str());
//        cout<<"A"<<endl;

        if(cli.kind_of_operation == "vrG")
            // GLOBAL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_GLOBAL(mesh,graph,cli.debug_mode);
        else if(cli.kind_of_operation == "vrL" || cli.kind_of_operation == "vrH")
            // LOCAL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex(tree,mesh,graph,cli.debug_mode);
        else if(cli.kind_of_operation == "vrP")
            // PARALLEL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_PARALLEL(tree,mesh,graph);
    }

    cerr << "RAM peak for generating a V-Rips complex: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    stringstream ss;
    ss<<"VRgen_eps_"<<cli.eps;

    if(!cli.debug_mode) //NOT IN DEBUG MODE
    {
        stringstream ss2;
        ss2 << get_path_without_file_extension(cli.mesh_path) << "_" << ss.str() << "_mesh_stats.txt";
        ofstream output(ss2.str());
        mesh.print_mesh_stats(output);
        Writer::write_mesh(get_path_without_file_extension(cli.mesh_path),ss.str(),mesh);
    }

//    if(cli.debug_mode)
//    {
//        cerr<<"computing the adjacency matrix of the GENERATED mesh"<<endl;
//        stringstream ss2;
//        ss2 << get_file_name(cli.mesh_path) << "_" << ss.str();
//        Adjacency_Graph g(mesh);
//        tree.visit(adj_graph_extractor::extract_adjacency_graph_Simplicial_wrapper,tree.get_root(),mesh,g);
//        g.write_graph_file(ss2.str(),mesh);
////        mesh.print_mesh_adj_graph(ss2.str());
//    }

    tree.add_and_resort_top_cells(mesh,reindexer);

    cerr << "RAM peak for add/resort the top cells in the tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    if (cli.get_index_stats)
        stats.get_index_statistics(tree,mesh,false);
    if (cli.count_sub_simplices)
        stats.get_unique_cp_cells_counters(tree,mesh);

    if(cli.debug_mode) //IN DEBUG MODE
    {
        tuple<IAstar,bool,IAstar_stats> params;
        get<0>(params) = IAstar(mesh);
        get<1>(params) = cli.debug_mode;
        get<2>(params) = IAstar_stats();

        Timer time;
        time.start();
        tree.visit(IAstar_Generator::global_generation_Simplicial_wrapper,tree.get_root(),mesh,params);
        time.stop();
        time.print_elapsed_time("Build global IA-star ");

        stringstream stream2;
        stream2<<get_path_without_file_extension(cli.mesh_path)<<"_"<<ss.str()<<".iastar";
        get<0>(params).save_IAstar(stream2.str(),reindexer.get_original_vertices_ordering(),reindexer.get_original_tops_ordering());
        /// *** ///
    }

    return (EXIT_SUCCESS);
}

// NOTA: old approach, the visit function is implemented into the Mesh_Repairer class
template<class T,class C, class Top> int repair_mesh(T& tree, Mesh<C,Top>& mesh, cli_parameters &cli)
{
    if(cli.kind_of_operation == "fix")
    {
        cerr<<"[WARNING] The mesh repair is a limited and experimental feature.."<<endl;
        Mesh_Repairer<C,Top> repairer;
        repairer.repair_mesh(tree.get_root(),mesh);
        repairer.print_repair_output(get_path_without_file_extension(cli.mesh_path),mesh);
    }
    return (EXIT_SUCCESS);
}
