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

#include "main_simplicial_functions.h"

void load_mesh(Simplicial_Mesh &mesh, cli_parameters &cli)
{
    if(get_file_extension(cli.mesh_path) == "tri")
    {
        if (!Reader::read_triangle_mesh(mesh, cli.mesh_path, cli.verbose_encoding))
        {
            cerr << "[load_mesh] Error Loading mesh file. Execution Stopped." << endl;
        }
    }
    else if(get_file_extension(cli.mesh_path) == "ts")
    {
        if (!Reader::read_tetra_mesh(mesh, cli.mesh_path, cli.verbose_encoding))
        {
            cerr << "[load_mesh] Error Loading mesh file. Execution Stopped." << endl;
        }
    }
    else if(get_file_extension(cli.mesh_path) == "gmv")
    {
        if (!Reader::read_gmv_mesh(mesh, cli.mesh_path))
        {
            cerr << "[load_mesh] Error Loading mesh file. Execution Stopped." << endl;
        }
    }
    else if(get_file_extension(cli.mesh_path) == "sierp")
    {
        if (!Reader::read_sierp_mesh(mesh, cli.mesh_path))
        {
            cerr << "[load_mesh] Error Loading mesh file. Execution Stopped." << endl;
        }
    }
    else if(get_file_extension(cli.mesh_path) == "off")
    {
        if (!Reader::read_simplicial_mesh(mesh, cli.mesh_path, cli.verbose_encoding,cli.force3D))
        {
            cerr << "[load_mesh] Error Loading mesh file. Execution Stopped." << endl;
        }
    }
}

int extract_VTop_relations(Stellar_Tree &tree, Simplicial_Mesh &mesh)
{
    Timer time;
    pair<bool,int> p; p.first = false; p.second = 0;

    time.start();
    tree.visit(topological_queries::extract_VTop_Simplicial_wrapper,tree.get_root(),mesh,p);
    time.stop();
    time.print_elapsed_time("Local VT computation ");

//    /// *** PARALLEL *** ///
//    time.start();
//    tree.parallel_visit(topological_queries::extract_VTop_Simplicial_wrapper,mesh,p);
//    time.stop();
//    time.print_elapsed_time("Local VT computation --- parallel --- ");
//    /// ****** ///

    /// this enables the statistic computation
    p.first = true;
    tree.visit(topological_queries::extract_VTop_Simplicial_wrapper,tree.get_root(),mesh,p);
    cerr<<"max_list_summation: "<<p.second<<endl;

    return (EXIT_SUCCESS);
}

int extract_VV_relations(Stellar_Tree &tree, Simplicial_Mesh &mesh)
{
    Timer time;
    pair<bool,int> p; p.first = false; p.second = 0;

    time.start();
    tree.visit(topological_queries::extract_VV_Simplicial_wrapper,tree.get_root(),mesh,p);
    time.stop();
    time.print_elapsed_time("Local VV computation ");

//    /// *** PARALLEL *** ///
//    time.start();
//    tree.parallel_visit(topological_queries::extract_VV_Simplicial_wrapper,mesh,p);
//    time.stop();
//    time.print_elapsed_time("Local VV computation --- parallel --- ");
//    /// ****** ///

    /// this enables the statistic computation
    p.first = true;
    tree.visit(topological_queries::extract_VV_Simplicial_wrapper,tree.get_root(),mesh,p);
    cerr<<"max_list_summation: "<<p.second<<endl;

    return (EXIT_SUCCESS);
}

int extract_links_relations(Stellar_Tree &tree, Simplicial_Mesh &mesh)
{
    Timer time;
    tuple<bool,int,int,long> p (false,0,0,0);

    time.start();
    tree.visit(topological_queries::extract_links_Simplicial_wrapper,tree.get_root(),mesh,p);
    time.stop();
    time.print_elapsed_time("Local LINKS computation ");

    /// *** PARALLEL *** ///
    time.start();
    tree.parallel_visit(topological_queries::extract_links_Simplicial_wrapper,mesh,p);
    time.stop();
    time.print_elapsed_time("Local LINKS computation --- parallel --- ");
    /// ****** ///

    /// this enables the statistic computation
    get<0>(p) = true;
    tree.visit(topological_queries::extract_links_Simplicial_wrapper,tree.get_root(),mesh,p);
    cerr<<"max_vertices: "<<get<1>(p)<<endl;
    cerr<<"max_edges: "<<get<2>(p)<<endl;
    cerr<<"max_references: "<<get<3>(p)<<endl;

    return (EXIT_SUCCESS);
}

int extract_p_faces(Stellar_Tree& tree, Simplicial_Mesh &mesh)
{
    Timer time;


    ///we exclude the vertices and the top simplices
    ///then we visit all the implicitly encoded simplexes
    cerr<<"[OPERATION] Extract p-faces"<<endl;

    double tot_timings = 0;

    stringstream ss;
    for(int d=1; d<mesh.get_implicitly_encoded_cells_num(); d++) // from triangles to d-2 cells
    {
        tuple<int,bool,unsigned> p(d,true,0);
        time.start();
        tree.visit(topological_queries::extract_p_faces_Simplicial_wrapper,tree.get_root(),mesh,p);
        time.stop();
        ss << "[TIME] Extracting "<<d<<"-faces ";
        time.print_elapsed_time(ss.str());
        tot_timings += time.get_elapsed_time();
        cerr<<"[STAT] max "<<d<<"-faces auxiliary structure size: "<<get<2>(p)<<endl;
        ss.str("");
    }

    cerr<<"[TIME][SEQUENTIAL] Total extraction timings "<<tot_timings<<endl;

    /// *** PARALLEL *** ///
    tot_timings = 0;

    for(int d=1; d<mesh.get_implicitly_encoded_cells_num(); d++) // from triangles to d-2 cells
    {
        tuple<int,bool,unsigned> p(d,true,0);
        time.start();
        tree.parallel_visit(topological_queries::extract_p_faces_Simplicial_wrapper,mesh,p);
        time.stop();
        ss << "[TIME] Extracting "<<d<<"-faces ";
        time.print_elapsed_time(ss.str());
        tot_timings += time.get_elapsed_time();
        cerr<<"[STAT] max "<<d<<"-faces auxiliary structure size: "<<get<2>(p)<<endl;
        ss.str("");
    }

    cerr<<"[TIME][OPENMP] Total extraction timings "<<tot_timings<<endl;

    return (EXIT_SUCCESS);
}

int extract_p_faces_v2(Stellar_Tree& tree, Simplicial_Mesh &mesh)
{
    Timer time;


    ///we exclude the vertices and the top simplices
    ///then we visit all the implicitly encoded simplexes
    cerr<<"[OPERATION] Extract p-faces - v2"<<endl;

    double tot_timings = 0;

    stringstream ss;
    for(int d=1; d<mesh.get_implicitly_encoded_cells_num(); d++) // from triangles to d-2 cells
    {
        tuple<int,bool,unsigned> p(d,true,0);
        time.start();
        tree.visit(topological_queries::extract_p_faces_Simplicial_wrapper_v2,tree.get_root(),mesh,p);
        time.stop();
        ss << "[TIME] Extracting "<<d<<"-faces ";
        time.print_elapsed_time(ss.str());
        tot_timings += time.get_elapsed_time();
        cerr<<"[STAT] max "<<d<<"-faces auxiliary structure size: "<<get<2>(p)<<endl;
        ss.str("");
    }

    cerr<<"[TIME][SEQUENTIAL] Total extraction timings "<<tot_timings<<endl;

    /// *** PARALLEL *** ///
    tot_timings = 0;

    for(int d=1; d<mesh.get_implicitly_encoded_cells_num(); d++) // from triangles to d-2 cells
    {
        tuple<int,bool,unsigned> p(d,true,0);
        time.start();
        tree.parallel_visit(topological_queries::extract_p_faces_Simplicial_wrapper_v2,mesh,p);
        time.stop();
        ss << "[TIME] Extracting "<<d<<"-faces ";
        time.print_elapsed_time(ss.str());
        tot_timings += time.get_elapsed_time();
        cerr<<"[STAT] max "<<d<<"-faces auxiliary structure size: "<<get<2>(p)<<endl;
        ss.str("");
    }

    cerr<<"[TIME][OPENMP] Total extraction timings "<<tot_timings<<endl;

    return (EXIT_SUCCESS);
}

int extract_adjacencies(Stellar_Tree &tree, Simplicial_Mesh &mesh, Reindexer &reindexer, cli_parameters &cli)
{
    Timer time;

    if(cli.kind_of_operation == "adjG")
    {
        tuple<IAstar,bool,IAstar_stats> params;
        get<0>(params) = IAstar(mesh);
        get<1>(params) = cli.debug_mode;
        get<2>(params) = IAstar_stats();

        time.start();
        tree.visit(IAstar_Generator::global_generation_Simplicial_wrapper,tree.get_root(),mesh,params);
        time.stop();
        time.print_elapsed_time("Build global IA-star ");
        get<2>(params).print_global_stats(get<0>(params));

        if(cli.debug_mode)
        {
            stringstream stream;
            stream<<get_path_without_file_extension(cli.mesh_path)<<".iastar";
//            stream<<get_path_without_file_extension(cli.mesh_path)<<".ia";
            get<0>(params).save_IAstar(stream.str(),reindexer.get_original_vertices_ordering(),reindexer.get_original_tops_ordering());
//            get<0>(params).save_IAstar_full(stream.str(),mesh);
            //        reindexer.reset_original_ordering_variables();
        }

        reindexer.reset_original_ordering_variables();

        if(cli.debug_mode)
        {
            get<2>(params).print_partial_timings();
            //        gia.print_non_manifold_adjacencies();

            /// query time!!
            /// (1) extract VT for all vertices
            time.start();
            get<0>(params).extract_all_VTop(mesh/*,variables.debug_mode*/);
            time.stop();
            time.print_elapsed_time("Extract VT relations on IA-star ");

            //        gia.print_non_manifold_adjacencies();

            /// (2) validate the model
            Connectedness_Validator validator;
            Connectedness_Validator_stats stats;
            validator.validate_connectedness_IAstar(get<0>(params),mesh,stats);
        }
    }
    else if(cli.kind_of_operation == "adjL")
    {
        pair<bool,IAstar_stats> params;
        params.first = cli.debug_mode;
        params.second = IAstar_stats(mesh.get_top_cells_types(),cli.debug_mode);

        time.start();
        tree.visit(IAstar_Generator::local_generation_Simplicial_wrapper,tree.get_root(),mesh,params);
        time.stop();
        time.print_elapsed_time("Build local IA-star ");

        if(cli.debug_mode)
        {
            params.second.print_local_stats();
            params.second.print_partial_timings();
            params.second.reset_stats(mesh.get_top_cells_types());
        }
    }
    else if(cli.kind_of_operation == "adjGP")
    {
        tuple<IAstar,bool,IAstar_stats> params;
        /// *** PARALLEL *** ///
        get<0>(params) = IAstar(mesh);
        get<1>(params) = cli.debug_mode;
        get<2>(params) = IAstar_stats();

        time.start();
        tree.parallel_visit(IAstar_Generator::global_generation_Simplicial_wrapper,mesh,params);
        time.stop();
        time.print_elapsed_time("Build global IA-star --- parallel --- ");
        get<2>(params).print_global_stats(get<0>(params));
        get<0>(params).compute_storage(mesh);
        /// *** ///

        if(cli.debug_mode)
            get<2>(params).print_partial_timings();
    }
    else if(cli.kind_of_operation == "adjLP")
    {
        pair<bool,IAstar_stats> params;
        params.first = cli.debug_mode;
        params.second = IAstar_stats(mesh.get_top_cells_types(),cli.debug_mode);
        /// *** PARALLEL *** ///
        time.start();
        tree.parallel_visit(IAstar_Generator::local_generation_Simplicial_wrapper,mesh,params);
        time.stop();
        time.print_elapsed_time("Build local IA-star --- parallel --- ");

        if(cli.debug_mode)
        {
            params.second.print_local_stats();
            params.second.print_partial_timings();
            params.second.reset_stats(mesh.get_top_cells_types());
        }
        /// *** ///
    }

    return (EXIT_SUCCESS);
}

int extract_global_skeleton(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli)
{
    if(cli.kind_of_operation == "skel")
    {
        skeleton_graph n_graph = skeleton_graph(mesh.get_vertices_num());

        Timer time;
        time.start();

        /// (1) extract skeleton ///
        tree.visit(skeleton_generator::extract_global_directed_1skeleton_Simplicial_wrapper,tree.get_root(),mesh,n_graph);

        time.stop();
        time.print_elapsed_time("1-skeleton extraction time: ");

        n_graph.save_skeleton_graph(get_path_without_file_extension(cli.mesh_path));
//        n_graph.write_gv_file(get_file_name(variables.mesh_path));
    }

    return (EXIT_SUCCESS);
}

int validate_connectivity(Stellar_Tree& tree, Simplicial_Mesh& mesh, cli_parameters &cli)
{
    if(cli.kind_of_operation == "validate")
    {
        pair<Connectedness_Validator_stats,bool> params;
        params.first = Connectedness_Validator_stats();
        params.second = cli.debug_mode;

        Timer t;

        /// check 0-connectedness
        params.first.cc.assign(mesh.get_vertices_num(),0);

        if(cli.debug_mode)
            params.first.reset_time_variables();///for debug

        t.start();
        tree.visit(Connectedness_Validator::validate_0_connectedness_Simplicial_wrapper,tree.get_root(),mesh,params);
        t.stop();
        t.print_elapsed_time("[validator] check_connectness: ");

        if(cli.debug_mode)
        {
            params.first.print_time_variables();
            params.first.reset_time_variables();
        }

        params.first.print_validation_result();
        params.first.print_skeleton_stats(params.second);
        params.first.reset();

        params.first.cc.assign(mesh.get_vertices_num(),0);
        t.start();
        tree.visit(Connectedness_Validator::validate_0_connectednessVVs_Simplicial_wrapper,tree.get_root(),mesh,params);
        t.stop();
        t.print_elapsed_time("[validator] check_connectness (with VVs): ");

        if(cli.debug_mode)
        {
            params.first.print_time_variables();
            params.first.reset_time_variables();
        }

        params.first.print_validation_result();
        params.first.reset();

        params.first.cc.assign(mesh.get_top_cells_num(mesh.get_top_cells_types()-1),0); //init an array for the highest dimensional top simplexes

        t.start();
        tree.visit(Connectedness_Validator::validate_d_connectedness_Simplicial_wrapper,tree.get_root(),mesh,params);
        t.stop();
        t.print_elapsed_time("[validator] check_d-1_connecteness: ");

        if(cli.debug_mode)
        {
            params.first.print_time_variables();
            params.first.reset_time_variables();
        }

        params.first.print_validation_result();
        params.first.print_pseudo_manifold_result();
        params.first.reset();
    }
    return (EXIT_SUCCESS);
}

int validate_link_conditions(Stellar_Tree& tree, Simplicial_Mesh& mesh, cli_parameters &cli)
{
    if(cli.kind_of_operation == "validate")
    {
        edge_link_map faulty_links;
        Timer t;
        t.start();
        tree.visit(Links_Validator::validate_links_Simplicial_wrapper,tree.get_root(),mesh,faulty_links);
        t.stop();
        t.print_elapsed_time("[validator] check link conditions: ");
        Links_Validator::print_validation_result(faulty_links,cli.debug_mode);

        faulty_links.clear();

        t.start();
        tree.parallel_visit(Links_Validator::validate_links_Simplicial_wrapper,mesh,faulty_links);
        t.stop();
        t.print_elapsed_time("[validator] check link conditions --- parallel ---: ");
        Links_Validator::print_validation_result(faulty_links,cli.debug_mode);
    }
    return (EXIT_SUCCESS);
}

int simplify(Stellar_Tree& tree, Reindexer &reindexer, Simplicial_Mesh& mesh, cli_parameters &cli, Statistics &stats)
{
    Timer time;

    if((cli.kind_of_operation == "contractW" || cli.kind_of_operation == "contractT") &&
            cli.exec_weighted_simplification()) ///--> AGGIUNTO PER IL LAVORO ACM.. PER IL LAVORO CGF NON SERVONO QUESTE STATISTICHE
    {
        cerr<<"[STAT] initial mesh "<<endl;
        cerr<<"------------"<<endl;
        mesh.print_mesh_stats(cerr);
        cerr<<"------------"<<endl;

        pair<int,int> tmp_p = make_pair(-1,-1);
        vector<pair<int,int> > v_stats;
        v_stats.assign(mesh.get_vertices_num(),tmp_p);
        tree.visit(topological_queries::extract_VTop_VE_sizes_Simplicial_wrapper,tree.get_root(),mesh,v_stats);

        /// prints the 1-skeleton of the complex
        skeleton_graph n_graph = skeleton_graph(mesh.get_vertices_num());
        time.start();
        tree.visit(skeleton_generator::extract_global_directed_1skeleton_Simplicial_wrapper,tree.get_root(),mesh,n_graph);
        time.stop();
        time.print_elapsed_time("[TIME] computing directed 1-skeleton: ");
//        cerr << "[RAM] peak for computing directed 1-skeleton: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//        n_graph.write_gv_file(get_file_name(cli.mesh_path),reindexer.get_original_vertices_ordering(),mesh);
        n_graph.write_cesium_file(get_path_without_file_extension(cli.mesh_path),reindexer.get_original_vertices_ordering(),v_stats,mesh);
        n_graph.compute_skeleton_stats();
        n_graph.reset(); /// clear the skeleton graph
        v_stats.clear();

        n_graph = skeleton_graph(mesh.get_vertices_num());
        /// (1) extract skeleton ///
        time.start();
        tree.visit(skeleton_generator::extract_global_undirected_1skeleton_Simplicial_wrapper,tree.get_root(),mesh,n_graph);
        time.stop();
        time.print_elapsed_time("[TIME] computing undirected 1-skeleton: ");
        n_graph.compute_skeleton_stats();
        n_graph.reset(); /// clear the skeleton graph

//        /// prints an adjacency graph for the complete complex
////        cerr<<"[OUTPUT] computing the adjacency graph of the INITIAL mesh"<<endl;
//        Adjacency_Graph gg(mesh);
//        time.start();
//        tree.visit(Adjacency_Graph_Extractor::extract_lite_adjacency_graph_Simplicial_wrapper,tree.get_root(),mesh,gg);
//        time.stop();
//        time.print_elapsed_time("[TIME] computing adjacency graph: ");
////        cerr << "[RAM] peak for computing adjacency graph: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////        gg.write_graph_file(get_file_name(cli.mesh_path),mesh);
////        gg.write_fixed_graph_file(get_file_name(cli.mesh_path),mesh);
//        gg.write_cesium_graph_file(get_file_name(cli.mesh_path),mesh);

//        /// compute graphs statistics
//        cout<<"[STATS] Initial Graph Statistics (d-1 / d-2 degrees)"<<endl;
//        gg.compute_adj_graph_stats();
//        gg.reset(); /// clear the adjacency graph

//        cerr<<"[STATS] Initial VTop degrees stats"<<endl;
        vector<variance> vtop_stats;
        vtop_stats.assign(mesh.get_top_cells_types()+1,variance());
        time.start();
        tree.visit(topological_queries::extract_VTop_stats_Simplicial_wrapper,tree.get_root(),mesh,vtop_stats);
        for(int d=0; d<mesh.get_top_cells_types(); d++)
        {
            if(vtop_stats[d].min!=INT_MAX && mesh.get_top_cells_num(d) > 0)
            {
                vtop_stats[d].avg /= mesh.get_top_cells_num(d);
                cerr<<"   "<<d+1<<"-simplices "<<vtop_stats[d].min<<" "<<vtop_stats[d].avg<<" "<<vtop_stats[d].max<<endl;
            }
        }
        vtop_stats[vtop_stats.size()-1].avg /= mesh.count_top_cells_num();        
        time.stop();
        time.print_elapsed_time("[TIME] computing VTop cardinality: ");
        cout<<"[STATS] INITIAL Facet-degree: "<<vtop_stats[vtop_stats.size()-1].min<<" "<<vtop_stats[vtop_stats.size()-1].avg<<" "<<vtop_stats[vtop_stats.size()-1].max<<endl;
//        cerr << "[RAM] peak for computing VTop cardinality: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

//        pair<Connectedness_Validator_stats,bool> params;
//        params.first = Connectedness_Validator_stats();
//        params.second = cli.debug_mode;

//        params.first.cc.assign(mesh.get_vertices_num(),0);
//        time.start();
//        tree.visit(Connectedness_Validator::validate_0_connectednessVVs_Simplicial_wrapper,tree.get_root(),mesh,params);
//        time.stop();
//        time.print_elapsed_time("[VALIDATOR] check_connectness (with VVs): ");

//        if(cli.debug_mode)
//        {
//            params.first.print_time_variables();
//            params.first.reset_time_variables();
//        }

//        params.first.print_validation_result();
//        params.first.reset();
    }
    /// *** ///

    /// === output file name === ///
    /// if kv=INT_MAX we do not want to have a hierarchical decomposition
    string kv = (cli.vertices_per_leaf == 0) ? "inf" : to_string(cli.vertices_per_leaf);
    string strategy;
    if(cli.exec_weighted_simplification()) /// weighted edges + vertices degrees
        strategy = (string_management::get_file_extension(cli.weights_file) == "friend_weights") ? "friendship_weight" : "distance_weight";
    else
        strategy = "no_weight";
    stringstream ss;
    ss << cli.kind_of_operation << "_" << strategy << "_kv_" << kv;

    stringstream ss2;
    ss2 << get_path_without_file_extension(cli.mesh_path) << "_" << ss.str();
    /// === === ///

    /// here we discretize between weighted and unweighted simplification process
    /// and between top-based and weak link conditions checks
    if(cli.kind_of_operation == "contractW")
    {
        Contraction_Simplifier_Weak simplifier;
        if(cli.len) /// triangle mesh + edge length metric
            simplifier.simplify_length(tree,mesh,cli);
        else if(!cli.exec_weighted_simplification()) ///not weighted
            simplifier.simplify(tree,mesh,cli);
    }
    else if(cli.kind_of_operation == "contractT")
    {
        Contraction_Simplifier_Top simplifier;
        if(cli.exec_weighted_simplification()) /// weighted edges + vertices degrees
            simplifier.simplify_weighted(tree,mesh,cli,reindexer.get_original_vertices_ordering(),ss2.str());
        else if(cli.len) /// triangle mesh + edge length metric
            simplifier.simplify_length(tree,mesh,cli);
        else ///not weighted
            simplifier.simplify(tree,mesh,cli,reindexer.get_original_vertices_ordering(),ss2.str());
    }

    if(cli.kind_of_operation == "contractW" || cli.kind_of_operation == "contractT")
    {
        /// NOTA ENABLE FOR DEBUG
//        mesh.print_mesh(cout);

//        time.start();
        cerr<<"[TREE] remove unnecessary splitting"<<endl;
        tree.visit_and_unify(tree.get_root(),mesh);
        cerr<<"[TREE] reindexing the index"<<endl;
        reindexer.reorganize_index_and_mesh(tree,mesh,false);
//        time.stop();
//        time.print_elapsed_time("[TIME] final tree updates: ");
//        cerr << "[RAM peak] for compacting and exploit again the spatial coherence after the simplification: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

//        pair<Connectedness_Validator_stats,bool> params;
//        params.first.cc.assign(mesh.get_vertices_num(),0);
//        time.start();
//        tree.visit(Connectedness_Validator::validate_0_connectednessVVs_Simplicial_wrapper,tree.get_root(),mesh,params);
//        time.stop();
//        time.print_elapsed_time("[VALIDATOR] check_connectness (with VVs): ");

//        if(cli.debug_mode)
//        {
//            params.first.print_time_variables();
//            params.first.reset_time_variables();
//        }

//        params.first.print_validation_result();
//        params.first.reset();

        if(cli.get_index_stats)
            stats.get_index_statistics(tree,mesh,false); /// for DEBUG

//        /// write the saved mesh in OFF format
//        Writer::write_mesh(get_path_without_file_extension(cli.mesh_path),ss.str(),mesh);
    }

    if((cli.kind_of_operation == "contractW" || cli.kind_of_operation == "contractT") &&
            cli.exec_weighted_simplification()) ///--> AGGIUNTO PER IL LAVORO ACM.. PER IL LAVORO CGF NON SERVONO QUESTE STATISTICHE
    {
        pair<int,int> tmp_p = make_pair(-1,-1);
        vector<pair<int,int> > v_stats;
        v_stats.assign(mesh.get_vertices_num(),tmp_p);
        tree.visit(topological_queries::extract_VTop_VE_sizes_Simplicial_wrapper,tree.get_root(),mesh,v_stats);

//        cerr<<"[OUTPUT] compute the 1-skeleton of the SIMPLIFIED mesh"<<endl;
        skeleton_graph n_graph_s = skeleton_graph(mesh.get_vertices_num());
        /// (1) extract skeleton ///
        time.start();
        tree.visit(skeleton_generator::extract_global_directed_1skeleton_Simplicial_wrapper,tree.get_root(),mesh,n_graph_s);
        time.stop();
        time.print_elapsed_time("[TIME] computing directed 1-skeleton: ");
//        cerr << "[RAM] peak for computing directed 1-skeleton: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//        n_graph_s.write_gv_file(ss2.str(),reindexer.get_original_vertices_ordering(),mesh);
        n_graph_s.write_cesium_file(ss2.str(),reindexer.get_original_vertices_ordering(),v_stats,mesh);
        n_graph_s.compute_skeleton_stats();

        n_graph_s.reset(); /// clear the skeleton graph
        v_stats.clear();

        n_graph_s = skeleton_graph(mesh.get_vertices_num());
        /// (1) extract skeleton ///
        time.start();
        tree.visit(skeleton_generator::extract_global_undirected_1skeleton_Simplicial_wrapper,tree.get_root(),mesh,n_graph_s);
        time.stop();
        time.print_elapsed_time("[TIME] computing undirected 1-skeleton: ");
        n_graph_s.compute_skeleton_stats();
        n_graph_s.reset(); /// clear the skeleton graph

//        /// compute graphs statistics
//        cout<<"[STATS] Final Graph Statistics (d-1 / d-2 degrees)"<<endl;
//        Adjacency_Graph gg_s = Adjacency_Graph(mesh);
//        tree.visit(Adjacency_Graph_Extractor::extract_lite_adjacency_graph_Simplicial_wrapper,tree.get_root(),mesh,gg_s);
//        gg_s.compute_adj_graph_stats();

//        cerr<<"[STATS] Final VTop degrees stats"<<endl;
        vector</*tuple<int,double,int>*/variance> vtop_stats;
        vtop_stats.assign(mesh.get_top_cells_types()+1,variance());
        time.start();
        tree.visit(topological_queries::extract_VTop_stats_Simplicial_wrapper,tree.get_root(),mesh,vtop_stats);
        for(int d=0; d<mesh.get_top_cells_types(); d++)
        {
            if(vtop_stats[d].min!=INT_MAX && mesh.get_top_cells_num(d) > 0)
            {
                vtop_stats[d].avg /= mesh.get_top_cells_num(d);
//                cerr<<"   "<<d+1<<"-simplices "<<vtop_stats[d].min<<" "<<vtop_stats[d].avg<<" "<<vtop_stats[d].max<<endl;
            }
        }
        vtop_stats[vtop_stats.size()-1].avg /= mesh.count_top_cells_num();
        time.stop();
        time.print_elapsed_time("[TIME] computing VTop cardinality: ");
//        cerr << "[RAM] peak for computing VTop cardinality: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        cout<<"[STATS] FINAL Facet-degree: "<<vtop_stats[vtop_stats.size()-1].min<<" "<<vtop_stats[vtop_stats.size()-1].avg<<" "<<vtop_stats[vtop_stats.size()-1].max<<endl;

        cerr << "[RAM] peak POST-simplification: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

        /// we compute the blockers only on demand
        if(cli.extract_blockers != DEFAULT)
        {
            /// WARNING -> moved the blockers computation as it is the most expensive operation
            cerr<<"[OUTPUT] Blocker extraction on the SIMPLIFIED mesh"<<endl;
            Blockers_Extractor be;
            blocker_set blockers = be.extract_blockers(tree,mesh,(cli.extract_blockers=="blockersC"),reindexer.get_original_vertices_ordering(),ss2.str());
            cerr << "[RAM] peak for computing blockers: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
            //        int a; cin>>a;

//            cerr<<"[OUTPUT] computing the adjacency graph of the SIMPLIFIED mesh"<<endl;
//            gg_s.reset(); /// clear the adjacency graph
//            /// prints an adjacency graph for the complete complex
//            gg_s = Adjacency_Graph(mesh,blockers);
//            time.start();
//            tree.visit(Adjacency_Graph_Extractor::extract_lite_adjacency_graph_Simplicial_wrapper,tree.get_root(),mesh,gg_s);
//            time.stop();
//            time.print_elapsed_time("[TIME] computing adjacency graph: ");
//            //        cerr << "[RAM] peak for computing adjacency graph: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//            //        gg_s.write_graph_file(ss2.str(),mesh);
//            //        gg_s.write_blocker_graph_file(ss2.str());
//            //        gg_s.write_fixed_graph_file(ss2.str(),mesh);
//            gg_s.write_cesium_graph_file(ss2.str(),mesh);
//            //        gg_s.write_fixed_blockers_file(ss2.str(),mesh);
//            //        gg_s.write_cesium_blockers_file(ss2.str(),mesh);

//            cerr << "[RAM] peak POST-simplification: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }
        /// DISABILITATO RAMO ELSE COMPUTAZIONE BLOCKERS
//        else
//        {
////            /// simply write the cesium-file without blockers
////            gg_s.write_cesium_graph_file(ss2.str(),mesh);
////            Blockers_Extractor be;
////            be.extract_aux_structures_statistics(tree,mesh);
//        }
    }

    return (EXIT_SUCCESS);
}
