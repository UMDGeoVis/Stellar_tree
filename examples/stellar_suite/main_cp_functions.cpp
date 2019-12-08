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

#include "main_cp_functions.h"

void load_mesh(CP_Mesh &mesh, cli_parameters &cli)
{
    if(get_file_extension(cli.mesh_path) == "off")
    {
        if (!Reader::read_cell_mesh(mesh, cli.mesh_path))
        {
            cerr << "[load_mesh] Error Loading mesh file. Execution Stopped." << endl;
        }
    }
    else if(get_file_extension(cli.mesh_path) == "hex")
    {
        if (!Reader::read_hexa_mesh(mesh, cli.mesh_path))
        {
            cerr << "[load_mesh] Error Loading mesh file. Execution Stopped." << endl;
        }
    }
}

int extract_VTop_relations(Stellar_Tree &tree, CP_Mesh &mesh)
{
    Timer time;
    pair<bool,int> p; p.first = false; p.second = 0;

    time.start();
    tree.visit(topological_queries::extract_VTop_CP_wrapper,tree.get_root(),mesh,p);
    time.stop();
    time.print_elapsed_time("Local VT computation ");

    /// *** PARALLEL *** ///
    time.start();
    tree.parallel_visit(topological_queries::extract_VTop_CP_wrapper,mesh,p);
    time.stop();
    time.print_elapsed_time("Local VT computation --- parallel --- ");
    /// ****** ///

    /// this enables the statistic computation
    p.first = true;
    tree.visit(topological_queries::extract_VTop_CP_wrapper,tree.get_root(),mesh,p);
    cerr<<"max_list_summation: "<<p.second<<endl;

    return (EXIT_SUCCESS);
}

int extract_VV_relations(Stellar_Tree &tree, CP_Mesh &mesh)
{
    Timer time;
    pair<bool,int> p; p.first = false; p.second = 0;

    time.start();
    tree.visit(topological_queries::extract_VV_CP_wrapper,tree.get_root(),mesh,p);
    time.stop();
    time.print_elapsed_time("Local VT computation ");

    /// *** PARALLEL *** ///
    time.start();
    tree.parallel_visit(topological_queries::extract_VV_CP_wrapper,mesh,p);
    time.stop();
    time.print_elapsed_time("Local VT computation --- parallel --- ");
    /// ****** ///

    /// this enables the statistic computation
    p.first = true;
    tree.visit(topological_queries::extract_VV_CP_wrapper,tree.get_root(),mesh,p);
    cerr<<"max_list_summation: "<<p.second<<endl;

    return (EXIT_SUCCESS);
}

int extract_links_relations(Stellar_Tree &tree, CP_Mesh &mesh)
{
    Timer time;
    tuple<bool,int,int,long> p (false,0,0,0);

    time.start();
    tree.visit(topological_queries::extract_links_CP_wrapper,tree.get_root(),mesh,p);
    time.stop();
    time.print_elapsed_time("Local LINKS computation ");

    /// *** PARALLEL *** ///
    time.start();
    tree.parallel_visit(topological_queries::extract_links_CP_wrapper,mesh,p);
    time.stop();
    time.print_elapsed_time("Local LINKS computation --- parallel --- ");
    /// ****** ///

    /// this enables the statistic computation
    get<0>(p) = true;
    tree.visit(topological_queries::extract_links_CP_wrapper,tree.get_root(),mesh,p);
    cerr<<"max_vertices: "<<get<1>(p)<<endl;
    cerr<<"max_edges: "<<get<2>(p)<<endl;
    cerr<<"max_references: "<<get<3>(p)<<endl;

    return (EXIT_SUCCESS);
}

int extract_p_faces(Stellar_Tree& tree, CP_Mesh &mesh)
{
    Timer time;
    tuple<int,bool,unsigned> p(0,true,0);

    ///we exclude the vertices and the top simplices
    ///then we visit all the implicitly encoded simplexes
    cerr<<"[OPERATION] Extract p-faces"<<endl;

    double tot_timings = 0;

//    cerr<<"0-simplices num: "<<mesh.get_vertices_num()<<endl;
//    tot_simplices += mesh.get_vertices_num();
//    cerr<<mesh.get_implicitly_encoded_cells_num()<<"-simplices num: "<<mesh.get_top_cells_num(mesh.get_top_cells_types()-1)<<endl;
//    tot_simplices += mesh.get_top_cells_num(mesh.get_top_cells_types()-1);
    stringstream ss;
    for(int d=1; d<mesh.get_implicitly_encoded_cells_num(); d++) // from triangles to d-2 cells
    {
        get<0>(p) = d;
        time.start();
        tree.visit(topological_queries::extract_p_faces_CP_wrapper,tree.get_root(),mesh,p);
        time.stop();
        ss << "[TIME] Extracting "<<d<<"-faces ";
        time.print_elapsed_time(ss.str());
        tot_timings += time.get_elapsed_time();
        cerr<<"[STAT] max "<<d<<"-faces auxiliary structure size: "<<get<2>(p)<<endl;
        get<2>(p) = 0;        
        ss.str("");
    }

    cerr<<"[TIME][SEQUENTIAL] Total extractiong timings "<<tot_timings<<endl;

    /// *** PARALLEL *** ///
    tot_timings = 0;

    for(int d=1; d<mesh.get_implicitly_encoded_cells_num(); d++) // from triangles to d-2 cells
    {
        tuple<int,bool,unsigned> p(d,true,0);
        time.start();
        tree.parallel_visit(topological_queries::extract_p_faces_CP_wrapper,mesh,p);
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

int extract_p_faces_v2(Stellar_Tree& tree, CP_Mesh &mesh)
{
    Timer time;
    tuple<int,bool,unsigned> p(0,true,0);

    ///we exclude the vertices and the top simplices
    ///then we visit all the implicitly encoded simplexes
    cerr<<"[OPERATION] Extract p-faces - v2"<<endl;

    double tot_timings = 0;

//    cerr<<"0-simplices num: "<<mesh.get_vertices_num()<<endl;
//    tot_simplices += mesh.get_vertices_num();
//    cerr<<mesh.get_implicitly_encoded_cells_num()<<"-simplices num: "<<mesh.get_top_cells_num(mesh.get_top_cells_types()-1)<<endl;
//    tot_simplices += mesh.get_top_cells_num(mesh.get_top_cells_types()-1);
    stringstream ss;
    for(int d=1; d<mesh.get_implicitly_encoded_cells_num(); d++) // from triangles to d-2 cells
    {
        get<0>(p) = d;
        time.start();
        tree.visit(topological_queries::extract_p_faces_CP_wrapper_v2,tree.get_root(),mesh,p);
        time.stop();
        ss << "[TIME] Extracting "<<d<<"-faces ";
        time.print_elapsed_time(ss.str());
        tot_timings += time.get_elapsed_time();
        cerr<<"[STAT] max "<<d<<"-faces auxiliary structure size: "<<get<2>(p)<<endl;
        get<2>(p) = 0;
        ss.str("");
    }

    cerr<<"[TIME][SEQUENTIAL] Total extractiong timings "<<tot_timings<<endl;

    /// *** PARALLEL *** ///
    tot_timings = 0;

    for(int d=1; d<mesh.get_implicitly_encoded_cells_num(); d++) // from triangles to d-2 cells
    {
        tuple<int,bool,unsigned> p(d,true,0);
        time.start();
        tree.parallel_visit(topological_queries::extract_p_faces_CP_wrapper_v2,mesh,p);
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

int extract_adjacencies(Stellar_Tree &tree, CP_Mesh &mesh, Reindexer &reindexer, cli_parameters &cli)
{
    Timer time;

    if(cli.kind_of_operation == "adjG")
    {
        tuple<IAstar,bool,IAstar_stats> params;
        get<0>(params) = IAstar(mesh);
        get<1>(params) = cli.debug_mode;
        get<2>(params) = IAstar_stats();

        time.start();
        tree.visit(IAstar_Generator::global_generation_CP_wrapper,tree.get_root(),mesh,params);
        time.stop();
        time.print_elapsed_time("Build global IA-star ");
        get<2>(params).print_global_stats(get<0>(params));

        if(cli.debug_mode)
        {
            stringstream stream;
            stream<<get_path_without_file_extension(cli.mesh_path)<<".iastar";
            get<0>(params).save_IAstar(stream.str(),reindexer.get_original_vertices_ordering(),reindexer.get_original_tops_ordering());
            reindexer.reset_original_ordering_variables();
        }

        if(cli.debug_mode)
        {
            get<2>(params).print_partial_timings();

            /// query time!!
            /// (1) extract VT for all vertices
            time.start();
            get<0>(params).extract_all_VTop(mesh/*,variables.debug_mode*/);
            time.stop();
            time.print_elapsed_time("Extract VT relations on IA-star ");

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
        tree.visit(IAstar_Generator::local_generation_CP_wrapper,tree.get_root(),mesh,params);
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
        tree.parallel_visit(IAstar_Generator::global_generation_CP_wrapper,mesh,params);
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
        tree.parallel_visit(IAstar_Generator::local_generation_CP_wrapper,mesh,params);
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

int validate_connectivity(Stellar_Tree& tree, CP_Mesh& mesh, cli_parameters &cli)
{
    if(cli.kind_of_operation == "validate")
    {
        pair<Connectedness_Validator_stats,bool> params;
        params.first = Connectedness_Validator_stats();
        params.second = cli.debug_mode;

        Timer t;

        /// check 0-connectedness
        params.first.cc.assign(mesh.get_vertices_num(),0);

        if(params.second)
            params.first.reset_time_variables();///for debug

        t.start();
        tree.visit(Connectedness_Validator::validate_0_connectedness_CP_wrapper,tree.get_root(),mesh,params);
        t.stop();
        t.print_elapsed_time("[validator] check_connectness: ");

        if(params.second)
        {
            params.first.print_time_variables();
            params.first.reset_time_variables();
        }

        params.first.print_validation_result();
        params.first.print_skeleton_stats(params.second);
        params.first.reset();

        params.first.cc.assign(mesh.get_vertices_num(),0);
        t.start();
        tree.visit(Connectedness_Validator::validate_0_connectednessVVs_CP_wrapper,tree.get_root(),mesh,params);
        t.stop();
        t.print_elapsed_time("[validator] check_connectness (with VVs): ");

        if(params.second)
        {
            params.first.print_time_variables();
            params.first.reset_time_variables();
        }

        params.first.print_validation_result();
        params.first.reset();

        params.first.cc.assign(mesh.get_top_cells_num(mesh.get_top_cells_types()-1),0); //init an array for the highest dimensional top simplexes

        t.start();
        tree.visit(Connectedness_Validator::validate_d_connectedness_CP_wrapper,tree.get_root(),mesh,params);
        t.stop();
        t.print_elapsed_time("[validator] check_d-1_connecteness: ");

        if(params.second)
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

int validate_link_conditions(Stellar_Tree& tree, CP_Mesh& mesh, cli_parameters &cli)
{
    if(cli.kind_of_operation == "validate")
    {
        edge_link_map faulty_links;
        Timer t;
        t.start();
        tree.visit(Links_Validator::validate_links_CP_wrapper,tree.get_root(),mesh,faulty_links);
        t.stop();
        t.print_elapsed_time("[validator] check link conditions: ");
        Links_Validator::print_validation_result(faulty_links,cli.debug_mode);

        faulty_links.clear();

        t.start();
        tree.parallel_visit(Links_Validator::validate_links_CP_wrapper,mesh,faulty_links);
        t.stop();
        t.print_elapsed_time("[validator] check link conditions --- parallel ---: ");
        Links_Validator::print_validation_result(faulty_links,cli.debug_mode);

    }
    return (EXIT_SUCCESS);
}
