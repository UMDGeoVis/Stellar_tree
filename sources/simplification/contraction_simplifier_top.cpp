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

#include "contraction_simplifier_top.h"

void Contraction_Simplifier_Top::simplify(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli, ivect &v_orig, string base_file_name)
{
    cerr<<"==Homology preserving simplification - top-link=="<<endl;

    cerr<<"[NOTICE] Cache size: "<<cli.cache_size<<endl;
    LRU_Cache<int,leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VTop relations
    contraction_parameters params(&tree,mesh.get_top_cells_types());

    if(cli.debug_prints_mode)
        params.enable_debug_prints();
    if(cli.debug_mode)
        params.enable_debug_mode();

    Timer time;
    int simplification_round;
    int round = 1;

    if(params.is_enable_debug_mode())
    {
        cerr<<"[DEBUG] Computing Initial VTop cardinality"<<endl;
        vector<variance> vtop_stats;
        vtop_stats.assign(mesh.get_top_cells_types()+1,variance());
        time.start();
        tree.visit(topological_queries::extract_VTop_stats_Simplicial_wrapper,tree.get_root(),mesh,vtop_stats);
        time.stop();
        time.print_elapsed_time("[TIME] Computing VTop cardinality: ");
        for(int d=0; d<mesh.get_top_cells_types(); d++)
        {
            if(vtop_stats[d].min!=INT_MAX && mesh.get_top_cells_num(d) > 0)
            {
//                vtop_stats[d].avg /= mesh.get_top_cells_num(d);
                vtop_stats[d].avg /= vtop_stats[d].v_counter; //mesh.get_vertices_num();
                cerr<<"   "<<d+1<<"-simplices "<<vtop_stats[d].min<<" "<<vtop_stats[d].avg<<" "<<vtop_stats[d].max;
                cerr<<" (number of vertices incident in a "<<d+1<<"-simplex: "<<vtop_stats[d].v_counter<<")"<<endl;
            }
        }
//        vtop_stats[vtop_stats.size()-1].avg /= mesh.count_top_cells_num();
        vtop_stats[vtop_stats.size()-1].avg /= mesh.get_vertices_num();
        cerr<<"[STATS] INITIAL VTop degree: "<<vtop_stats[vtop_stats.size()-1].min<<" "<<vtop_stats[vtop_stats.size()-1].avg<<" "<<vtop_stats[vtop_stats.size()-1].max<<endl;
    }

    time.start();
    while(1)
    {
        simplification_round = params.get_contracted_edges_num();
        tree.visit_with_cache(Contraction_Simplifier_Top::simplify_leaf,tree.get_root(),mesh,cache,params);

        // PARTIAL SIMPLIFICATION STATS
        cerr<<"=== end-of-round "<<round<<") --> contracted edges: ";
        cerr<<params.get_contracted_edges_num()-simplification_round<<endl;
        round++;

        if(cli.debug_mode)
        {
            time.stop();
            time.print_elapsed_time("   [TIME] executing a simplification round: ");
            time.start();
            cerr << "   [RAM] peak for executing a simplification round: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }

        if(simplification_round == params.get_contracted_edges_num())
            break;

        cache.reset();
    }
    time.stop();
    if(!cli.debug_mode)
        time.print_elapsed_time("[TIME] Edge contraction simplification: ");
    else
        params.print_simplification_partial_timings();
    params.print_simplification_counters();

    cerr << "[RAM peak] for contracting a simplicial complex: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    /// finally we have to update/compress the mesh and the tree
    Contraction_Simplifier_Weak::update_mesh_and_tree(tree,mesh,params);

    /// we output the list of the checked and contracted edges during the procedure
    if(params.is_enable_debug_mode())
    {
        cerr<<"[DEBUG] Writing checked-edges order (for Skeleton-Blockers)"<<endl;
        params.write_checked_edges(base_file_name,v_orig);

        cerr<<"[DEBUG] Writing simplified mesh"<<endl;
        stringstream ss; ss<<"_kv_"<<cli.vertices_per_leaf<<"_top_link_simplification_";
        Writer::write_mesh(string_management::get_path_without_file_extension(cli.mesh_path),ss.str(),mesh);

        cerr<<"[DEBUG] Computing Final VTop cardinality"<<endl;
        vector<variance> vtop_stats;
        vtop_stats.assign(mesh.get_top_cells_types()+1,variance());
        time.start();
        tree.visit(topological_queries::extract_VTop_stats_Simplicial_wrapper,tree.get_root(),mesh,vtop_stats);
        time.stop();
        time.print_elapsed_time("[TIME] Computing VTop cardinality: ");
        for(int d=0; d<mesh.get_top_cells_types(); d++)
        {
            if(vtop_stats[d].min!=INT_MAX && mesh.get_top_cells_num(d) > 0)
            {
                vtop_stats[d].avg /= vtop_stats[d].v_counter; //mesh.get_vertices_num();
                cerr<<"   "<<d+1<<"-simplices "<<vtop_stats[d].min<<" "<<vtop_stats[d].avg<<" "<<vtop_stats[d].max;
                cerr<<" (number of vertices incident in a "<<d+1<<"-simplex: "<<vtop_stats[d].v_counter<<")"<<endl;
            }
        }
        vtop_stats[vtop_stats.size()-1].avg /= mesh.count_top_cells_num();
        cerr<<"[STATS] FINAL VTop degree: "<<vtop_stats[vtop_stats.size()-1].min<<" "<<vtop_stats[vtop_stats.size()-1].avg<<" "<<vtop_stats[vtop_stats.size()-1].max<<endl;
    }
}

void Contraction_Simplifier_Top::simplify_weighted(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli, ivect &v_orig, string base_file_name)
{
    cerr<<"==Weighted Homology preserving simplification - top-link=="<<endl;

    cerr<<"[NOTICED] Cache size: "<<cli.cache_size<<endl;
    LRU_Cache<int,leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VTop relations
    w_contraction_parameters params(&tree,mesh.get_top_cells_types()/*,cli.no_homology_preservation*/);
    params.init_v_communities_array(mesh.get_vertices_num());
    params.init_v_degrees_array(mesh.get_vertices_num());

    if(cli.debug_prints_mode)
        params.enable_debug_prints();
    if(cli.debug_mode)
        params.enable_debug_mode();

    /// (1) extract the vertices degree
//    cout<<"extract the vertices degree"<<endl;    
    Timer time;
    time.start();
    tree.visit(topological_queries::extract_v_degrees_Simplicial_wrapper,tree.get_root(),mesh,params.get_v_degrees_array());
    time.stop();
    time.print_elapsed_time("[TIME] computing vertices degrees: ");
//    cerr << "[RAM] peak for computing vertices degrees: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
    params.compute_v_degree_stats(mesh);

    /// (2) get (again) the original-2-resorted vertices array
    /// + read the edges weight file giving to the edge the new vertices indices
//    cout<<"get (again) the original-2-resorted vertices array"<<endl;
    ivect v_orig_to_new;
    v_orig_to_new.assign(v_orig.size(),0);
    for(unsigned v=0; v<v_orig.size(); v++)
        v_orig_to_new[v_orig[v] - 1] = v + 1;
    params.init_history_map(cli.weights_file,v_orig_to_new);
    v_orig_to_new.clear();

    cerr << "[RAM] peak PRE-simplification: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;


    int simplification_round;
    int round = 1;

    cerr<<"==Simplification=="<<endl;
    time.start();
    while(1)
    {
        simplification_round = params.get_contracted_edges_num();
        tree.visit_with_cache(Contraction_Simplifier_Top::simplify_weighted_leaf,tree.get_root(),mesh,cache,params);

        // PARTIAL SIMPLIFICATION STATS
        cerr<<"=== end-of-round "<<round<<") --> contracted edges: ";
        cerr<<params.get_contracted_edges_num()-simplification_round<<endl;
        round++;

        if(cli.debug_mode)
        {
            time.stop();
            time.print_elapsed_time("   [TIME] executing a simplification round: ");
            time.start();
            cerr << "   [RAM] peak for executing a simplification round: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }

        if(simplification_round == params.get_contracted_edges_num())
            break;

        cache.reset();
    }
    time.stop();
    if(!cli.debug_mode)
        time.print_elapsed_time("[TIME] Edge contraction simplification: ");
    else
        params.print_simplification_partial_timings();
    params.print_simplification_counters();

    cerr << "[RAM peak] for CONTRACTING a simplicial complex: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

//    int a; cin>>a;

    params.compute_v_community_stats(mesh);
    params.write_v_community(base_file_name,v_orig,mesh);

//    if(params.is_do_not_preserve_homology())
//    {
//        cerr<<"[TAPULLO] extra check only for non-homology preserving simplification"<<endl;
//        /// TAPULLO TIME ///
//        /// we force a check for deleted tops that are not actually catched during the simplification
//        for(int d=0; d<mesh.get_top_cells_types(); d++)
//        {
//            for(int i=1; i<mesh.get_top_cells_num(d); i++)
//            {
//                if(!mesh.is_top_cell_removed(d,i))
//                {
//                    Top_Simplex &top = mesh.get_top_cell(d,i);
//                    for(int v=0; v<top.get_vertices_num(); v++)
//                    {
//                        if(mesh.is_vertex_removed(abs(top.TV(v))))
//                        {
//                            cerr<<"   [REMOVED] "<<i<<"] "<<top<<"has vertex "<<abs(top.TV(v))<<" deleted"<<endl;
//                            mesh.remove_top_cell(d,i);
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//    }

    /// finally we have to update/compress the mesh and the tree
    Contraction_Simplifier_Weak::update_mesh_and_tree(tree,mesh,params);

    /// buggy ---- TO DO -> CHECK WHY!
//    params.init_v_degrees_array(mesh.get_vertices_num());
//    time.start();
//    tree.visit(topological_queries::extract_v_degrees_Simplicial_wrapper,tree.get_root(),mesh,params.get_v_degrees_array());
//    time.stop();
//    time.print_elapsed_time("[TIME] computing vertices degrees: ");
////    cerr << "[RAM] peak for computing vertices degrees: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//    params.compute_v_degree_stats(mesh);

    params.compute_summation_history_map_stats();

    /// we output the list of the checked and contracted edges during the procedure
//    if(params.is_enable_debug_mode())
    params.write_checked_edges(base_file_name,v_orig);
}

void Contraction_Simplifier_Top::simplify_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params)
{
    Timer time;

    if(params.is_enable_debug_mode())
    {
        time.start();
    }

    leaf_VT vtops;
    n.extract_local_VTop(mesh,vtops);

//    if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//        cerr << "[RAM peak] after extracting VTops: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////        int a; cin>>a; /// <====== WARNING!!!!
//    }

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.vt_gen += time.get_elapsed_time();

        /// get local VT statistics
        long vt_refs = vtops.size();
        for(auto it=vtops.begin(); it!=vtops.end(); ++it)
        {
            vt_refs += it->size();
        }
        if(params.local_vt_max_references < vt_refs)
            params.local_vt_max_references = vt_refs;

        time.start();
    }

    edge_queue edges;
    Contraction_Simplifier_Weak::extract_target_edges(n,mesh,edges);

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.get_edges += time.get_elapsed_time();

        /// get edge queue statistics
        if(params.edge_q_max_size < (int)edges.size())
            params.edge_q_max_size = edges.size();        
    }

//    if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//        cerr << "[RAM peak] for gathering initial target edges: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////        int a; cin>>a; /// <====== WARNING!!!!
//    }

//    bool stepbystep = false;

    while(!edges.empty())
    {
        ivect e = edges.front();
        edges.pop();

        /// for debug only
        if(params.is_enable_debug_mode())
        {
            if(Contraction_Simplifier_Weak::is_degenerate_edge(e,n,mesh))
            {
                int a; cin>>a;
            }
        }

        if(Contraction_Simplifier_Weak::skip_edge(e,mesh))
        {
            continue;
        }

        /// get VT0 -- VT1 -- ET and (if exists) the leaf block indexing the e[0] (by definition e[1] is inside b)
        Node_Stellar *outer_v_block;
        VT *vt1, *vt0;
        ET et;
        Contraction_Simplifier_Weak::get_edge_relations(e,et,vt0,vt1,outer_v_block,n,mesh,vtops,cache,params);

//        if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//        {
//            params.ram_peak = MemoryUsage().getValue_in_MB(false);
//            cerr << "[RAM peak] after get_edge_relations: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////            int a; cin>>a; /// <====== WARNING!!!!
//        }

//        print_container_of_containers_content("ET of e: ",et);

        if(params.is_enable_debug_prints())
        {
            cout<<n<<endl;
            print_container_content("current edge: ",e);
//            if(e[0]==164 && e[1]==584)
//            if(stepbystep)
//            {
//                cout<<"-- VT0 "<<e[0]<<endl;
//                Contraction_Simplifier_Weak::print_verbose_coboundary_relation(*vt0,mesh);
//                cout<<"-- VT1 "<<e[1]<<endl;
//                Contraction_Simplifier_Weak::print_verbose_coboundary_relation(*vt1,mesh);
//                cout<<"-- ET"<<endl;
//                Contraction_Simplifier_Weak::print_verbose_coboundary_relation(et,mesh);
//            }
        }

        bool is_contracted = false;

////        if(e[0]==1 && e[1]==604)
////        {
//            if(Contraction_Simplifier_Top::is_tops_link_condition_valid(e,et,*vt0,*vt1,mesh,params) !=
//                    Contraction_Simplifier_Top::is_tops_link_condition_valid_v2(e,et,*vt0,*vt1,mesh,params))
//            {
//                params.enable_debug_prints();

//                cout<<"result v1 method (current): "<<Contraction_Simplifier_Top::is_tops_link_condition_valid(e,et,*vt0,*vt1,mesh,params)<<endl;
//                cout<<"result v2 method (adherent to algorithm): "<<Contraction_Simplifier_Top::is_tops_link_condition_valid_v2(e,et,*vt0,*vt1,mesh,params)<<endl;
//                params.disable_debug_prints();
//                print_container_content("edge e: ",e);
////                print_container_of_containers_content("ET of e: ",et);
////                print_container_of_containers_content("VT of v0: ",*vt0);
////                print_container_of_containers_content("VT of v1: ",*vt1);
////                cout<<"-- VT0 "<<e[0]<<endl;
////                Contraction_Simplifier_Weak::print_verbose_coboundary_relation(*vt0,mesh);
////                cout<<"-- VT1 "<<e[1]<<endl;
////                Contraction_Simplifier_Weak::print_verbose_coboundary_relation(*vt1,mesh);
////                cout<<"-- ET"<<endl;
////                Contraction_Simplifier_Weak::print_verbose_coboundary_relation(et,mesh);
////                cout<<"result weak link condition: "<<Contraction_Simplifier_Weak::is_weak_link_condition_valid_top_down_version(e,et,*vt0,*vt1,mesh,params)<<endl;
//                int a; cin>>a;
//            }
////        }

//        if(Contraction_Simplifier_Top::is_tops_link_condition_valid(e,et,*vt0,*vt1,mesh,params))
        if(Contraction_Simplifier_Top::is_tops_link_condition_valid_v2(e,et,*vt0,*vt1,mesh,params))
        {
//            if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//            {
//                params.ram_peak = MemoryUsage().getValue_in_MB(false);
//                cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//                cerr << "[RAM peak] after checking top link condition (true): " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////                int a; cin>>a; /// <====== WARNING!!!!
//            }

            if(params.is_enable_debug_mode())
                time.start();

            Contraction_Simplifier_Weak::contract_edge(e,et,*vt0,*vt1,vtops,*outer_v_block,edges,n,mesh,cache,params);

            if(params.is_enable_debug_mode())
            {
                time.stop();
                params.contract_edge += time.get_elapsed_time();

                if(Contraction_Simplifier_Weak::exists_degenerate_tops(*vt0,mesh))
                {
                    int a; cin>>a;
                }
            }

            is_contracted = true;

//            if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//            {
//                params.ram_peak = MemoryUsage().getValue_in_MB(false);
//                cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//                cerr << "[RAM peak] after contracting edge: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////                int a; cin>>a; /// <====== WARNING!!!!
//            }
        }

        if(params.is_enable_debug_mode())
        {
            /// we extract the complete list of edge contraction for the skeleton blockers
            params.add_checked_edge(e,is_contracted); /// JUST FOR DEBUG

            /// get edge queue statistics
            if(params.edge_q_max_size < (int)edges.size())
                params.edge_q_max_size = edges.size();
        }

//        if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//        {
//            params.ram_peak = MemoryUsage().getValue_in_MB(false);
//            cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//            cerr << "[RAM peak] after checking top link condition (false): " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////            int a; cin>>a; /// <====== WARNING!!!!
//        }
    }

    if(params.is_enable_debug_mode())
    {
        time.start();
    }

    cache.insert(n.get_v_start(),vtops);

//    if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//        cerr << "[RAM peak] after inserting in cache: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////        int a; cin>>a; /// <====== WARNING!!!!
//    }

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.cache_insert += time.get_elapsed_time();

        /// get cache statistics
        Contraction_Simplifier_Weak::compute_cache_stats(cache,params);        

        time.start();
    }

    // we remove the top simplices from the current leaf block
    n.compact_top_cell_arrays(mesh);

//    if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//        cerr << "[RAM peak] after compacting current leaf block: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////        int a; cin>>a; /// <====== WARNING!!!!
//    }

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.clean_t_lists += time.get_elapsed_time();

        if(!Contraction_Simplifier_Weak::cache_is_synced(n,mesh,params.get_tree(),cache))
        {
            int a; cin>>a;
        }
    }


}

//void Contraction_Simplifier_Top::simplify_weighted_with_rejected_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int,leaf_VT> &cache, w_contraction_parameters &params)
//{
////    if(params.is_enable_debug_mode()) /// WARNING NOT PRINT MODE BUT DEBUG!!!!
////        cout<<n<<endl;

//    Timer time;

//    if(params.is_enable_debug_mode())
//    {
//        time.start();
//    }

//    leaf_VT vtops;
//    n.extract_local_VTop(mesh,vtops);

//    if(params.is_enable_debug_mode())
//    {
//        time.stop();
//        params.vt_gen += time.get_elapsed_time();

//        /// get local VT statistics
//        long vt_refs = vtops.size();
//        for(auto it=vtops.begin(); it!=vtops.end(); ++it)
//        {
//            vt_refs += it->size();
//        }
//        if(params.local_vt_max_references < vt_refs)
//            params.local_vt_max_references = vt_refs;

//        time.start();
//    }

//    edge_pqueue edges;
//    Contraction_Simplifier_Weak::extract_target_edges(n,mesh,edges,params);

////    edge_list recheck_edges;
//    rejected_edge_map r_edges;

////    int max_pqueue = edges.size();
////    int max_r_list = 0;
////    int num_change_state = 0;
////    int num_contraction = 0;

//    if(params.is_enable_debug_mode())
//    {
//        time.stop();
//        params.get_edges += time.get_elapsed_time();

//        /// get edge queue statistics
//        if(params.edge_q_max_size < (int)edges.size())
//            params.edge_q_max_size = edges.size();
//    }

////    cout<<"before while loop"<<endl;

//    /// disabilitato controllo per vedere se abbiamo edge che cambiano stato
//    /// perche' da studi preliminari e' pressoche' nulla la possibilita' che questo avvenga
//    /// ~~ 1 su 1000 cambia stato
//    /// computazionalmente, processare un edge una volta sola e' molto rilevante

//    while(!edges.empty())
//    {
//        edge_weight ew = edges.top();
//        edges.pop();

////        cout<<"before checked edge"<<endl;
//        int contracted = check_edge(ew,edges,vtops,n,mesh,cache,params);
////        cout<<"checked edge"<<endl;
//        if(contracted > 0) /// contraction executed
//        {
////            if(params.is_enable_debug_mode())
////                num_contraction++;

//            int target, deleted;
//            if(mesh.is_vertex_removed(ew.get_edge()[0]))
//            {
//                target = ew.get_edge()[1];
//                deleted = ew.get_edge()[0];
//            }
//            else
//            {
//                target = ew.get_edge()[0];
//                deleted = ew.get_edge()[1];
//            }

//            r_edges.erase(deleted);

//            rejected_edge_map::iterator it=r_edges.find(target);
//            if(it!=r_edges.end())
//            {
//                int_vect e_r; e_r.assign(2,0);
//                set<pair<int,short> > &s = it->second;
//                for(auto itv=s.begin(); itv!=s.end();)
//                {
//                    e_r[0] = target;
//                    e_r[1] = itv->first;
//                    sort(e_r.begin(),e_r.end());

//                    edge_weight ew_r(e_r,itv->second);

//                    if(check_edge(ew_r,edges,vtops,n,mesh,cache,params) > 0) /// contraction executed
//                    {
//                        itv = s.erase(itv);

//                        if(params.is_enable_debug_mode())
//                            params.increment_change_state_edges();
//                    }
//                    else
//                        ++itv;
//                }

//            }

////            for(auto it=recheck_edges.begin(); it!=recheck_edges.end(); )
////            {
////                if(check_edge(*it,edges,vtops,n,mesh,cache,params) > 0) /// contraction executed
////                {
////                    it = recheck_edges.erase(it);

////                    if(params.is_enable_debug_mode())
////                        num_change_state++;
////                }
////                else
////                    ++it;
////            }
//        }
//        else if(contracted < 0) /// link condition fail
//        {
////            recheck_edges.push_back(ew);
//            r_edges[ew.get_edge()[0]].insert(make_pair(ew.get_edge()[1],ew.get_weight()));
//            r_edges[ew.get_edge()[1]].insert(make_pair(ew.get_edge()[0],ew.get_weight()));
//        }

////        if(params.is_enable_debug_mode())
////        {
////            if(edges.size() > max_pqueue)
////                max_pqueue = edges.size();
////            if(r_edges.size() > max_r_list)
////                max_r_list = r_edges.size();
////        }
//    }

//    if(params.is_enable_debug_mode())
//    {
//        time.start();
//    }

//    cache.insert(n.get_v_start(),vtops);

//    if(params.is_enable_debug_mode())
//    {
//        time.stop();
//        params.cache_insert += time.get_elapsed_time();

//        /// get cache statistics
//        Contraction_Simplifier_Weak::compute_cache_stats(cache,params);

////        time.start();
//    }

//    /// non ha piu' senso compattare qui.. perche' se cambia qualcosa nell'intorno della foglia lo faccio direttamente quando contraggo
////    // we remove the top simplices from the current leaf block
////    n.compact_top_cell_arrays(mesh);

//    /// TO ENABLE ONLY FOR DEBUGGING THE PROCEDURE..
//    /// THE OPERATION IS COSTLY
////    if(params.is_enable_debug_mode())
////    {
//////        time.stop();
//////        params.clean_t_lists += time.get_elapsed_time();

////        if(!Contraction_Simplifier_Weak::cache_is_synced(n,mesh,params.get_tree(),cache))
////        {
////            int a; cin>>a;
////        }

//////        cout<<"max pqueue size: "<<max_pqueue<<endl;
//////        cout<<"max recheck list size: "<<max_r_list<<endl;
//////        cout<<"contraction number: "<<num_contraction<<endl;
//////        cout<<"edges that change state: "<<num_change_state<<endl;
////    }
//}

void Contraction_Simplifier_Top::simplify_weighted_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int,leaf_VT> &cache, w_contraction_parameters &params)
{
//    if(params.is_enable_debug_mode()) /// WARNING NOT PRINT MODE BUT DEBUG!!!!
//    {
//        cerr<<n<<endl;
//        n.print_node_stats();
//    }

    Timer time;

    if(params.is_enable_debug_mode())
    {
        time.start();
    }

    leaf_VT vtops;
    n.extract_local_VTop(mesh,vtops);

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.vt_gen += time.get_elapsed_time();

        /// get local VT statistics
        long vt_refs = vtops.size();
        for(auto it=vtops.begin(); it!=vtops.end(); ++it)
        {
            vt_refs += it->size();
        }
        if(params.local_vt_max_references < vt_refs)
            params.local_vt_max_references = vt_refs;

        time.start();
    }

    edge_pqueue edges;
    Contraction_Simplifier_Weak::extract_target_edges(n,mesh,edges,params);

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.get_edges += time.get_elapsed_time();

        /// get edge queue statistics
        if(params.edge_q_max_size < (int)edges.size())
            params.edge_q_max_size = edges.size();
    }

    /// disabilitato controllo per vedere se abbiamo edge che cambiano stato
    /// perche' da studi preliminari e' pressoche' nulla la possibilita' che questo avvenga
    /// ~~ 1 su 1000 cambia stato
    /// computazionalmente, processare un edge una volta sola e' molto rilevante
    while(!edges.empty())
    {
        edge_weight ew = edges.top();
        edges.pop();

        check_weighted_edge(ew,edges,vtops,n,mesh,cache,params);
    }

    if(params.is_enable_debug_mode())
    {
        time.start();
    }

    cache.insert(n.get_v_start(),vtops);

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.cache_insert += time.get_elapsed_time();

        /// get cache statistics
        Contraction_Simplifier_Weak::compute_cache_stats(cache,params);

//        params.print_simplification_partial_timings();
    }

    /// non ha piu' senso compattare qui.. perche' se cambia qualcosa nell'intorno della foglia lo faccio direttamente quando contraggo
//    // we remove the top simplices from the current leaf block
//    n.compact_top_cell_arrays(mesh);

    /// TO ENABLE ONLY FOR DEBUGGING THE PROCEDURE..
    /// THE OPERATION IS COSTLY
//    if(params.is_enable_debug_mode())
//    {
////        time.stop();
////        params.clean_t_lists += time.get_elapsed_time();

//        if(!Contraction_Simplifier_Weak::cache_is_synced(n,mesh,params.get_tree(),cache))
//        {
//            int a; cin>>a;
//        }

////        cout<<"max pqueue size: "<<max_pqueue<<endl;
////        cout<<"max recheck list size: "<<max_r_list<<endl;
////        cout<<"contraction number: "<<num_contraction<<endl;
////        cout<<"edges that change state: "<<num_change_state<<endl;
//    }
}

int Contraction_Simplifier_Top::check_weighted_edge(edge_weight &ew, edge_pqueue &edges, leaf_VT &vtops, Node_Stellar &n,
                                           Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, w_contraction_parameters &params)
{
    Timer time;

    ivect &e = ew.get_edge();

    /// for debug only
    if(params.is_enable_debug_mode() /*|| params.is_do_not_preserve_homology()*/)
    {
        if(Contraction_Simplifier_Weak::is_degenerate_edge(e,n,mesh))
        {
//            if(params.is_do_not_preserve_homology())
//            {
//                /// we have to remove from the history map the current edge
//                params.remove_edge(ew.get_edge());
//                return 0;
//            }

            print_container_content("degenerate edge: ",e);
            Node_Stellar *outer_v_block;
            VT *vt1, *vt0;
            ET et;
            Contraction_Simplifier_Weak::get_edge_relations(e,et,vt0,vt1,outer_v_block,n,mesh,vtops,cache,params);
            cout<<"-- ET"<<endl;
            Contraction_Simplifier_Weak::print_verbose_coboundary_relation(et,mesh);
            int a; cin>>a;
        }
    }

    if(Contraction_Simplifier_Weak::skip_edge(e,mesh))
    {
        /// we have to remove from the history map the current edge
        params.remove_edge(ew.get_edge());
        return 0;
    }

    /// get VT0 -- VT1 -- ET and (if exists) the leaf block indexing the outer extreme)
    Node_Stellar *outer_v_block;
    VT *vt1, *vt0;
    ET et;
    Contraction_Simplifier_Weak::get_edge_relations(e,et,vt0,vt1,outer_v_block,n,mesh,vtops,cache,params);

    if(params.is_enable_debug_prints())
    {
        params.enable_debug_prints();
        cout<<n<<endl;
        print_container_content("current edge: ",e);
        cout<<"weight: "<<ew.get_weight()<<endl;
    }

    int contracted = -1;

    /// nota: we contract always to the vertex with the initial highest degree
    params.set_contraction_order( false );


    ///////// NOTA ////////////
    /// we contract if we are not preserving the homology and the current edge has a weight below the average
    /// OR if the top-link condition is verified
    if(/*(params.is_do_not_preserve_homology() && params.is_weight_below_avg(e)) ||*/
//            Contraction_Simplifier_Top::is_tops_link_condition_valid(e,et,*vt0,*vt1,mesh,params))
             Contraction_Simplifier_Top::is_tops_link_condition_valid_v2(e,et,*vt0,*vt1,mesh,params))
    {
        if(params.is_enable_debug_mode())
            time.start();

        /// nota: we contract always to the vertex with the initial highest degree
        params.set_contraction_order( (params.get_v_degree(e[0]) < params.get_v_degree(e[1])) );

        if(params.is_enable_debug_prints())
            cout<<"CONTRACT edge: "<<e[0]<<" "<<e[1]<<" swapped contraction order? "<<params.swap_order()<<endl;

        Contraction_Simplifier_Weak::contract_weighted_edge(e,ew.get_weight(),et,*vt0,*vt1,vtops,*outer_v_block,edges,n,mesh,cache,params);

        contracted = 1;

        /// we keep track of the vertices communities
        if(params.swap_order())
        {
//            if(params.is_do_not_preserve_homology())
//            {
//                Contraction_Simplifier_Weak::remove_duplicated_tops(*vt1,mesh);
//                Contraction_Simplifier_Weak::remove_degenerate_tops(*vt1,mesh,false);
//            }
            if(params.is_enable_debug_mode())
            {
                time.stop();
                params.contract_edge += time.get_elapsed_time();

                if(Contraction_Simplifier_Weak::exists_duplicated_tops(*vt1,mesh))
                {
                    int a; cin>>a;
                }
                if(Contraction_Simplifier_Weak::exists_degenerate_tops(*vt1,mesh))
                {
                    int a; cin>>a;
                }
            }
            params.add_v_to_community(e[1],e[0]);
        }
        else
        {
//            if(params.is_do_not_preserve_homology())
//            {
//                Contraction_Simplifier_Weak::remove_duplicated_tops(*vt0,mesh);
//                Contraction_Simplifier_Weak::remove_degenerate_tops(*vt0,mesh,false);
//            }
            if(params.is_enable_debug_mode())
            {
                time.stop();
                params.contract_edge += time.get_elapsed_time();

                if(Contraction_Simplifier_Weak::exists_duplicated_tops(*vt0,mesh))
                {
                    int a; cin>>a;
                }
                if(Contraction_Simplifier_Weak::exists_degenerate_tops(*vt0,mesh))
                {
                    int a; cin>>a;
                }
            }
            params.add_v_to_community(e[0],e[1]);
        }
    }

    if(params.is_enable_debug_mode())
    {
        /// we extract the complete list of edge contraction for the skeleton blockers
        if(params.swap_order())
        {
            ivect eswapped = e;
            swap(eswapped[0],eswapped[1]);
            params.add_checked_edge(eswapped,(contracted==1));
        }
        else
            params.add_checked_edge(e,(contracted==1)); /// JUST FOR DEBUG

        /// get edge queue statistics
        if(params.edge_q_max_size < (int)edges.size())
            params.edge_q_max_size = edges.size();
    }

    return contracted;
}

bool Contraction_Simplifier_Top::is_tops_link_condition_valid(ivect &e, ET &et, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, contraction_parameters &params)
{
    Timer time;
    iset v_in_link0, v_in_link1, v_in_linke;

    if(params.is_enable_debug_mode())
        time.start();

    /// we have to extract just the vertices in the link of e and its extrema
    Contraction_Simplifier_Top::get_vertices_in_link(e[0],vt0,v_in_link0,mesh);
    Contraction_Simplifier_Top::get_vertices_in_link(e[1],vt1,v_in_link1,mesh);
    Contraction_Simplifier_Top::get_vertices_in_link(e,et,v_in_linke,mesh);

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.get_links += time.get_elapsed_time();
    }

    if(params.is_enable_debug_prints())
    {
        cout<<"v0 - num vertices in link: "<<v_in_link0.size()<<endl;
        cout<<"v1 - num vertices in link: "<<v_in_link1.size()<<endl;
        cout<<"e - num vertices in link: "<<v_in_linke.size()<<endl;
    }

    if(params.is_enable_debug_mode())
        time.start();

    /// (1) we check first the weak link condition on the vertices
    iset v_intersection = v_in_link0;
    intersect_sets(v_intersection,v_in_link1);
    if(v_intersection != v_in_linke)
    {
        if(params.is_enable_debug_mode())
        {
            time.stop();
            params.check_link += time.get_elapsed_time();

            params.failed_vertices_lc_num++;
        }
        if(params.is_enable_debug_prints())
            cout<<"[FALSE] weak link condition violated by vertices"<<endl;
        return false;
    }

//    if(params.is_enable_debug_prints())
//        print_container_content("=== v_intersection: ",v_intersection);

    set<ivect> shared_v0, shared_v1;
    Contraction_Simplifier_Top::get_shared_simplices(vt0,v_intersection,1,mesh,shared_v0); /// WE SKIP THE EDGES, AS WE WANT SUB-SIMPLICES WITH SIZE > 1
    Contraction_Simplifier_Top::get_shared_simplices(vt1,v_intersection,1,mesh,shared_v1); /// WE SKIP THE EDGES, AS WE WANT SUB-SIMPLICES WITH SIZE > 1

//    for(unsigned i=1; i<vt0.size(); i++) /// WE CAN SKIP THE EDGES, AS WE WANT SUB-SIMPLICES WITH SIZE > 1
//    {
//        for(auto it=vt0[i].begin(); it!=vt0[i].end(); ++it)
//        {
//            Top_Simplex &top = mesh.get_top_cell(i,*it);
//            int_vect shared;
//            top.get_shared_cell(v_intersection,shared);
//            shared_v0.insert(shared);
//        }
//    }
//    for(unsigned i=1; i<vt1.size(); i++) /// WE CAN SKIP THE EDGES, AS WE WANT SUB-SIMPLICES WITH SIZE > 1
//    {
//        for(auto it=vt1[i].begin(); it!=vt1[i].end(); ++it)
//        {
//            Top_Simplex &top = mesh.get_top_cell(i,*it);
//            int_vect shared;
//            top.get_shared_cell(v_intersection,shared);
//            shared_v1.insert(shared);
//        }
//    }

    for(auto it0=shared_v0.begin(); it0!=shared_v0.end(); ++it0)
    {
        for(auto it1=shared_v1.begin(); it1!=shared_v1.end(); ++it1)
        {
            ivect intersect = *it0;
            intersect_containers(intersect,*it1);
            if(intersect.size() > 0)
            {
//                if(params.is_enable_debug_prints())
//                    print_container_content("intersect: ",intersect);
                if(!Contraction_Simplifier_Weak::exists_simplex_in_coboundary(intersect.size()-1,intersect,et,mesh))
                {
                    if(params.is_enable_debug_prints())
                    {
                        cout<<"[FALSE] identified shared sub-simplex that is not in the ET"<<endl;
                    }
                    if(params.is_enable_debug_mode())
                    {
                        time.stop();
                        params.check_link += time.get_elapsed_time();

                        params.failed_top_lc_num++;
                    }
                    return false;
                }
//                else if(params.is_enable_debug_prints())
//                    cout<<"= contained in ET (intersect) = "<<endl;
            }
        }
    }

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.check_link += time.get_elapsed_time();

        params.ok_lc++;
    }

    /// if I arrive here means that I do not have found any "blocker"
    if(params.is_enable_debug_prints())
        cout<<"[TRUE] nothing identified.. contraction valid"<<endl;
    return true;
}

bool Contraction_Simplifier_Top::is_tops_link_condition_valid_v2(ivect &e, ET &et, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, contraction_parameters &params)
{
    Timer time;

    if(params.is_enable_debug_mode())
        time.start();

    VT T0 = vt0; difference_of_container_of_containers(T0,et);
    VT T1 = vt1; difference_of_container_of_containers(T1,et);    

//    vector< vector<ivect> > explicit_etop; explicit_etop.assign(et.size(),vector<ivect>());
//    Contraction_Simplifier_Weak::explicit_top_coboundary(et,explicit_etop,mesh);

    /// special-case: e is a top edge
    ///    we remove e from VT(T0) and VT(T1)
    bool is_top_edge = false;
    int id = -1;
    if(/*id != -1*/get_num_elements_in_container_of_containers(et)==0) //if ET is empty then e is a top edge
    {
        id = Contraction_Simplifier_Weak::get_top_edge_id(e,T0,mesh);
        if(id==-1)
        {
            cout << "e is NOT a top edge with id: "<<id<<endl;
            int a; cin>>a;
        }
//        et[0].push_back(id);
        ivect::iterator position = find(T0[0].begin(), T0[0].end(), id);
        if (position != T0[0].end())
            T0[0].erase(position);
        position = find(T1[0].begin(), T1[0].end(), id);
        if (position != T1[0].end())
            T1[0].erase(position);
    }
    /// Notice this removal is only on the local auxiliary data structures..
    ///    if the contraction is successfull we have add the e to ET,
    ///    in order to successfully remove it from the corresponding VT during the update phase    

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.get_links += time.get_elapsed_time();
    }

    if(params.is_enable_debug_mode())
        time.start();

    for (int d0=0;d0<T0.size();d0++)
    {
        for(auto t0 : T0[d0])
        {
            ivect t0_vect;
            mesh.get_top_cell(d0,t0).get_positive_sorted_vertices(t0_vect);

            for (int d1=0;d1<T1.size();d1++)
            {
                for (auto t1 : T1[d1])
                {
//                    Top_Simplex &top1 = mesh.get_top_cell(d1,t1);
                    ivect intersect = t0_vect;
                    ivect t1_vect;
                    mesh.get_top_cell(d1,t1).get_positive_sorted_vertices(t1_vect);

                    intersect_containers(intersect,t1_vect);
                    if(intersect.size() > 0)
                    {
        //                if(params.is_enable_debug_prints())
        //                    print_container_content("intersect: ",intersect);
                        if(params.is_enable_debug_mode())
                        {
                            time.stop();
                            params.check_link += time.get_elapsed_time();
                            time.start();
                        }

                        if(!Contraction_Simplifier_Weak::exists_simplex_in_coboundary(intersect.size()-1,intersect,et,mesh))
//                        if(!Contraction_Simplifier_Weak::exists_simplex_in_coboundary_v2(intersect.size()-1,intersect,explicit_etop))
                        {
                            if(params.is_enable_debug_prints())
                            {
                                print_container_content("t0_vect: ",t0_vect);
                                print_container_content("t1_vect: ",t1_vect);
                                print_container_content("intersect: ",intersect);
                                cout<<"[FALSE] identified shared sub-simplex that is not in the ET"<<endl;
                            }
                            if(params.is_enable_debug_mode())
                            {
                                time.stop();
                                params.check_simplex_in_cob += time.get_elapsed_time();

                                params.failed_top_lc_num++;
                            }
                            return false;
                        }
        //                else if(params.is_enable_debug_prints())
        //                    cout<<"= contained in ET (intersect) = "<<endl;
                        if(params.is_enable_debug_mode())
                        {
                            time.stop();
                            params.check_simplex_in_cob += time.get_elapsed_time();
                            time.start();
                        }
                    }
                }
            }
        }
    }

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.check_link += time.get_elapsed_time();

        params.ok_lc++;
    }

    /// if I arrive here means that I do not have found any "blocker"
    if(params.is_enable_debug_prints())
        cout<<"[TRUE] nothing identified.. contraction valid"<<endl;

    /// I commented this out as in contract_edge we check if the ET is empty --> then e is a top-edge
    /// that contradicts what I do here
//    /// If the contraction is successfull and e is a top edge,
//    ///    we have add e to its ET, in order to successfully remove it from the corresponding VT during the update phase
//    if(is_top_edge)
//        et[0].push_back(id);

    return true;
}

void Contraction_Simplifier_Top::get_vertices_in_link(int v_id, VT &vt, iset &v_in_link, Simplicial_Mesh &mesh)
{
    for(unsigned d=0; d<vt.size();d++)
    {
        for(unsigned i=0; i<vt[d].size();i++)
        {
            Top_Simplex &top = mesh.get_top_cell(d,vt[d][i]);
            ivect &v_arr = top.get_vertices_array();

            ///(1) we extract the vertices of t that are in the link of v_id
            ivect v_in_v_link(top.get_vertices_num()-1);
            copy_if(v_arr.begin(),v_arr.end(),v_in_v_link.begin(), [v_id](int i){return (i!=v_id);} );
            sort(v_in_v_link.begin(),v_in_v_link.end());

            v_in_link.insert(v_in_v_link.begin(),v_in_v_link.end());
        }
    }
}

void Contraction_Simplifier_Top::get_vertices_in_link(const ivect &e, ET &et, iset &v_in_link, Simplicial_Mesh &mesh)
{
    for(unsigned d=0; d<et.size();d++)
    {
        for(unsigned i=0; i<et[d].size();i++)
        {
            Top_Simplex &top = mesh.get_top_cell(d,et[d][i]);
            ivect &v_arr = top.get_vertices_array();

            ///(1) we extract the vertices of t that are in the link of e
            int v_id1 = e[0], v_id2 = e[1];
            ivect v_in_e_link(top.get_vertices_num()-2);
            copy_if(v_arr.begin(),v_arr.end(),v_in_e_link.begin(), [v_id1,v_id2](int i){return (i!=v_id1 && i!=v_id2);} );
            sort(v_in_e_link.begin(),v_in_e_link.end());

            v_in_link.insert(v_in_e_link.begin(),v_in_e_link.end());
        }
    }
}

void Contraction_Simplifier_Top::get_shared_simplices(VT &vt, iset &v_intersection, int d, Simplicial_Mesh &mesh, set<ivect> &shared_s)
{
    for(unsigned i=d; i<vt.size(); i++) /// WE CAN SKIP THE EDGES, AS WE WANT SUB-SIMPLICES WITH SIZE > 1
    {
        for(auto it=vt[i].begin(); it!=vt[i].end(); ++it)
        {
            Top_Simplex &top = mesh.get_top_cell(i,*it);
//            cout<<top<<endl;
            ivect shared;
            top.get_shared_cell(v_intersection,shared);
            if(shared.size() > 0)
                shared_s.insert(shared);
        }
    }
}


