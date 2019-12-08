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

#include "contraction_simplifier_weak.h"
#include "stellar_decomposition/reindexer.h"
#include "utilities/container_utilities.h"

void Contraction_Simplifier_Weak::simplify(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli)
{
    cerr<<"==Homology preserving simplification - weak-link condition=="<<endl;

    cerr<<"[NOTICED] Cache size: "<<cli.cache_size<<endl;
    LRU_Cache<int,leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VTop relations
    contraction_parameters params(&tree,mesh.get_top_cells_types());

    if(cli.debug_prints_mode)
        params.enable_debug_prints();
    if(cli.debug_mode)
        params.enable_debug_mode();

    Timer time;
    int simplification_round;
    int round = 1;

    time.start();
    while(1)
    {
        simplification_round = params.get_contracted_edges_num();
        tree.visit_with_cache(Contraction_Simplifier_Weak::simplify_leaf,tree.get_root(),mesh,cache,params);

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
}


void Contraction_Simplifier_Weak::simplify_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params)
{
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

        /// get VT0 -- VT1 -- ET and (if exists) the leaf block indexing the outer extreme)
        Node_Stellar *outer_v_block;
        VT *vt1, *vt0;
        ET et;
        Contraction_Simplifier_Weak::get_edge_relations(e,et,vt0,vt1,outer_v_block,n,mesh,vtops,cache,params);

        if(params.is_enable_debug_prints() /*&& params.ram_peak != MemoryUsage().getValue_in_MB(false)*/)
        {
            cout<<n<<endl;
            print_container_content("current edge: ",e);
//            params.ram_peak = MemoryUsage().getValue_in_MB(false);
//            cerr << "[RAM peak] for gathering the topological relations of an edge: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }

        if(Contraction_Simplifier_Weak::is_weak_link_condition_valid_top_down_version(e,et,*vt0,*vt1,mesh,params))
//        if(Contraction_Simplifier::is_weak_link_condition_valid_bottom_up_version(e,et,*vt0,*vt1,mesh,params)) /// MUCH SLOWER --> DISABLED
        {
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

//            if(params.is_enable_debug_prints() && params.ram_peak != MemoryUsage().getValue_in_MB(false))
//            {
//                cout<<n<<endl;
//                print_container_content("current edge: ",e);
//                params.ram_peak = MemoryUsage().getValue_in_MB(false);
//                cerr << "[RAM peak] after contracting an edge: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//            }
        }

//        if(params.is_enable_debug_prints() && params.ram_peak != MemoryUsage().getValue_in_MB(false))
//        {
//            cout<<n<<endl;
//            print_container_content("current edge: ",e);
//            params.ram_peak = MemoryUsage().getValue_in_MB(false);
//            cerr << "[RAM peak] after checking/contracting an edge: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//        }

        /// get edge queue statistics
        if(params.is_enable_debug_mode())
        {
            if(params.edge_q_max_size < (int)edges.size())
                params.edge_q_max_size = edges.size();
        }
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

        time.start();
    }

    // we remove the top simplices from the current leaf block
    n.compact_top_cell_arrays(mesh);

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.clean_t_lists += time.get_elapsed_time();
    }
}





void Contraction_Simplifier_Weak::extract_target_edges(Node_Stellar &n, Simplicial_Mesh &mesh, edge_queue &edges)
{
    set<ivect> e_set;
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                Top_Simplex& t = mesh.get_top_cell(d,*t_id);
                ivect e; e.assign(2,0);

                for(int i=0; i<t.get_vertices_num(); i++)
                {
                    for(int j=i+1; j<t.get_vertices_num(); j++)
                    {
                        t.TE(e,i,j);
                        if(n.indexes_vertex(e[1])) /// we process an edge only if it has all the extrema already processed
                        {
                            auto r = e_set.insert(e);
                            if(r.second)
                                edges.push(e);
                        }
                    }
                }
            }
        }
    }
}

void Contraction_Simplifier_Weak::extract_target_edges(Node_Stellar &n, Simplicial_Mesh &mesh, edge_pqueue &edges, w_contraction_parameters &params)
{
    set<ivect> e_set;
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                Top_Simplex& t = mesh.get_top_cell(d,*t_id);
                ivect e; e.assign(2,0);

                for(int i=0; i<t.get_vertices_num(); i++)
                {
                    for(int j=i+1; j<t.get_vertices_num(); j++)
                    {
                        t.TE(e,i,j);
                        if(n.indexes_vertex(e[1])) /// we process an edge only if it has all the extrema already processed
                        {
                            auto r = e_set.insert(e);
                            if(r.second)
                            {
                                edge_history &h = params.get_edge_history(e);
                                edge_weight ew = edge_weight(e,h.initial_weight);
                                edges.push(ew);
                            }
                        }
                    }
                }
            }
        }
    }
}

bool Contraction_Simplifier_Weak::is_weak_link_condition_valid_top_down_version(ivect &e, ET &et, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, contraction_parameters &params)
{
    Timer time;
    s_link link0, link1, linke;

    if(params.is_enable_debug_mode())
        time.start();

    Contraction_Simplifier_Weak::get_link_lite(e[0],vt0,link0,mesh);
    Contraction_Simplifier_Weak::get_link_lite(e[1],vt1,link1,mesh);
    Contraction_Simplifier_Weak::get_link_lite(e,et,linke,mesh);

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.get_links += time.get_elapsed_time();
    }

    if(params.is_enable_debug_mode())
        time.start();

    /// (1) we check first the weak link condition on the vertices
    iset v_intersection = link0.get_v_in_link();
    intersect_sets(v_intersection,link1.get_v_in_link());
    if(v_intersection != linke.get_v_in_link())
    {
        if(params.is_enable_debug_prints() /*&& params.ram_peak != MemoryUsage().getValue_in_MB(false)*/)
        {
//            print_container_content("current edge: ",e);
            cout<<"[FALSE] weak link condition violated by vertices"<<endl;
//            params.ram_peak = MemoryUsage().getValue_in_MB(false);
//            cerr << "[RAM peak] for checking weak link condition: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }

        if(params.is_enable_debug_mode())
        {
            time.stop();
            params.check_link += time.get_elapsed_time();

            params.failed_vertices_lc_num++;
        }
        return false;
    }

    /// (2) then, we check the weak link condition on the d-1 faces (as we have already extracted all of them
    int d = linke.get_s_in_link_size()-1;
    set<ivect > s_intersection = link0.get_s_in_link(d);
    intersect_sets(s_intersection,link1.get_s_in_link(d));
    if(s_intersection != linke.get_s_in_link(d))
    {
        if(params.is_enable_debug_prints() /*&& params.ram_peak != MemoryUsage().getValue_in_MB(false)*/)
        {
            print_container_content("current edge: ",e);
            cout<<"[FALSE] weak link condition violated by "<<d<<"-simplices"<<endl;
//            params.ram_peak = MemoryUsage().getValue_in_MB(false);
//            cerr << "[RAM peak] for checking weak link condition: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }

        if(params.is_enable_debug_mode())
        {
            time.stop();
            params.check_link += time.get_elapsed_time();

            params.failed_d_simpl_lc_num++;
        }
        return false;
    }

    /// (3) finally we iteratively check the weak link condition on the d-2 simplices ... up to 1-simplices
    for(d=linke.get_s_in_link_size()-2; d>=0; d--)
    {
        /// ad ogni giro occorre estrarre dai simplessi di maggior dimensione i sub-simplessi di dimensione d
        /// visitiamo i d+1 simplessi nel link per estrarre i d-simplessi
        Contraction_Simplifier_Weak::extract_sub_d_simplices(d+1,link0);
        Contraction_Simplifier_Weak::extract_sub_d_simplices(d+1,link1);
        Contraction_Simplifier_Weak::extract_sub_d_simplices(d+1,linke);

        /// gli d+1 simplessi non sono più necessari
        /// quindi cancelliamo il loro link
        link0.clear_s_link(d+1);
        link1.clear_s_link(d+1);
        linke.clear_s_link(d+1);

        s_intersection = link0.get_s_in_link(d);
        intersect_sets(s_intersection,link1.get_s_in_link(d));
        if(s_intersection != linke.get_s_in_link(d))
        {
            if(params.is_enable_debug_prints() /*&& params.ram_peak != MemoryUsage().getValue_in_MB(false)*/)
            {
                print_container_content("current edge: ",e);
                cout<<"[FALSE] weak link condition violated by "<<d<<"-simplices"<<endl;
//                params.ram_peak = MemoryUsage().getValue_in_MB(false);
//                cerr << "[RAM peak] for checking weak link condition: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
            }

            if(params.is_enable_debug_mode())
            {
                time.stop();
                params.check_link += time.get_elapsed_time();

                params.failed_sub_simpl_lc_num[d]++;
            }
            return false;
        }
    }

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.check_link += time.get_elapsed_time();

        params.ok_lc++;
    }

    if(params.is_enable_debug_prints() /*&& params.ram_peak != MemoryUsage().getValue_in_MB(false)*/)
    {
        print_container_content("current edge: ",e);
        cout<<"[TRUE] weak link condition verified"<<endl;
//        print_container_of_containers_content("VT(v0): ",vt0);
//        print_container_of_containers_content("VT(v1): ",vt1);
//        print_container_of_containers_content("ET(e): ",et);
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[RAM peak] for checking weak link condition: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
    }

    return true;
}

//bool Contraction_Simplifier_Weak::is_weak_link_condition_valid_bottom_up_version(ivect &e, ET &et, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, contraction_parameters &params)
//{
//    Timer time;
//    s_link link0, link1, linke;

//    if(params.is_enable_debug_mode())
//        time.start();

//    Contraction_Simplifier_Weak::get_link_lite(e[0],vt0,link0,mesh);
//    Contraction_Simplifier_Weak::get_link_lite(e[1],vt1,link1,mesh);
//    Contraction_Simplifier_Weak::get_link_lite(e,et,linke,mesh);

//    if(params.is_enable_debug_mode())
//    {
//        time.stop();
//        params.get_links += time.get_elapsed_time();
//    }

//    if(params.is_enable_debug_mode())
//        time.start();

//    /// (1) we check first the weak link condition on the vertices
//    iset v_intersection = link0.get_v_in_link();
//    intersect_sets(v_intersection,link1.get_v_in_link());
//    if(v_intersection != linke.get_v_in_link())
//    {
//        if(params.is_enable_debug_mode())
//        {
//            time.stop();
//            params.check_link += time.get_elapsed_time();

//            params.failed_vertices_lc_num++;
//        }
//        return false;
//    }

//    /// (2) then, we check the weak link condition on the d-1 faces (as we have already extracted all of them
//    int d = linke.get_s_in_link_size()-1;
//    set<ivect> s_intersection = link0.get_s_in_link(d);
//    intersect_sets(s_intersection,link1.get_s_in_link(d));
//    if(s_intersection != linke.get_s_in_link(d))
//    {
//        if(params.is_enable_debug_mode())
//        {
//            time.stop();
//            params.check_link += time.get_elapsed_time();

//            params.failed_d_simpl_lc_num++;
//        }
//        return false;
//    }

//    /// i d-1 simplessi non sono più necessari
//    /// quindi cancelliamo il loro link
//    link0.clear_s_link(d);
//    link1.clear_s_link(d);
//    linke.clear_s_link(d);

//    /// (3) finally we iteratively check the weak link condition on the 1-simplices ... up to d-2-simplices
//    for(d=0; d<linke.get_s_in_link_size()-2; d++)
//    {
//        /// ad ogni livello devo estrarre il d-simplessi
//        /// dato che devo lavorare a livello di top
//        /// ad ogni iterazione, nelle VT, parto dai top di dimensione d e salgo
//        Contraction_Simplifier_Weak::extract_d_simplices(d,vt0,link0,mesh);
//        Contraction_Simplifier_Weak::extract_d_simplices(d,vt1,link1,mesh);
//        Contraction_Simplifier_Weak::extract_d_simplices(d,et,linke,mesh);

//        s_intersection = link0.get_s_in_link(d);
//        intersect_sets(s_intersection,link1.get_s_in_link(d));
//        if(s_intersection != linke.get_s_in_link(d))
//        {
//            if(params.is_enable_debug_mode())
//            {
//                time.stop();
//                params.check_link += time.get_elapsed_time();

//                params.failed_sub_simpl_lc_num[d]++;
//            }
//            return false;
//        }

//        /// i d simplessi non sono più necessari
//        /// quindi cancelliamo il loro link
//        link0.clear_s_link(d);
//        link1.clear_s_link(d);
//        linke.clear_s_link(d);
//    }

//    if(params.is_enable_debug_mode())
//    {
//        time.stop();
//        params.check_link += time.get_elapsed_time();

//        params.ok_lc++;
//    }

//    return true;
//}


void Contraction_Simplifier_Weak::contract_edge(ivect &e, ET &et, VT &vt0, VT &vt1, leaf_VT &vtops, Node_Stellar &outer_v_block, edge_queue &edges,
                                           Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params)
{
    difference_of_container_of_containers(vt0,et); // vt0 now contains the difference VT0 - ET
    difference_of_container_of_containers(vt1,et); // vt1 now contains the difference VT1 - ET

//    if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[==INSIDE==] contract_edge"<<endl;
//        cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//        cerr << "[RAM peak] after difference_of_c_of_c: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////        int a; cin>>a; /// <====== WARNING!!!!
//    }

    // we have to check if edge is a top edge
    // this happens if #[ETop(e)]=0 and e \in VTop(v0)
    int e_top_id = -1;
    if(get_num_elements_in_container_of_containers(et) == 0)
    {
        Contraction_Simplifier_Weak::check_and_remove_top_edge(e,e_top_id,vt0,vt1,mesh,params);
    }

//    if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[==INSIDE==] contract_edge"<<endl;
//        cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//        cerr << "[RAM peak] after remove top-edge: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////        int a; cin>>a; /// <====== WARNING!!!!
//    }

    /// prior checking the d-1 faces we update
    /// (1) the corresponding vt0 relation (by adding the tops in vt1-et
    /// (2) then the outer_v_block (if e is a cross edge)
    /// (3) and, finally, we add the new edges crossing b to the edge queue
    Contraction_Simplifier_Weak::update(e,vt0,vt1,n,outer_v_block,edges,mesh,params);

//    if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[==INSIDE==] contract_edge"<<endl;
//        cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//        cerr << "[RAM peak] after update: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////        int a; cin>>a; /// <====== WARNING!!!!
//    }

    // we check if the boundary d-1 faces have an adjacent in the corresponding vt, if not we add the face to the top d-1 simplices list
    // of the leaf block indexing it and to the corresponding VTop (if local or cached)
    Contraction_Simplifier_Weak::check_for_new_tops(e[1],vt0,et,n,mesh,vtops,cache,params.get_tree(),params);

//    if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[==INSIDE==] contract_edge"<<endl;
//        cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//        cerr << "[RAM peak] after check_for_new_tops: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////        int a; cin>>a; /// <====== WARNING!!!!
//    }

    // we remove v2 and the top simplices in et
    Contraction_Simplifier_Weak::remove_from_mesh(e[1],e_top_id,et,mesh,params);
    // finally we clear the VTop(v2)
    vt1.clear();
    et.clear();

//    if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[==INSIDE==] contract_edge"<<endl;
//        cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//        cerr << "[RAM peak] after remove_from_mesh: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
////        int a; cin>>a; /// <====== WARNING!!!!
//    }

//    if(params.is_enable_debug_prints() && params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        print_container_content("current edge: ",e);
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[RAM peak] end of contract_edge procedure: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//    }
}

void Contraction_Simplifier_Weak::contract_weighted_edge(ivect &e, double e_weight, ET &et, VT &vt0, VT &vt1, leaf_VT &vtops, Node_Stellar &outer_v_block, edge_pqueue &edges,
                                           Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, w_contraction_parameters &params)
{
    bool top_outside_origin = false;

    difference_of_container_of_containers(vt0,et); // vt0 now contains the difference VT0 - ET
    difference_of_container_of_containers(vt1,et); // vt1 now contains the difference VT1 - ET

    // we have to check if edge is a top edge
    // this happens if #[ETop(e)]=0 and e \in VTop(v0)
    int e_top_id = -1;
    if(get_num_elements_in_container_of_containers(et) == 0)
    {
        Contraction_Simplifier_Weak::check_and_remove_top_edge(e,e_top_id,vt0,vt1,mesh,params);
    }

    if(params.swap_order())
    {
        if(params.is_enable_debug_prints())
            cout<<"SWAP order"<<endl;

        /// prior checking the d-1 faces we have to update the corresponding vtops relations
        top_outside_origin = Contraction_Simplifier_Weak::update(e[1],e[0],e_weight,vt1,vt0,n,outer_v_block,edges,mesh,params);
//        cout<<"update_mesh_leaf_and_VTop"<<endl;

        // we check if the boundary d-1 faces have an adjacent in the corresponding vt, if not we add the face to the top d-1 simplices list
        // of the leaf block indexing it and to the corresponding VTop (if local or cached)
        Contraction_Simplifier_Weak::check_for_new_tops(e[0],vt1,et,n,mesh,vtops,cache,params.get_tree(),params);
//        cout<<"check_for_new_tops"<<endl;

        // we remove v2 and the top simplices in et
        Contraction_Simplifier_Weak::remove_from_mesh(e[0],e_top_id,et,mesh,params);
        // finally we clear the VTop(v2)
        vt0.clear();

//        if(params.is_enable_debug_prints())
//            print_container_of_containers_content("updated-VT1: ",vt1);
        /// we remove the top simplices from the origin leaf block
        if(top_outside_origin)
            outer_v_block.compact_top_cell_arrays(mesh);
    }
    else
    {
        if(params.is_enable_debug_prints())
            cout<<"REGULAR order"<<endl;

        /// prior checking the d-1 faces we have to update the corresponding vtops relations
        top_outside_origin = Contraction_Simplifier_Weak::update(e[0],e[1],e_weight,vt0,vt1,outer_v_block,n,edges,mesh,params);
//        cout<<"update_mesh_leaf_and_VTop"<<endl;

        // we check if the boundary d-1 faces have an adjacent in the corresponding vt, if not we add the face to the top d-1 simplices list
        // of the leaf block indexing it and to the corresponding VTop (if local or cached)
        Contraction_Simplifier_Weak::check_for_new_tops(e[1],vt0,et,n,mesh,vtops,cache,params.get_tree(),params);
//        cout<<"check_for_new_tops"<<endl;

        // we remove v2 and the top simplices in et
        Contraction_Simplifier_Weak::remove_from_mesh(e[1],e_top_id,et,mesh,params);
        // finally we clear the VTop(v2)
        vt1.clear();

//        if(params.is_enable_debug_prints())
//            print_container_of_containers_content("updated-VT1: ",vt0);
        if(top_outside_origin)
            n.compact_top_cell_arrays(mesh);
    }

    et.clear();

//    if(params.is_enable_debug_prints() && params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        print_container_content("current edge: ",e);
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[RAM peak] end of contract_edge procedure: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//    }

//    return top_outside_origin;
}

void Contraction_Simplifier_Weak::get_edge_relations(ivect &e, ET &et, VT *&vt0, VT *&vt1, Node_Stellar *&outer_v_block,
                                                Node_Stellar &n, Simplicial_Mesh &mesh, leaf_VT &vtops, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params)
{
    Timer time;

    if(params.is_enable_debug_mode())
        time.start();

    outer_v_block = NULL;
    /// inverted order as I only need the block indexing v1
    vt1 = Contraction_Simplifier_Weak::get_VTop(e[1],n,mesh,vtops,cache,params.get_tree(),outer_v_block,params);
    vt0 = Contraction_Simplifier_Weak::get_VTop(e[0],n,mesh,vtops,cache,params.get_tree(),outer_v_block,params);

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.get_vts += time.get_elapsed_time();
        time.start();
    }

    Contraction_Simplifier_Weak::get_ETop(e,et,n,mesh,vtops);

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.get_et += time.get_elapsed_time();
    }
}

VT* Contraction_Simplifier_Weak::get_VTop(int v_id, Node_Stellar &n, Simplicial_Mesh &mesh, leaf_VT &vtops, LRU_Cache<int, leaf_VT> &cache,
                                          Stellar_Tree &tree, Node_Stellar *&v_block, contraction_parameters &params)
{
    int local_index;

    if(n.indexes_vertex(v_id))
    {
        if(params.is_enable_debug_prints()/* || v_id ==2355*/)
            cout<<"[get_VTop] "<<v_id<<" -> LOCAL VERTEX "<<n<<endl;

        local_index = v_id - n.get_v_start();
        VT* vt = &(vtops[local_index]);
        Contraction_Simplifier_Weak::clean_coboundary(*vt,mesh);
        v_block = &n;
        return vt;
    }
    else
    {
        if(params.is_enable_debug_prints()/* || v_id ==2355*/)
            cout<<"[get_VTop] "<<v_id<<" -> EXTERNAL VERTEX "<<n<<endl;

        tree.get_leaf_indexing_vertex(tree.get_root(),v_id,v_block);
        local_index = v_id - v_block->get_v_start();

        LRU_Cache<int,leaf_VT>::mapIt it_c = cache.find(v_block->get_v_start());
        if(it_c == cache.end())
        {
            if(params.is_enable_debug_prints()/* || v_id ==2355*/)
                cout<<"    -> LEAF BLOCK OUTSIDE CACHE - REGEN "<<*v_block<<endl;

            leaf_VT lVT;
            v_block->extract_local_VTop(mesh,lVT);
            it_c = cache.insert(v_block->get_v_start(),lVT);
        }
        else
        {
            if(params.is_enable_debug_prints()/* || v_id ==2355*/)
            {
                cout<<"    -> LEAF BLOCK IN CACHE - CLEAN "<<*v_block<<endl;
//                if(v_id==164)
//                {
////                    v_block->print_top_cells_lists(0,mesh);
////                    print_container_of_containers_content((it_c->second)[local_index]);
//                    Contraction_Simplifier_Weak::print_verbose_coboundary_relation((it_c->second)[local_index],mesh);
//                    leaf_VT lVT;
//                    v_block->extract_local_VTop(mesh,lVT);
//                    Contraction_Simplifier_Weak::print_verbose_coboundary_relation(lVT[local_index],mesh);
//                }
            }

            if(params.is_enable_debug_prints()/* || v_id ==2355*/)
                cout<<"num_elem_in_vtop: "<<get_num_elements_in_container_of_containers((it_c->second)[local_index])<<endl;
//                print_container_of_containers_content("VTop(2355) ",(it_c->second)[local_index]);

            Contraction_Simplifier_Weak::clean_coboundary((it_c->second)[local_index],mesh);

            if(params.is_enable_debug_prints()/* || v_id ==2355*/)
                cout<<"num_elem_in_vtop: "<<get_num_elements_in_container_of_containers((it_c->second)[local_index])<<endl;
//                print_container_of_containers_content("CleanedVTop(2355) ",(it_c->second)[local_index]);
        }

        return &(it_c->second)[local_index];
    }
}

void Contraction_Simplifier_Weak::get_ETop(ivect &e, ET &et, Node_Stellar &n, Simplicial_Mesh &mesh, leaf_VT &vtops)
{
    int other_v, local_v_id;
    if(n.indexes_vertex(e[0]))
    {
        other_v = e[1];
        local_v_id = e[0] - n.get_v_start();
    }
    else
    {
        other_v = e[0];
        local_v_id = e[1] - n.get_v_start();
    }

    VT &vt = vtops[local_v_id];

    et.assign(vt.size(), ivect());

    for(unsigned d=1; d<vt.size();d++) /// we skip the edges
    {
        for(unsigned i=0; i<vt[d].size();i++)
        {
            Top_Simplex &top = mesh.get_top_cell(d,vt[d][i]);
            if(top.has_vertex(other_v))
                et[d].push_back(vt[d][i]);
        }
    }
}

void Contraction_Simplifier_Weak::get_link_lite(int v_id, VT &vt, s_link &link, Simplicial_Mesh &mesh)
{
    link = s_link(vt.size()-1);

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
            ///(2) we add the vertices to the link of v_id
            link.insert_vertices(v_in_v_link);

            if(v_in_v_link.size() > 1)
            {
                ///(3) we add the simplex composed by those vertices to the link of v_id
                link.insert_simplex(v_in_v_link);
            }
        }
    }
}

void Contraction_Simplifier_Weak::get_link_lite(const ivect &e, ET &et, s_link &link, Simplicial_Mesh &mesh)
{
    link= s_link(et.size()-1);

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

            ///(2) we add the vertices to the link of v_id
            link.insert_vertices(v_in_e_link);

            if(v_in_e_link.size() > 1)
            {
                ///(3) we add the simplex composed by those vertices to the link of v_id
                link.insert_simplex(v_in_e_link);
            }
        }
    }
}

void Contraction_Simplifier_Weak::extract_sub_d_simplices(int d, s_link &link)
{
    set<ivect> &d_link = link.get_s_in_link(d);

    for(set<ivect>::iterator it=d_link.begin(); it!=d_link.end(); ++it)
    {
        Top_Simplex s = Top_Simplex(*it);
        for(int f=0; f<s.get_dfaces_num(); f++)
        {
            ivect d_face;
            s.TF(d_face,f);
            link.insert_simplex(d_face);
        }
    }
}

void Contraction_Simplifier_Weak::extract_d_simplices(int d, VT &cob, s_link &link, Simplicial_Mesh &mesh)
{
    for(unsigned i=d; i<cob.size(); i++)
    {
        for(ivect_iter it=cob[i].begin(); it!=cob[i].end(); ++it)
        {
            Top_Simplex &top = mesh.get_top_cell(i,*it);
            for(int s=0; s<top.get_sub_types_num(d+1); s++) /// +1 perche' la posizione nel cobounday e' differente rispetto a quella nel top
            {
                ivect sub;
                top.get_d_cell(sub,d+1,s); /// +1 (nel cob gli edge sono alla posizione 0, mentre nel top alla pos 1)
                link.insert_simplex(sub);
            }
        }
    }
}

void Contraction_Simplifier_Weak::clean_coboundary(VT &cob, Simplicial_Mesh &mesh)
{
    for(unsigned d=0; d<cob.size(); d++)
    {
        for(ivect_iter it=cob[d].begin(); it!=cob[d].end(); )
        {
            if(mesh.is_top_cell_removed(d,*it))
                cob[d].erase(it);
            else
                ++it;
        }
    }
}

void Contraction_Simplifier_Weak::check_and_remove_top_edge(const ivect &e, int &e_top_id, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, contraction_parameters &params)
{
    ivect &edges = vt0[0];

    for(ivect_iter it=edges.begin(); it!=edges.end(); ++it)
    {
        Top_Simplex &edge = mesh.get_top_cell(0,/*edges[i]*/*it);
        if(edge.are_equal(e))
        {
            if(params.is_enable_debug_prints())
                cout<<"[top-edge] removal: "<<e[0]<<" "<<e[1]<<endl;

            e_top_id = *it;
            vt0[0].erase(it);

            /// we have to remove the top edge also from the other vt in order to avoid misleading set unions
            ivect &edges1 = vt1[0];
            for(ivect_iter it1=edges1.begin(); it1!=edges1.end(); ++it1)
            {
                if(e_top_id == *it1)
                {
                    vt1[0].erase(it1);
                    break;
                }
            }

            return;
        }
    }
}

void Contraction_Simplifier_Weak::check_for_new_tops(int v, VT &vt, ET &et, Node_Stellar &n, Simplicial_Mesh &mesh,
                                               leaf_VT &vtops, LRU_Cache<int, leaf_VT> &cache, Stellar_Tree &tree, contraction_parameters &params)
{
//    vector< vector<ivect> > explicit_vtop; explicit_vtop.assign(vt.size(),vector<ivect>());
//    Contraction_Simplifier_Weak::explicit_top_coboundary(vt,explicit_vtop,mesh);

    for(unsigned d=1; d<et.size(); d++) //we skip the edges
    {
        ivect f;
        for(ivect_iter it=et[d].begin(); it!=et[d].end(); ++it)
        {
            Top_Simplex &top = mesh.get_top_cell(d,*it);
            int pos = top.vertex_index(v);
            top.TF(f,pos);

//            if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//            {
//                params.ram_peak = MemoryUsage().getValue_in_MB(false);
//                cerr << "[==INSIDE==] check_for_new_tops"<<endl;
//                cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//                cerr << "[RAM peak] after TF: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//        //        int a; cin>>a; /// <====== WARNING!!!!
//            }

            if(!Contraction_Simplifier_Weak::exists_simplex_in_coboundary(d,f,vt,mesh))
//            if(!Contraction_Simplifier_Weak::exists_simplex_in_coboundary_v2(d,f,explicit_vtop))
            {
//                if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//                {
//                    params.ram_peak = MemoryUsage().getValue_in_MB(false);
//                    cerr << "[==INSIDE==] check_for_new_tops"<<endl;
//                    cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//                    cerr << "[RAM peak] after exists_simplex_in_coboundary (true): " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//            //        int a; cin>>a; /// <====== WARNING!!!!
//                }

                int t_id = mesh.add_top_cell(d-1,f);
                if(params.is_enable_debug_prints())
                {
                    cout<<"[add top] "<<t_id<<"] ";
                    print_container_content(f);
                    cout<<endl;
                }


//                tree.insert_top_cell(tree.get_root(),mesh.get_top_cell(d-1,t_id),t_id,d-1);
                vector<Node_Stellar*> leaves;
                tree.get_leaves_indexing_cell(tree.get_root(),f,leaves);
                for(vector<Node_Stellar*>::iterator l_it=leaves.begin(); l_it!=leaves.end(); ++l_it)
                {
                    if(params.is_enable_debug_prints())
                        cout<<**l_it<<" | ";
                    (*l_it)->add_top_cell(d-1,t_id);
                }
                if(params.is_enable_debug_prints())
                    cout<<endl;

//                if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//                {
//                    params.ram_peak = MemoryUsage().getValue_in_MB(false);
//                    cerr << "[==INSIDE==] check_for_new_tops"<<endl;
//                    cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//                    cerr << "[RAM peak] after get_leaves_indexing_cell: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//                    cerr << "[STAT] leaves pointer array size: "<<leaves.size()<<endl;
//            //        int a; cin>>a; /// <====== WARNING!!!!
//                }

                /// finally, we add the new top simplex in the VTop of its vertices
                for(ivect_iter it_v= f.begin(); it_v!= f.end(); ++it_v)
                {
                    Contraction_Simplifier_Weak::add_new_top_in_VTop(*it_v,d-1,t_id,n,leaves,vtops,cache,params);
                }

//                if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//                {
//                    params.ram_peak = MemoryUsage().getValue_in_MB(false);
//                    cerr << "[==INSIDE==] check_for_new_tops"<<endl;
//                    cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//                    cerr << "[RAM peak] after add_new_top_in_VTop: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//            //        int a; cin>>a; /// <====== WARNING!!!!
//                }

//                if(params.is_enable_debug_prints() && f.size() == 2 && f[0]==171 && f[1]==493)
//                {
//                    cout<<"===TARGET EDGE==="<<endl;
//                    int a; cin>>a;
//                }
            }

//            if(params.ram_peak != MemoryUsage().getValue_in_MB(false))
//            {
//                params.ram_peak = MemoryUsage().getValue_in_MB(false);
//                cerr << "[==INSIDE==] check_for_new_tops"<<endl;
//                cerr << "[STAT] current leaf: "<< n.get_v_start() << " " << n.get_v_end() << " -- tot vertices: " << mesh.get_vertices_num() << endl;
//                cerr << "[RAM peak] after exists_simplex_in_coboundary (false): " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//        //        int a; cin>>a; /// <====== WARNING!!!!
//            }
        }
    }
}

bool Contraction_Simplifier_Weak::exists_simplex_in_coboundary(int start_d, const ivect &s, VT &vt, Simplicial_Mesh &mesh)
{
    for(unsigned top_dim = start_d; top_dim < vt.size(); top_dim++)
    {
        for(ivect_iter it_t=vt[top_dim].begin(); it_t!=vt[top_dim].end(); ++it_t)
        {
            Top_Simplex &top = mesh.get_top_cell(top_dim,*it_t);

            if(top.has_cell(s))
            {
                return true;
            }
        }
    }
    return false;
}

bool Contraction_Simplifier_Weak::exists_simplex_in_coboundary_v2(int start_d, ivect &s, vector<vector<ivect> > &explicit_vt)
{
    for(unsigned top_dim = start_d; top_dim < explicit_vt.size(); top_dim++)
    {
        for (auto top : explicit_vt[top_dim])
        {
            if(contains_container(top,s))
//            intersect_containers(top,s);
//            if(top == s)
//            auto res = search(begin(top), end(top), begin(s), end(s));
//            if(res != end(top))
                return true;
        }
    }
    return false;
}

void Contraction_Simplifier_Weak::add_new_top_in_VTop(int v_id, int d, int t_id, Node_Stellar &n, vector<Node_Stellar *> &leaves, leaf_VT &vtops,
                                                      LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params)
{
    int local_index;

    if(n.indexes_vertex(v_id))
    {
        if(params.is_enable_debug_prints())
            cout<<"[add_new_top_in_VTop] - local vtop update - "<<v_id<<" -> "<<t_id<<endl;

        local_index = v_id - n.get_v_start();
        vtops[local_index][d].push_back(t_id);
    }
    else
    {
        for(vector<Node_Stellar*>::iterator l_it=leaves.begin(); l_it!=leaves.end(); ++l_it)
        {
            if((*l_it)->indexes_vertex(v_id))
            {
                local_index = v_id - (*l_it)->get_v_start();
                LRU_Cache<int,leaf_VT>::mapIt it_c = cache.find((*l_it)->get_v_start());
                if(it_c != cache.end())
                {
                    if(params.is_enable_debug_prints())
//                        cout<<"[add_new_top_in_VTop] - cache update - "<<v_id<<" -> "<<t_id<<endl;
                        cout<<"[add_new_top_in_VTop] - cache update - v: "<<v_id<<" -> t: "<<t_id<<" -- d:"<<d<<" -- leaf: "<<(*l_it)->get_v_start()<<" -- lvid: "<<local_index<<endl;

                    leaf_VT &lvt = it_c->second;
                    lvt[local_index][d].push_back(t_id);

                    if(params.is_enable_debug_prints())
                        cout<<" ---- updated"<<endl;
                }
            }
        }
    }
}

void Contraction_Simplifier_Weak::update(const ivect &e, VT& vt, VT& difference, Node_Stellar &n, Node_Stellar &v_block,
                                                            edge_queue &edges, Simplicial_Mesh &mesh, contraction_parameters &params)
{
    set<ivect> e_set; /// we insert the new edges first in this set to avoid duplicate insertions in the queue

    for(unsigned d=0; d<difference.size(); d++)
    {
        for(ivect_iter it=difference[d].begin(); it!=difference[d].end(); ++it)
        {
            Top_Simplex &t = mesh.get_top_cell(d,*it);

            /// before updating the top simplex, we check
            /// if the leaf block indexing e[0] does not contain the current top d-simplex we have to add it
            /// NOTA: there is one possible case.. as leaf block n already indexes the top simplices in e[1]
            /// NOTA2: we have just to check that n and v_block are different (as if they are equal the edge is internal in n)
            if(n.get_v_start() != v_block.get_v_start() && !v_block.indexes_top(t))
            {
//                if(params.is_enable_debug_prints())
//                    cout<<"[add top to leaf] "<<t<<" -> "<<*v_block<<endl;
                v_block.add_top_cell(d,*it);
            }

            /// then we update the top simplex changing e[1] with e[0]
            int pos = t.vertex_index(e[1]);
            t.setTV(pos,e[0]);

            /// we have to add the new edges in the queue
            ivect new_e; new_e.assign(2,0);
            for(int i=0; i<t.get_vertices_num(); i++)
            {
                if(i!=pos)
                {
                    t.TE(new_e,pos,i);

                    if(n.indexes_vertex(new_e[1])) /// we process an edge only if it has alle the extrema already processed
                        e_set.insert(new_e);
                }
            }            
        }
    }

    /// we push the new "unique" edges in the queue
    for(auto it=e_set.begin(); it!=e_set.end(); ++it)
    {
        edges.push(*it);
    }

    /// finally we update the VTop relation of e[0]
    unify_container_of_containers(vt,difference);
}

bool Contraction_Simplifier_Weak::update(int v1, int v2, double e_weight, VT& vt1, VT &vt2, Node_Stellar &target, Node_Stellar &origin,
                                                            edge_pqueue &edges, Simplicial_Mesh &mesh, w_contraction_parameters &params)
{
    set<edge_weight> e_set; /// we insert the new edges first in this set to avoid duplicate insertions in the queue

//    if(target.get_v_start()==2353 && target.get_v_end()==2484)
//    {
//        cout<<target<<endl;
//        cout<<"contracting edge: "<<v1<<" "<<v2<<endl;
//        target.print_top_cells_list(2);
//        cout<<mesh.get_top_cell(2,6550)<<endl;
//        cout<<mesh.get_top_cell(2,6551)<<endl;
//        int a; cin>>a;
//    }
    bool top_outside_origin = false;

    for(unsigned d=0; d<vt2.size(); d++)
    {
        for(ivect_iter it=vt2[d].begin(); it!=vt2[d].end(); ++it)
        {
            Top_Simplex &t = mesh.get_top_cell(d,*it);

//            if(params.is_enable_debug_prints())
//            {
//                cout<<*it<<"] "<<t<<endl;
//                cout<<"origin: "<<origin<<endl;
//                cout<<"target: "<<target<<endl;
//            }

            /// before updating the top simplex, we check
            /// if the leaf block indexing e[0] does not contain the current top d-simplex we have to add it
            /// NOTA: there is one possible case.. as leaf block n already indexes the top simplices in v2
            /// NOTA2: we have just to check that n and v_block are different (as if they are equal the edge is internal in n)
            if(origin.get_v_start() != target.get_v_start() && !target.indexes_top(t))
            {
                if(params.is_enable_debug_mode())
                {
                    if(target.check_duplicate_top_cells_array(d,*it))
                    {
                        cout<<target<<endl;
                        target.print_top_cells_array(d);
                        cout<<"[add top to leaf] -> "<<*it<<"] "<<t<<endl;
                        int a; cin>>a;
                    }
                }
//                if(params.is_enable_debug_prints())
////                if(target.get_v_start()==2353 && target.get_v_end()==2484 && d==2 && (*it==6550 || *it==6551))
//                {
//                    cout<<"[add top to leaf] -> "<<*it<<"] "<<t<<endl;
//                    int a; cin>>a;
//                }
                target.add_top_cell(d,*it);
            }

            /// then we update the top simplex changing e[1] with e[0]
            int pos = t.vertex_index(v2);
            t.setTV(pos,v1);

//            if(params.is_enable_debug_prints())
//                cout<<"UPDATED "<<*it<<"] "<<t<<endl;
            /// if after updating the top simplex it goes out the origin block
            if(!top_outside_origin && origin.get_v_start() != target.get_v_start() && !origin.indexes_top(t))
            {
                /// we activate the flag to update the top simplices array
                top_outside_origin = true;
            }

            /// in the next loop we have also to update the entries in the edge history map

            /// we have to add the new edges in the queue
            ivect new_e; new_e.assign(2,0);
            ivect old_e; old_e.assign(2,0);
            for(int i=0; i<t.get_vertices_num(); i++)
            {
                if(i!=pos)
                {
                    /// first we recover from the history map of the deleted
                    old_e[0] = v2;
                    old_e[1] = abs(t.TV(i));
                    sort(old_e.begin(),old_e.end());

                    /// controllo che non ha piu' senso da quando rimuovo la entry dell'edge vecchio
                    ///
//                    if(params.find(old_e)==params.end())
//                    {
//                        cout<<"NON HO EDGE IN CACHE STORICO: "<<old_e[0]<<" "<<old_e[1]<<endl;
//                        int a; cin>>a;
//                    }

                    /// 1) check if the old-edge has been already deleted from the map history
                    /// if so, do nothing
                    edge_history_map::iterator it_old = params.find(old_e);
                    if(it_old!=params.end())
                    {
                        edge_history old_h = it_old->second;
//                                            params.get_edge_history(old_e);
//                        if(params.is_enable_debug_prints())
//                            cout<<"old_h: "<<old_h<<endl;

                        t.TE(new_e,pos,i);

                        /// the new edge gets an entry in the edge history map
                        /// we have to check if we have already inserted the edge in the map
                        edge_history_map::iterator it = params.find(new_e);
                        if(it == params.end())
                        {
                            edge_history eh = old_h;
                            eh.weight_summation += old_h.weight_summation;
                            eh.last_added_weight = e_weight;

//                            if(params.is_enable_debug_prints())
//                                cout<<"new-eh: "<<eh<<endl;

                            params.add_edge(new_e,eh);
                        }

                        /// we process an edge only if it has alle the extrema already processed
                        if((params.swap_order() && target.indexes_vertex(new_e[1])) ||
                                origin.indexes_vertex(new_e[1]))
                        {
                            edge_weight ew_add(new_e,old_h.initial_weight);
                            e_set.insert(ew_add);
                        }

                         /// not using the it_old iterator as this can be become invalid after the insertion of new_e
                        params.remove_edge(old_e);
                    }
                }
            }

            /// NOTA: we do not have to check if a top d-simplex goes outside the current leaf block n
            ///       as this is done when entering/exiting
        }
    }

    /// we push the new "unique" edges in the queue
    for(auto it=e_set.begin(); it!=e_set.end(); ++it)
    {
        edges.push(*it);
    }

//    if(params.is_enable_debug_prints())
//    {
//        print_container_of_containers_content("VT1: ",vt1);
//        print_container_of_containers_content("VT2: ",vt2);
//    }

    /// finally we update the VTop relation of e[0]
    unify_container_of_containers(vt1,vt2);

//    if(params.is_enable_debug_prints())
//        print_container_of_containers_content("updated-VT1: ",vt1);

//    if(target.get_v_start()==2353 && target.get_v_end()==2484)
//    {
//        int a; cin>>a;
//    }
    return top_outside_origin;
}



void Contraction_Simplifier_Weak::remove_from_mesh(int to_delete_v, int e_top_id, ET &et, Simplicial_Mesh &mesh, contraction_parameters &params)
{
    for(unsigned d=0; d<et.size(); d++)
    {
        for(ivect_iter it=et[d].begin(); it!=et[d].end(); ++it)
        {
            mesh.remove_top_cell(d,*it);
            params.increment_counter(d);
        }
    }
    // if the current edge is a top-edge, then, we have to delete it from the mesh
    if(e_top_id != -1)
    {
        mesh.remove_top_cell(0,e_top_id);
        params.increment_counter(0);
    }

    mesh.remove_vertex(to_delete_v);
    params.increment_contracted_edges_counter();
}

void Contraction_Simplifier_Weak::update_mesh_and_tree(Stellar_Tree &tree, Simplicial_Mesh &mesh, contraction_parameters &params)
{
    Timer time;

    ///  UPDATE OF MESH AND TREE
    ivect new_v_positions;
    vector<ivect> new_top_positions;
    ivect surviving_vertices;

    time.start();
//    cerr<<"[TREE] compact vertices lists"<<endl;
    tree.compact_vertices_lists(tree.get_root(),mesh,surviving_vertices);
    time.stop();
    time.print_elapsed_time("[TIME] Compact tree vertices lists: ");

//    print_container_content("surviving vertices: ",surviving_vertices);
//    mesh.print_mesh(cout);
//    int a; cin>>a;

    time.start();
//    cerr<<"[MESH] compact"<<endl;
    Mesh_Updater mu;
    mu.clean_vertices_array(mesh,new_v_positions,surviving_vertices);
    /// NEW: the update_and_compact procedure check internally if we have removed all the top d-simplices
    boost::dynamic_bitset<> all_deleted = mu.update_and_clean_top_cells_arrays(mesh,new_v_positions,new_top_positions,params.get_counters());
    time.stop();
    time.print_elapsed_time("[TIME] Compact and update mesh: ");

    cerr<<"[STAT] mesh "<<endl;
    cerr<<"------------"<<endl;
    mesh.print_mesh_stats(cerr);
    cerr<<"------------"<<endl;

    time.start();
//    cerr<<"[TREE] update indices in the tree"<<endl;
    tree.update_tree(tree.get_root(),new_v_positions,new_top_positions,all_deleted);
    time.stop();
    time.print_elapsed_time("[TIME] Update tree (top-simplices): ");

//    Reindexer r;
//    r.reorganize_index_and_mesh(tree,mesh,false);

    cerr << "[RAM peak] for updating the mesh and the tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
}

///////////////////////// FOR DEBUG ONLY ///////////////////////////////////////////
int Contraction_Simplifier_Weak::get_top_edge_id(const ivect &e, VT &vt0, Simplicial_Mesh &mesh)
{
    ivect &edges = vt0[0];

    for(ivect_iter it=edges.begin(); it!=edges.end(); ++it)
    {
        Top_Simplex &edge = mesh.get_top_cell(0,*it);
        if(edge.are_equal(e))
        {
            return *it;
        }
    }
    return -1;
}

bool Contraction_Simplifier_Weak::exists_degenerate_tops(Node_Stellar &n, Simplicial_Mesh &mesh)
{
    bool deg = false;
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;
            Top_Simplex &top = mesh.get_top_cell(d,*t_id);
            if(!mesh.is_top_cell_removed(top) && top.is_degenerate())
            {
                deg = true;
                cout<<"[degenerate top-simplex] "<<top<<endl;
            }
        }
    }
    return deg;
}

bool Contraction_Simplifier_Weak::exists_degenerate_tops(VT& vt, Simplicial_Mesh &mesh)
{
    bool deg = false;
    for(unsigned d=0; d<vt.size(); d++)
    {
        for(ivect_iter it=vt[d].begin(); it!=vt[d].end(); ++it)
        {
            Top_Simplex &top = mesh.get_top_cell(d,*it);
            if(!mesh.is_top_cell_removed(top) && top.is_degenerate())
            {
                deg = true;
                cout<<"[degenerate top-simplex] "<<top<<endl;
            }
        }
    }
    return deg;
}

bool Contraction_Simplifier_Weak::exists_duplicated_tops(VT& vt, Simplicial_Mesh &mesh)
{
    bool deg = false;
    for(unsigned d=0; d<vt.size(); d++)
    {
        for(unsigned i=0; i<vt[d].size(); i++)
        {
            if(!mesh.is_top_cell_removed(d,vt[d][i]))
            {
                Top_Simplex &top = mesh.get_top_cell(d,vt[d][i]);

                for(unsigned j=i+1; j<vt[d].size(); j++)
                {
                    Top_Simplex &top2 = mesh.get_top_cell(d,vt[d][j]);

                    if(!mesh.is_top_cell_removed(top2) && top == top2)
                    {
                        deg = true;
                        cout<<"[duplicated top-simplex] "<<top<<endl;
                    }
                }
            }
        }
    }
    return deg;
}

//void Contraction_Simplifier_Weak::remove_duplicated_tops(VT& vt, Simplicial_Mesh &mesh)
//{
//    for(unsigned d=0; d<vt.size(); d++)
//    {
//        for(unsigned i=0; i<vt[d].size(); i++)
//        {
//            if(!mesh.is_top_cell_removed(d,vt[d][i]))
//            {
//                Top_Simplex &top = mesh.get_top_cell(d,vt[d][i]);

//                for(unsigned j=i+1; j<vt[d].size(); j++)
//                {
//                    Top_Simplex &top2 = mesh.get_top_cell(d,vt[d][j]);

//                    if(!mesh.is_top_cell_removed(top2) && top == top2)
//                    {
//                        mesh.remove_top_cell(d,vt[d][j]);
//                    }
//                }
//            }
//        }
//    }
//}

//void Contraction_Simplifier_Weak::remove_degenerate_tops(VT& vt, Simplicial_Mesh &mesh, bool print)
//{
//    ///
//    for(unsigned d=0; d<vt.size(); d++)
//    {
//        for(ivect_iter it=vt[d].begin(); it!=vt[d].end(); ++it)
//        {
//            if(!mesh.is_top_cell_removed(d,*it))
//            {
//                Top_Simplex &top = mesh.get_top_cell(d,*it);

//                if(print)
//                    cout<<*it<<"] "<<top<<endl;

//                if(top.is_degenerate())
//                {
//                    if(print)
//                        cout<<" is degenerate"<<endl;
//                    /// if a top is degenerate we remove it
//                    mesh.remove_top_cell(d,*it);
//                }
//                else
//                {
//                    ivect vids;
//                    top.get_positive_sorted_vertices(vids);
//                    /// we check if a top is not a real top..
//                    /// i.e. is in the boundary of another top
//                    /// then we remove it
//                    bool found = false;
//                    for(unsigned d2=d+1; d2<vt.size(); d2++)
//                    {
//                        for(ivect_iter it2=vt[d2].begin(); it2!=vt[d2].end(); ++it2)
//                        {
//                            if(!mesh.is_top_cell_removed(d2,*it2))
//                            {
//                                Top_Simplex &top2 = mesh.get_top_cell(d2,*it2);
//                                if(top2.has_cell(vids))
//                                {
//                                    mesh.remove_top_cell(d,*it);
//                                    found = true;
//                                }
//                            }
//                            if(found)
//                                break;
//                        }
//                        if(found)
//                            break;
//                    }
//                }
//            }
//        }
//    }

//}

bool Contraction_Simplifier_Weak::everything_is_synced(const ivect &e, VT &vt0, Node_Stellar &n, Node_Stellar *&v_block, Stellar_Tree &tree, LRU_Cache<int,leaf_VT> &cache)
{
    /// we have to sync only for cross edges
    if(n.get_v_start() != v_block->get_v_start())
    {
        /// (1) check if the v_block representation is synced with the leaf block of the tree
        Node_Stellar *target_leaf;
        tree.get_leaf_indexing_vertex(tree.get_root(),e[0],target_leaf);

        vector<ivect> v_block_expanded_lists;
        v_block_expanded_lists.assign(v_block->get_num_top_cells_encoded(),ivect());
        for(int d=0; d<v_block->get_num_top_cells_encoded();d++)
        {
            for(RunIteratorPair itPair = v_block->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
                v_block_expanded_lists[d].push_back(*itPair.first);
        }

        vector<ivect> target_leaf_expanded_lists;
        target_leaf_expanded_lists.assign(target_leaf->get_num_top_cells_encoded(),ivect());
        for(int d=0; d<target_leaf->get_num_top_cells_encoded();d++)
        {
            for(RunIteratorPair itPair = target_leaf->make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
                target_leaf_expanded_lists[d].push_back(*itPair.first);
        }

        sort_container_of_containers(v_block_expanded_lists);
        sort_container_of_containers(target_leaf_expanded_lists);

        for(unsigned d=0; d<v_block_expanded_lists.size(); d++)
        {
            if(v_block_expanded_lists[d] != target_leaf_expanded_lists[d])
            {
                cout<<"[error] the leaf block representation is not synced"<<endl;
                v_block->print_top_cells_arrays();
                target_leaf->print_top_cells_arrays();

                return false;
            }
        }

        /// (2) check if the local VT is synced with its cache representation
        LRU_Cache<int,leaf_VT>::mapIt it = cache.find(v_block->get_v_start());
        VT &vt_c = (it->second)[e[0]-v_block->get_v_start()];

        sort_container_of_containers(vt0);
        sort_container_of_containers(vt_c);

        for(unsigned d=0; d<vt0.size(); d++)
        {
            if(vt0[d] != vt_c[d])
            {
                cout<<"[error] the VT representation is not synced"<<endl;
                print_container_of_containers_content(vt0);
                print_container_of_containers_content(vt_c);

                return false;
            }
        }
    }
    return true;
}

bool Contraction_Simplifier_Weak::cache_is_synced(Node_Stellar &n, Simplicial_Mesh &mesh, Stellar_Tree &tree, LRU_Cache<int,leaf_VT> &cache)
{
    for(LRU_Cache<int,leaf_VT>::mapIt it=cache.begin(); it!=cache.end(); ++it)
    {
        Node_Stellar *target_leaf;
        tree.get_leaf_indexing_vertex(tree.get_root(),it->first,target_leaf);
//        leaf_VT &cached_vt = (it->second);

        leaf_VT vt;
        target_leaf->extract_local_VTop(mesh,vt);

//        cout<<"target vt"<<endl;
//        for(auto it=vt.begin(); it!=vt.end(); ++it)
//        {
//            print_container_of_containers_content(*it);
//        }
//        cout<<"cached vt"<<endl;
//        for(auto it=cached_vt.begin(); it!=cached_vt.end(); ++it)
//        {
//            print_container_of_containers_content(*it);
//        }

        for(unsigned v=0; v!=vt.size(); v++)
        {
            if(!mesh.is_vertex_removed(target_leaf->get_v_start()+v))
            {
                sort_container_of_containers(vt[v]);
                sort_container_of_containers((it->second)[v]);

                if(vt[v] != (it->second)[v])
                {
                    VT vl = vt[v];
                    VT vc = (it->second)[v];


                    int num_elem_leaf = get_num_elements_in_container_of_containers(vt[v]);
                    int num_elem_cache = get_num_elements_in_container_of_containers((it->second)[v]);

                    if(num_elem_cache < num_elem_leaf)
                    {
                        difference_of_container_of_containers(vl,vc);

                        Contraction_Simplifier_Weak::clean_coboundary(vl,mesh);
                        if(get_num_elements_in_container_of_containers(vl)>0)
                        {
                            cout<<"[error] the VT representation is not synced"<<endl;
                            cout<<"checked in leaf block: "<<n<<endl;
                            cout<<"cache size: "<<cache.get_actual_size()<<endl;
                            cout<<"leaf block: "<<*target_leaf<<endl;

                            cout<<"v("<<target_leaf->get_v_start()+v<<") leaf->"<<num_elem_leaf<<" cache->"<<num_elem_cache<<endl;
                            cout<<"==leaf representation=="<<endl;
                            print_container_of_containers_content(vt[v]);
                            cout<<"==cache representation=="<<endl;
                            print_container_of_containers_content((it->second)[v]);
                            cout<<"==missing top=="<<endl;
                            Contraction_Simplifier_Weak::print_verbose_coboundary_relation(vl,mesh);
                            return false;
                        }
                    }
                    else
                    {
                        difference_of_container_of_containers(vc,vl);

                        Contraction_Simplifier_Weak::clean_coboundary(vc,mesh);
                        if(get_num_elements_in_container_of_containers(vc)>0)
                        {
                            cout<<"[error] the VT representation is not synced"<<endl;
                            cout<<"cache size: "<<cache.get_actual_size()<<endl;
                            cout<<"leaf block: "<<*target_leaf<<endl;

                            cout<<"v("<<target_leaf->get_v_start()+v<<") leaf->"<<num_elem_leaf<<" cache->"<<num_elem_cache<<endl;
                            cout<<"==leaf representation=="<<endl;
                            print_container_of_containers_content(vt[v]);
                            cout<<"==cache representation=="<<endl;
                            print_container_of_containers_content((it->second)[v]);
                            cout<<"==missing top=="<<endl;
                            Contraction_Simplifier_Weak::print_verbose_coboundary_relation(vc,mesh);
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}

void Contraction_Simplifier_Weak::compute_cache_stats(LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params)
{
    if(params.cache_max_size < cache.get_actual_size())
        params.cache_max_size = cache.get_actual_size();
}

void Contraction_Simplifier_Weak::print_verbose_coboundary_relation(VT &vt, Simplicial_Mesh &mesh)
{
    int d=0;
    for(VT::iterator dit=vt.begin(); dit!=vt.end(); ++dit)
    {
        if(dit->size() > 0)
            cout<<"D"<<d<<endl;
        for(ivect_iter tit=dit->begin(); tit!=dit->end(); ++tit)
            cout<<"   "<<*tit<<"] "<<mesh.get_top_cell(d,*tit)<<endl;
        d++;
//        int a; cin>>a;
    }
}

bool Contraction_Simplifier_Weak::is_degenerate_edge(const ivect &e,Node_Stellar &, Simplicial_Mesh &)
{
    if(e[0] == e[1])
    {
//        cerr<<"[is_degenerate_edge] YES"<<endl;
//        cerr<<n<<endl;
//        if(n.indexes_vertex(e[0]))
//            cerr<<" INDEXED vertex"<<endl;
//        else
//            cerr<<" NOT INDEXED vertex"<<endl;
////        n.print_top_cells_arrays(mesh);
        return true;
    }
    return false;
}

/// just for debug
//void Contraction_Simplifier_Weak::explicit_top_coboundary(VT &vt, vector< vector<ivect> > &explicit_top_cob, Simplicial_Mesh &mesh)
//{
//    ivect v;
//    for(unsigned i=0; i<vt.size(); i++)
//    {
//        for(auto it=vt[i].begin(); it!=vt[i].end(); ++it)
//        {
//            Top_Simplex &top = mesh.get_top_cell(i,*it);
////            cout<<top<<endl;
//            top.get_positive_sorted_vertices(v);
//            explicit_top_cob[i].push_back(v);
//        }
//    }
//}
