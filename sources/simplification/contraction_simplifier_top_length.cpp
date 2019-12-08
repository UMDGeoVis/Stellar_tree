#include "contraction_simplifier_top.h"

void Contraction_Simplifier_Top::simplify_length(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli)
{
    cerr<<"==Length-based Homology preserving simplification - top-link condition=="<<endl;
    mesh.print_mesh_stats(cerr);
//    cout<<mesh.get_top_cells_types()<<endl;
//    cout<<mesh.get_top_cell(1,1).get_vertices_num()<<endl;
    Writer::write_mesh_VTK(string_management::get_path_without_file_extension(cli.mesh_path),mesh);

    /// (COMMENTED FEB, 26TH, 2019)
    /// THERE IS NO QUADRIC COMPUTATION INVOLVED... SO THIS SHOULD BE NOT NECESSARY
    /// FOR VISUALIZATION PURPOSES THIS SHOULD WORK ON BOTH TRI AND TETRA MESH
//    if(mesh.get_top_cells_types() > 2 ||
//            (mesh.get_top_cells_types() == 2 && (mesh.get_top_cells_num(0) > 1 || mesh.get_top_cell(1,1).get_vertices_num() != 3)))
//    {
//        cerr << "[ERROR] The simplification based on the quadrics is defined only for pure triangular meshes." << endl;
//        return;
//    }

    LRU_Cache<int,leaf_VT> cache(cli.cache_size); // the key is v_start while the value are the VTop relations
    contraction_parameters params(&tree,mesh.get_top_cells_types());

    if(cli.debug_prints_mode)
        params.enable_debug_prints();
    if(cli.debug_mode)
        params.enable_debug_mode();

    Timer time;
    int simplification_round;
    int round = 1;

    tree.visit(Contraction_Simplifier_Weak::extract_edges_lengths,tree.get_root(),mesh,params);
    params.avg_err /= params.num_edges;
    cout<<"min_edge_l: "<<params.min_err<<" - avg_edge_l: "<<params.avg_err<<" - max_edge_l: "<<params.max_err<<endl;
    tree.visit(Contraction_Simplifier_Weak::extract_edge_length_stats,tree.get_root(),mesh,params);
    cout<<"exactly_min: "<<params.min_c<<" - exactly_max: "<<params.max_c<<endl;
    cout<<"distribution: "<<params.histogram[0]<<" "<<params.histogram[1]<<" "<<params.histogram[2]<<" "<<params.histogram[3]<<" "<<endl;

//    cout<<"range: " <<fabs(q_max-q_min)<<" "<<fabs(q_max-q_min)*cli.eps<<endl;
//    params.threshold = (params.min_err < 0) ? params.min_err + fabs(params.max_err-params.min_err) * cli.eps :
//                                              params.min_err - fabs(params.max_err-params.min_err) * cli.eps;
//    params.threshold = params.avg_err * 1.75;
//    params.threshold = params.avg_err / 1.25;
    if(cli.eps == -1)
        params.threshold = INT_MAX;//params.max_err;
    else
        params.threshold = cli.eps;
    cout<<"simplificafion threshold: "<<params.threshold<<endl;
//    int a; cin>>a;

//    cout<<"HERE-B"<<endl;

    time.start();
    while(1)
    {
        simplification_round = params.get_contracted_edges_num();
        tree.visit_with_cache(Contraction_Simplifier_Top::simplify_tri_length_leaf,tree.get_root(),mesh,cache,params);

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

        // for debug -- re-compute the min-max qem of the edges and set accordingly the simplification threshold
//        params.qem_min = INT_MAX; params.qem_max = INT_MIN;
//        tree.visit(Contraction_Simplifier_Weak::extract_edge_qem_stats,tree.get_root(),mesh,params);
//        cout<<"min_edge_qem: "<<params.qem_min<<" - max_edge_qem: "<<params.qem_max<<endl;
//        params.threshold = (params.qem_min < 0) ? params.qem_min + fabs(params.qem_max-params.qem_min) * cli.eps :
//                                                  params.qem_min - fabs(params.qem_max-params.qem_min) * cli.eps;
//        cout<<"updated-simplificafion threshold: "<<params.threshold<<endl;

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

    stringstream ss;
    if(cli.eps == -1)
        ss<<"length_Top_simplification_max";
    else
        ss<<"length_Top_simplification_"<<params.threshold;
    Writer::write_mesh_VTK(string_management::get_path_without_file_extension(cli.mesh_path),ss.str(),cli.vertices_per_leaf,mesh);

    double initial_min = params.min_err;
    params.min_err = INT_MAX; params.max_err = INT_MIN; params.avg_err = 0;
    tree.visit(Contraction_Simplifier_Weak::extract_edges_lengths,tree.get_root(),mesh,params);
    params.avg_err /= params.num_edges;
    cout<<"min_edge_l: "<<params.min_err<<" - avg_edge_l: "<<params.avg_err<<" - max_edge_l: "<<params.max_err<<endl;
    tree.visit(Contraction_Simplifier_Weak::extract_edge_length_stats,tree.get_root(),mesh,params);
    cout<<"exactly_min: "<<params.min_c<<" - exactly_max: "<<params.max_c<<endl;
    cout<<"distribution: "<<params.histogram[0]<<" "<<params.histogram[1]<<" "<<params.histogram[2]<<" "<<params.histogram[3]<<" "<<endl;

    if(cli.eps == -1)
    {
        cout<<"thresholds for the next rounds"<<endl;
        cout<<"range: " <<fabs(params.max_err-initial_min)<<endl;
        params.threshold = //(params.min_err < 0) ? params.min_err + fabs(params.max_err-params.min_err) * 0.01 :
                                                  initial_min + fabs(params.max_err-initial_min) * 0.01;
        cout<<"1% threshold: "<<params.threshold<<endl;
        params.threshold = //(params.min_err < 0) ? params.min_err + fabs(params.max_err-params.min_err) * 0.02 :
                                                  initial_min + fabs(params.max_err-initial_min) * 0.02;
        cout<<"2% threshold: "<<params.threshold<<endl;
        params.threshold = //(params.min_err < 0) ? params.min_err + fabs(params.max_err-params.min_err) * 0.03 :
                                                  initial_min + fabs(params.max_err-initial_min) * 0.03;
        cout<<"3% threshold: "<<params.threshold<<endl;
        params.threshold = //(params.min_err < 0) ? params.min_err + fabs(params.max_err-params.min_err) * 0.05 :
                                                  initial_min + fabs(params.max_err-initial_min) * 0.05;
        cout<<"5% threshold: "<<params.threshold<<endl;
        params.threshold = //(params.min_err < 0) ? params.min_err + fabs(params.max_err-params.min_err) * 0.1 :
                                                  initial_min + fabs(params.max_err-initial_min) * 0.1;
        cout<<"10% threshold: "<<params.threshold<<endl;
        params.threshold = //(params.min_err < 0) ? params.min_err + fabs(params.max_err-params.min_err) * 0.15 :
                                                  initial_min + fabs(params.max_err-initial_min) * 0.15;
        cout<<"15% threshold: "<<params.threshold<<endl;
        params.threshold = //(params.min_err < 0) ? params.min_err + fabs(params.max_err-params.min_err) * 0.2 :
                                                  initial_min + fabs(params.max_err-initial_min) * 0.2;
        cout<<"20% threshold: "<<params.threshold<<endl;
    }

}

void Contraction_Simplifier_Top::simplify_tri_length_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params)
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

    edge_pqueue edges;
    Contraction_Simplifier_Weak::extract_target_edges_length(n,mesh,edges,params);

    if(params.is_enable_debug_mode())
    {
        time.stop();
        params.get_edges += time.get_elapsed_time();

        /// get edge queue statistics
        if(params.edge_q_max_size < (int)edges.size())
            params.edge_q_max_size = edges.size();
    }

//    bool simplified = false;

//    cout<<n<<endl;
//    if(edges.empty())
//        cout<<"--> NO EDGES"<<endl;

    while(!edges.empty())
    {
        edge_weight ew = edges.top();
        edges.pop();

//        cout<<ew<<endl;

        /// for debug only
        if(params.is_enable_debug_mode())
        {
            if(Contraction_Simplifier_Weak::is_degenerate_edge(ew.get_edge(),n,mesh))
            {
                int a; cin>>a;
            }
        }

        if(Contraction_Simplifier_Weak::skip_edge(ew.get_edge(),mesh))
        {
            continue;
        }

        /// get VT0 -- VT1 -- ET and (if exists) the leaf block indexing the outer extreme)
        Node_Stellar *outer_v_block;
        VT *vt1, *vt0;
        ET et;
        Contraction_Simplifier_Weak::get_edge_relations(ew.get_edge(),et,vt0,vt1,outer_v_block,n,mesh,vtops,cache,params);

        if(params.is_enable_debug_prints() /*&& params.ram_peak != MemoryUsage().getValue_in_MB(false)*/)
        {
            cout<<n<<endl;
            print_container_content("current edge: ",ew.get_edge());
//            params.ram_peak = MemoryUsage().getValue_in_MB(false);
//            cerr << "[RAM peak] for gathering the topological relations of an edge: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        }

//        if(Contraction_Simplifier_Top::is_tops_link_condition_valid(ew.get_edge(),et,*vt0,*vt1,mesh,params))
        if(Contraction_Simplifier_Top::is_tops_link_condition_valid_v2(ew.get_edge(),et,*vt0,*vt1,mesh,params))
        {
//            simplified = true;
            if(params.is_enable_debug_mode())
                time.start();

//            cout<<"   CONTRACTED"<<endl;

            Contraction_Simplifier_Weak::contract_edge_length(ew.get_edge(),et,*vt0,*vt1,vtops,*outer_v_block,edges,n,mesh,cache,params);

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

//    int a; cin>>a;

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

////    /// DEBUGGING
//    if(simplified)
//    {
//        cout<<n<<endl;
//        cout<<"SIMPLIFICATION EXECUTED.. WAITING FOR NEXT STEP"<<endl;
////        Contraction_Simplifier_Weak::update_mesh_and_tree(params.get_tree(),mesh,params);
////        cout<<"PRINTING MESH"<<endl;
////        Writer::write_2D_3D_mesh_VTK("one_leaf_simplification",mesh);
//        int a; cin>>a;
//    }
}
