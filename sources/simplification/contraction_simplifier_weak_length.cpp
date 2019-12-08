#include "contraction_simplifier_weak.h"

void Contraction_Simplifier_Weak::simplify_length(Stellar_Tree &tree, Simplicial_Mesh &mesh, cli_parameters &cli)
{
    cerr<<"==Length-based Homology preserving simplification - weak-link condition=="<<endl;
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
        tree.visit_with_cache(Contraction_Simplifier_Weak::simplify_tri_length_leaf,tree.get_root(),mesh,cache,params);

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
        ss<<"length_Weak_simplification_max";
    else
        ss<<"length_Weak_simplification_"<<params.threshold;
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

void Contraction_Simplifier_Weak::simplify_tri_length_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int,leaf_VT> &cache, contraction_parameters &params)
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

        if(Contraction_Simplifier_Weak::is_weak_link_condition_valid_top_down_version(ew.get_edge(),et,*vt0,*vt1,mesh,params))
//        if(Contraction_Simplifier::is_weak_link_condition_valid_bottom_up_version(e,et,*vt0,*vt1,mesh,params))
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

void Contraction_Simplifier_Weak::extract_target_edges_length(Node_Stellar &n, Simplicial_Mesh &mesh, edge_pqueue &edges, contraction_parameters &params)
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
//                        cout<<params.initialQuadric[e[0]-1][0]<<" "<<params.initialQuadric[e[1]-1][0]<<" : "<<params.threshold<<endl;
                        if(n.indexes_vertex(e[1]))  /// we process an edge only if it has all the extrema already processed
                        {
//                            cout<<"inserting edge: "<<e[0]<<" "<<e[1]<<endl;
                            auto r = e_set.insert(e);
                            if(r.second)
                            {
                                double error = mesh.get_vertex(e[0]).distance(mesh.get_vertex(e[1]));
                                if(error <= params.threshold)
//                                if(error >= params.threshold)
                                {
                                    edge_weight ew = edge_weight(e,error);
                                    edges.push(ew);
                                }
//                                cout<<e[0]<<" "<<e[1]<<" -length: "<<error<<" -threshold: "<<params.threshold<<endl;
//                                int a; cin>>a;
                            }
                        }
                    }
                }
            }
        }
    }
}

void Contraction_Simplifier_Weak::contract_edge_length(ivect &e, ET &et, VT &vt0, VT &vt1, leaf_VT &vtops, Node_Stellar &outer_v_block, edge_pqueue &edges,
                                           Node_Stellar &n, Simplicial_Mesh &mesh, LRU_Cache<int, leaf_VT> &cache, contraction_parameters &params)
{
    difference_of_container_of_containers(vt0,et); // vt0 now contains the difference VT0 - ET
    difference_of_container_of_containers(vt1,et); // vt1 now contains the difference VT1 - ET

    // we have to check if edge is a top edge
    // this happens if #[ETop(e)]=0 and e \in VTop(v0)
    int e_top_id = -1;
    if(get_num_elements_in_container_of_containers(et) == 0)
    {
        Contraction_Simplifier_Weak::check_and_remove_top_edge(e,e_top_id,vt0,vt1,mesh,params);
    }

    /// prior checking the d-1 faces we have to update the corresponding vtops relations
    Contraction_Simplifier_Weak::update_length(e,vt0,vt1,n,outer_v_block,edges,mesh,params);

    // we check if the boundary d-1 faces have an adjacent in the corresponding vt, if not we add the face to the top d-1 simplices list
    // of the leaf block indexing it and to the corresponding VTop (if local or cached)
    Contraction_Simplifier_Weak::check_for_new_tops(e[1],vt0,et,n,mesh,vtops,cache,params.get_tree(),params);

    // extra: update the coordinate position of v1 (set it as the mid point of edge e
    Vertex<COORDBASETYPE> &v0 = mesh.get_vertex(e[0]);
    Vertex<COORDBASETYPE> &v1 = mesh.get_vertex(e[1]);
    for(int c=0; c<v0.get_dimension(); c++)
    {
        v0.setC(c,(v0.getC(c)+v1.getC(c))/2.0);
    }

    // we remove v2 and the top simplices in et
    Contraction_Simplifier_Weak::remove_from_mesh(e[1],e_top_id,et,mesh,params);
    // finally we clear the VTop(v2)
    vt1.clear();
    et.clear();

//    if(params.is_enable_debug_prints() && params.ram_peak != MemoryUsage().getValue_in_MB(false))
//    {
//        print_container_content("current edge: ",e);
//        params.ram_peak = MemoryUsage().getValue_in_MB(false);
//        cerr << "[RAM peak] end of contract_edge procedure: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//    }
}

void Contraction_Simplifier_Weak::update_length(const ivect &e, VT& vt, VT& difference, Node_Stellar &n, Node_Stellar &v_block,
                                                edge_pqueue &edges, Simplicial_Mesh &mesh, contraction_parameters &params)
{
//    cout<<e[0]<<" "<<e[1]<<endl;

    set<edge_weight> e_set; /// we insert the new edges first in this set to avoid duplicate insertions in the queue

//    cout<<"REMOVING "<<e[1]<<endl;

    for(unsigned d=0; d<difference.size(); d++)
    {
        for(ivect_iter it=difference[d].begin(); it!=difference[d].end(); ++it)
        {
            Top_Simplex &t = mesh.get_top_cell(d,*it);
//            cout<<*it<<"] old: "<<t<<endl;

            /// before updating the top simplex, we check
            /// if the leaf block indexing e[0] does not contain the current top d-simplex we have to add it
            /// NOTA: there is one possible case.. as leaf block n already indexes the top simplices in e[1]
            /// NOTA2: we have just to check that n and v_block are different (as if they are equal the edge is internal in n)
            if(n.get_v_start() != v_block.get_v_start() && !v_block.indexes_top(t))
            {
//                if(params.is_enable_debug_prints())
//                cout<<"[add top to leaf] "<<t<<" -> "<<v_block<<endl;
                v_block.add_top_cell(d,*it);
            }

            /// then we update the top simplex changing e[1] with e[0]
            int pos = t.vertex_index(e[1]);
            t.setTV(pos,e[0]);

//            cout<<e[1]<<" pos: "<<pos<<endl;
//            cout<<*it<<"] new: "<<t<<endl;

            /// we have to add the new edges in the queue
            ivect new_e; new_e.assign(2,0);
            for(int i=0; i<t.get_vertices_num(); i++)
            {
                if(i!=pos)
                {
                    t.TE(new_e,pos,i);
//                    cout<<"newe: "<<new_e[0]<<" "<<new_e[1]<<endl;

                    if(n.indexes_vertex(new_e[1])) /// we process an edge only if it has alle the extrema already processed
                    {
                        double error = mesh.get_vertex(e[0]).distance(mesh.get_vertex(e[1]));
                        if(error <= params.threshold)
//                        if(error >= params.threshold)
                        {
                            edge_weight ew_add(new_e,error);
                            e_set.insert(ew_add);
                        }
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

    /// finally we update the VTop relation of e[0]
    unify_container_of_containers(vt,difference);

//    /// DEBUGGING PRINT
//    print_verbose_coboundary_relation(vt,mesh);
//    int a; cin>>a;
}

void Contraction_Simplifier_Weak::extract_edges_lengths(Node_Stellar &n, Simplicial_Mesh &mesh, contraction_parameters &params)
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
//                        cout<<params.initialQuadric[e[0]-1][0]<<" "<<params.initialQuadric[e[1]-1][0]<<" : "<<params.threshold<<endl;
                        if(n.indexes_vertex(e[1])) /// we process an edge only if it has all the extrema already processed
                        {
                            auto r = e_set.insert(e);
                            if(r.second)
                            {
                                double error = mesh.get_vertex(e[0]).distance(mesh.get_vertex(e[1]));

                                if(error < params.min_err)
                                    params.min_err = error;
                                if(error > params.max_err)
                                    params.max_err = error;
                                params.avg_err += error;
                            }
                        }
                    }
                }
            }
        }
    }
    params.num_edges += e_set.size();
}

void Contraction_Simplifier_Weak::extract_edge_length_stats(Node_Stellar &n, Simplicial_Mesh &mesh, contraction_parameters &params)
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
//                        cout<<params.initialQuadric[e[0]-1][0]<<" "<<params.initialQuadric[e[1]-1][0]<<" : "<<params.threshold<<endl;
                        if(n.indexes_vertex(e[1])) /// we process an edge only if it has all the extrema already processed
                        {
                            auto r = e_set.insert(e);
                            if(r.second)
                            {
                                double error = mesh.get_vertex(e[0]).distance(mesh.get_vertex(e[1]));

                                if(error == params.min_err)
                                    params.min_c++;
                                if(error == params.max_err)
                                    params.max_c++;

                                if(params.min_err <= error && error < (params.min_err+params.avg_err)/2.0)
                                    params.histogram[0]++;
                                else if((params.min_err+params.avg_err)/2.0 <= error && error < params.avg_err)
                                    params.histogram[1]++;
                                else if(params.avg_err <= error && error < (params.avg_err+params.max_err)/2.0)
                                    params.histogram[2]++;
                                else
                                    params.histogram[3]++;
                            }
                        }
                    }
                }
            }
        }
    }
}
