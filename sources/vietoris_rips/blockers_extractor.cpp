#include "blockers_extractor.h"

blocker_set Blockers_Extractor::extract_blockers(Stellar_Tree &tree, Simplicial_Mesh &mesh, bool extract_from_maximal_cliques, ivect &orig_v_pos, string blockers_file_name)
{
//    cerr<<"[Blocker identification]"<<endl;
    Timer time;
//    pair<blocker_set,leaf_VT> p = make_pair(blocker_set(),leaf_VT());
    pair<int,int> fs_stats; fs_stats.first = 0; fs_stats.second = 0;
    tuple<blocker_set,leaf_VT,pair<int,int> > tupla = make_tuple(blocker_set(),leaf_VT(),fs_stats);

    double tot_time = 0;


//    leaf_VT all_vt;

    time.start();
    Blockers_Extractor::extract_global_VTops(get<1>(tupla),mesh);
    time.stop();
    time.print_elapsed_time("[P-TIME] Global VTop extraction: ");
    tot_time += time.get_elapsed_time();


    /// hierarchical visit - BUGGY -
//    time.start();
//    tree.visit(Blockers_Extractor::extract_leaf_blockers,tree.get_root(),mesh,tupla);
//    time.stop();
//    time.print_elapsed_time("[P-TIME] (local)MissingS and Blockers extraction: ");
//    tot_time += time.get_elapsed_time();
////    tree.parallel_visit(Blockers_Extractor::extract_leaf_blockers,mesh,p);
//    cerr<<"[TOT-TIME] "<<tot_time<<endl;
//    cerr<<"[STATS] (local-missing) max_d: "<<get<2>(tupla).first<<" max-entities: "<<get<2>(tupla).second<<endl;
//    Blockers_Extractor::get_blockers_stats(get<0>(tupla));
//    Blockers_Extractor::write_blockers_file(get<0>(tupla),blockers_file_name,orig_v_pos);

    if(extract_from_maximal_cliques)
    {
        /// === CLIQUES EXTRACTION FOR STATISTICAL PURPOSE === ///
        cerr << "BLOCKER EXTRACTION THROUGH MAXIMAL CLIQUES" << endl;
        top_simplices_map tops;
        //// reset the variables ////
        skeleton_graph graph = skeleton_graph(mesh.get_vertices_num());
        /// global visit
        time.start();
        Blockers_Extractor::extract_undirected_1skeleton(mesh,graph);
        time.stop();
        time.print_elapsed_time("undirected 1-skeleton extraction: ");
        tot_time += time.get_elapsed_time();
        time.start();
        Blockers_Extractor::get_maximal_cliques(graph,tops,mesh);
        time.stop();
        time.print_elapsed_time("maximal-cliques extraction: ");
        tot_time += time.get_elapsed_time();
        /// get cliques stats
        Blockers_Extractor::get_maximal_cliques_stats(tops);
        /// === === ///

        get<0>(tupla).clear();
        time.start();
        Blockers_Extractor::extract_blockers_from_maximal_cliques(tops,get<1>(tupla),mesh,get<0>(tupla));
        time.stop();
        time.print_elapsed_time("clique-blockers extraction: ");
        tot_time += time.get_elapsed_time();
        cerr<<"[TOT-TIME-CLIQUES] "<<tot_time<<endl;
//        Blockers_Extractor::get_blockers_stats(get<0>(tupla));
    }
    else
    {
        /// === failing simplices [global visit] === ///
        cerr << "BLOCKER EXTRACTION THROUGH MISSING SIMPLICES" << endl;
        skeleton_graph graph = skeleton_graph(mesh.get_vertices_num());
        time.start();
        tree.visit(skeleton_generator::extract_global_directed_1skeleton_Simplicial_wrapper,tree.get_root(),mesh,graph);
        time.stop();
        time.print_elapsed_time("directed 1-skeleton extraction: ");
        tot_time += time.get_elapsed_time();
        time.start();
        vector<set<ivect> > failing_simplices;
        failing_simplices.assign(mesh.get_top_cells_types(),set<ivect>());
        /// extract the failing simplices
        Blockers_Extractor::extract_failing_simplices(graph,mesh,get<1>(tupla),failing_simplices);
        time.stop();
        time.print_elapsed_time("[P-TIME] MissingS extraction: ");
        tot_time += time.get_elapsed_time();
        Blockers_Extractor::get_failing_simplices_stats(failing_simplices);

        /// BLOCKER EXTRACTIN BASED ON THE FAILING SIMPLICES
        time.start();
        Blockers_Extractor::extract_blockers(failing_simplices,get<1>(tupla),mesh,get<0>(tupla));
        time.stop();
        time.print_elapsed_time("[P-TIME] Blockers extraction: ");
        tot_time += time.get_elapsed_time();
        cerr<<"[TOT-TIME] "<<tot_time<<endl;
    //    /// FOR DEBUG ONLY - it checks if a blocker actually is a top simplex
    //    Blockers_Extractor::check_blockers(p.first,mesh);
        /// === === ///
    }

    Blockers_Extractor::get_blockers_stats(get<0>(tupla));
    Blockers_Extractor::write_blockers_file(get<0>(tupla),blockers_file_name,orig_v_pos);

    return get<0>(tupla);
}

void Blockers_Extractor::extract_aux_structures_statistics(Stellar_Tree &tree, Simplicial_Mesh &mesh)
{
        Timer time;
        pair<int,int> fs_stats; fs_stats.first = 0; fs_stats.second = 0;
        tuple<blocker_set,leaf_VT,pair<int,int> > tupla = make_tuple(blocker_set(),leaf_VT(),fs_stats);

        double tot_time = 0;

        time.start();
        Blockers_Extractor::extract_global_VTops(get<1>(tupla),mesh);
        time.stop();
        time.print_elapsed_time("[P-TIME] Global VTop extraction: ");
        tot_time += time.get_elapsed_time();

        /// === CLIQUES EXTRACTION FOR STATISTICAL PURPOSE === ///
        top_simplices_map tops;
        //// reset the variables ////
        skeleton_graph graph = skeleton_graph(mesh.get_vertices_num());
        /// global visit
        time.start();
        Blockers_Extractor::extract_undirected_1skeleton(mesh,graph);
        time.stop();
        time.print_elapsed_time("[P-TIME] undirected 1-skeleton extraction: ");
        tot_time += time.get_elapsed_time();
        time.start();
        Blockers_Extractor::get_maximal_cliques(graph,tops,mesh);
        time.stop();
        time.print_elapsed_time("[P-TIME] maximal-cliques extraction: ");
        tot_time += time.get_elapsed_time();
        /// get cliques stats
        Blockers_Extractor::get_maximal_cliques_stats(tops);
        /// === === ///
        cerr<<"[TOT-TIME-CLIQUES] "<<tot_time<<endl;
        tot_time = 0;
        /// === failing simplices [global visit] === ///
        cerr << "BLOCKER EXTRACTION THROUGH MISSING SIMPLICES" << endl;
        graph.reset();
        graph = skeleton_graph(mesh.get_vertices_num());
        time.start();
        tree.visit(skeleton_generator::extract_global_directed_1skeleton_Simplicial_wrapper,tree.get_root(),mesh,graph);
        time.stop();
        time.print_elapsed_time("[P-TIME] directed 1-skeleton extraction: ");
        tot_time += time.get_elapsed_time();
        time.start();
        vector<set<ivect> > failing_simplices;
        failing_simplices.assign(mesh.get_top_cells_types(),set<ivect>());
        /// extract the failing simplices
        Blockers_Extractor::extract_failing_simplices(graph,mesh,get<1>(tupla),failing_simplices);
        time.stop();
        time.print_elapsed_time("[P-TIME] MissingS extraction: ");
        tot_time += time.get_elapsed_time();
        Blockers_Extractor::get_failing_simplices_stats(failing_simplices);
        cerr<<"[TOT-TIME-MISSING] "<<tot_time<<endl;
}

void Blockers_Extractor::add_blockers(Stellar_Tree &tree, Simplicial_Mesh &mesh, vector<blocker_array> &b_arr)
{
    /// (1) we add the blockers as top simplices in the mesh
    vector<ivect> blocker_to_id;
    blocker_to_id.assign(b_arr.size(),ivect());
    for(unsigned d=0; d<b_arr.size(); d++)
    {
        blocker_to_id[d].assign(b_arr[d].size(),-1);
        for(unsigned i=0; i<b_arr[d].size(); i++)
        {
            int b_top_id = mesh.add_top_cell(d,b_arr[d][i]);
            blocker_to_id[d][i] = b_top_id;
        }
    }

    /// (2) then we add them in the hierarchy
    tree.visit(Blockers_Extractor::add_blockers_to_leaf,tree.get_root(),mesh,blocker_to_id);
}

void Blockers_Extractor::extract_blockers_from_maximal_cliques(top_simplices_map &tops, leaf_VT &all_vt, Simplicial_Mesh &mesh, blocker_set &bs)
{
    int current_d = 0;
    for(top_simplices_map::iterator it=tops.begin(); it!=tops.end(); ++it)
    {
//        list<int_vect> &tops_d = it->second;
        vector<ivect> tops_d(it->second.begin(),it->second.end());
#pragma omp parallel for
//        for(list<int_vect>::iterator itt=tops_d.begin(); itt!=tops_d.end(); ++itt)
        for(unsigned i=0; i<tops_d.size(); i++)
        {
//            Top_Simplex blocker_candidate(*itt);
            Top_Simplex blocker_candidate(tops_d[i]);
            bool is_top = false;

            current_d = blocker_candidate.get_vertices_num() - 2;

//            cout<<blocker_candidate<<flush;

            /// we check for a top simplex only if it could be represented
            if(current_d < mesh.get_top_cells_types())
            {
                int first_v = blocker_candidate.TV(0);
                VT &vt = all_vt[first_v-1];
                for(auto it=vt[current_d].begin(); it!=vt[current_d].end(); ++it)
                {
                    Top_Simplex &top = mesh.get_top_cell(current_d,*it);
                    if(top == blocker_candidate)
                    {
                        is_top = true;
                        break;
                    }
                }
            }

//            cout<<" A "<<flush;

            if(!is_top && blocker_candidate.get_vertices_num() > 2) /// an edge cannot be a blocker
            {
                blocker_set::iterator iter;
                #pragma omp critical
                {
                    iter = bs.find(/**itt*/tops_d[i]);
                }
                /// we have not inserted it yet
                if(iter==bs.end())
                {
//                    if(Blockers_Extractor::is_blocker(blocker_candidate,current_d,mesh,bs,p.second))
                    if(Blockers_Extractor::is_blocker(blocker_candidate,current_d,all_vt,mesh,bs))
                    {
                        #pragma omp critical
                        {
                            bs.insert(/**itt*/tops_d[i]);
                        }
                    }
                }

//                cout<<"A "<<flush;

            }

//            cout<<"A "<<endl;
        }
//        current_d++;
    }
}

void Blockers_Extractor::extract_blockers(failing_array &failing_simplices, leaf_VT &all_vt, Simplicial_Mesh &mesh, blocker_set &p)
{
//    Timer time;
//    time.start();

//#pragma omp parallel for
    for(unsigned d=0; d<failing_simplices.size(); d++)
    {
//        set<int_vect> &d_set = failing_simplices[d];
        vector<ivect> d_vect(failing_simplices[d].begin(),failing_simplices[d].end());

#pragma omp parallel for
//        for(set<int_vect>::iterator it=d_set.begin(); it!=d_set.end(); ++it)
        for(unsigned i=0; i<d_vect.size(); i++)
        {
            Top_Simplex b = Top_Simplex(/**it*/d_vect[i]);

            if(d==1) /// the target triangles are by definition blocker (as a target simplex is not in Sigma)
#pragma omp critical
            {
                p.insert(b);
            }
            else if(Blockers_Extractor::is_blocker(b,b.get_vertices_num()-2,all_vt,mesh,failing_simplices,p))
            {
//                cout<<" is blocker"<<endl;
#pragma omp critical
                {
                    p.insert(b);
                }
            }
//            int a; cin>>a;
        }
    }
//    time.stop();
//    time.print_elapsed_time("[P-TIME] blocker extraction: ");

//    return time.get_elapsed_time();
}

void Blockers_Extractor::get_blockers_stats(blocker_set &blockers)
{
    cout<<"[STATS] Blockers number: "<<blockers.size()<<endl;
    int max_blocker = 0;
    for(auto it=blockers.begin(); it!=blockers.end(); ++it)
    {
        if(max_blocker < it->get_vertices_num())
            max_blocker = it->get_vertices_num();
    }
    cout<<"max_blocker dim: "<<max_blocker-1<<endl;
    ivect bs; bs.assign(max_blocker-1,0);
    for(auto it=blockers.begin(); it!=blockers.end(); ++it)
         bs[it->get_vertices_num()-2]++;
    for(unsigned i=0; i<bs.size(); i++)
        cout<<i+1<<"-blockers: "<<bs[i]<<endl;
}

void Blockers_Extractor::write_blockers_file(blocker_set &blockers, string file_name, ivect &orig_v_pos)
{
    stringstream ss;
    ss << file_name << ".blockers";
    ofstream output(ss.str());

    output << blockers.size() << endl;
    for(auto it=blockers.begin(); it!=blockers.end(); ++it)
    {
        const Top_Simplex &top = *it;
        for(int v=0; v<top.get_vertices_num(); v++)
        {
            output << orig_v_pos[abs(top.TV(v))-1] << " ";
        }
        output << endl;
    }
    output.close();
}

void Blockers_Extractor::extract_leaf_blockers(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<blocker_set,leaf_VT,pair<int,int> > &tupla)
{
//    cout<<n<<endl;
//    n.print_top_cells_arrays(mesh);
//    leaf_VT vtops;
//    n.extract_local_VTop(mesh,vtops);

    vector<set<ivect> > failing_simplices;
    failing_simplices.assign(mesh.get_top_cells_types(),set<ivect>());

    skeleton_graph sg(n.get_v_end()-n.get_v_start());
    skeleton_generator::extract_local_directed_1skeleton(n,mesh,sg);

    Blockers_Extractor::extract_failing_simplices(sg,mesh,get<1>(tupla),failing_simplices);

    Blockers_Extractor::extract_blockers(failing_simplices,get<1>(tupla),mesh,get<0>(tupla));

    ///get failing simplices stats
    Blockers_Extractor::get_failing_simplices_stats(failing_simplices,get<2>(tupla));

//    unordered_map<int,int_set> ov;
//    top_simplices_map tops;

//    cout<<"A"<<endl;
//    Blockers_Extractor::extract_undirected_1skeleton(n,mesh,sg,ov);
//    cout<<"B"<<endl;
//    Blockers_Extractor::get_local_maximal_cliques(sg,ov,tops,n);
//    cout<<"C"<<endl;

//    int current_d = 0;
//    for(top_simplices_map::iterator it=tops.begin(); it!=tops.end(); ++it)
//    {
//        list<int_vect> &tops_d = it->second;
//        for(list<int_vect>::iterator itt=tops_d.begin(); itt!=tops_d.end(); ++itt)
//        {
//            Top_Simplex blocker_candidate(*itt);
//            current_d = blocker_candidate.get_vertices_num();
//            bool is_top = false;
//            for(RunIteratorPair itPair = n.make_t_array_iterator_pair(current_d); itPair.first != itPair.second; ++itPair.first)
//            {
//                Top_Simplex &top = mesh.get_top_cell(current_d,*itPair.first);
//                if(top == blocker_candidate)
//                {
//                    if(p.second)
//                        cout<<"===[TOP]=== "<<blocker_candidate<<endl;
//                    is_top = true;
//                    break;
//                }
//            }

//            if(!is_top && blocker_candidate.get_vertices_num() > 2) /// an edge cannot be a blocker
//            {
//                if(p.second)
//                    cout<<"possible blocker: "<<blocker_candidate<<endl;

//                blocker_set::iterator iter;
//                #pragma omp critical
//                {
//                    iter = p.first.find(*itt);
//                }
//                /// we have not inserted it yet
//                if(iter==p.first.end())
//                {
////                    if(p.second)
////                        cout<<"possible blocker: "<<blocker_candidate<<endl;

//                    if(Blockers_Extractor::is_blocker(blocker_candidate,current_d,n,mesh,p.first))
//                    {
//                        if(p.second)
//                            cout<<"====="<<blocker_candidate<<" is a blocker"<<endl;

//                        #pragma omp critical
//                        {
//                            p.first.insert(*itt);
//                        }
//                    }
//                    //                int a; cin>>a;
//                }
//                else
//                    cout << "     already into the blocker set"<<endl;

//            }

//        }
//        current_d++;
//    }

//    int a; cin>>a;
}

void Blockers_Extractor::add_blockers_to_leaf(Node_Stellar &n, Simplicial_Mesh &mesh, vector<ivect> &blockers_id)
{
    for(unsigned d=0; d<blockers_id.size(); d++)
    {
        for(unsigned i=0; i<blockers_id[d].size(); i++)
        {
            Top_Simplex &blocker_top = mesh.get_top_cell(d,blockers_id[d][i]);
            if(n.indexes_top(blocker_top))
            {
                /// (1) if a leaf block indexes a blockers.. then we push it into the leaf block
                n.add_top_cell(d,blockers_id[d][i]);

                /// (2) we check if some of the top simplex indexed by n are d-1 faces of the blocker
                ///     -> we flag for removal..
                ivect f;
                for(int v=0; v<blocker_top.get_vertices_num(); v++)
                {
                    blocker_top.TF(f,v);
                    if(n.indexes_cell(f))
                    {
                        int target_d = f.size() -2;
                        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(target_d); itPair.first != itPair.second; ++itPair.first)
                        {
                            RunIterator const& t_id = itPair.first;

                            if(!mesh.is_top_cell_removed(target_d,*t_id))
                            {
                                Top_Simplex& t = mesh.get_top_cell(target_d,*t_id);
                                if(t.are_equal(f))
                                {
                                    mesh.remove_top_cell(target_d,*t_id);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void Blockers_Extractor::extract_undirected_1skeleton(Node_Stellar &n, Simplicial_Mesh &mesh, skeleton_graph &graph, unordered_map<int, iset> &ov)
{
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                Top_Simplex& t = mesh.get_top_cell(d,*t_id);

                for(int v=0; v<t.get_vertices_num(); v++)
                {
                    int v_id = abs(t.TV(v));

                    for(int v2=v+1; v2<t.get_vertices_num(); v2++)
                    {
                        int v2_id = abs(t.TV(v2));

                        if(n.indexes_vertex(v_id))
                            graph.insert(v_id-n.get_v_start(),v2_id);
                        else
                            ov[v_id].insert(v2_id);

                        if(n.indexes_vertex(v2_id))
                            graph.insert(v2_id-n.get_v_start(),v_id);
                        else
                            ov[v2_id].insert(v_id);
                    }
                }
            }
        }
    }
}

void Blockers_Extractor::extract_undirected_1skeleton(Simplicial_Mesh &mesh, skeleton_graph &graph)
{
    for(int d=0; d<mesh.get_top_cells_types();d++)
    {
        for(int i=1; i<=mesh.get_top_cells_num(d); i++)
        {
            if(!mesh.is_top_cell_removed(d,i))
            {
                Top_Simplex& t = mesh.get_top_cell(d,i);

                for(int v=0; v<t.get_vertices_num(); v++)
                {
                    int v_id = abs(t.TV(v));

                    for(int v2=v+1; v2<t.get_vertices_num(); v2++)
                    {
                        int v2_id = abs(t.TV(v2));

                        graph.insert(v_id-1,v2_id);
                        graph.insert(v2_id-1,v_id);
                    }
                }
            }
        }
    }
}

void Blockers_Extractor::get_local_maximal_cliques(skeleton_graph &graph, unordered_map<int, iset> &outer_vertices, top_simplices_map &tops, Node_Stellar &n)
{
    iset points;
    for(int v_id=n.get_v_start(); v_id<n.get_v_end(); v_id++)
        points.insert(v_id);
    for(unordered_map<int, iset>::iterator it=outer_vertices.begin(); it!=outer_vertices.end(); ++it)
        points.insert(it->first);

    local_generation_parameters p(0,NULL,false); /// dummy variable

//    top_simplices_map tops;
    iset setR, setX;
    VietorisRips_Generator::get_local_maximal_cliques_rec(graph,outer_vertices,setR, points, setX, tops, n, n.get_v_start(),p);
}

void Blockers_Extractor::get_maximal_cliques(skeleton_graph &graph, top_simplices_map &tops, Simplicial_Mesh &mesh)
{
    iset points;
    for(int i=1; i<=mesh.get_vertices_num(); i++)
    {
        points.insert(i);
    }

    iset setR, setX;
    unsigned maxX = 0;
    VietorisRips_Generator::get_maximal_cliques(graph,1,setR, points, setX, tops, false, maxX);
}

void Blockers_Extractor::get_maximal_cliques_stats(top_simplices_map &tops)
{
    int max_clique = 0;
    for(top_simplices_map::iterator it=tops.begin(); it!=tops.end(); ++it)
    {
        const list<ivect> &clist = it->second;
        for(list<ivect>::const_iterator cit=clist.begin(); cit!=clist.end(); ++cit)
        {
            if((int)cit->size() > max_clique)
                max_clique = cit->size();
        }
    }
    cout<<"[STATS] max clique DIM: "<<max_clique-1<<endl;
    ivect csum; csum.assign(max_clique-1,0);
    for(top_simplices_map::iterator it=tops.begin(); it!=tops.end(); ++it)
    {
        const list<ivect> &clist = it->second;
        for(list<ivect>::const_iterator cit=clist.begin(); cit!=clist.end(); ++cit)
        {
            csum[cit->size()-2]++;
        }
    }
    int tot=0;
    for(unsigned i=0; i<csum.size(); i++)
    {
        cout<<i+1<<"-cliques: "<<csum[i]<<endl;
        tot += csum[i];
    }
    cout<<"TOT "<<tot<<endl;
//    int a; cin>>a;
}

bool Blockers_Extractor::is_blocker(Top_Simplex &blocker_candidate, int current_d, Node_Stellar &n, Simplicial_Mesh &mesh, blocker_set &bs)
{
    ivect f;
    int num_represented_d_1faces = 0;
    for(int v=0; v<blocker_candidate.get_vertices_num(); v++)
    {
        blocker_candidate.TF(f,v);
        /// for debug only
        if(!n.indexes_cell(f))
        {
            cout<<n<<" does not index ";
            print_container_content(f);
            cout<<endl;
        }
//        cout<<current_d<<endl;
        if(Blockers_Extractor::is_in_Sigma(f,f.size()-2,n,mesh))
        {
            num_represented_d_1faces++;
        }
        else
        {
            if(f.size() > 2) /// we cannot have an edge as a blocker
            {
                Top_Simplex bc2(f);
                blocker_set::iterator iter;
                cout<<" sub "<<bc2<<" candidate blocker"<<endl;
                #pragma omp critical
                {
                    iter = bs.find(f);
                }
                if(iter==bs.end())
                {
                    if(Blockers_Extractor::is_blocker(bc2,current_d-1,n,mesh,bs)) /// NOTA: the current dimension is reduced by 1
                    {
                        cout<<" --> "<<bc2<<" is a blocker"<<endl;
    //                    int a; cin>>a;
                        #pragma omp critical
                        {
                            bs.insert(f);
                        }
                    }
                }
                else
                    cout << "     already into the blocker set"<<endl;
            }
        }
    }

    return (num_represented_d_1faces == blocker_candidate.get_vertices_num());
}

bool Blockers_Extractor::is_in_Sigma(ivect &f, int current_d, Node_Stellar &n, Simplicial_Mesh &mesh)
{
    for(int d=current_d; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                Top_Simplex& t = mesh.get_top_cell(d,*t_id);
                if(t.has_cell(f))
                {
                    cout<<"   ";
                    print_container_content(f);
                    cout<<" is in the boundary of "<<t<<endl;
//                    int a; cin>>a;
                    return true;
                }
            }
        }
    }

    return false;
}

bool Blockers_Extractor::is_blocker(Top_Simplex &blocker_candidate, int current_d, leaf_VT &all_vt, Simplicial_Mesh &mesh, blocker_set &bs)
{
    ivect f;
    int num_represented_d_1faces = 0;
    for(int v=0; v<blocker_candidate.get_vertices_num(); v++)
    {
        blocker_candidate.TF(f,v);

        if(Blockers_Extractor::is_in_Sigma(f,f.size()-2,all_vt,mesh,false))
           num_represented_d_1faces++;
        else
        {
            if(f.size() > 2) /// we cannot have an edge as a blocker
            {
                Top_Simplex bc2(f);
                blocker_set::iterator iter;

                iter = bs.find(f);

                if(iter==bs.end())
                {
                    if(Blockers_Extractor::is_blocker(bc2,current_d-1,all_vt,mesh,bs)) /// NOTA: the current dimension is reduced by 1
                    {
                        #pragma omp critical
                        {
                            bs.insert(f);
                        }
                    }
                }
            }
        }
    }

    return (num_represented_d_1faces == blocker_candidate.get_vertices_num());
}

bool Blockers_Extractor::is_blocker(Top_Simplex &blocker_candidate, int current_d, leaf_VT &all_vt, Simplicial_Mesh &mesh, failing_array &failing_simplices, blocker_set &bs/*, bool debug*/)
{
    ivect f;
    int num_represented_d_1faces = 0;
    for(int v=0; v<blocker_candidate.get_vertices_num(); v++)
    {
        blocker_candidate.TF(f,v);

        set<ivect>::iterator fiter = failing_simplices[f.size()-2].find(f);
        if(fiter!=failing_simplices[f.size()-2].end())
            continue;

        blocker_set::iterator biter = bs.find(f);
        if(biter!=bs.end())
            continue;

//        if(debug)
//            print_container_content("current d-1 face: ",f);
        if(Blockers_Extractor::is_in_Sigma(f,f.size()-2,all_vt,mesh,false))
           num_represented_d_1faces++;
        else
        {
            if(f.size() > 2) /// we cannot have an edge as a blocker
            {
                Top_Simplex bc2(f);
                if(Blockers_Extractor::is_blocker(bc2,current_d-1,all_vt,mesh,failing_simplices,bs)) /// NOTA: the current dimension is reduced by 1
                {
//                    if(debug)
//                        cout<<" =====> "<<bc2<<" is a blocker"<<endl;

#pragma omp critical
                    {
                        bs.insert(f);
                    }
                }
            }
        }
    }

    return (num_represented_d_1faces == blocker_candidate.get_vertices_num());
}

//bool Blockers_Extractor::is_in_Sigma(int_vect &f, int current_d, Simplicial_Mesh &mesh, bool debug)
//{
//    for(int d=current_d; d<mesh.get_top_cells_types(); d++)
//    {
//        for(int i=1; i<=mesh.get_top_cells_num(d); i++)
//        {
//            if(!mesh.is_top_cell_removed(d,i))
//            {
//                Top_Simplex& t = mesh.get_top_cell(d,i);

//                if(t.has_cell(f))
//                {
//                    if(debug)
//                        cout<<" ---> is in the boundary of "<<t<<endl;
//                    return true;
//                }
//            }
//        }
//    }

//    return false;
//}

bool Blockers_Extractor::is_in_Sigma(ivect &f, int current_d, leaf_VT &all_vt, Simplicial_Mesh &mesh, bool debug)
{
    for(auto v=f.begin(); v!=f.end(); ++v)
    {
        VT &vt = all_vt[*v-1];
        for(unsigned d=current_d; d<vt.size(); d++)
        {
            for(unsigned i=0; i<vt[d].size(); i++)
            {
                if(!mesh.is_top_cell_removed(d,vt[d][i]))
                {
                    Top_Simplex& t = mesh.get_top_cell(d,vt[d][i]);
                    if(t.has_cell(f))
                    {
                        if(debug)
                            cout<<" ---> is in the boundary of "<<t<<endl;
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

void Blockers_Extractor::check_blockers(blocker_set &bs, Simplicial_Mesh &mesh)
{
    for(blocker_set::iterator it=bs.begin(); it!=bs.end(); ++it)
    {
        Top_Simplex blocker = *it;
        int current_d = blocker.get_vertices_num() -2;

        /// first we check if a blocker is a top simplex
        for(int i=1; i<=mesh.get_top_cells_num(current_d); i++)
        {
            if(!mesh.is_top_cell_removed(current_d,i))
            {
                Top_Simplex& t = mesh.get_top_cell(current_d,i);
                if(t == blocker)
                {
                    cerr<<"the blocker is a top simplex"<<endl;
                    cerr<<blocker<<" -- "<<t<<endl;
                }
            }
        }
    }
}

void Blockers_Extractor::extract_global_VTops(leaf_VT &all_vt, Simplicial_Mesh &mesh)
{
    VT tmp; tmp.assign(mesh.get_top_cells_types(),ivect());
    all_vt.assign(mesh.get_vertices_num(),tmp);

    for(int d=0; d<mesh.get_top_cells_types(); d++)
    {
        for(int i=1; i<=mesh.get_top_cells_num(d); i++)
        {
            if(!mesh.is_top_cell_removed(d,i))
            {
                Top_Simplex &top = mesh.get_top_cell(d,i);
                for(int v=0; v<top.get_vertices_num(); v++)
                {
                    int v_pos = abs(top.TV(v)) -1;
                    all_vt[v_pos][d].push_back(i);
                }
            }
        }
    }
}

void Blockers_Extractor::extract_failing_simplices(skeleton_graph &graph, Simplicial_Mesh &mesh, leaf_VT &all_vt, vector<set<ivect> > &failing_simplices)
{
//    Timer time;

//    double t = 0;

//    time.start();
//    Blockers_Extractor::extract_global_VTops(all_vt,mesh);
//    time.stop();
//    time.print_elapsed_time("Global VTop extraction: ");
//    t += time.get_elapsed_time();

//    time.start();
#pragma omp parallel for
    for(int v=0; v<graph.graph_size(); v++)
    {
        ivect e; e.assign(2,-1);
        iset &adjs = graph.get_graph(v);
        int real_id = v+1;
        e[0] = real_id;

        for(iset_iter it=adjs.begin(); it!=adjs.end(); ++it)
        {
            e[1] = *it;
            VT &vt0 = all_vt[real_id-1];
            VT &vt1 = all_vt[*it-1];
            ET et; et.assign(vt0.size(),ivect());
            intersect_container_of_containers(vt0,vt1,et,true);

            Blockers_Extractor::extract_failing_simplices(e,et,vt0,vt1,mesh,failing_simplices);
        }
    }
//    time.stop();
//    time.print_elapsed_time("Failing-Simplices extraction: ");
//    t += time.get_elapsed_time();

//    return t;
}

//void Blockers_Extractor::extract_failing_simplices(Node_Stellar &n, skeleton_graph &graph, Simplicial_Mesh &mesh, leaf_VT &all_vt, vector<set<int_vect> > &failing_simplices)
//{
//    int_vect e; e.assign(2,-1);

//    for(int v=0; v<graph.graph_size(); v++)
//    {
//        int_set &adjs = graph.get_graph(v);
//        int real_id = v+1;
//        e[0] = real_id;

//        for(int_set_iter it=adjs.begin(); it!=adjs.end(); ++it)
//        {
//            e[1] = *it;
//            VT &vt0 = all_vt[real_id-1];
//            VT &vt1 = all_vt[*it-1];
//            ET et; et.assign(vt0.size(),int_vect());
//            intersect_container_of_containers(vt0,vt1,et,true);

//            Blockers_Extractor::extract_failing_simplices(e,et,vt0,vt1,mesh,failing_simplices);
//        }
//    }
//}

void Blockers_Extractor::extract_failing_simplices(ivect &e, ET &et, VT &vt0, VT &vt1, Simplicial_Mesh &mesh, vector<set<ivect> > &failing_simplices)
{
    iset v_in_link0, v_in_link1, v_in_linke;

    /// we have to extract just the vertices in the link of e and its extrema
    Contraction_Simplifier_Top::get_vertices_in_link(e[0],vt0,v_in_link0,mesh);
    Contraction_Simplifier_Top::get_vertices_in_link(e[1],vt1,v_in_link1,mesh);
    Contraction_Simplifier_Top::get_vertices_in_link(e,et,v_in_linke,mesh);

    /// (1) we intersect the vertices links
    iset v_intersection = v_in_link0;
    intersect_sets(v_intersection,v_in_link1);

    set<ivect> shared_v0, shared_v1;
    Contraction_Simplifier_Top::get_shared_simplices(vt0,v_intersection,0,mesh,shared_v0);
    Contraction_Simplifier_Top::get_shared_simplices(vt1,v_intersection,0,mesh,shared_v1);

    for(auto it0=shared_v0.begin(); it0!=shared_v0.end(); ++it0)
    {
        for(auto it1=shared_v1.begin(); it1!=shared_v1.end(); ++it1)
        {
            ivect intersect = *it0;
            intersect_containers(intersect,*it1);
            if(intersect.size() > 0)
            {
                if(!Contraction_Simplifier_Weak::exists_simplex_in_coboundary(intersect.size()-1,intersect,et,mesh))
                {
                    /// we get the full-fail-simplex
                    intersect.push_back(e[0]);
                    intersect.push_back(e[1]);
                    sort(intersect.begin(),intersect.end());
#pragma omp critical
                    {
                        failing_simplices[intersect.size()-2].insert(intersect);
                    }
                }
            }
        }
    }
}

void Blockers_Extractor::get_failing_simplices_stats(vector<set<ivect> > &failing_simplices)
{
    int max_fail_simplex_dim = 0;
    ivect fsum; fsum.assign(failing_simplices.size(),0);
    for(unsigned d=0; d<failing_simplices.size(); d++)
    {
        const set<ivect> &d_set = failing_simplices[d];
        if(d_set.size() > 0)
            max_fail_simplex_dim = d+1;
        fsum[d] = d_set.size();
    }
    int tot=0;
    cout<<"[STATS] max failing simplex dim: "<<max_fail_simplex_dim<<endl;
    for(unsigned i=0; i<fsum.size(); i++)
    {
        cout<<i+1<<"-failing-simplices: "<<fsum[i]<<endl;
        tot += fsum[i];
    }
    cout<<"TOT "<<tot<<endl;
}

void Blockers_Extractor::get_failing_simplices_stats(vector<set<ivect> > &failing_simplices, pair<int,int> &s)
{
    int max_fail_simplex_dim = 0;
    int tot=0;
//    int_vect fsum; fsum.assign(failing_simplices.size(),0);
    for(unsigned d=0; d<failing_simplices.size(); d++)
    {
        const set<ivect> &d_set = failing_simplices[d];
        if(d_set.size() > 0)
            max_fail_simplex_dim = d+1;
        /*fsum[d]*/tot += d_set.size();
    }

    if(max_fail_simplex_dim > s.first)
        s.first = max_fail_simplex_dim;
    if(tot > s.second)
        s.second = tot;

//    cerr<<"max failing simplex size: "<<max_fail_simplex_dim<<endl;
//    for(unsigned i=0; i<fsum.size(); i++)
//    {
//        cerr<<i+1<<"-failing-simplices: "<<fsum[i]<<endl;
//        tot += fsum[i];
//    }
//    cerr<<"TOT "<<tot<<endl;
}
