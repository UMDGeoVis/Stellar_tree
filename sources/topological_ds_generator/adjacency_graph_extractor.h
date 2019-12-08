#ifndef ADJACENCY_GRAPH_EXTRACTOR_H
#define ADJACENCY_GRAPH_EXTRACTOR_H

#include "utilities/container_utilities.h"
#include "stellar_decomposition/run_iterator.h"
#include "stellar_decomposition/node_stellar.h"

#include "topological_ds/adjacency_graph.h"

class Adjacency_Graph_Extractor
{
 public:
    Adjacency_Graph_Extractor() {}

    static void extract_lite_adjacency_graph_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, Adjacency_Graph& graph);
    static void extract_lite_adjacency_graph_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, Adjacency_Graph& graph);

    static void extract_full_adjacency_graph_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, Adjacency_Graph& graph);
    static void extract_full_adjacency_graph_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, Adjacency_Graph& graph);

    static void extract_local_lite_adjacency_graph_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, Local_Adjacency_Graph& graph);
    static void extract_local_lite_adjacency_graph_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, Local_Adjacency_Graph& graph);

    static void extract_local_full_adjacency_graph_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, Local_Adjacency_Graph& graph);
    static void extract_local_full_adjacency_graph_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, Local_Adjacency_Graph& graph);

    static void extract_local_1skeleton_graph_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, string mesh_name);
    static void extract_local_1skeleton_graph_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, string mesh_name);

private:
    template<class C, class T> static void extract_adjacency_graph(Node_Stellar &n, Mesh<C,T> &mesh, Adjacency_Graph& graph);
    template<class C, class T> static void extract_adjacency_graph_lite(Node_Stellar &n, Mesh<C,T> &mesh, Adjacency_Graph& graph);
    template<class C, class T> static void extract_adjacency_graph_blockers_lite(Node_Stellar &n, Mesh<C,T> &mesh, Adjacency_Graph& graph);
    template<class C, class T> static void extract_local_adjacency_graph(Node_Stellar &n, Mesh<C,T> &mesh, Local_Adjacency_Graph &graph);
    template<class C, class T> static void extract_local_adjacency_graph_lite(Node_Stellar &n, Mesh<C,T> &mesh, Local_Adjacency_Graph& graph);
    template<class C, class T> static void extract_local_1skeleton_graph(Node_Stellar &n, Mesh<C,T> &mesh, string mesh_name);

    template<class C, class T> static void extract_d1_d2_adjacents(T &t, int t_id, int d, bool is_blocker, Node_Stellar &n, Mesh<C,T> &mesh, Adjacency_Graph &graph);
//    template<class C, class T> static void extract_d1_d2_adjacent_blockers(T &t, int t_id, int d, Adjacency_Graph &graph);
};

template<class C, class T> void Adjacency_Graph_Extractor::extract_adjacency_graph(Node_Stellar &n, Mesh<C,T> &mesh, Adjacency_Graph &graph)
{
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIterator runIt = n.t_array_begin_iterator(d), runEnd = n.t_array_end_iterator(d); runIt != runEnd; ++runIt)
//        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
//            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*runIt))
            {
                T& top = mesh.get_top_cell(d,*runIt);

                ivect top_v; top.get_positive_sorted_vertices(top_v);
                int top_id = graph.get_global_id(d,*runIt);/*global_ids[i][j-1];*/
                ivect top2_v;

                // loop on all the successive top simplices
                for(unsigned d2=d; d2<n.get_num_top_cells_encoded(); d2++)
                {
//                    for(unsigned j2=j+1; j2<=tops[d2].size(); j2++)
                    for(RunIterator runIt2 = n.t_array_begin_iterator(d2), runEnd2 = n.t_array_end_iterator(d2); runIt2 != runEnd2; ++runIt2)
                    {
                        if(d!=d2 || (d==d2 && runIt!=runIt2)) /// we are not processing a top-d-simplex against itself
                        {
                            T& top2 = mesh.get_top_cell(d2,*runIt2);
                            int top2_id = graph.get_global_id(d2,*runIt2);/*global_ids[d2][j2-1];*/
                            top2.get_positive_sorted_vertices(top2_v);
                            intersect_containers(top2_v,top_v);

                            if(top2_v.size() > 1) // we skip those top simplices that are incident just in a vertex
                            {
                                // we just push each edge once
                                pair<int,short> p;
                                p.first = top2_id;
                                p.second = top2_v.size();
                                graph.add_pair(top_id,p);
//                                adj_matrix[top_id].push_back(p);
    //                                cout << top_id << " -- " << top2_id << " color = " << palette[p.second].a << " " << palette[p.second].b << " " << palette[p.second].c << endl;
    //                                adj_matrix[top_id][top2_id] = top2_v.size()-1;
    //                                adj_matrix[top2_id][top_id] = top2_v.size()-1;
                            }
    //                            if(top2_v.size() == top_v.size()-1) // d-1 adjacency
    //                            {
    //                                d_1_adjs[top_id].push_back(top2_id);
    //                            }
                        }
                    }
                }
            }
        }
    }
}

template<class C, class T> void Adjacency_Graph_Extractor::extract_adjacency_graph_lite(Node_Stellar &n, Mesh<C,T> &mesh, Adjacency_Graph &graph)
{
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIterator runIt = n.t_array_begin_iterator(d), runEnd = n.t_array_end_iterator(d); runIt != runEnd; ++runIt)
        {
            if(!mesh.is_top_cell_removed(d,*runIt))
            {
                T& top = mesh.get_top_cell(d,*runIt);
                Adjacency_Graph_Extractor::extract_d1_d2_adjacents(top,*runIt,d,false,n,mesh,graph);

//                int_vect top_v; top.get_positive_sorted_vertices(top_v);
//                int top_id = graph.get_global_id(d,*runIt);
//                int_vect top2_v;

//                /// first check d-1 adjacencies
//                for(RunIterator runIt2 = n.t_array_begin_iterator(d), runEnd2 = n.t_array_end_iterator(d); runIt2 != runEnd2; ++runIt2)
//                {
//                    if(runIt != runIt2) ///
//                    {
//                        T& top2 = mesh.get_top_cell(d,*runIt2);
//                        int top2_id = graph.get_global_id(d,*runIt2);
//                        top2.get_positive_sorted_vertices(top2_v);
//                        intersect_containers(top2_v,top_v);

//                        if(top2_v.size() == top.get_vertices_num()-1 ) /// we have a d-1 adjacency
//                        {
//                            // we just push each edge once
//                            pair<int,short> p;
//                            p.first = top2_id;
//                            p.second = top2_v.size();
//                            graph.add_pair(top_id,p);
//                        }
//                    }
//                }

//                if(d>0) /// at least top 2-simplices
//                {
//                    int d2 = d-1;
//                    for(RunIterator runIt2 = n.t_array_begin_iterator(d2), runEnd2 = n.t_array_end_iterator(d2); runIt2 != runEnd2; ++runIt2)
//                    {
//                        T& top2 = mesh.get_top_cell(d2,*runIt2);
//                        int top2_id = graph.get_global_id(d2,*runIt2);
//                        top2.get_positive_sorted_vertices(top2_v);
//                        intersect_containers(top2_v,top_v);

//                        if(top2_v.size() == top2.get_vertices_num()-1 ) /// we have a d-2 adjacency
//                        {
//                            // we just push each edge once
//                            pair<int,short> p;
//                            p.first = top2_id;
//                            p.second = top2_v.size();
//                            graph.add_pair(top_id,p);
//                        }
//                    }
//                }
            }
        }
    }
}

template<class C, class T> void Adjacency_Graph_Extractor::extract_adjacency_graph_blockers_lite(Node_Stellar &n, Mesh<C,T> &mesh, Adjacency_Graph &graph)
{
//    cout<<n<<endl;

    ivect f;
    for(int d=0; d<graph.get_blockers_num();d++)
    {
        for(int i=1; i<=graph.get_blockers_num(d);i++)
        {
            Top_Simplex &blocker = graph.get_blocker(d,i);
//            cout<<"=====blocker: "<<blocker<<endl;
            if(n.indexes_top(blocker))
            {
                int blocker_global_id = graph.get_blocker_global_id(d,i);

                /// (1) first check if the d-1 faces are top simplices
                ///  then, remove these from the graph..
                ///  and veiculate the arcs to the blockers
                for(int v=0; v<blocker.get_vertices_num(); v++)
                {
                    blocker.TF(f,v);
                    int d = f.size() -2;
                    for(RunIterator runIt = n.t_array_begin_iterator(d), runEnd = n.t_array_end_iterator(d); runIt != runEnd; ++runIt)
                    {
                        /// if the top-simplex has not been deleted during the simplification
                        /// or if it the immediate boundary of a blocker
                        if(!mesh.is_top_cell_removed(d,*runIt) && !graph.is_top_a_blocker_d1face(d,*runIt))
                        {
                            T &top = mesh.get_top_cell(d,*runIt);

                            if(top.are_equal(f))
                            {
//                                cout<<"top: "<<top<<endl;
//                                cout<<"clear top -> push it into blocker"<<endl;
                                int top_global_id = graph.get_global_id(d,*runIt);
                                set<pair<int,short> > &adjs = graph.get_adj_set(top_global_id);
                                graph.add_pairs_set(blocker_global_id,adjs);
                                graph.reset_adj_set(top_global_id);
                                graph.set_top_as_blocker_d1face(d,*runIt);
//                                int a; cin>>a;
                            }
                        }
                    }
                }
                /// (2) check if exists a d-1 or d-2 adjacency with the indexed top simplices of the leaf block
                /// NOTA: how to flag those alredy considered in the previous step? <<=====
                ///
                Adjacency_Graph_Extractor::extract_d1_d2_adjacents(blocker,i,d,true,n,mesh,graph);
            }
        }
    }

    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIterator runIt = n.t_array_begin_iterator(d), runEnd = n.t_array_end_iterator(d); runIt != runEnd; ++runIt)
        {
            if(!mesh.is_top_cell_removed(d,*runIt) && !graph.is_top_a_blocker_d1face(d,*runIt))
            {
                T& top = mesh.get_top_cell(d,*runIt);
//                cout<<"==target-top: "<<top<<endl;
                Adjacency_Graph_Extractor::extract_d1_d2_adjacents(top,*runIt,d,false,n,mesh,graph);
            }
        }
    }
}



template<class C, class T> void Adjacency_Graph_Extractor::extract_d1_d2_adjacents(T &t, int t_id, int d, bool is_blocker, Node_Stellar &n, Mesh<C,T> &mesh, Adjacency_Graph &graph)
{
//    cout<<"extract_d1_d2_adjacents"<<endl;

    ivect top_v; t.get_positive_sorted_vertices(top_v);
//    print_container_content("sorted_vs: ",top_v);

//    cout<<"A "<<is_blocker<<" d: "<<d<<" id: "<<t_id<<endl;
//    vector<int_vect> &gi = graph.get_global_ids();
//    cout<<"graph_types: "<<gi.size()<<endl;
//    cout<<"graph_d_types: "<<gi[d].size()<<endl;

    int global_t_id;
    if(is_blocker)
        global_t_id = graph.get_blocker_global_id(d,t_id);
    else
        global_t_id = graph.get_global_id(d,t_id);
    ivect top2_v;

//    cout<<"B "<<is_blocker<<endl;
//    cout<<" d: "<<d<<" id: "<<t_id<<" global id: "<<global_t_id<<endl;

    /// first check d-1 adjacencies
    for(RunIterator runIt2 = n.t_array_begin_iterator(d), runEnd2 = n.t_array_end_iterator(d); runIt2 != runEnd2; ++runIt2)
    {
        /// if we are processing a blocker we analyze all the top simplices
        /// nota: we skip a top simplex if we have identified it as a d-1 face of a blocker
        if(!graph.is_top_a_blocker_d1face(d,*runIt2) && (is_blocker || t_id != *runIt2))
        {
            T& top2 = mesh.get_top_cell(d,*runIt2);

            int top2_id = graph.get_global_id(d,*runIt2);
            top2.get_positive_sorted_vertices(top2_v);
            intersect_containers(top2_v,top_v);

//            cout<<" same-d: "<<d<<" id: "<<*runIt2<<" global id: "<<top2_id<<endl;

            if(top2_v.size() == t.get_vertices_num()-1 ) /// we have a d-1 adjacency
            {
//                cout<<"top2: "<<top2<<endl;
//                cout<<"     d-1 adjacency"<<endl;
                // we just push each edge once
                pair<int,short> p;
                p.first = top2_id;
                p.second = top2_v.size(); /// we write the correct label (if we refer to the adj it should be -1)
                graph.add_pair(global_t_id,p);

//                if(p.first==-1)
//                {
//                    int a; cin>>a;
//                }
            }
        }
    }

    if(d>0) /// at least top 2-simplices
    {
        int d2 = d-1;
        for(RunIterator runIt2 = n.t_array_begin_iterator(d2), runEnd2 = n.t_array_end_iterator(d2); runIt2 != runEnd2; ++runIt2)
        {
            if(!graph.is_top_a_blocker_d1face(d2,*runIt2))
            {
                T& top2 = mesh.get_top_cell(d2,*runIt2);

                int top2_id = graph.get_global_id(d2,*runIt2);
                top2.get_positive_sorted_vertices(top2_v);
                intersect_containers(top2_v,top_v);

//                cout<<" d-1: "<<d2<<" id: "<<*runIt2<<" global id: "<<top2_id<<endl;

                if(top2_v.size() == top2.get_vertices_num()-1 ) /// we have a d-2 adjacency with respect to the higher dimensional simplex
                {
//                    cout<<"top2: "<<top2<<endl;
//                    cout<<"     d-2 adjacency"<<endl;
                    // we just push each edge once
                    pair<int,short> p;
                    p.first = top2_id;
                    p.second = top2_v.size(); /// we write the correct label (if we refer to the adj it should be -1)
                    graph.add_pair(global_t_id,p);

//                    if(p.first==-1)
//                    {
//                        int a; cin>>a;
//                    }
                }
            }
        }
    }
}

template<class C, class T> void Adjacency_Graph_Extractor::extract_local_adjacency_graph(Node_Stellar &n, Mesh<C,T> &mesh, Local_Adjacency_Graph &graph)
{
    graph.set_ids(n);

    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIterator runIt = n.t_array_begin_iterator(d), runEnd = n.t_array_end_iterator(d); runIt != runEnd; ++runIt)
        {
            if(!mesh.is_top_cell_removed(d,*runIt))
            {
                T& top = mesh.get_top_cell(d,*runIt);

                ivect top_v; top.get_positive_sorted_vertices(top_v);
                pair<int,int> key = make_pair(d,*runIt);
                int top_id = graph.get_global_id(key);
                ivect top2_v;

                // loop on all the successive top simplices
                for(unsigned d2=d; d2<n.get_num_top_cells_encoded(); d2++)
                {
                    for(RunIterator runIt2 = n.t_array_begin_iterator(d2), runEnd2 = n.t_array_end_iterator(d2); runIt2 != runEnd2; ++runIt2)
                    {
                        if(d!=d2 || (d==d2 && runIt!=runIt2)) /// we are not processing a top-d-simplex against itself
                        {
                            T& top2 = mesh.get_top_cell(d2,*runIt2);
                            pair<int,int> key2 = make_pair(d2,*runIt2);
                            int top2_id = graph.get_global_id(key2);
                            top2.get_positive_sorted_vertices(top2_v);
                            intersect_containers(top2_v,top_v);

                            if(top2_v.size() > 1) // we skip those top simplices that are incident just in a vertex
                            {
                                // we just push each edge once
                                pair<int,short> p;
                                p.first = top2_id;
                                p.second = top2_v.size();
                                graph.add_pair(top_id,p);
//                                adj_matrix[top_id].push_back(p);
    //                                cout << top_id << " -- " << top2_id << " color = " << palette[p.second].a << " " << palette[p.second].b << " " << palette[p.second].c << endl;
    //                                adj_matrix[top_id][top2_id] = top2_v.size()-1;
    //                                adj_matrix[top2_id][top_id] = top2_v.size()-1;
                            }
    //                            if(top2_v.size() == top_v.size()-1) // d-1 adjacency
    //                            {
    //                                d_1_adjs[top_id].push_back(top2_id);
    //                            }
                        }
                    }
                }
            }
        }
    }

    graph.write_graph_file(n,mesh);
}



template<class C, class T> void Adjacency_Graph_Extractor::extract_local_adjacency_graph_lite(Node_Stellar &n, Mesh<C,T> &mesh, Local_Adjacency_Graph &graph)
{
    graph.set_ids(n);

    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIterator runIt = n.t_array_begin_iterator(d), runEnd = n.t_array_end_iterator(d); runIt != runEnd; ++runIt)
        {
            if(!mesh.is_top_cell_removed(d,*runIt))
            {
                T& top = mesh.get_top_cell(d,*runIt);

                ivect top_v; top.get_positive_sorted_vertices(top_v);
                pair<int,int> key = make_pair(d,*runIt);
                int top_id = graph.get_global_id(key);
                ivect top2_v;

                /// first check d-1 adjacencies
                for(RunIterator runIt2 = n.t_array_begin_iterator(d), runEnd2 = n.t_array_end_iterator(d); runIt2 != runEnd2; ++runIt2)
                {
                    if(runIt != runIt2) ///
                    {
                        T& top2 = mesh.get_top_cell(d,*runIt2);
                        pair<int,int> key2 = make_pair(d,*runIt2);
                        int top2_id = graph.get_global_id(key2);
                        top2.get_positive_sorted_vertices(top2_v);
                        intersect_containers(top2_v,top_v);

                        if(top2_v.size() == top.get_vertices_num()-1 ) /// we have a d-1 adjacency
                        {
                            // we just push each edge once
                            pair<int,short> p;
                            p.first = top2_id;
                            p.second = top2_v.size();
                            graph.add_pair(top_id,p);
                        }
                    }
                }

                if(d>0) /// at least top 2-simplices
                {
                    int d2 = d-1;
                    for(RunIterator runIt2 = n.t_array_begin_iterator(d2), runEnd2 = n.t_array_end_iterator(d2); runIt2 != runEnd2; ++runIt2)
                    {
                        T& top2 = mesh.get_top_cell(d2,*runIt2);
                        pair<int,int> key2 = make_pair(d2,*runIt2);
                        int top2_id = graph.get_global_id(key2);
                        top2.get_positive_sorted_vertices(top2_v);
                        intersect_containers(top2_v,top_v);

                        if(top2_v.size() == top2.get_vertices_num()-1 ) /// we have a d-2 adjacency
                        {
                            // we just push each edge once
                            pair<int,short> p;
                            p.first = top2_id;
                            p.second = top2_v.size();
                            graph.add_pair(top_id,p);
                        }
                    }
                }
            }
        }
    }

    graph.write_graph_file(n,mesh);
}


template<class C, class T> void Adjacency_Graph_Extractor::extract_local_1skeleton_graph(Node_Stellar &n, Mesh<C,T> &mesh, string mesh_name)
{
    leaf_VV vvs;
    /// 1) extract local VV relations
    n.extract_local_VV(mesh,vvs);

    stringstream ss;
    ss << mesh_name << "_L_" << n.get_v_start() << "_" << n.get_v_end() << "_1skel_graph.gv";

    ofstream out_graph(ss.str());
    out_graph << "graph {" <<endl;

    for(unsigned i=0; i<vvs.size(); i++)
    {
        int real_id = i+n.get_v_start();
        set<int> &vv = vvs[i];
        for(auto it=vv.begin(); it!=vv.end(); ++it)
        {
            if(real_id < *it)
                out_graph << real_id << "--" << *it << endl;
        }
    }

    out_graph << "}" <<endl;
    out_graph.close();
}

#endif // ADJACENCY_GRAPH_EXTRACTOR_H
