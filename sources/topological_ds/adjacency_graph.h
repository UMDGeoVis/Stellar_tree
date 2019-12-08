#ifndef ADJACENCY_GRAPH_H
#define ADJACENCY_GRAPH_H

#include <sstream>
#include <fstream>
#include <cmath>
#define PI 3.14159265358979323846

#include "stellar_tree/mesh.h"
#include "utilities/basic_wrappers.h"
#include "stellar_decomposition/node_stellar.h"
#include "utilities/color.h"

typedef set<Top_Simplex> blocker_set;
typedef vector<Top_Simplex> blocker_array;

using namespace std;

class Adjacency_Graph
{
public:
    template<class C, class T> Adjacency_Graph(Mesh<C,T> &mesh);
    template<class C, class T> Adjacency_Graph(Mesh<C,T> &mesh, blocker_set &b_set);

    void reset()
    {
        palette.clear();
        global_ids.clear();
        adj_sets.clear();
        blockers.clear();
        blockers_global_ids.clear();
        hided_tops.clear();
    }

    inline vector<color>& get_palette() { return palette; }

    inline int get_global_id(int d, int pos) { return global_ids[d][pos-1]; } /// nota
    inline bool is_top_a_blocker_d1face(int d, int pos) { return (global_ids[d][pos-1]==-1); }
    inline void set_top_as_blocker_d1face(int d, int pos) { global_ids[d][pos-1] = -1; }
    inline vector<ivect>& get_global_ids() { return global_ids; }

    inline void add_pair(int top_id, pair<int,short> &p) { adj_sets[top_id].insert(p); }
    inline void add_pairs_set(int top_id, set<pair<int,short> > &p) { adj_sets[top_id].insert(p.begin(),p.end()); }
    inline vector<set<pair<int,short> > >& get_adj_sets() { return adj_sets; }
    inline set<pair<int,short> >& get_adj_set(int i) { return adj_sets[i]; }
    inline void reset_adj_set(int top_id)
    {
        adj_sets[top_id].clear();
        hided_tops.insert(top_id); /// if we reset a top adj set it means that the top is a d-1 face of a blocker.. thus in practice we hide it in the adj graph
    }

    inline Top_Simplex& get_blocker(int d, int pos) { return blockers[d][pos-1]; } /// nota
    inline int get_blockers_num() { return blockers.size(); }
    inline int get_blockers_num(int d) { return blockers[d].size(); }
    inline vector<blocker_array>& get_blocker_array() { return blockers; }
    inline vector<ivect>& get_blockers_global_ids() { return blockers_global_ids; }
    inline int get_blocker_global_id(int d, int pos) { return blockers_global_ids[d][pos-1]; } /// nota

    template<class C, class T> void write_graph_file(string mesh_name, Mesh<C,T> &mesh);
    template<class C, class T> void write_fixed_graph_file(string mesh_name, Mesh<C,T> &mesh);
    template<class C, class T> void write_cesium_graph_file(string mesh_name, Mesh<C,T> &mesh);
    void write_blockers_file(string mesh_name);
    template<class C, class T> void write_fixed_blockers_file(string mesh_name, Mesh<C,T> &mesh);
    template<class C, class T> void write_cesium_blockers_file(string mesh_name, Mesh<C,T> &mesh);

    void compute_adj_graph_stats();

private:
    vector<color> palette;
    vector<ivect> global_ids;
    vector<set<pair<int,short> > > adj_sets;

    vector<blocker_array> blockers;
    vector<ivect> blockers_global_ids;
    int blocker_range_start, blocker_range_end;
    iset hided_tops;

    void set_blockers_adjacencies();
    void extract_d1_d2_adjacent_blockers(Top_Simplex &t, int t_id, int d);

    template<class C, class T> Point<C> get_reference_point(T& t, Mesh<C,T> &mesh, boost::dynamic_bitset<> &unused_vertices);
    dvect cartesian(double lat, double lon);
    dvect spherical(dvect &c);
};

template<class C, class T> Adjacency_Graph::Adjacency_Graph(Mesh<C,T> &mesh)
{
    get_hsv_palette(palette,mesh.get_implicitly_and_explicitly_encoded_cells_num());
//    print_palette(palette); /// for debug

    //give new 'unique' position indices for all the top simplices
    global_ids.assign(mesh.get_top_cells_types(),ivect());

//    cout<<"     generating global position indices"<<endl;
    int tot_simpl_num = 0, id=0;
    for(unsigned i=0; i<mesh.get_top_cells_types(); i++)
    {
        global_ids[i].assign(mesh.get_top_cells_num(i),0);
        tot_simpl_num += mesh.get_top_cells_num(i);
        for(unsigned j=0; j<mesh.get_top_cells_num(i); j++)
        {
            global_ids[i][j] = id;
            id++;
        }
    }

    blocker_range_start = blocker_range_end = -1;

    adj_sets.assign(tot_simpl_num,set<pair<int,short> >());
}

template<class C, class T> Adjacency_Graph::Adjacency_Graph(Mesh<C,T> &mesh, blocker_set &b_set)
{
    get_hsv_palette(palette,mesh.get_implicitly_and_explicitly_encoded_cells_num());
//    print_palette(palette); /// for debug

    //give new 'unique' position indices for all the top simplices
    global_ids.assign(mesh.get_top_cells_types(),ivect());

//    cout<<"     generating global position indices"<<endl;
    int tot_simpl_num = 0, id=0;
    for(unsigned i=0; i<mesh.get_top_cells_types(); i++)
    {
        global_ids[i].assign(mesh.get_top_cells_num(i),0);
        tot_simpl_num += mesh.get_top_cells_num(i);
        for(unsigned j=0; j<mesh.get_top_cells_num(i); j++)
        {
            global_ids[i][j] = id;
            id++;
        }
    }

    if(b_set.size() == 0) /// no blockers
        blocker_range_start = blocker_range_end = -1;
    else
    {
        blocker_range_start = id; /// first entry inside the range

        blockers.assign(mesh.get_top_cells_types(),blocker_array());
        for(auto it=b_set.begin(); it!=b_set.end(); ++it)
        {
            const Top_Simplex &blocker = *it;
            blockers[blocker.get_vertices_num()-2].push_back(blocker);
        }
        blockers_global_ids.assign(mesh.get_top_cells_types(),ivect());
        for(unsigned i=0; i<blockers.size(); i++)
        {
            blockers_global_ids[i].assign(blockers[i].size(),0);
            tot_simpl_num += blockers[i].size();
            for(unsigned j=0; j<blockers[i].size(); j++)
            {
                blockers_global_ids[i][j] = id;
                id++;
            }
        }

        blocker_range_end = id; /// nota: first entry outside the range
    }

    adj_sets.assign(tot_simpl_num,set<pair<int,short> >()); /// in the adj list we have room also for the blockers!

    this->set_blockers_adjacencies();

    b_set.clear(); /// erase the blocker set
}

template<class C, class T> void Adjacency_Graph::write_graph_file(string mesh_name, Mesh<C,T> &mesh)
{
    stringstream ss/*, ss2*//*, ss3, ss4*/;
    ss << mesh_name << "_graph.gv";
//        ss2 << mesh_name << "_d_1_graph.gv";
//        ss3 << mesh_name << "_adj_lists.txt";
//        ss4 << mesh_name << "_adj_lists_01.txt";

    // graphviz compatible graph format
    ofstream out_graph(ss.str());
    out_graph << "graph {" <<endl;
    out_graph << "node [style=bold];" << endl;
    // (1) write top nodes attributes
    for(unsigned i=0; i<global_ids.size(); i++)
    {
        for(unsigned j=0; j<global_ids[i].size(); j++)
        {
            /// we skip those tops that a d-1 face of a blocker
            if(global_ids[i][j] != -1)
            {
                T &t = mesh.get_top_cell(i,j+1);
                color &c = palette[t.get_vertices_num()-1];
                out_graph << global_ids[i][j] << " [fillcolor=\"" << c.a << "," << c.b << "," << c.c << "\""
                          <<",color=\"" << c.a << "," << c.b << "," << c.c << "\""
                          <<",label=\""<<t.get_vertices_num()-1<<"\"];" <<endl;
            }
        }
    }
    // (2) write blockers
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
        for(unsigned j=0; j<blockers_global_ids[i].size(); j++)
        {
            T &t = blockers[i][j];
            color &c = palette[t.get_vertices_num()-1];
            out_graph << blockers_global_ids[i][j] << " [fillcolor=\"" << c.a << "," << c.b << "," << c.c << "\""
                      <<",color=\"" << c.a << "," << c.b << "," << c.c << "\""
                      <<",label=\""<<t.get_vertices_num()-1<<"\",shape=box];" <<endl;
        }
    }
    // (2) write arcs attribute (contains both the arcs to blockers)
    for(unsigned i=0; i<adj_sets.size(); i++)
    {
//            out_graph <<i<< " -- { ";
//            for(unsigned j=0; j<adj_lists[i].size(); j++)
        for(set<pair<int,short> >::iterator it=adj_sets[i].begin(); it!=adj_sets[i].end(); ++it)
        {
            if(hided_tops.find(it->first) == hided_tops.end())
            {
                color &c = palette[it->second];
                out_graph << i << "--" << it->first << "[color=\"" << c.a << "," << c.b << "," << c.c << "\","
                          <<"label=\""<< it->second-1 <<"\"];"
                          <<endl;
            }
        }
//            out_graph << "};" <<endl;
    }
    out_graph << "}" <<endl;
    out_graph.close();
}

template<class C, class T> void Adjacency_Graph::write_fixed_graph_file(string mesh_name, Mesh<C,T> &mesh)
{
    stringstream ss/*, ss2*//*, ss3, ss4*/;
    ss << mesh_name << "_fixed_graph.gv";
//        ss2 << mesh_name << "_d_1_graph.gv";
//        ss3 << mesh_name << "_adj_lists.txt";
//        ss4 << mesh_name << "_adj_lists_01.txt";
    boost::dynamic_bitset<> unused_vertices(mesh.get_vertices_num());

    // graphviz compatible graph format
    ofstream out_graph(ss.str());
    out_graph << "graph {" <<endl;
    out_graph << "node [style=bold];" << endl;
    // (1) write top nodes attributes
    for(unsigned i=0; i<global_ids.size(); i++)
    {
        for(unsigned j=0; j<global_ids[i].size(); j++)
        {
            /// we skip those tops that a d-1 face of a blocker
            if(global_ids[i][j] != -1)
            {
                T &t = mesh.get_top_cell(i,j+1);
                color &c = palette[t.get_vertices_num()-1];
                Point<C> ref = get_reference_point(t,mesh,unused_vertices);

                out_graph << global_ids[i][j] << " [fillcolor=\"" << c.a << "," << c.b << "," << c.c << "\""
                          <<",color=\"" << c.a << "," << c.b << "," << c.c << "\""
                         <<",label=\""<<t.get_vertices_num()-1<<"\""
                        <<",pos=\""<<ref.getC(0)<<","<<ref.getC(1)<<"!\"];"<<endl;
            }
        }
    }
    // (2) write blockers
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
        for(unsigned j=0; j<blockers_global_ids[i].size(); j++)
        {
            T &t = blockers[i][j];
            color &c = palette[t.get_vertices_num()-1];
            Point<C> ref = get_reference_point(t,mesh,unused_vertices);

            out_graph << blockers_global_ids[i][j] << " [fillcolor=\"" << c.a << "," << c.b << "," << c.c << "\""
                      <<",color=\"" << c.a << "," << c.b << "," << c.c << "\""
                      <<",label=\""<<t.get_vertices_num()-1<<"\""
                      <<",pos=\""<<ref.getC(0)<<","<<ref.getC(1)<<"!\",shape=box];" <<endl;
        }
    }
    // (2) write arcs attribute (contains both the arcs to blockers)
    for(unsigned i=0; i<adj_sets.size(); i++)
    {
        for(set<pair<int,short> >::iterator it=adj_sets[i].begin(); it!=adj_sets[i].end(); ++it)
        {
            if(hided_tops.find(it->first) == hided_tops.end())
            {
                color &c = palette[it->second];
                out_graph << i << "--" << it->first << "[color=\"" << c.a << "," << c.b << "," << c.c << "\","
                          <<"label=\""<< it->second-1 <<"\"];"
                         <<endl;
            }
        }
    }
    out_graph << "}" <<endl;
    out_graph.close();
}

template<class C, class T> void Adjacency_Graph::write_cesium_graph_file(string mesh_name, Mesh<C,T> &mesh)
{
    stringstream ss, ssb;
    ss << mesh_name << "_graph.json";
    ssb << mesh_name << "_blockers_graph.json";
    boost::dynamic_bitset<> unused_vertices(mesh.get_vertices_num());

    map<int,short> top_dim_map; /// auxiliary map that keeps track of the dimension of each top/blocker
    int valid_entries = 0;
    for(unsigned i=0; i<global_ids.size(); i++)
    {
        /// first we count the effectively initialized
        for(unsigned j=0; j<global_ids[i].size(); j++)
            if(global_ids[i][j] != -1)
                valid_entries++;
    }

    int top_s_counter_valid = 0;

    ofstream out_graph(ss.str());
    out_graph << "{ \"nodes\":[" <<endl;
    ofstream out_bgraph(ssb.str());
    out_bgraph << "{ \"nodes\":[" <<endl;
    // (1) write top nodes attributes
    for(unsigned i=0; i<global_ids.size(); i++)
    {
        /// first we count the effectively initialized
        int c = 0;
        for(unsigned j=0; j<global_ids[i].size(); j++)
            if(global_ids[i][j] != -1)
                c++;

        int c2 = 1;
        for(unsigned j=0; j<global_ids[i].size(); j++)
        {
            /// we skip those tops that a d-1 face of a blocker
            if(global_ids[i][j] != -1)
            {
                T &t = mesh.get_top_cell(i,j+1);
//                color &c = palette[t.get_vertices_num()-1];
                Point<C> ref = get_reference_point(t,mesh,unused_vertices);

                out_graph << "{ \"id\":" << global_ids[i][j]
                          << ",\"dim\":" << t.get_vertices_num()-1
                          << ",\"long\":" << ref.getC(0)
                          << ",\"lat\":" << ref.getC(1)
                          << ",\"isBlocker\":0}";
                if(c2 < c)
                    out_graph << "," << endl;
//                else
//                    out_graph << endl;
                c2++;

                top_dim_map[global_ids[i][j]] = t.get_vertices_num()-1;

                top_s_counter_valid++;
            }
        }

        if(c > 0)
        {
            if(top_s_counter_valid < valid_entries)
                out_graph << "," << endl;
//            else
//                out_graph << endl;
        }
    }
    if(blockers_global_ids.size() > 0) /// if we have blockers
        out_graph << "," << endl;
    else
        out_graph << endl;
    // (2) write blockers
    int blocker_counter = 0;
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
        for(unsigned j=0; j<blockers_global_ids[i].size(); j++)
            blocker_counter++;
    }
    int new_counter_b = 0;
//    out_graph << "]," << endl << "\"blockers\":[" << endl;
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
        for(unsigned j=0; j<blockers_global_ids[i].size(); j++)
        {
            T &t = blockers[i][j];
            Point<C> ref = get_reference_point(t,mesh,unused_vertices);

            out_graph << "{ \"id\":" << blockers_global_ids[i][j]
                      << ",\"dim\":" << t.get_vertices_num()-1
                      << ",\"long\":" << ref.getC(0)
                      << ",\"lat\":" << ref.getC(1)
                      << ",\"isBlocker\":1}";
            out_bgraph << "{ \"id\":" << blockers_global_ids[i][j]
                      << ",\"dim\":" << t.get_vertices_num()-1
                      << ",\"long\":" << ref.getC(0)
                      << ",\"lat\":" << ref.getC(1)
                      << ",\"isBlocker\":1}";

            if(j < blockers_global_ids[i].size()-1)
            {
                out_graph << "," << endl;
                out_bgraph << "," << endl;
            }
//            else
//                out_graph << endl;
            top_dim_map[blockers_global_ids[i][j]] = t.get_vertices_num()-1;
            new_counter_b++;
        }
        if(blockers_global_ids[i].size() > 0)
        {
            if(new_counter_b < blocker_counter)
            {
                out_graph << "," << endl;
                out_bgraph << "," << endl;
            }
            else
            {
                out_graph << endl;
                out_bgraph << endl;
            }
        }
    }
//    valid_entries += blocker_counter;
//    cout<<"valid entries: "<<valid_entries<<endl;
    // (2) write arcs attribute (contains both the arcs to blockers)
    int a_c = 0;
    for(unsigned i=0; i<adj_sets.size(); i++)
    {
        if(adj_sets[i].size() > 0)
            a_c++;
    }
    int local_counter = 0;

    out_graph << "]," << endl << "\"arcs\":[" << endl;
    for(unsigned i=0; i<adj_sets.size(); i++)
    {
        if(adj_sets[i].size() > 0)
        {
            local_counter++;

            /// first we have to unde
            short d = top_dim_map[i];
            short d1 = d, d2 = d-1;

//            cout<<d<<" "<<d1<<" "<<d2<<endl;

            int numd1 = 0, numd2 = 0;
            for(set<pair<int,short> >::iterator it=adj_sets[i].begin(); it!=adj_sets[i].end(); ++it)
            {
                if(hided_tops.find(it->first) == hided_tops.end())
                {
//                    cout<<"label: "<<it->second<<endl;
                    if(it->second == d1)
                        numd1++;
                    else if(it->second == d2)
                        numd2++;
                }
            }

//            cout<<adj_sets[i].size()<<" "<<numd1<<" "<<numd2<<endl;
//            int a; cin>>a;

            out_graph << "{ \"id\":" << i <<",\"d1len\":"<<numd1<<",\"d1adjs\": [";
            int c=0;
            for(set<pair<int,short> >::iterator it=adj_sets[i].begin(); it!=adj_sets[i].end(); ++it)
            {
                if(hided_tops.find(it->first) == hided_tops.end())
                {
                    if(it->second == d1)
                    {
                        out_graph<<it->first;
                        if(c < numd1-1)
                            out_graph << ",";
                        c++;
                    }
                }
            }
            c=0;
            out_graph << "],\"d2len\":"<<numd2<<",\"d2adjs\": [";
            for(set<pair<int,short> >::iterator it=adj_sets[i].begin(); it!=adj_sets[i].end(); ++it)
            {
                if(hided_tops.find(it->first) == hided_tops.end())
                {
                    if(it->second == d2)
                    {
                        out_graph<<it->first;
                        if(c < numd2-1)
                            out_graph << ",";
                        c++;
                    }
                }
            }
            out_graph << "]}";
//            cout<<local_counter<<" "<<a_c<<endl;
            if(local_counter < a_c)
                out_graph << "," << endl;
            else
                out_graph << endl;
        }

    }
    out_graph << "]}" <<endl;
    out_graph.close();

    // (2) write arcs attribute (contains both the arcs to blockers)
    a_c = 0;
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
//        if(blockers_global_ids[i].size() > 0)
        a_c += blockers_global_ids[i].size();
    }
    local_counter = 0;
    out_bgraph << "]," << endl << "\"arcs\":[" << endl;
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
        for(unsigned j=0; j<blockers_global_ids[i].size(); j++)
        {
            local_counter++;

            int blocker_pos = blockers_global_ids[i][j];

            /// first we have to unde
            short d = top_dim_map[blocker_pos];
            short d1 = d, d2 = d-1;

//            cout<<d<<" "<<d1<<" "<<d2<<endl;

            int numd1 = 0, numd2 = 0;
            for(set<pair<int,short> >::iterator it=adj_sets[blocker_pos].begin(); it!=adj_sets[blocker_pos].end(); ++it)
            {
                if(blocker_range_start <= it->first && it->first < blocker_range_end)
                {
                    if(it->second == d1)
                        numd1++;
                    else if(it->second == d2)
                        numd2++;
                }
            }

//            cout<<numd1<<" "<<numd2<<endl;

            out_bgraph << "{ \"id\":" << blocker_pos <<",\"d1len\":"<<numd1<<",\"d1adjs\": [";
            int c=0;
            for(set<pair<int,short> >::iterator it=adj_sets[blocker_pos].begin(); it!=adj_sets[blocker_pos].end(); ++it)
            {
                /// we output just the blockers
                if(blocker_range_start <= it->first && it->first < blocker_range_end)
                {
                    if(it->second == d1)
                    {
                        out_bgraph<<it->first;
                        if(c < numd1-1)
                            out_bgraph << ",";
                        c++;
                    }
                }
            }
            c=0;
            out_bgraph << "],\"d2len\":"<<numd2<<",\"d2adjs\": [";
            for(set<pair<int,short> >::iterator it=adj_sets[blocker_pos].begin(); it!=adj_sets[blocker_pos].end(); ++it)
            {
                /// we output just the blockers
                if(blocker_range_start <= it->first && it->first < blocker_range_end)
                {
                    if(it->second == d2)
                    {
                        out_bgraph<<it->first;
                        if(c < numd2-1)
                            out_bgraph << ",";
                        c++;
                    }
                }
            }
            out_bgraph << "]}";
//            cout<<local_counter<<" "<<a_c<<endl;
            if(local_counter < a_c)
                out_bgraph << "," << endl;
            else
                out_bgraph << endl;
        }
    }
    out_bgraph << "]}" <<endl;
    out_bgraph.close();
}

template<class C, class T> void Adjacency_Graph::write_cesium_blockers_file(string mesh_name, Mesh<C,T> &mesh)
{
    stringstream ss;
    ss << mesh_name << "_blockers_graph.json";
    boost::dynamic_bitset<> unused_vertices(mesh.get_vertices_num());

    map<int,short> top_dim_map; /// auxiliary map that keeps track of the dimension of each top/blocker
    int valid_entries = 0;
    for(unsigned i=0; i<global_ids.size(); i++)
    {
        /// first we count the effectively initialized
        for(unsigned j=0; j<global_ids[i].size(); j++)
            if(global_ids[i][j] != -1)
                valid_entries++;
    }

    ofstream out_graph(ss.str());
    out_graph << "{ \"nodes\":[" <<endl;
    // (1) write blockers
    int blocker_counter = 0;
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
        for(unsigned j=0; j<blockers_global_ids[i].size(); j++)
            blocker_counter++;
    }
    int new_counter_b = 0;

    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
        for(unsigned j=0; j<blockers_global_ids[i].size(); j++)
        {
            T &t = blockers[i][j];
            Point<C> ref = get_reference_point(t,mesh,unused_vertices);

            out_graph << "{ \"id\":" << blockers_global_ids[i][j]
                      << ",\"dim\":" << t.get_vertices_num()-1
                      << ",\"long\":" << ref.getC(0)
                      << ",\"lat\":" << ref.getC(1)
                      << ",\"isBlocker\":1}";

            if(j < blockers_global_ids[i].size()-1)
                out_graph << "," << endl;

            top_dim_map[blockers_global_ids[i][j]] = t.get_vertices_num()-1;
            new_counter_b++;
        }
        if(blockers_global_ids[i].size() > 0)
        {
//            cout<<new_counter_b<<" "<<blocker_counter<<endl;
            if(new_counter_b < blocker_counter)
                out_graph << "," << endl;
            else
                out_graph << endl;
        }
    }
    // (2) write arcs attribute (contains both the arcs to blockers)
    int a_c = 0;
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
//        if(blockers_global_ids[i].size() > 0)
        a_c += blockers_global_ids[i].size();
    }
    int local_counter = 0;
    out_graph << "]," << endl << "\"arcs\":[" << endl;
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
        for(unsigned j=0; j<blockers_global_ids[i].size(); j++)
        {
            local_counter++;

            int blocker_pos = blockers_global_ids[i][j];

            /// first we have to unde
            short d = top_dim_map[blocker_pos];
            short d1 = d, d2 = d-1;

//            cout<<d<<" "<<d1<<" "<<d2<<endl;

            int numd1 = 0, numd2 = 0;
            for(set<pair<int,short> >::iterator it=adj_sets[blocker_pos].begin(); it!=adj_sets[blocker_pos].end(); ++it)
            {
                if(blocker_range_start <= it->first && it->first < blocker_range_end)
                {
                    if(it->second == d1)
                        numd1++;
                    else if(it->second == d2)
                        numd2++;
                }
            }

//            cout<<numd1<<" "<<numd2<<endl;

            out_graph << "{ \"id\":" << blocker_pos <<",\"d1len\":"<<numd1<<",\"d1adjs\": [";
            int c=0;
            for(set<pair<int,short> >::iterator it=adj_sets[blocker_pos].begin(); it!=adj_sets[blocker_pos].end(); ++it)
            {
                /// we output just the blockers
                if(blocker_range_start <= it->first && it->first < blocker_range_end)
                {
                    if(it->second == d1)
                    {
                        out_graph<<it->first;
                        if(c < numd1-1)
                            out_graph << ",";
                        c++;
                    }
                }
            }
            c=0;
            out_graph << "],\"d2len\":"<<numd2<<",\"d2adjs\": [";
            for(set<pair<int,short> >::iterator it=adj_sets[blocker_pos].begin(); it!=adj_sets[blocker_pos].end(); ++it)
            {
                /// we output just the blockers
                if(blocker_range_start <= it->first && it->first < blocker_range_end)
                {
                    if(it->second == d2)
                    {
                        out_graph<<it->first;
                        if(c < numd2-1)
                            out_graph << ",";
                        c++;
                    }
                }
            }
            out_graph << "]}";
//            cout<<local_counter<<" "<<a_c<<endl;
            if(local_counter < a_c)
                out_graph << "," << endl;
            else
                out_graph << endl;
        }
    }
    out_graph << "]}" <<endl;
    out_graph.close();
}

template<class C, class T> void Adjacency_Graph::write_fixed_blockers_file(string mesh_name, Mesh<C,T> &mesh)
{
    stringstream ss/*, ss2*//*, ss3, ss4*/;
    ss << mesh_name << "_fixed_blocker_graph.gv";
    boost::dynamic_bitset<> unused_vertices(mesh.get_vertices_num());
//        ss2 << mesh_name << "_d_1_graph.gv";
//        ss3 << mesh_name << "_adj_lists.txt";
//        ss4 << mesh_name << "_adj_lists_01.txt";

    // graphviz compatible graph format
    ofstream out_graph(ss.str());
    out_graph << "graph {" <<endl;
    out_graph << "node [style=bold];" << endl;
    // (2) write blockers
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
        for(unsigned j=0; j<blockers_global_ids[i].size(); j++)
        {
            T &t = blockers[i][j];
            color &c = palette[t.get_vertices_num()-1];
            Point<C> ref = get_reference_point(t,mesh,unused_vertices);

            out_graph << blockers_global_ids[i][j] << " [fillcolor=\"" << c.a << "," << c.b << "," << c.c << "\""
                      <<",color=\"" << c.a << "," << c.b << "," << c.c << "\""
                      <<",label=\""<<t.get_vertices_num()-1<<"\""
                      <<",pos=\""<<ref.getC(0)<<","<<ref.getC(1)<<"!\",shape=box];" <<endl;
        }
    }
    // (2) write arcs attribute (connecting blockers)
    for(unsigned i=0; i<blockers_global_ids.size(); i++)
    {
        for(unsigned j=0; j<blockers_global_ids[i].size(); j++)
        {
            int blocker_pos = blockers_global_ids[i][j];
            for(set<pair<int,short> >::iterator it=adj_sets[blocker_pos].begin(); it!=adj_sets[blocker_pos].end(); ++it)
            {
                /// we output just the blockers
                if(blocker_range_start <= it->first && it->first < blocker_range_end)
                {
                    color &c = palette[it->second];
                    out_graph << blocker_pos << "--" << it->first << "[color=\"" << c.a << "," << c.b << "," << c.c << "\","
                              <<"label=\""<< it->second-1 <<"\"];"
                              <<endl;
                }
            }
        }

    }
    out_graph << "}" <<endl;
    out_graph.close();
}

template<class C, class T> Point<C> Adjacency_Graph::get_reference_point(T& t, Mesh<C,T> &mesh, boost::dynamic_bitset<> &unused_vertices)
{
//    Point<C> centroid;
    Box<C> b;
    /// first compute the bounding box of the
    for(int v=0; v<t.get_vertices_num(); v++)
    {
        Vertex<C> &v0 = mesh.get_vertex(abs(t.TV(v)));
        if(!unused_vertices[abs(t.TV(v))-1])
        {
            unused_vertices.set(abs(t.TV(v))-1);
            return v0;
        }
//        cout<<v0<<endl;
        if (v == 0) {
//            centroid = Point<C>(v0.get_dimension());
            b = Box<C>(v0, v0);
        } else {
            b.resize_to_contain(v0);
        }
    }

    dvect min_c = cartesian(b.getMinPoint().getC(1),b.getMinPoint().getC(0));
//    cout<<"min: ";
//    for(unsigned i=0; i<min_c.size(); i++)
//        cout<<min_c[i]<<" ";
//    cout<<endl;
    dvect max_c = cartesian(b.getMaxPoint().getC(1),b.getMaxPoint().getC(0));
//    cout<<"max: ";
//    for(unsigned i=0; i<max_c.size(); i++)
//        cout<<max_c[i]<<" ";
//    cout<<endl;
    dvect midpoint_c; midpoint_c.assign(min_c.size(),0);
    for(unsigned i=0; i<min_c.size(); i++)
    {
        midpoint_c[i] = (min_c[i] + max_c[i]) / 2.0;
    }
//    cout<<"mid: ";
//    for(unsigned i=0; i<midpoint_c.size(); i++)
//        cout<<midpoint_c[i]<<" ";
//    cout<<endl;
    dvect midpoint_g = spherical(midpoint_c);
//    cout<<"mid_deg: ";
//    for(unsigned i=0; i<midpoint_g.size(); i++)
//        cout<<midpoint_g[i]<<" ";
//    cout<<endl;
    Point<C> centroid = Point<C>(midpoint_g);
//    cout<<centroid<<endl;
//    int a; cin>>a;
    return Point<C>(midpoint_g);

//    b.getCenter(centroid);

//    return centroid;
}

/// ===================== ///
/// LOCAL ADJACENCY GRAPH ///
class Local_Adjacency_Graph
{
public:
    template<class C, class T> Local_Adjacency_Graph(string mesh_name, Mesh<C,T> &mesh)
    {
        this->mesh_name = mesh_name;
        get_hsv_palette(palette,mesh.get_implicitly_and_explicitly_encoded_cells_num());
    }

    inline void set_ids(Node_Stellar &n)
    {
        adj_lists.clear();
        global_ids.clear();
        int id=0;
        for(int d=0; d<n.get_num_top_cells_encoded();d++)
        {
            for(RunIterator runIt = n.t_array_begin_iterator(d), runEnd = n.t_array_end_iterator(d); runIt != runEnd; ++runIt)
            {
                pair<int,int> key = make_pair(d,*runIt);
                global_ids.insert(make_pair(key,id));
                id++;
            }
        }
    }

    inline vector<color>& get_palette() { return palette; }

    inline int get_global_id(pair<int,int> pos) { return global_ids[pos]; }
    inline map<pair<int,int>,int> & get_global_ids() { return global_ids; }

    inline void add_pair(int top_id, pair<int,short> &p) { adj_lists[top_id].insert(p); }
    inline map<int,set<pair<int,short> > >& get_adj_lists() { return adj_lists; }

    template<class C, class T> void write_graph_file(Node_Stellar &n, Mesh<C,T> &mesh);

private:
    vector<color> palette;
    map<pair<int,int>,int> global_ids;
    map<int,set<pair<int,short> > > adj_lists;
    string mesh_name;
};

template<class C, class T> void Local_Adjacency_Graph::write_graph_file(Node_Stellar &n, Mesh<C,T> &mesh)
{
    stringstream ss/*, ss2*//*, ss3, ss4*/;
    ss << mesh_name << "_L_" << n.get_v_start() << "_" << n.get_v_end() << "_graph.gv";
//        ss2 << mesh_name << "_d_1_graph.gv";
//        ss3 << mesh_name << "_adj_lists.txt";
//        ss4 << mesh_name << "_adj_lists_01.txt";

    // graphviz compatible graph format
    ofstream out_graph(ss.str());
    out_graph << "graph {" <<endl;
    out_graph << "node [style=bold];" << endl;
    // (1) write node attributes
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIterator runIt = n.t_array_begin_iterator(d), runEnd = n.t_array_end_iterator(d); runIt != runEnd; ++runIt)
        {
            pair<int,int> key = make_pair(d,*runIt);
            T &t = mesh.get_top_cell(d,*runIt);
            color &c = palette[t.get_vertices_num()-1];
            out_graph << global_ids[key] << " [fillcolor=\"" << c.a << "," << c.b << "," << c.c << "\""
                      <<",color=\"" << c.a << "," << c.b << "," << c.c << "\""
                      <<",label=\""<<t.get_vertices_num()-1<<"\"];" <<endl;
        }
    }
    // (2) write edge attribute
    for(auto itt=adj_lists.begin(); itt!=adj_lists.end(); ++itt)
    {
//            out_graph <<i<< " -- { ";
//            for(unsigned j=0; j<adj_lists[i].size(); j++)
        for(set<pair<int,short> >::iterator it=itt->second.begin(); it!=itt->second.end(); ++it)
        {
            color &c = palette[it->second];
            out_graph << itt->first << "--" << it->first << "[color=\"" << c.a << "," << c.b << "," << c.c << "\","
                      <<"label=\""<< it->second-1 <<"\"];"
                      <<endl;
        }
//            out_graph << "};" <<endl;
    }
    out_graph << "}" <<endl;
    out_graph.close();
}

#endif // ADJACENCY_GRAPH_H
