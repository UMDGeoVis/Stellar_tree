#include "adjacency_graph.h"
#include "utilities/container_utilities.h"

#include <tuple>

void Adjacency_Graph::set_blockers_adjacencies()
{
    for(int d=0; d<get_blockers_num();d++)
    {
        for(int i=1; i<=get_blockers_num(d);i++)
        {
            Top_Simplex &blocker = get_blocker(d,i);
            extract_d1_d2_adjacent_blockers(blocker,i,d);
        }
    }
}

void Adjacency_Graph::extract_d1_d2_adjacent_blockers(Top_Simplex &t, int t_id, int d)
{
    ivect b_v; t.get_positive_sorted_vertices(b_v);
    int global_t_id = get_blocker_global_id(d,t_id);
    ivect b2_v;

    /// first check d-1 adjacencies
    for(int i=t_id+1; i<=get_blockers_num(d);i++)
    {
        Top_Simplex& b2 = get_blocker(d,i);
        int b2_id = get_blocker_global_id(d,i);
        b2.get_positive_sorted_vertices(b2_v);
        intersect_containers(b2_v,b_v);

        if((int)b2_v.size() == t.get_vertices_num()-1 ) /// we have a d-1 adjacency
        {
            // we just push each edge once
            pair<int,short> p;
            p.first = b2_id;
            p.second = b2_v.size();
            add_pair(global_t_id,p);
        }
    }

    if(d>0) /// at least top 2-simplices
    {
        int d2 = d-1;
        for(int i=1; i<=get_blockers_num(d2);i++)
        {
            Top_Simplex& b2 = get_blocker(d2,i);
            int b2_id = get_blocker_global_id(d2,i);
            b2.get_positive_sorted_vertices(b2_v);
            intersect_containers(b2_v,b_v);

            if((int)b2_v.size() == b2.get_vertices_num()-1 ) /// we have a d-2 adjacency with respect to the higher dimensional blocker
            {
                // we just push each edge once
                pair<int,short> p;
                p.first = b2_id;
                p.second = b2_v.size();
                add_pair(global_t_id,p);
            }
        }
    }
}

void Adjacency_Graph::write_blockers_file(string mesh_name)
{
    stringstream ss/*, ss2*//*, ss3, ss4*/;
    ss << mesh_name << "_blocker_graph.gv";
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
            Top_Simplex &t = blockers[i][j];
            color &c = palette[t.get_vertices_num()-1];
            out_graph << blockers_global_ids[i][j] << " [fillcolor=\"" << c.a << "," << c.b << "," << c.c << "\""
                      <<",color=\"" << c.a << "," << c.b << "," << c.c << "\""
                      <<",label=\""<<t.get_vertices_num()-1<<"\",shape=box];" <<endl;
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

void Adjacency_Graph::compute_adj_graph_stats()
{
//    vector<tuple<int,double,int> > tmp;
//    tmp.assign(2,make_tuple(INT_MAX,0,0)); /// lite adjacency (i.e. d-1 and d-2 adjacencies)
//    vector<tuple<int,double,int> > graph_stats;
//    graph_stats.assign(this->global_ids.size(),make_tuple(INT_MAX,0,0));

    int entries = 0;
    for(unsigned t=0; t<this->adj_sets.size(); t++)
    {
        entries += adj_sets[t].size();
    }
    cout<<"  nodes: "<<adj_sets.size()<<" -- entries-summation: "<<entries<<endl;

    for(unsigned d=0; d<this->global_ids.size(); d++)
    {
        if(this->global_ids[d].size() > 0)
            cout<<"  top-"<<d+1<<"-simpl "/*<<endl*/;

        map<int,int> min;
        map<int,double> avg;
        map<int,int> max;

        for(unsigned i=0; i<this->global_ids[d].size(); i++)
        {
            map<int,int> counters;
//            cout<<i<<endl;
            int t_pos = global_ids[d][i];

            if(hided_tops.find(t_pos)!=hided_tops.end())
                continue;

            set<pair<int,short> > &adjs = adj_sets[t_pos];
//            cout<<"size "<<adjs.size()<<endl;
            for(set<pair<int,short> >::iterator it=adjs.begin(); it!=adjs.end(); ++it)
            {
                const pair<int,short> &p = *it;
                counters[p.second]++;
//                cout<<p.first<<" "<<p.second<<endl;
//                int a; cin>>a;
            }

            for(map<int,int>::iterator it=counters.begin(); it!=counters.end(); ++it)
            {
                if(min.find(it->first)!=min.end())
                {
                    if(min[it->first] > it->second)
                        min[it->first] = it->second;
                    if(max[it->first] < it->second)
                        max[it->first] = it->second;
                    avg[it->first] += it->second;
                }
                else
                {
                    min[it->first] = it->second;
                    avg[it->first] = it->second;
                    max[it->first] = it->second;
                }
            }
        }


        for(map<int,double>::iterator it=avg.begin(); it!=avg.end(); ++it)
        {
            it->second /= (double)global_ids[d].size();
            cout<<"[ "<<it->first-1<<"-adj "<<min[it->first]<<" "<<it->second<<" "<<max[it->first]<<" ] "/*<<endl*/;
        }
        if(this->global_ids[d].size() > 0)
            cout << endl;
    }
}

dvect Adjacency_Graph::cartesian(double lat, double lon)
{
    dvect ret;
    double c = cos(lat*PI/180.0) * cos(lon*PI/180.0);
    ret.push_back(c);
    c = cos(lat*PI/180.0) * sin(lon*PI/180.0);
    ret.push_back(c);
    c = sin(lat*PI/180.0);
    ret.push_back(c);
    return ret;
}

dvect Adjacency_Graph::spherical(dvect &c)
{
    dvect ret; ret.assign(2,0);
    double r = sqrt(c[0]*c[0] + c[1]*c[1]);
//    cout<<r<<endl;

    if (r == 0)
    {
        if (c[2] > 0)
        {
            ret[1] = 90;
            ret[2] = 0;
        }
        else if (c[2] < 0)
        {
            ret[1] = -90;
            ret[0] = 0;
        }
        // else return Undefined // (x,y,z) == (0,0,0)
        return ret;
    }
    else
    {
        ret[1] = atan2(c[2], r)*180.0/PI; // latitude
        ret[0] = atan2(c[1], c[0])*180.0/PI; // longitude
        return ret;
//        return (, ) // atan2 must return *degrees*
    }
}
