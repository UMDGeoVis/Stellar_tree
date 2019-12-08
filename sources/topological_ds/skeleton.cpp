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

#include "skeleton.h"
#include <boost/algorithm/minmax.hpp>
#include "connectedness_validator_stats.h"

void skeleton_graph::print_skeleton()
{
    for(unsigned v=0; v<graphs.size(); v++)
    {
        if(graphs[v].size() > 0)
        {
            cout<<v+1<<"] ";
            for(iset_iter it=graphs[v].begin(); it!=graphs[v].end(); ++it)
                cout<<*it<<" ";
            cout<<endl;
        }
    }
}

void skeleton_graph::print_skeleton(int pos)
{
    if(graphs[pos].size() > 0)
    {
        cout<<pos+1<<"] ";
        for(iset_iter it=graphs[pos].begin(); it!=graphs[pos].end(); ++it)
            cout<<*it<<" ";
        cout<<endl;
    }
}

void skeleton_graph::compute_skeleton_stats()
{
    long long entries = 0;
    int min_entries = INT_MAX;
    int max_entries = INT_MIN;
    for(unsigned v=0; v<graphs.size(); v++)
    {
        int l_size = adjacent_vertices_num(v);
        if(min_entries > l_size)
            min_entries = l_size;
        if(max_entries < l_size)
            max_entries = l_size;
        entries += l_size;
    }
    cout<<"[STATS] skeleton nodes: "<<graphs.size()<<" -- entries-summation: "<<entries<<endl;
    cout<<"  neighbors -- min: "<<min_entries<<" -- avg: "<<entries/(double)graphs.size()<<" -- max: "<<max_entries<<endl;
}

void skeleton_graph::set_skeleton_stats(unsigned &maxV, long long &maxE)
{
    long long entries = 0;
    for(unsigned v=0; v<graphs.size(); v++)
        entries += adjacent_vertices_num(v);
    if(maxV < graphs.size())
        maxV = graphs.size();
    if(maxE < entries)
        maxE = entries;
}

void skeleton_graph::save_skeleton_graph(string filename)
{
    stringstream ss;
    ss << filename << "_1skeleton_edges.txt";
    ofstream output(ss.str().c_str());

    /// count the edges
    int num_edges = 0;
    for(vector<iset>::iterator it=graphs.begin(); it!=graphs.end(); ++it)
    {
        num_edges += it->size();
    }

    output<<graphs.size()<<" "<<num_edges<<endl;

    for(unsigned v=0; v<graphs.size(); v++)
    {
        for(iset_iter its=begin(v); its!=end(v); ++its)
        {
            //                output<<v+1<<" "<<*its<<endl;
            output<<v<<" "<<*its-1<<endl; /// per vis di ciccio
        }
    }
    output.close();
}

void skeleton_graph::write_gv_file(string filename, ivect &v_orig)
{
    stringstream ss;
    ss << filename << "_global_skeleton.gv";

    ofstream output(ss.str().c_str());
    output << "graph {" <<endl;

    for(unsigned v=0; v<graphs.size(); v++)
    {
        for(iset_iter its=begin(v); its!=end(v); ++its)
        {
            output<<v_orig[v]<<"--"<<v_orig[*its-1]<<endl;
        }
    }

    output << "}" <<endl;
    output.close();
}

void skeleton_graph::check_skeleton_graph()
{
    int empty_lists = 0;
    for(vector<iset>::iterator it=graphs.begin(); it!=graphs.end(); ++it)
    {
        if(it->size() == 0)
            empty_lists++;
    }
    cout<<"[STATS] EMPTY ADJACENCY LISTS: "<<empty_lists<<endl;
    iqueue q;
    int flag;
    Connectedness_Validator_stats stats; stats.cc.assign(graphs.size(),0);
    boost::dynamic_bitset<> is_visited = boost::dynamic_bitset<>(graphs.size(),0);
    for(int i = 0; i<graphs.size(); i++)
    {
//        cout<<"ci entro"<<endl;
        if(!is_visited[i])
        {
            if(stats.cc[i]==0)
            {
                flag = i;
                stats.cc[i]=flag;
                stats.connected_components.insert(make_pair(flag,deque<int>()));
            }
            else
            {
                flag = stats.cc[i];
            }

            is_visited.set(i);

            for(iset_iter it=graphs[i].begin(); it!=graphs[i].end(); ++it)
            {
                if(stats.cc[*it -1] == 0)
                {
                    stats.cc[*it -1] = flag;
                    stats.connected_components[flag].push_back(*it);
                }
                else if(stats.cc[*it -1] != flag) // found two path that in reality are the same
                {
                    // If I have two path..
                    // I take the smaller id and I unify the two in the global array
                    //
                    // this approach seems to be safer than fixing this intermediate connected components in the end
                    //

                    boost::tuples::tuple<int,int> minmax = boost::minmax(stats.cc[i], flag);

                    deque<int> &mind = stats.connected_components[minmax.get<0>()];
                    deque<int> &maxd = stats.connected_components[minmax.get<1>()];
                    mind.insert(mind.end(),maxd.begin(),maxd.end());

                    for(deque<int>::iterator it=maxd.begin(); it!=maxd.end(); ++it)
                    {
                        stats.cc[*it -1] = minmax.get<0>();
                    }

                    stats.connected_components.erase(minmax.get<1>());
                    flag = minmax.get<0>();
                }

                if(!is_visited[*it-1])
                {
                    q.push(*it);
                    is_visited.set(*it-1);
                }
            }

            while(!q.empty())
            {
                int v = q.front();
                q.pop();

                for(iset_iter it2=graphs[v-1].begin(); it2!=graphs[v-1].end(); ++it2)
                {
//                    Connectedness_Validator::set_flag(*it2,flag,stats);
                    if(stats.cc[*it2 -1] == 0)
                    {
                        stats.cc[*it2 -1] = flag;
                        stats.connected_components[flag].push_back(*it2);
                    }
                    else if(stats.cc[*it2 -1] != flag) // found two path that in reality are the same
                    {
                        // If I have two path..
                        // I take the smaller id and I unify the two in the global array
                        //
                        // this approach seems to be safer than fixing this intermediate connected components in the end
                        //

                        boost::tuples::tuple<int,int> minmax = boost::minmax(stats.cc[i], flag);

                        deque<int> &mind = stats.connected_components[minmax.get<0>()];
                        deque<int> &maxd = stats.connected_components[minmax.get<1>()];
                        mind.insert(mind.end(),maxd.begin(),maxd.end());

                        for(deque<int>::iterator it=maxd.begin(); it!=maxd.end(); ++it)
                        {
                            stats.cc[*it -1] = minmax.get<0>();
                        }

                        stats.connected_components.erase(minmax.get<1>());
                        flag = minmax.get<0>();
                    }

                    if(!is_visited[*it2-1])
                    {
                        q.push(*it2);
                        is_visited.set(*it2-1);
                    }
                }
            }

            if(is_visited.count() == is_visited.size())
                break;
        }
    }
    stats.print_validation_result();
}
