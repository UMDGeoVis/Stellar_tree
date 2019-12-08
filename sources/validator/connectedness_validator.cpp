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

#include "connectedness_validator.h"

void Connectedness_Validator::set_flag(int id, int &flag, Connectedness_Validator_stats &stats)
{
    if(stats.cc[id -1] == 0)
    {
        stats.cc[id -1] = flag;
        stats.connected_components[flag].push_back(id);
    }
    else if(stats.cc[id -1] != flag) // found two path that in reality are the same
    {
        // If I have two path..
        // I take the smaller id and I unify the two in the global array
        //
        // this approach seems to be safer than fixing this intermediate connected components in the end
        //

        boost::tuples::tuple<int,int> minmax = boost::minmax(stats.cc[id -1], flag);

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
}

void Connectedness_Validator::visit_skeleton(Node_Stellar &n, connVal_skeleton& skeleton, Connectedness_Validator_stats &stats)
{
    iqueue q;
    int flag;

    // start new paths from the current leaf
    // (if every edge is linked with something external we should do nothing in this step)
    for(int i=0; i<skeleton.get_edges_num(); i++)
    {
        if(!skeleton.is_visited(i))
        {
            if(stats.cc[skeleton.get_edge(i)[0] -1] == 0 && stats.cc[skeleton.get_edge(i)[1] -1] == 0)
            {
                flag = skeleton.get_edge(i)[0];

                // every time we start a graph visit we have a new a connected components,
                // we flag each connected components with an integer, i.e. the index of the first vertex
                stats.connected_components.insert(make_pair(flag,deque<int>()));
            }
            else
            {
                if(stats.cc[skeleton.get_edge(i)[0] -1] != 0)
                    flag = stats.cc[skeleton.get_edge(i)[0] -1];
                else if(stats.cc[skeleton.get_edge(i)[1] -1] != 0)
                    flag = stats.cc[skeleton.get_edge(i)[1] -1];
            }

            skeleton.set_visited(i);
            q.push(i);

            // visit the queue
            visit_edge_queue(n,skeleton,q,flag,stats);

            // if we visit all the edges then we have done within this leaf
            if(skeleton.all_visited())
            {
                return;
            }
        }
    }
}

void Connectedness_Validator::visit_edge_queue(Node_Stellar &n, connVal_skeleton &skeleton, iqueue &q, int &flag, Connectedness_Validator_stats &stats)
{
    int edge_id;
    ivect edge;
    while(!q.empty())
    {
        edge_id = q.front();
        q.pop();
        edge = skeleton.get_edge(edge_id);

        for(unsigned j=0; j<edge.size(); j++)
        {
            Connectedness_Validator::set_flag(edge[j],flag,stats);

            if(n.indexes_vertex(edge[j]))
            {
                ivect &ve = skeleton.get_ve(edge[j]-n.get_v_start());
                for(ivect_iter it=ve.begin(); it!=ve.end(); ++it)
                {
                    if(!skeleton.is_visited(*it))
                    {
                        q.push(*it);
                        skeleton.set_visited(*it);
                    }
                }
            }
        }
    }
}

void Connectedness_Validator::visit_VVs(Node_Stellar &n, leaf_VV &vvs, Connectedness_Validator_stats &stats)
{
    iqueue q;
    int flag;
    boost::dynamic_bitset<> is_visited = boost::dynamic_bitset<>(n.get_real_v_array_size(),0);

    for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& v_id = itPair.first;

        if(!is_visited[*v_id-n.get_v_start()])
        {
            if(stats.cc[*v_id-1]==0)
            {
                flag = *v_id;
                stats.cc[*v_id-1]=flag;
                stats.connected_components.insert(make_pair(flag,deque<int>()));
            }
            else
            {
                flag = stats.cc[*v_id-1];
            }

            is_visited.set(*v_id-n.get_v_start());

            for(iset_iter it=vvs[*v_id-n.get_v_start()].begin(); it!=vvs[*v_id-n.get_v_start()].end(); ++it)
            {
                Connectedness_Validator::set_flag(*it,flag,stats);

                if(n.indexes_vertex(*it) && !is_visited[*it-n.get_v_start()])
                {
                    q.push(*it);
                    is_visited.set(*it-n.get_v_start());
                }
            }

            while(!q.empty())
            {
                int v = q.front();
                q.pop();

                for(iset_iter it2=vvs[v-n.get_v_start()].begin(); it2!=vvs[v-n.get_v_start()].end(); ++it2)
                {
                    Connectedness_Validator::set_flag(*it2,flag,stats);

                    if(n.indexes_vertex(*it2) && !is_visited[*it2-n.get_v_start()])
                    {
                        q.push(*it2);
                        is_visited.set(*it2-n.get_v_start());
                    }
                }
            }

            if(is_visited.count() == is_visited.size())
                return;
        }
    }
}

void Connectedness_Validator::visit_local_adjacencies(connVal_d_1_adj& local, Connectedness_Validator_stats &stats)
{
    int flag;
    iqueue top_queue;
    Timer time;

    for(int i=1; i<=local.get_index_array_size(); i++) // <========= PAY ATTENTION
    {
        if(!local.is_visited(i))
        {
            if(stats.cc[local.get_index(i) -1] == 0)
            {
                flag = local.get_index(i);

                /// every time we start a graph visit we have a new a connected components,
                /// we flag each connected components with an integer, i.e. the index of the first vertex
                stats.connected_components.insert(make_pair(flag,deque<int>()));
            }
            else
            {
                flag = stats.cc[local.get_index(i) -1];
            }

            local.set_visited( i );
            top_queue.push( i ); // here we need the local indexing

            time.start();

            /// visit the queue
            visit_top_queue(local,top_queue,flag,stats);

            time.stop();
            stats.time_set_flags += time.get_elapsed_time();

            /// if we visit all the edges then we have done within this leaf
            if(local.all_visited())
            {
                return;
            }
        }
    }
}

void Connectedness_Validator::visit_top_queue(connVal_d_1_adj& local, iqueue &q, int &flag, Connectedness_Validator_stats &stats)
{
    while(!q.empty())
    {
        int t = q.front();
        set_flag(local.get_index(t),flag,stats); // here we have to
        q.pop();

        Top_nmAdj &t_adj = local.get_adj(t);

        for(int j=0; j<t_adj.get_adj_num(); j++)
        {
            if(t_adj.get_adj_num(j) > 1)
                stats.is_pseudo_manifold = false;
            for(int w=0; w<t_adj.get_adj_num(j); w++)
            {
                if(!local.is_visited(t_adj.get_adj(j,w)))
                {
                    q.push(t_adj.get_adj(j,w));
                    local.set_visited(t_adj.get_adj(j,w));
                }
            }
        }
    }
}

// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// IA* code ////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////


void Connectedness_Validator::visit_skeleton(connVal_skeleton& skeleton, Connectedness_Validator_stats &stats)
{
    iqueue q;
    int flag;

    // start new paths from the current leaf
    // (if every edge is linked with something external we should do nothing in this step)
    for(int i=0; i<skeleton.get_edges_num(); i++)
    {
        if(!skeleton.is_visited(i))
        {
            if(stats.cc[skeleton.get_edge(i)[0] -1] == 0 && stats.cc[skeleton.get_edge(i)[1] -1] == 0)
            {
                flag = skeleton.get_edge(i)[0];

                // every time we start a graph visit we have a new a connected components,
                // we flag each connected components with an integer, i.e. the index of the first vertex
                stats.connected_components.insert(make_pair(flag,deque<int>()));
            }
            else
            {
                if(stats.cc[skeleton.get_edge(i)[0] -1] != 0)
                    flag = stats.cc[skeleton.get_edge(i)[0] -1];
                else if(stats.cc[skeleton.get_edge(i)[1] -1] != 0)
                    flag = stats.cc[skeleton.get_edge(i)[1] -1];
            }

            skeleton.set_visited(i);
            q.push(i);

            // visit the queue
            visit_edge_queue(skeleton,q,flag,stats);

            // if we visit all the edges then we have done within this leaf
            if(skeleton.all_visited())
            {
                return;
            }
        }
    }
}

void Connectedness_Validator::visit_edge_queue(connVal_skeleton &skeleton, iqueue &q, int &flag, Connectedness_Validator_stats &stats)
{
    int edge_id;
    ivect edge;
    while(!q.empty())
    {
        edge_id = q.front();
        q.pop();
        edge = skeleton.get_edge(edge_id);

        for(unsigned j=0; j<edge.size(); j++)
        {
            Connectedness_Validator::set_flag(edge[j],flag,stats);

            ivect &ve = skeleton.get_ve(edge[j]-1);
            for(ivect_iter it=ve.begin(); it!=ve.end(); ++it)
            {
                if(!skeleton.is_visited(*it))
                {
                    q.push(*it);
                    skeleton.set_visited(*it);
                }
            }
        }
    }
}
