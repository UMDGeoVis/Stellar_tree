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

#include "halfedge_generator.h"

void HalfEdge_Generator::print_half_edges()
{
    for(unsigned i=0; i<half_edges.size(); i++)
    {
        cerr<<i<<"]";
        half_edges[i].print();
    }
}

void HalfEdge_Generator::print_stats()
{
    cerr<<"max num half edges: "<<max_half_edges<<endl;
    cerr<<"max num half vertices: "<<max_half_vertices<<endl;
    cerr<<"max num half faces: "<<max_half_faces<<endl;
    cerr<<"max edges list size: "<<max_edges_list_size<<endl;
}

void HalfEdge_Generator::fix_local_half_edges_ordering()
{
    boost::dynamic_bitset<> visited(half_edges.size());

    iqueue he_queue;
    he_queue.push(0);
    visited.set(0);

    unsigned count = 1;

    /// the first half edge gives the convention on the ordering

    while(!he_queue.empty() || count < half_edges.size())
    {
        /// case [0] not all the simplices are in the same connected component..
        /// this should happen only in the local approach!
        if(he_queue.empty() && count < half_edges.size())
        {
            for (boost::dynamic_bitset<>::size_type i = 0; i < visited.size(); ++i)
            {
                if(!visited.test(i))
                {
                    he_queue.push(i);
                    visited.set(i);
                    count++;
                    break;
                }
            }
        }

        fix_half_edges_ordering(he_queue,visited,count);
    }
}

void HalfEdge_Generator::fix_global_half_edges_ordering()
{
    boost::dynamic_bitset<> visited(half_edges.size());

    iqueue he_queue;
    he_queue.push(0);
    visited.set(0);

    unsigned count = 1; /// unused.. added for compatibility

    /// the first half edge gives the convention on the ordering
    while(!he_queue.empty())
    {
        fix_half_edges_ordering(he_queue,visited,count);
    }
}

void HalfEdge_Generator::fix_half_edges_ordering(iqueue &he_queue, boost::dynamic_bitset<> &visited, unsigned &count)
{
    int current = he_queue.front();
    he_queue.pop();
    HalfEdge &he = half_edges[current];

    if(he.get_opposite_half_edge() != -1 && !visited[he.get_opposite_half_edge()])
    {
        HalfEdge &opposite = half_edges[he.get_opposite_half_edge()];

        if(he.get_pointed_vertex() == opposite.get_pointed_vertex())
        {
            /// fix the opposite half edge ordering
            /// reverse the ordering of all the half edges with the same reference face
            int local_current_he = he.get_opposite_half_edge();
            int local_next_he = opposite.get_next_half_edge();

            do
            {
                HalfEdge &next_he = half_edges[local_next_he];

                int tmp_next = next_he.get_next_half_edge();
                next_he.set_next_half_edge(local_current_he);

                local_current_he = local_next_he;
                local_next_he = tmp_next;
            }
            while(he.get_opposite_half_edge() != local_next_he);

            opposite.set_next_half_edge(local_current_he);
        }

        visited.set(he.get_opposite_half_edge());
        he_queue.push(he.get_opposite_half_edge());
        count++;
    }

    if(!visited[he.get_next_half_edge()])
    {
        visited.set(he.get_next_half_edge());
        he_queue.push(he.get_next_half_edge());
        count++;
    }
}

void HalfEdge_Generator::test_local_half_edge_structure()
{
    boost::dynamic_bitset<> visited(half_edges.size());

    iqueue he_queue;
    he_queue.push(0);
    visited.set(0);

    unsigned count = 0;
    int conn_components = 1;

    /// the first half edge gives the convention on the ordering

    while(!he_queue.empty() || count < half_edges.size())
    {
        /// case [0] not all the simplices are in the same connected component..
        /// this should happen only in the local approach!
        if(he_queue.empty() && count < half_edges.size())
        {
            for (boost::dynamic_bitset<>::size_type i = 0; i < visited.size(); ++i)
                if(!visited.test(i))
                {
                    he_queue.push(i);
                    conn_components++;
                    break;
                }
        }

        test_half_edge_structure(he_queue,visited);

        count = visited.count();
    }

    print_test_result(conn_components,visited);
}

void HalfEdge_Generator::test_global_half_edge_structure()
{
    boost::dynamic_bitset<> visited(half_edges.size());

    iqueue he_queue;
    he_queue.push(0);
    visited.set(0);

    /// the first half edge gives the convention on the ordering

    while(!he_queue.empty())
    {
        test_half_edge_structure(he_queue,visited);
    }

    print_test_result(1,visited);
}

void HalfEdge_Generator::test_half_edge_structure(iqueue &he_queue, boost::dynamic_bitset<> &visited)
{
    int current = he_queue.front();
    he_queue.pop();
    HalfEdge &he = half_edges[current];

    if(he.get_opposite_half_edge() != -1 && !visited[he.get_opposite_half_edge()])
    {
        visited.set(he.get_opposite_half_edge());
        he_queue.push(he.get_opposite_half_edge());
    }

    if(!visited[he.get_next_half_edge()])
    {
        visited.set(he.get_next_half_edge());
        he_queue.push(he.get_next_half_edge());
    }
}

void HalfEdge_Generator::print_test_result(int conn_components, boost::dynamic_bitset<> &visited)
{
    if(visited.count() == visited.size())
    {
        cerr<<"[test_half_edge_structure] PASSED"<<endl;
        if(conn_components > 1)
        {
            cerr<<"[log] number of connected components: "<<conn_components<<endl;
        }
    }
    else
    {
        cerr<<"[test_half_edge_structure] ERROR"<<endl;
        print_half_edges();
        int counter = 0;
        for(size_t i=0; i<visited.size(); i++)
            if(!visited[i])
            {
                counter++;
                cerr<<i<<" "; half_edges[i].print();
            }
        cout<<"not visited "<<counter<<endl;
        int a; cin>>a;
    }
}
