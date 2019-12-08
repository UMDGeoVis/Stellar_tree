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

#ifndef IG_H
#define IG_H

#include <vector>

#include "ig_node.h"

using namespace std;

class IG
{
public:
    IG(int num_top_type)
    {
        nodes.assign(num_top_type,set<IG_node*,IG_node::IG_node_compare>());
    }

    /// when I add I have not the co-boundary set, so I do not need to update the co-boundary of the element in the boundary of node
    inline void add_node(int d, IG_node *&node)
    {
        pair<IG_node::IG_nodes_iterator,bool> ret = nodes[d].insert(node);
        if(!ret.second) /// already into
        {
            delete node;
            node = *ret.first;
        }
    }

    inline void print_ig()
    {
        for(unsigned d=0; d<nodes.size(); d++)
        {
            cout<<"--------------"<<d<<"-dimensional cells--------------"<<endl;
            for(IG_node::IG_nodes_iterator it=nodes[d].begin(); it!=nodes[d].end(); ++it)
                cout<<**it<<endl;
        }
    }

    inline void print_stats(IG_stats &stats) /// for global approach
    {
        for(unsigned d=0; d<nodes.size(); d++)
        {
            for(IG_node::IG_nodes_iterator it=nodes[d].begin(); it!=nodes[d].end(); ++it)
                (*it)->get_stats(d,stats);
        }
        stats.print_stats();
    }

    inline void update_stats(IG_stats &stats) /// for local approach
    {
        for(unsigned d=0; d<nodes.size(); d++)
        {
            long long n = 0, b = 0, c = 0;
            for(IG_node::IG_nodes_iterator it=nodes[d].begin(); it!=nodes[d].end(); ++it)
                (*it)->update_stats(n,b,c);
            stats.set_max_stats(d,n,b,c);
        }
    }

    inline void clear_IG()
    {
        for(unsigned i=0; i<nodes.size(); i++)
        {
            for(IG_node::IG_nodes_iterator it=nodes[i].begin(); it!=nodes[i].end(); )
            {
                IG_node* n = *it;
                nodes[i].erase(it++);
                delete n;
            }
        }
    }

private:
    vector<IG_node::IG_nodes_set> nodes;
};



#endif // IG_H
