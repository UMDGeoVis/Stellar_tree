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

#ifndef STATISTICS_IG_H
#define STATISTICS_IG_H

#include <vector>
#include <iostream>

using namespace std;

class IG_stats
{
public:
    IG_stats() {}
    IG_stats(int num_simpl_types)
    {
        num_d_nodes.assign(num_simpl_types,0);
        num_d_boundary_arcs.assign(num_simpl_types,0);
        num_d_coboundary_arcs.assign(num_simpl_types,0);
    }

    inline void set_max_stats(int d, long long n, long long b, long long c)
    {
        if(num_d_nodes[d] < n)
            num_d_nodes[d] = n;
        if(num_d_boundary_arcs[d] < b)
            num_d_boundary_arcs[d] = b;
        if(num_d_coboundary_arcs[d] < c)
            num_d_coboundary_arcs[d] = c;
    }

    inline void print_stats() /// for global approach
    {
        for(unsigned d=0; d<num_d_nodes.size(); d++)
        {
            cerr<<"-----"<<d<<"-dimensional stats-----"<<endl;
            cerr<<"num_nodes: "<<num_d_nodes[d]<<endl;
            cerr<<"num_boundary_arcs: "<<num_d_boundary_arcs[d]<<endl;
            cerr<<"num_coboundary_arcs: "<<num_d_coboundary_arcs[d]<<endl;
        }
        cerr<<"-----TOT stats-----"<<endl;
        long long tot_nodes = 0, tot_b = 0, tot_c = 0;
        for(unsigned d=0; d<num_d_nodes.size(); d++)
        {
            tot_nodes += num_d_nodes[d];
            tot_b += num_d_boundary_arcs[d];
            tot_c += num_d_coboundary_arcs[d];
        }
        cerr<<"tot_nodes: "<<tot_nodes<<endl;
        cerr<<"tot_boundary_arcs: "<<tot_b<<endl;
        cerr<<"tot_coboundary_arcs: "<<tot_c<<endl;
    }

    vector<long long> num_d_nodes;
    vector<long long> num_d_boundary_arcs;
    vector<long long> num_d_coboundary_arcs;
};

#endif // STATISTICS_IG_H
