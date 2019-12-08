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

#ifndef IG_NODE_H
#define IG_NODE_H

#include "utilities/basic_wrappers.h"

#include "topological_ds_generator/ig_statistics.h"

using namespace std;

class IG_node
{
public:
    struct IG_node_compare
    {
      bool operator() (const IG_node* lhs, const IG_node* rhs) const
      {return *lhs<*rhs;}
    };

    typedef set<IG_node*,IG_node_compare> IG_nodes_set;
    typedef IG_nodes_set::iterator IG_nodes_iterator;

    IG_node() { this->vertex = -1;}
    virtual ~IG_node()
    {
        this->boundary.clear();
        this->coboundary.clear();
    }

    inline void set_vertex(int v) { this->vertex = v; }
    inline int get_vertex() const { return this->vertex; }

    inline void set_boundary(IG_node *node) { this->boundary.insert(node); }
    inline void set_coboundary(IG_node *node) { this->coboundary.insert(node); }

    inline IG_nodes_set& get_boundary() { return this->boundary; }
    inline IG_nodes_set& get_coboundary() { return this->coboundary; }

    inline void get_vertices(iset &vertices) const
    {
        if(this->vertex != -1)
            vertices.insert(this->vertex);
        else
        {
            for(IG_nodes_iterator it=this->boundary.begin(); it!=this->boundary.end(); ++it)
                (*it)->get_vertices(vertices);
        }
    }

    inline friend bool operator<(const IG_node &q, const IG_node &p)
    {
        iset qvert, pvert;
        q.get_vertices(qvert);
        p.get_vertices(pvert);
        return qvert < pvert;
    }
    inline bool operator==(const IG_node &s) const
    {
        iset qvert, pvert;
        get_vertices(qvert);
        s.get_vertices(pvert);
        return qvert == pvert;
    }
    inline bool operator!=(const IG_node &s) const { return !(*this==s); }
    inline friend std::ostream& operator<<(std::ostream& out, const IG_node& p)
    {
        p.print_vertices(out); out<<endl;

        if(p.boundary.size() > 0)
        {
            out<<"boundary --> ";
            for(IG_nodes_iterator it=p.boundary.begin(); it!=p.boundary.end(); ++it)
                (*it)->print_vertices(out);
            out<<endl;
        }
        if(p.coboundary.size() > 0)
        {
            out<<"coboundary --> ";
            for(IG_nodes_iterator it=p.coboundary.begin(); it!=p.coboundary.end(); ++it)
            {
                (*it)->print_vertices(out);
            }
            out<<endl;
        }
        return out;
    }

    inline void print_vertices(std::ostream& out) const
    {
        if(this->vertex != -1)
            out<< "V[ "<<this->vertex<<" ] ";
        else
        {
            iset vertices;
            this->get_vertices(vertices);
            out<<"[ ";
            for(iset_iter it=vertices.begin(); it!=vertices.end(); ++it)
                out<<*it<<" ";
            out<<"] ";
            /// for debug
            if(vertices.size() == 0)
            {
                out<<"EMPTY "<<this<<endl;
                int a; cin>>a;
            }
        }
    }

    inline void get_stats(int d, IG_stats &stats)
    {
        stats.num_d_nodes[d]++;
        stats.num_d_boundary_arcs[d] += this->boundary.size();
        stats.num_d_coboundary_arcs[d] += this->coboundary.size();
    }

    inline void update_stats(long long &n, long long &b, long long &c)
    {
        n++;
        b += this->boundary.size();
        c += this->coboundary.size();
    }

private:
    int vertex; /// initialized only on the 0-cells
    IG_nodes_set boundary;
    IG_nodes_set coboundary;
};

#endif // IG_NODE_H
