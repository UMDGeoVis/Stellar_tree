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

#ifndef SKELETON_H
#define SKELETON_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <climits>
#include "utilities/basic_wrappers.h"
#include "stellar_tree/mesh.h"

using namespace std;

class skeleton_graph
{
public:
    skeleton_graph() {}
    skeleton_graph(double eps, int num_vertices)
    {
        this->epsilon = eps;
        graphs.assign(num_vertices,iset());
    }
    skeleton_graph(int num_vertices)
    {
        this->epsilon = 0;
        graphs.assign(num_vertices,iset());
    }

    inline double get_epsilon() { return epsilon; }

    inline void insert(int pos, int pivot) { graphs[pos].insert(pivot); }
    template<class C> inline void insert(int pos, C &c) { graphs[pos].insert(c.begin(),c.end()); }
//    template<class C> inline void conditional_insert(int pos, C &c)
//    {
//        ///
//        for(typename C::iterator it=c.begin(); it!=c.end(); ++it)
//        {
//            if(pos < *it)
//                graphs[pos].insert(*it);
//        }
//    }
    template<class C> inline void conditional_insert(int pos, C &c, int offset)
    {
        ///
        for(typename C::iterator it=c.begin(); it!=c.end(); ++it)
        {
            if(pos < *it)
                graphs[pos-offset].insert(*it);
        }
    }
    inline iset& get_graph(int pos) { return graphs[pos]; }
    inline iset_iter begin(int pos) { return graphs[pos].begin(); }
    inline iset_iter end(int pos) { return graphs[pos].end(); }
    inline iset_iter find(int pos, int pivot) { return graphs[pos].find(pivot); }
    inline unsigned adjacent_vertices_num(int pos) { return graphs[pos].size(); }

    inline vector<iset>::iterator begin() { return graphs.begin(); }
    inline vector<iset>::iterator end() { return graphs.end(); }
    inline int graph_size() { return graphs.size(); }
    inline int graph_size(int pos) { return graphs[pos].size(); }

    inline void reset() { graphs.clear(); }

    void print_skeleton();
    void print_skeleton(int pos);
    template<class N> void print_skeleton(N &n);

    void compute_skeleton_stats();
    void set_skeleton_stats(unsigned &maxV, long long &maxE);

    void save_skeleton_graph(string filename);
    void write_gv_file(string filename, ivect &v_orig);
    template<class C, class T> void write_gv_file(string filename, ivect &v_orig, Mesh<C,T> &mesh);
    template<class C, class T> void write_cesium_file(string filename, ivect &v_orig, vector<pair<int,int> > &v_stats, Mesh<C,T> &mesh);

    void check_skeleton_graph();

private:
    vector<iset> graphs;
    double epsilon; /// refers to the maximum distance between two vertices
};

template<class C, class T> void skeleton_graph::write_gv_file(string filename, ivect &v_orig, Mesh<C,T> &mesh)
{
    stringstream ss;
    ss << filename << "_skeleton.gv";

    ofstream output(ss.str().c_str());
    output << "graph {" <<endl;

    for(unsigned v=0; v<graph_size(); v++)
    {
        output << v_orig[v] <<" [pos=\""<<mesh.get_vertex(v+1).getC(0)<<","<<mesh.get_vertex(v+1).getC(1)<<"!\"];" <<endl;
    }

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

template<class C, class T> void skeleton_graph::write_cesium_file(string filename, ivect &v_orig, vector<pair<int, int> > &v_stats, Mesh<C,T> &mesh)
{
    stringstream ss;
    ss << filename << "_skeleton.json";

    ofstream output(ss.str().c_str());
    output << "{ \"nodes\":[" <<endl;
    for(int v=0; v<graph_size(); v++)
    {
        output << "{ \"id\":" << v_orig[v] <<",\"long\":"<<mesh.get_vertex(v+1).getC(0)<<",\"lat\":"<<mesh.get_vertex(v+1).getC(1)
               <<",\"sizeVTop\":"<<v_stats[v].first<<",\"sizeVE\":"<<v_stats[v].second<<"}";
        if(v < graph_size()-1)
            output<<","<<endl;
        else
            output<<endl;
    }
    output << "]," << endl << "\"arcs\":[" << endl;
    for(unsigned v=0; v<graphs.size(); v++)
    {
//        if(this->graph_size(v) > 0) /// we write only the adjacency lists greater than 0
        {
            output << "{ \"id\":" << v_orig[v] <<",\"len\":"<<this->graph_size(v)<<",\"adjs\": [";
            int counter = 0;
            for(iset_iter its=begin(v); its!=end(v); ++its)
            {
                output<<v_orig[*its-1];
                if(counter < graph_size(v)-1)
                    output << ",";
                counter++;
            }
            output << "]}";
            if(v < graph_size()-1)
                output<<","<<endl;
            else
                output<<endl;
        }
    }

    output << "]}" <<endl;
    output.close();
}

//class local_skeleton_graph
//{
//public:
//    local_skeleton_graph() {}

//    inline void insert(int pos, int pivot) { graphs[pos].insert(pivot); }
//    template<class C> inline void insert(int pos, C c) { graphs[pos].insert(c.begin(),c.end()); }
//    inline iset& get_graph(int pos) { return graphs[pos]; }
//    inline iset_iter begin(int pos) { return graphs[pos].begin(); }
//    inline iset_iter end(int pos) { return graphs[pos].end(); }
//    inline iset_iter find(int pos, int pivot) { return graphs[pos].find(pivot); }
//    inline unsigned adjacent_vertices_num(int pos) { return graphs[pos].size(); }
//    inline void reset() { graphs.clear(); }

//private:
//    map<int,iset> graphs;
//};

template<class N> void skeleton_graph::print_skeleton(N &n)
{
    for(unsigned v=0; v<graphs.size(); v++)
    {
        cout<<v+n.get_v_start()<<"] ";
        for(iset_iter it=graphs[v].begin(); it!=graphs[v].end(); ++it)
            cout<<*it<<" ";
        cout<<endl;
    }
}

#endif // SKELETON_H
