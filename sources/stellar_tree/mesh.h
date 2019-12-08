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

#ifndef MESH_H
#define MESH_H

#include <vector>
#include <cstdlib>
#include <climits>
#include <list>
#include <algorithm>
#include <numeric>

#include "vertex.h"
#include "box.h"
#include "stellar_decomposition/top_cp_cell.h"
#include "stellar_decomposition/top_simplex.h"
#include "utilities/basic_wrappers.h"

#define COORDBASETYPE double

using namespace std;

/**
 * @brief An enum type defining the top cells types that can be represented by a this->
 *
 */
enum CellType {LINE=0, TRIANGLE=1, QUAD=2, TETRA=3, HEX=4, nDIM_SIMPL=5, UNSET=-1};
static const char *cell_names[] = { "LINE", "TRIANGLE", "QUAD", "TETRA", "HEX", "SIMPL" }; /**< TODO */

/**
 * @brief A class representing a generic indexed mesh parametrized by the type of top simplexes.
 * Two templates parameters must be specified, (1) the coordinate precision of the vertices and (2) the type of top cells represented by the mesh
 * NOTA: the range of vertices and top cells indexes is defined from 1 to #num_entities (included)
 * 0 is a special value not used to index any element.
 */
template<class C, class T> class Mesh
{
public:
    /**
     * @brief A public alias for using the iterator associated to the vertices vector
     *
     */
    typedef typename vector<Vertex<C> >::iterator v_array_iter;
    /**
     * @brief A public alias for using the iterator associated to the top d-cells vector
     */
    typedef typename vector<T>::iterator top_d_array_iter;

    /**
     * @brief A constructor method
     *
     */
    Mesh()
    {
        domain = Box<C>();
        vertices = vector<Vertex<C> >();
    }
    /**
     * @brief A constructor method
     *
     * @param orig
     */
    Mesh(const Mesh& orig)
    {
        this->vertices = orig.vertices;
        this->domain = orig.domain;
        this->tops = orig.tops;
    }
    /**
     * @brief A destructor method
     *
     */
    virtual ~Mesh()
    {
        vertices.clear();
        tops.clear();
    }

    // ADD PROCEDURES //
    /**
     * @brief A public method that adds a vertex to the vertices list
     *
     * @param v a Vertex<C>& argument, representing the vertex to add
     * @return an integer representing the position index of the vertex (and the current size of the vertices array)
     *
     */
    inline int add_vertex(Vertex<C> &v)
    {
        this->vertices.push_back(v);
        return this->vertices.size();
    }
    /**
     * @brief Public method that adds a top cell to the specific array
     *
     * @param pos, the position index of the array
     * @param t, the top cell T& to add
     * @return an integer representing the position index of the top cell (and the current size of the top d-cells array)
     *
     */
    inline int add_top_cell(int pos, T& t)
    {
        this->tops[pos].push_back(t);
        return tops[pos].size();
    }
    /**
     * @brief Public method that adds a top cell to the specific array
     *
     * @param pos, the position index of the array
     * @param v_ids, an integer vector containing the position indexes of the vertices forming the top cell to add
     * @return an integer representing the position index of the top cell (and the current size of the top d-cells array)
     *
     */
    inline int add_top_cell(int pos, const ivect &v_ids)
    {
        T top = T(v_ids);
        tops[pos].push_back(top);
        return tops[pos].size();
    }
    /**
     * @brief A public method that adds an array of top cells in the respective array
     *
     * @param pos, the position index of the array encoding the top cells type
     * @param l, the array containing the top cells to add
     */
    template<class L> inline void add_top_cells(int pos, L &l)
    {
        for(typename L::iterator it=l.begin(); it!=l.end(); ++it)
        {
            this->add_top_cell(pos,*it);
        }
    }
    /**
     * @brief A public method that adds a the top cells array
     *
     */
    inline void add_top_cells_array() { tops.push_back(vector<T>()); }
    /**
     * @brief A public method that adds 'num' top cells arrays
     * @param num, an integer argument representing the number of arrays to add
     */
    inline void add_top_cells_arrays(int num)
    {
        for(int i=0; i<num; ++i)
            tops.push_back(vector<T>());
    }
    /**
     * @brief A public method that adds a top cells type
     *
     * @param type, the top cells type to add
     */
    inline void add_type(CellType type) { this->top_type.push_back(type); }
    /**
     * @brief A public method that adds 'num' top cells types
     * @param num, an integer argument representing the number of types to add
     * @param type, the top cells type to add
     */
    inline void add_types(int num, CellType type)
    {
        for(int i=0; i<num; ++i)
            top_type.push_back(type);
    }

    // SET PROCEDURES //
    /**
     * @brief A public method that init the vertex position pos of the vertices array with vertex v
     * NOTA: the index of the vertex is +1 respect to the real vertex position in the array
     * @param pos, the position index in the vertices array
     * @param v, the Vertex<C>& to set at position pos
     */
    inline void set_vertex(int pos, Vertex<C>& v) { this->vertices[pos-1] = v; }
    /**
     * @brief A public method that init the top cell position pos of the top cells array with top cell t
     * NOTA: the index of the top cell is +1 respect to the real top cell position in the array
     * @param dim, point to the corresponding top cell array
     * @param pos, the position index in the top cell
     * @param t, the top cell T& to set at position pos
     */
    inline void set_top_cell(int dim, int pos, T& t) { this->tops[dim][pos-1] = t; }
    /**
     * @brief A public method that init the top cell position pos of the top cells array with an array of indices v_ids, representing the boundary of a top cell
     * NOTA: the index of the top cell is +1 respect to the real top cell position in the array
     * @param dim, point to the corresponding top cell array
     * @param pos, the position index in the top cell
     * @param v_ids, the array of vertices indices forming the boundary of the top cell
     */
    inline void set_top_cell(int dim, int pos, ivect& v_ids) { this->tops[dim][pos-1] = T(v_ids); }
    /**
     * @brief A public method that sets the mesh domain
     *
     * @param d a Box<C>& argument, representing the mesh domain
     */
    inline void set_domain(Box<C>& d) { this->domain = d; }

    // GET PROCEDURES //
    /**
     * @brief A public method that returns the vertex at the i-th position in the mesh list
     * NOTA: the index of the vertex is +1 respect to the real vertex position in the array
     * @param id an integer argument, representing the position in the list
     * @return Vertex<C>&, the vertex at the id-th position in the list
     */
    inline Vertex<C>& get_vertex(int id) { return this->vertices[id-1]; }
    /**
     * @brief A public method that returns the number of mesh vertices
     * @return an integer, representing the number of vertices
     */
    inline int get_vertices_num() { return this->vertices.size(); }
    /**
     * @brief A public method that returns a reference to the vertices array
     * @return vector<Vertex<C> >&, the vertices array
     */
    inline vector<Vertex<C> >& get_vertices_array() { return this->vertices; }
    /**
     * @brief A public method that returns a reference to the top cells array
     * @return vector<vector<T> >&, the top cells array
     */
    inline vector<vector<T> >& get_top_cells_array() { return this->tops; }
    /**
     * @brief A public method that returns the mesh domain
     * @return Box<C>&, the mesh domain
     */
    inline Box<C>& get_domain() { return this->domain; }
    /**
     * @brief A public method that returns the top cell at the i-th position in a given top cells array
     * @param pos, the position index of the array
     * @param id, the top cell position index
     * @return T&, the reference of the top cell
     */
    inline T& get_top_cell(int pos, int id) { return this->tops[pos][id-1]; }
    /**
     * @brief A public method that returns the number of top cells types encoded by the mesh
     * @return int, the number of top cells types
     */
    inline int get_top_cells_types() { return this->tops.size(); }
    /**
     * @brief A public method that return the number of the top cells of a specific type
     *
     * @param pos, the position index of the array containing the specific top cells
     * @return int, the number of top cells
     */
    inline int get_top_cells_num(int pos) { return this->tops[pos].size(); }
    /**
     * @brief A public method counting the top cells encoded by the mesh
     *
     * @return int, the number of all top cells encoded by the mesh
     */
    inline int count_top_cells_num()
    {
        int counter = 0;
        for(unsigned i=0; i<tops.size(); i++)
            counter += this->tops[i].size();
        return counter;
    }
    /**
     * @brief A public method that returns the number of sub-cells implicitly encoded by the mesh
     *
     * @return int, the number of these sub-cells
     */
    inline int get_implicitly_encoded_cells_num()
    {
        return tops[tops.size()-1][0].get_sub_types_num();
    }
    /**
     * @brief A public method that returns the number of cells implicitly and explicitly encoded by the mesh
     *
     * @return int, the number of the cells implicitly and explicitly encoded
     */
    inline int get_implicitly_and_explicitly_encoded_cells_num()
    {
        return tops[tops.size()-1][0].get_sub_types_num() + 1; // +1 because we have to count also the top cells in the last level!
    }

    /**
     * @brief A public method that returns the begin iterator of the vertices array
     * @return v_array_iter, the begin iterator of the array
     */
    inline v_array_iter v_array_begin() { return this->vertices.begin(); }
    /**
     * @brief A public method that returns the end iterator of the vertices array
     * @return v_array_iter, the end iterator of the array
     */
    inline v_array_iter v_array_end() { return this->vertices.end(); }
    /**
     * @brief A public method that returns the begin iterator of the top d-cells array
     * @param d, the position index of the array
     * @return top_d_array_iter, the begin iterator of the array
     */
    inline top_d_array_iter top_d_array_begin(int d) { return this->tops[d].begin(); }
    /**
     * @brief A public method that returns the end iterator of the top d-cells array
     * @param d, the position index of the array
     * @return top_d_array_iter, the end iterator of the array
     */
    inline top_d_array_iter top_d_array_end(int d) { return this->tops[d].end(); }
    /**
     * @brief A public method that returns the type of top cells encoded at a given array position
     *
     * @param d, the position index of the array
     * @return CellType, the type of the top cells
     */
    inline CellType get_type(int d) const { return this->top_type[d]; }

    /**
     * @brief A public method that initializes the top cells array, adding 'dim' arrays.
     * The number of arrays added is coherent with the number of top cells types encoded by the mesh
     *
     * @param dim, the number of top cells types
     */
    inline void init_top_cells_array(int dim) { tops.assign(dim,vector<T>()); }
    /**
     * @brief A public method that reserve the space needed by the top d-cells array
     *
     * @param d an integer representing the position of the top d-cells in the top cells arrays
     * @param numT an integer representing the space required by the array of top d-cells
     */
    inline void reserve_top_cells_space(int d, int numT) { this->tops[d].reserve(numT); }
    /**
     * @brief A public method that reserve the space needed by the vertices array
     *
     * @param numV an integer representing the space required by the array of vertices
     */
    inline void reserve_vertices_space(int numV) { this->vertices.reserve(numV); }
    /**
     * @brief A public method that resets the vertices array
     *
     */
    inline void reset_vertices_vector() { this->vertices.clear(); }
    /**
     * @brief A public method that resets the top cells arrays
     *
     */
    inline void reset_top_cells_vector() { this->tops.clear(); }
    /**
     * @brief A public method that resizes the top d-cells array.
     * If new_size is smaller than the current size of the array, the position encoded after size are effectively deleted from the array.
     * @param d an integer representing the position of the top d-cells in the top cells arrays
     * @param new_size an integer representing the new size of the array
     */
    inline void resize_top_d_cell_array(int d, int new_size) { this->tops[d].resize(new_size); }

    /**
     * @brief A public method that checks if the vertex at position v is flagged as deleted
     *
     * @param v the position of the vertex in the vertices array
     * @return bool, true if the vertex is deleted, false otherwise
     */
    inline bool is_vertex_removed(int v) { return (get_vertex(v).getC(0)==INFINITY); }
    /**
     * @brief A public method that checks if the vertex v is flagged as deleted
     *
     * @param v the vertex to check
     * @return bool, true if the vertex is deleted, false otherwise
     */
    inline bool is_vertex_removed(Vertex<C> &v) { return (v.getC(0)==INFINITY); }
    /**
     * @brief A public method that checks if the top cell of dimension dim at position t is flagged as deleted
     *
     * @param dim position of the top cell in the array of arrays
     * @param t the top cell position index in that array
     * @return bool, true if the top cell is deleted, false otherwise
     */
    inline bool is_top_cell_removed(int dim, int t) { return (get_top_cell(dim,t).TV(0) == INT_MAX); }
    /**
     * @brief A public method that checks if the top cell t is flagged as deleted
     *
     * @param t the top cell to check
     * @return bool, true if the top cell is deleted, false otherwise
     */
    inline bool is_top_cell_removed(T &t) { return (t.TV(0) == INT_MAX); }
    /**
     * @brief A public method that flags as removed the vertex with position index v_id
     *
     * @param v_id the position of the vertex to delete
     */
    inline void remove_vertex(int v_id)
    {
        Vertex<C>& v = get_vertex(v_id);
        v.setC(0,INFINITY);
    }
    /**
     * @brief A public method that flags as removed the top cell in the dim-th array and with position index t
     *
     * @param dim, the array position
     * @param t the position index of the top cell to delete
     */
    inline void remove_top_cell(int dim, int t)
    {
        T& top = get_top_cell(dim,t);
        top.setTV(0,INT_MAX);
    }
    /**
     * @brief A public method that erase a vertex from the vertices array
     * NOTA: after the erase the position indices of the vertices following the iterator become invalid.
     * @param it, an iterator pointing to the position to erase
     * @return the iterator to the next valid entry of the array
     */
    inline v_array_iter erase_vertex(v_array_iter &it) { return this->vertices.erase(it); }
    /**
     * @brief A public method that erase a vertex from the top d-cells array
     * NOTA: after the erase the position indices of the top d-cells following the iterator become invalid.
     * @param d, an integer representing the array position
     * @param it, an iterator pointing to the position to erase
     * @return the iterator to the next valid entry of the array
     */
    inline top_d_array_iter erase_top_cell(int d, top_d_array_iter &it) { return this->tops[d].erase(it); }

    /**
     * @brief A public method that deletes the top d-cells array of the mesh
     * @param d, an integer representing the dimension of the top d-cells in the mesh arrays
     */
    inline void clear_top_d_cell_list(int d) { this->tops[d].clear(); }

    /**
     * @brief A public method that swaps the position of two vertices
     * @param v_id1, the position index of the first vertex
     * @param v_id2, the position index of the secton vertex
     */
    inline void vertices_swap(int v_id1, int v_id2) { swap(vertices[v_id1-1],vertices[v_id2-1]); }
    /**
     * @brief A public method that swaps the position of two top cells
     * @param pos, the position index of the array containing the two top cells
     * @param t_id1, the position index of the first top cell
     * @param t_id2, the position index of the second top cell
     */
    inline void tops_swap(int pos, int t_id1, int t_id2) { swap(this->tops[pos][t_id1-1],this->tops[pos][t_id2-1]); }

    /**
     * @brief A public method that returns a string containing the type of top cells encoded at a given array position
     *
     * @param d, the position index of the array
     * @return string containing the top cells type
     */
    inline string print_type(int d) { return (this->top_type[d] != UNSET? cell_names[this->top_type[d]] : "UNSET"); }
    /**
     * @brief A public method that checks if the the top cells in a given array are simplices
     *
     * @param d, the position index of the array
     * @return bool, true if the top cells in the d-array are simplices, false otherwise
     */
    inline bool is_simplex(int d) { return (this->top_type[d] == LINE || this->top_type[d] == TRIANGLE || this->top_type[d] == TETRA || this->top_type[d] == nDIM_SIMPL); }

    /**
     * @brief A public method that print to a specified output some statistics concering the indexed mesh
     * These statistics refer to the number of encoded vertices/top cells and the number of non-deleted entries.
     */
    template<class O> void print_mesh_stats(O &output);

    template<class O> void print_mesh(O &output);
//    //
//    void print_mesh_adj_graph(string mesh_name);
    /**
     * @brief A public method that computes some statistics concerning the vertices distance of the indexed mesh
     * These statistics refer to the minimum, average and maximum distance between two vertices
     * and the distribution of the vertices distances.
     */
    void get_vertices_distance_stats();

protected:
    ///A private varible representing the mesh domain
    Box<C> domain;
    ///A private varible representing the vertices list of the mesh
    vector<Vertex<C> > vertices;
    ///A private varible representing the top cells list of the mesh
    vector<vector<T> > tops;
    ///A private variable representing the top cells types
    vector<CellType> top_type;
};

template<class C, class T> template<class O> void Mesh<C,T>::print_mesh_stats(O &output)
{
    int real_vertices_counter = 0;
    for(unsigned i=1; i<=this->vertices.size(); i++)
    {
        if(!this->is_vertex_removed(i))
            real_vertices_counter++;
    }
    output<<"  Vertices: "<<real_vertices_counter<<" (with "<<this->vertices.size()<<" encoded entries) "<<endl;

    int tot_simpl_num = 0;
    for(unsigned i=0; i<tops.size(); i++)
        tot_simpl_num += tops[i].size();
    output<<"  Top-cells: "<<tot_simpl_num<<endl;
    if(tot_simpl_num > 0)
    {
        for(unsigned i=0; i<tops.size(); i++)
        {
            if(tops[i].size() > 0)
            {
                int counter = 0;
                for(unsigned j=1; j<=tops[i].size(); j++)
                    if(!this->is_top_cell_removed(i,j))
                        counter++;
                output<<"    "<<i+1<<"-simplices: "<<counter<<" (with "<<tops[i].size()<<" encoded entries) "<<endl;
            }
        }
    }
}

template<class C, class T> template<class O> void Mesh<C,T>::print_mesh(O &output)
{
    for(int d=0; d<get_top_cells_types(); d++)
        for(int i=1; i<=get_top_cells_num(d); i++)
            if(!is_top_cell_removed(d,i))
                output<<d<<"-"<<i<<"] "<<get_top_cell(d,i)<<endl;
}

template<class C, class T> void Mesh<C,T>::get_vertices_distance_stats()
{
    cerr<<"[DEBUG] Getting vertices distance stats"<<endl;

    long long num_test = 0;
    C min=INT_MAX, avg=0, max=INT_MIN;

#pragma omp parallel for
    for(unsigned i=0; i<vertices.size(); i++)
    {
        for(unsigned j=i+1; j<vertices.size(); j++)
        {
            C distance = vertices[i].distance(vertices[j]);

#pragma omp critical
            {
                if(min > distance)
                    min = distance;
                if(max < distance)
                    max = distance;
                avg += distance;
                num_test++;
            }
        }
    }

    avg = avg / num_test;
    cerr << "[Edges number] " << num_test << endl;
    cerr<<"[Vertices distance] min: "<<min<<" avg: "<<avg<<" max: "<<max<<endl;
    ivect edges_distribution;
    edges_distribution.assign(10,0);
    for(unsigned i=0; i<vertices.size(); i++)
    {
        for(unsigned j=i+1; j<vertices.size(); j++)
        {
            C distance = vertices[i].distance(vertices[j]);
            C ratio = distance * 100 / max;
            if(ratio > 0 && ratio <= 10)
                edges_distribution[0]++;
            else if(ratio > 10 && ratio <= 20)
                edges_distribution[1]++;
            else if(ratio > 20 && ratio <= 30)
                edges_distribution[2]++;
            else if(ratio > 30 && ratio <= 40)
                edges_distribution[3]++;
            else if(ratio > 40 && ratio <= 50)
                edges_distribution[4]++;
            else if(ratio > 50 && ratio <= 60)
                edges_distribution[5]++;
            else if(ratio > 60 && ratio <= 70)
                edges_distribution[6]++;
            else if(ratio > 70 && ratio <= 80)
                edges_distribution[7]++;
            else if(ratio > 80 && ratio <= 90)
                edges_distribution[8]++;
            else if(ratio > 90 && ratio <= 100)
                edges_distribution[9]++;
        }
    }
    for(unsigned i=0; i<edges_distribution.size(); i++)
    {
        C dmin = max / 100.0 * (i * 10.0), dmax = max / 100.0 * ((i+1) * 10.0);
        cerr << "num_edge between  " << dmin << " and " << dmax << " -> " << edges_distribution[i] << endl;
    }
}

/** @typedef CP_Mesh
 * @brief A Mesh<C,T> specialization, specific for Constant-Polytope (CP) complexes
 *
 */
typedef Mesh<COORDBASETYPE,Top_CP_Cell> CP_Mesh;
/** @typedef Simplicial_Mesh
 * @brief A Mesh<C,T> specialization, specific for Simplicial complexes
 *
 */
typedef Mesh<COORDBASETYPE,Top_Simplex> Simplicial_Mesh;

#endif // MESH_H
