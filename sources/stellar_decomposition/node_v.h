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

#ifndef NODE_V_H
#define NODE_V_H

#include <utilities/basic_wrappers.h>

#include "node.h"
#include "run_iterator.h"

/**
 * @brief A non-instantiable, extending Node, that encodes a list vertices contained by the node
 *
 */
template<class N> class Node_V : public Node<N>
{
public:
    /**
     * @brief A destructor method
     *
     */
    virtual ~Node_V() {}    
    /**
     * @brief A public method that adds a vertex position index to the corresponding list
     *
     * @param ind an integer argument, representing the vertex index to add
     */
    inline void add_vertex(int ind) { this->vertices.push_back(ind); }
    /**
     * @brief A public method that checks if a node contains the maximum number of vertices admitted
     *
     * @param maxVertex an integer argument, representing the maximum number of vertices admitted into the node
     * @return a boolean value, true if the limit is exceeded, false otherwise
     */
    inline bool is_full(int maxVertex) { return (this->vertices.size() > (unsigned)maxVertex); }
    /**
     * @brief A public method that clears space occupied by the two lists
     *
     */
    inline void clear_v_list() { this->vertices.clear(); }
    /**
     * @brief A public method that returns the starting vertex position index of the node
     *
     * @return an integer
     *
     * NOTA: this method can be used only after exploiting the vertices spatial coherence
     */
    inline int get_v_start() const { return abs(vertices[0]); }
    /**
     * @brief A public method that returns the first vertex position index not contained by the node
     *
     * @return an integer
     *
     * NOTA: this method can be used only after exploiting the vertices spatial coherence
     */
    inline int get_v_end() const { return (abs(vertices[0])+vertices[1])+1; }
    /**
     * @brief A public method that set the vertices range after exploiting the vertices spatial coherence
     *
     * @param start an integer with the first vertex position index of the node
     * @param end an integer with the first vertex position index not contained by the node
     */
    inline void set_v_range(int start, int end) { vertices.push_back(-start); vertices.push_back(end-start-1); }
    /**
     * @brief A public method that checks if a vertex is indexed by the node
     *
     * @param v_id is position index of the vertex to check
     * @return a boolean if the vertex is indexed by the node, false otherwise
     *
     * NOTA: this method can be used only after exploiting the vertices spatial coherence
     */
    inline bool indexes_vertex(int v_id) { return v_id >= get_v_start() && v_id < get_v_end(); }
    /**
     * @brief A public method that checks if the node is before the vertex with position index v_id
     *
     * @param v_id an integer representing the position index in the mesh of the target vertex
     * @return a boolean if the node is before the vertex, false otherwise
     *
     * NOTA: this method can be used only after exploiting the vertices spatial coherence
     */
    inline bool is_before_vertex(int v_id) { return v_id >= get_v_end(); } // the node is before the vertex v_id
    /**
     * @brief A public method that checks if the node is after the vertex with position index v_id
     *
     * @param v_id is position index of the vertex to check
     * @return a boolean if the node is after the vertex, false otherwise
     *
     * NOTA: this method can be used only after exploiting the vertices spatial coherence
     */
    inline bool is_after_vertex(int v_id) { return v_id < get_v_start(); }
    /**
     * @brief A public method that checks if a cell is indexed by the node
     *
     * @param s is an integer vector containing the cell to check
     * @return a boolean if the cell is indexed by the node, false otherwise
     *
     * NOTA: this method can be used only after exploiting the vertices spatial coherence
     */
    inline bool indexes_cell(ivect &s)
    {
        for(unsigned i=0; i<s.size(); i++)
        {
            if(this->indexes_vertex(s[i]))
            {
                return true;
            }
        }
        return false;
    }
    /**
     * @brief A public method that checks if a cell is indexed by the node and saves the first vertex indexed by the node
     *
     * @param s is an integer vector containing the cell to check
     * @param v_id an integer that saves (if any) the position index of a vertex of the cell indexed by the node
     * @return a boolean if the cell is indexed by the node, false otherwise
     *
     * NOTA: this method can be used only after exploiting the vertices spatial coherence
     */
    inline bool indexes_cell(ivect &s, int &v_id)
    {
        for(unsigned i=0; i<s.size(); i++)
        {
            if(this->indexes_vertex(s[i]))
            {
                v_id = s[i];
                return true;
            }
        }
        return false;
    }

    /**
     * @brief A public method that checks if a top-cell is indexed by the node
     *
     * @param t is the top-cell to check
     * @return a boolean if the cell is indexed by the node, false otherwise
     *
     * NOTA: this method can be used only after exploiting the vertices spatial coherence
     */
    template<class T> inline bool indexes_top(T &t)
    {
        for(int i=0; i<t.get_vertices_num(); i++)
            if(this->indexes_vertex(abs(t.TV(i))))
                return true;
        return false;
    }
    /**
     * @brief A public method that checks if a top-cell is completely indexed by the node
     *
     * @param t is the top-cell to check
     * @return a boolean if the cell is completely indexed by the node, false otherwise
     */
    template<class T> inline bool completely_indexes_top(T &t)
    {
        for(int i=0; i<t.get_vertices_num(); i++)
            if(!this->indexes_vertex(abs(t.TV(i))))
                return false;
        return true;
    }
    /**
     * @brief A public method that checks if a top-cell has the vertices indexed at most by the current node
     *
     * @param t is the top-cell to check
     * @return a boolean if the condition is verified, false otherwise
     */
    template<class T> inline bool has_all_vertices_visited(T &t)
    {
        for(int i=0; i<t.get_vertices_num(); i++)
            if(abs(t.TV(i)) >= get_v_end())
                return false;
        return true;
    }
    /**
     * @brief A public method that checks if a top-cell has the vertices indexed at most by the current node
     *
     * @param t is the top-cell to check
     * @return a boolean if the condition is verified, false otherwise
     */
    template<class T> inline bool has_lower_index_inside(T &t) { return indexes_vertex(t.min_v_index()); }

    /**
     * @brief A public method that save in a integer vector the vertices that indexed outside the current node
     *
     * @param t is the top-cell to check
     * @param outer_v_ids a integer vector that saves (if any) the vertices indexed outside the current node
     */
    template<class T> inline void get_outer_v_ids(T &t, ivect &outer_v_ids)
    {
        for(int i=0; i<t.get_vertices_num(); i++)
        {
            if(!this->indexes_vertex(abs(t.TV(i))))
                outer_v_ids.push_back(abs(t.TV(i)));
        }
    }
    ///
    /**
     * @brief
     *
     * @param out
     * @param p
     * @return std::ostream &operator
     */
    friend std::ostream& operator<<(std::ostream& out, const Node_V<N>& p)
    {
        if(p.is_leaf())
        {
            if(p.is_v_array_run_compressed())
                out <<"Leaf["<< p.get_v_start() << " " << p.get_v_end()<<"]";
            else
                out <<"Leaf[ v_size "<< p.get_real_v_array_size()<<"]";
        }
        else
            out <<"Node["<< p.get_v_start() << " " << p.get_v_end()<<"]";
        return out;
    }
    /**
     * @brief A public method that returns a pair of iterator to the vertices array
     *
     * @return RunIteratorPair
     */
    inline RunIteratorPair make_v_array_iterator_pair() { return run_iterator<int>::make_run_iterator_pair(vertices); }
    /**
     * @brief A public method that returns the iterator to the begin of the vertices array
     *
     * @return RunIterator
     */
    inline RunIterator v_array_begin_iterator() { return run_iterator<int>(vertices.begin(),vertices.end()); }
    /**
     * @brief A public method that returns the iterator to the end of the vertices array
     *
     * @return RunIterator
     */
    inline RunIterator v_array_end_iterator() { return run_iterator<int>(vertices.end()); }
    /**
     * @brief A public method that returns the number of indexed vertices
     *
     * @return int
     */
    inline int get_real_v_array_size() const { return run_iterator<int>(vertices.begin(),vertices.end()).elementCountFast(vertices); }
    /**
     * @brief A public method that returns the size of the vertices array
     *
     * @return int
     *
     * NOTA: after the exploit of the spatial coherence this procedure returns always 2, as each node encode a run for the vertices
     */
    inline int get_v_array_size() const { return this->vertices.size(); }

    /**
     * @brief A public method that update the list of vertices encoded by the leaf block
     *
     * @return vl contains the updated vertices list
     *
     * NOTA: this update can invalidate the run compression for the vertices
     */
    inline void set_v_array(ivect &vl) { this->vertices = vl; }

    /**
     * @brief A public method that check if the vertices list is encoded using the run compression encoding
     *
     * @return true if the run compression is used, false otherwise
     */
    inline bool is_v_array_run_compressed() const { return (this->vertices.size()==2 && this->vertices[0]<0); }

protected:
    /**
     * @brief A constructor method
     *
     */
    Node_V() : Node<N>() { }
    /**
     * @brief A constructor method
     *
     * @param orig
     */
    Node_V(const Node_V& orig) : Node<N>(orig) { this->vertices = orig.vertices; }

    /// A protected variable represeting the list of vertices indexed by the node
    ivect vertices;
};

#endif // NODE_V_H
