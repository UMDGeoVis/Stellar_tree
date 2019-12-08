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

#ifndef TOP_CELL_H
#define TOP_CELL_H

#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>

#include <utilities/basic_wrappers.h>

using namespace std;

/**
 * @brief An abstract-class representing a generic top cell.
 *
 * The class can represent top cells in arbitrary dimension and encodes just the boundary relation of a top cell.
 *
 */
class Top_Cell
{
public:
    /**
     * @brief A destructor method
     *
     */
    virtual ~Top_Cell() {}
    /**
     * @brief A public method that returns a vertex in the boundary of the cell
     * @param pos an integer representing the position index of the vertex in the boundary array
     * @return int
     */
    inline int TV(int pos) const { return this->vertices[pos]; }
    /**
     * @brief A public method that returns a vertex in the boundary of the cell saved into a vector
     * @param v an integer vector representing the vertex
     * @param pos an integer representing the position index of the vertex in the boundary array
     */
    inline void TV_vector(ivect &v, int pos) { v.clear(); v.push_back(this->vertices[pos]); }
    /**
     * @brief A public method that sets position pos in the boundary of the vertex
     * @param pos an integer representing the position index in the boundary array to set
     * @param newId an integer representing the position index of the vertex
     */
    inline void setTV(int pos, int newId) { this->vertices[pos] = newId; }
    /**
     * @brief A public methods that returns the number of vertices indices in the boundary of the cell
     * @return int
     */
    inline int get_vertices_num() const { return vertices.size(); }
    /**
     * @brief A public method that returns the vertices in the boundary of the top cell
     * @return a reference to the integer array
     */
    inline ivect& get_vertices_array() { return this->vertices; }
    /**
     * @brief A public method that returns a vector of integers containing the sorted (positive) vertices indices of the top cell
     *
     * @param v
     */
    inline void get_positive_sorted_vertices(ivect &v) const
    {
        v.assign(vertices.size(),-1);
        for(unsigned i=0; i<vertices.size(); i++)
            v[i] = abs(this->vertices[i]);
        sort(v.begin(),v.end());
    }
    /**
     * @brief A public method that checks if a top d-cell degenerates to a d-1 cell by having at least two identical vertices
     *
     * @return bool is true if the cell is degenerate, false otherwise
     */
    inline bool is_degenerate()
    {
        for(unsigned i=0; i<this->vertices.size(); i++)
        {
            for(unsigned j=i+1; j<this->vertices.size(); j++)
            {
                if(abs(this->vertices[i]) == abs(this->vertices[j]))
                    return true;
            }
        }
        return false;
    }    
    /**
     * @brief A public method that returns the position index of a vertex (if exists) in the local boundary array of the cell
     *
     * @param v_id an integer representing the position index of the vertex in the mesh
     * @return int an integer containing the position index of the vertex in the boundary array or -1 if the vertex is not in the boundary of the cell
     */
    inline int vertex_index(int v_id)
    {
        for(int j=0; j<get_vertices_num(); j++)
        {
            if(abs(this->TV(j)) == v_id)
                return j;
        }
        return -1;
    }
    /**
     * @brief A public method that checks if a vertex is on the boundary of the cell
     *
     * @param v_id an integer representing the position index of the vertex in the mesh
     * @return bool is true if the vertex is on the boundary, false otherwise
     */
    inline bool has_vertex(int v_id) { return (this->vertex_index(v_id)!=-1); }
    /**
     * @brief A public method that checks if cell c is a sub-cell of the top cell
     * @param c, the cell to test
     * @return true if c is a sub-cell of the top cell, false otherwise
     */
    bool has_cell(const ivect &c);
    /**
     * @brief A public method that returns the (possibly empty) shared sub-cell between two cells
     * @param cell, a generic integer container, that represents the other cell to intersect with the top cell
     * @param shared an integer vector, that contains the intersection between the two cells
     */
    template<class IntContainer> void get_shared_cell(IntContainer &cell, ivect &shared);
    /**
     * @brief A public method that returns the vertex with the lowest index
     *
     * @return int
     */
    int min_v_index();
    /**
     * @brief A public method that returns the vertex with the highest index
     *
     * @return int
     */
    int max_v_index();
    /**
     * @brief A public method that checks if two cells are equal
     * @param c, an integer vector representing the other cell to test
     * @return true, if the two cells are equal, false otherwise
     */
    bool are_equal(const ivect &c);

    /**
     * @brief
     *
     * @param p
     * @param q
     * @return bool operator
     */
    inline friend bool operator< (const Top_Cell &p, const Top_Cell &q)
    {
        if(p.get_vertices_num() < q.get_vertices_num())
            return true;
        else if(p.get_vertices_num() > q.get_vertices_num())
            return false;

        // here p and q have the same number of vertices
        ivect pv;
        p.get_positive_sorted_vertices(pv);
        ivect qv;
        q.get_positive_sorted_vertices(qv);

        return pv < qv;
    }
    /**
     * @brief
     *
     * @param p
     * @param q
     * @return bool operator
     */
    inline friend bool operator== (const Top_Cell &p, const Top_Cell &q)
    {
        if(p.get_vertices_num() != q.get_vertices_num())
            return false;

        boost::dynamic_bitset<> b(p.get_vertices_num());

        for(int i=0;i<p.get_vertices_num();i++)
        {
            for(int j=0;j<q.get_vertices_num();j++)
            {
                if(!b[j] && abs(p.vertices[i])==abs(q.vertices[j]))
                {
                    b[j] = 1;
                    break;
                }
            }
        }

        return b.count() == b.size();
    }
    /**
     * @brief
     *
     * @param p
     * @param q
     * @return bool operator
     */
    inline friend bool operator!= (const Top_Cell &p, const Top_Cell &q) { return !(p==q); }
    /**
     * @brief
     *
     * @param out
     * @param p
     * @return std::ostream &operator
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Top_Cell& p)
    {
        out << "[";
        for(int i=0; i<p.get_vertices_num(); i++)
        {
            out << p.TV(i) << " ";
        }
        out << "] ";
        return out;
    }
    /**
     * @brief A public method that saves the list of vertices in the boundary in a FILE* stream
     *
     * @param output
     */
    inline void save_to_file(FILE* output)
    {
        for(int i=0; i< get_vertices_num(); i++)
            fprintf(output,"%d ",abs(TV(i))-1);
    }

protected:
    /// A protected variable representing the vertices array in the boundary of the cell
    ivect vertices;

    /**
     * @brief A constructor method
     *
     */
    Top_Cell() { }
    /**
     * @brief A copy-constructor method
     *
     * @param orig
     */
    Top_Cell(const Top_Cell& orig) { this->vertices = orig.vertices; }
    /**
     * @brief A constructor method
     *
     * @param vertices, an integer vector containing the position indexes of vertices forming the cell
     */
    Top_Cell(const ivect &vertices) { this->vertices = vertices; }
};

template<class IntContainer> void Top_Cell::get_shared_cell(IntContainer &cell, ivect &shared)
{
    for(auto itv=cell.begin(); itv!=cell.end(); ++itv)
    {
        for(int v=0; v<get_vertices_num(); v++)
        {
            if(*itv == abs(TV(v)))
            {
                shared.push_back(*itv);
                break;
            }
        }
    }
}

#endif // TOP_CELL_H
