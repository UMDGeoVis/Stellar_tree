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

#ifndef GENERICTOPSIMPLEX_H
#define GENERICTOPSIMPLEX_H

#include "top_cell.h"
#include "../utilities/combinatorial_utilities.h"

/**
 * @brief A class, extending Top_CP_Cell, that represents a top simplex
 *
 * The class can represent top simplices in arbitrary dimension, but the algorithms are also optimized for lower dimensional simplices
 * as triangles (in 2D) and tetrahedra (in 3D)
 *
 */
class Top_Simplex: public Top_Cell
{
public:
    /**
     * @brief A constructor method
     *
     */
    Top_Simplex() : Top_Cell() { }
    /**
     * @brief A copy-constructor method
     *
     * @param orig
     */
    Top_Simplex(const Top_Simplex& orig) : Top_Cell(orig) { }
    /**
     * @brief A constructor method
     *
     * @param vertices an integer vector containing the position indexes of vertices forming the cell
     */
    Top_Simplex(const ivect &vertices) : Top_Cell(vertices) { }
    /**
     * @brief A destructor method
     *
     */
    virtual ~Top_Simplex() {}
    /**
     * @brief A public method that returns an edge in the boundary, given the position indices of its two (sorted) extrema
     *
     * @param e, an integer vector that in the end of the procedure contains the edge
     * @param i, the position index of the an extr
     * @param j
     */
    inline void TE(ivect &e, int i, int j)
    {
        e[0] = abs(this->vertices[i]);
        e[1] = abs(this->vertices[j]);
        sort(e.begin(),e.end());
    }
    /**
     * @brief A public method that returns the edge at position pos in the boundary of the simplex
     * The procedure returns a sorted vector representing the edge
     * The procedure defines optimizations for triangles and tetrahedra
     *
     * @param e a vector of integers representing the edge
     * @param pos an integer representing the position index of the edge in the boundary of the simplex
     */
    inline void TE(ivect &e, int pos)
    {
        if(this->get_vertices_num() == 3) // TRIANGLE
            this->TF(e,pos);
        else if(this->get_vertices_num() == 4) // TETRAHEDRON
            this->TetE(e,pos);
        else /// GENERIC SIMPLEX
            this->compute_d_simplex(e,1,pos);
    }
    /**
     * @brief A public method that returns the d-1 face at position pos in the boundary of a d-simplex
     * The procedure returns a sorted vector representing the face
     *
     * @param f a vector of integers representing the d-1 face
     * @param pos an integer representing the position index of the d-1 face in the boundary of the d-simplex
     */
    inline void TF(ivect &f, int pos)
    {
        TF_unsorted(f,pos);
        sort(f.begin(),f.end());
    }
    /**
     * @brief A public method that returns the d-1 face at position pos in the boundary of a d-simplex
     *
     * The procedure returns an unsorted vector representing the face
     *
     * @param f a vector of integers representing the d-1 face
     * @param pos an integer representing the position index of the d-1 face in the boundary of the d-simplex
     */
    inline void TF_unsorted(ivect &f, int pos)
    {
        f.clear();
        f.assign(vertices.size()-1,-1);
        for(unsigned i=1; i<vertices.size(); i++)
        {
            f[i-1] = abs(this->vertices[(pos+i)%vertices.size()]);
        }
    }
    /**
     * @brief A public methods that returns the number of edges in the boundary of the simplex
     * The procedure defines optimizations for triangle and tetrahedra
     *
     * @return int
     */
    inline int get_edges_num()
    {
        if(this->get_vertices_num() == 3) // TRIANGLE
            return this->get_dfaces_num();
        else if(this->get_vertices_num() == 4) // TETRAHEDRON
            return 6;
        else /// GENERIC SIMPLEX
            return compute_sub_types_num(1);
    }
    /**
     * @brief A public methods that returns the number of d-1 faces in the boundary of the d-simplex
     *
     * @return int
     */
    inline int get_dfaces_num() { return vertices.size(); }
    /**
     * @brief A public methods that returns the number of sub-simplices types in the boundary of the simplex
     *
     * @return int
     */
    inline int get_sub_types_num() { return vertices.size() - 1; }
    /**
     * @brief A public methods that returns the number of sub-simplices of a given dimension, in the boundary of the simplex
     *
     * The procedure defines optimizations for triangles and tetrahedra
     *
     * @param d an integer representing the dimension of sub-simplices
     * @return int
     */
    int get_sub_types_num(int d);
    /**
     * @brief A public method that returns a specific sub-simplex of a dimension d at a given position pos
     *
     * The procedure defines optimizations for triangles and tetrahedra
     *
     * @param s a vector of integers encoding the final sub-simplces
     * @param d an integer representing the dimension of the sub-simplex
     * @param pos an integer representing the position index of the sub-simplex in the boundary of the top simplex
     */
    void get_d_cell(ivect &s, int d, int pos);

private:

    /**
     * @brief A protected method that computes the number of sub-cells at a given dimension d
     *
     * @param d an integer representing the dimension of sub-cells
     * @return int
     */
    int compute_sub_types_num(int d);
    /**
     * @brief A protected method that computes the vertices in the boundary of a sub cell of dimension d and position index pos in the boundary of the top cell
     *
     * The output is a sorted vector of integers.
     *
     * @param s a vector of integers encoding the vertices of the sub-cell
     * @param d an integer representing the dimension of the sub-cell
     * @param pos an integer representing the position index of the sub-cell in the boundary of the top cell
     */
    void compute_d_simplex(ivect &s, int d, int pos);
    /**
     * @brief A protected method that computes the vertices in the boundary of a sub cell of dimension d and position index pos in the boundary of the top cell
     *
     * The output is a unsorted vector of integers.
     *
     * @param s a vector of integers encoding the vertices of the sub-cell
     * @param d an integer representing the dimension of the sub-cell
     * @param pos an integer representing the position index of the sub-cell in the boundary of the top cell
     */
    void compute_d_simplex_unsorted(ivect &s, int dim, int pos);
    /**
     * @brief A private method that computes the edge at position pos.
     * The output is a sorted vector
     * WARNING: The procedure is specific for tetrahedra
     *
     * @param e a vector of int representing the edge
     * @param pos the position of the edge in the boundary of the tetrahedron
     */
    void TetE(ivect &e, int pos);
};

#endif // GENERICTOPSIMPLEX_H
