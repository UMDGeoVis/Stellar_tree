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

#ifndef TOP_CP_CELL_H
#define TOP_CP_CELL_H

#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>

#include "top_cell.h"

using namespace std;

/**
 * @brief A class representing a generic top Constant-Polytope cell.
 *
 * The class can represent top cells in arbitrary dimension, but the algorithms are defined for
 * lower dimensional CP cells (quads in 2D and hexahedra in 3D)
 *
 */
class Top_CP_Cell : public Top_Cell
{
public:
    /**
     * @brief A constructor method
     *
     */
    Top_CP_Cell() : Top_Cell() { }
    /**
     * @brief A copy-constructor method
     *
     * @param orig
     */
    Top_CP_Cell(const Top_CP_Cell& orig) : Top_Cell(orig) { }
    /**
     * @brief A constructor method
     *
     * @param vertices, an integer vector containing the position indexes of vertices forming the CP-cell
     */
    Top_CP_Cell(const ivect &vertices) : Top_Cell(vertices) { }
    /**
     * @brief A destructor method
     *
     */
    virtual ~Top_CP_Cell() {}
    /**
     * @brief A public method that returns the edge at position pos in the boundary of the CP-cell
     *
     * The procedure returns a sorted vector representing the edge
     *
     * @param e a vector of integers representing the edge
     * @param pos an integer representing the position index of the edge in the boundary of the cell
     */
    inline void TE(ivect &e, int pos)
    {
        if(this->get_vertices_num() == 4) // QUAD
            this->QE(e,pos);
        else if(this->get_vertices_num() == 8) // HEXAHEDRON
            this->HE(e,pos);
    }
    /**
     * @brief A public method that returns the d-1 face at position pos in the boundary of a CP-d-cell
     * The procedure returns a sorted vector representing the face
     *
     * @param f a vector of integers representing the d-1 face
     * @param pos an integer representing the position index of the d-1 face in the boundary of the CP-d-cell
     */
    inline void TF(ivect &f, int pos)
    {
        if(this->get_vertices_num() == 4) // QUAD
            this->QE(f,pos);
        else if(this->get_vertices_num() == 8) // HEXAHEDRON
            this->HF(f,pos);
    }
    /**
     * @brief A public method that returns the d-1 face at position pos in the boundary of a CP-d-cell
     * The procedure returns a unsorted vector representing the face
     *
     * @param f a vector of integers representing the d-1 face
     * @param pos an integer representing the position index of the d-1 face in the boundary of the CP-d-cell
     */
    inline void TF_unsorted(ivect &f, int pos)
    {
        if(this->get_vertices_num() == 4) // QUAD
            this->QE_unsorted(f,pos);
        else if(this->get_vertices_num() == 8) // HEXAHEDRON
            this->HF_unsorted(f,pos);
    }
    /**
     * @brief A public methods that returns the number of edges in the boundary of the CP-cell
     * @return int
     */
    int get_edges_num();
    /**
     * @brief A public methods that returns the number of d-1 faces in the boundary of the d-cell
     * @return int
     */
    int get_dfaces_num();
    /**
     * @brief A public methods that returns the number of sub-cells types in the boundary of the cell
     * @return int
     */
    int get_sub_types_num();
    /**
     * @brief A public methods that returns the number of sub-cells of a given dimension, in the boundary of the cell
     * @param d an integer representing the dimension of sub-cells
     * @return int
     */
    int get_sub_types_num(int d);
    /**
     * @brief A public method that returns a specific sub-cell of a dimension d at a given position pos
     *
     * @param s a vector of integers encoding the final sub-cell
     * @param d an integer representing the dimension of the sub-cell
     * @param pos an integer representing the position index of the sub-cell in the boundary of the top cell
     */
    void get_d_cell(ivect &s, int d, int pos);
    /**
     * @brief A public method that returns the position indexes of the d-1 faces incident in a vertex
     *
     * @param v_pos an integer representing the position index of the vertex in the boundary array
     * @param fids a vector of integers that saves the position indexes of the faces incident in the vertex
     */
    inline void faces_in_vertex(int v_pos, ivect &fids)
    {
        if(this->get_vertices_num()==4)
            this->faces_in_vertex_quad(v_pos,fids);
        else if(this->get_vertices_num()==8)
            this->faces_in_vertex_hex(v_pos,fids);
    }
    /**
     * @brief A public method that returns the position indexes of the d-1 faces NOT incident in a vertex
     *
     * @param v_pos an integer representing the position index of the vertex in the boundary array
     * @param fids a vector of integers that saves the position indexes of the faces incident that are not incident in the vertex
     */
    inline void faces_NOT_in_vertex(int v_pos, ivect &fids)
    {
        if(this->get_vertices_num()==4)
            this->faces_NOT_in_vertex_quad(v_pos,fids);
        else if(this->get_vertices_num()==8)
            this->faces_NOT_in_vertex_hex(v_pos,fids);
    }

    /**
      * @brief A public method that returns the position indexes of the d-1 faces incident in a vertex
      * The procedure is defined only for the hexahedra! Undefined behavior otherwise.
      *
      * @param v_pos an integer representing the position index of the vertex in the boundary array
      * @param fids a vector of integers that saves the position indexes of the faces incident in the vertex
      */
    void edges_in_face(int f_pos, ivect &eids);

private:
    //////////////// QUAD SPECIFIC PROCEDURES ////////////////

    /**
     * @brief A private method that computes the edge at position pos.
     * The output is a sorted vector
     * WARNING: The procedure is specific for quads
     *
     * @param e a vector of int representing the edge
     * @param pos the position of the edge in the boundary of the quad
     */
    inline void QE(ivect &e, int pos)
    {
        QE_unsorted(e,pos);
        sort(e.begin(),e.end());
    }
    /**
     * @brief A private method that computes the edge at position pos.
     * The output is a unsorted vector
     * WARNING: The procedure is specific for quads
     *
     * @param e a vector of int representing the edge
     * @param pos the position of the edge in the boundary of the quad
     */
    inline void QE_unsorted(ivect &v, int pos)
    {
        v.clear();
        v.push_back(abs(this->vertices[(pos)%this->vertices.size()]));
        v.push_back(abs(this->vertices[(pos+1)%this->vertices.size()]));
    }
    /**
     * @brief A private method that returns the edges position indexes incident in a vertex
     * WARNING: The procedure is specific for quads
     *
     * @param v_pos an integer representing the position index of the vertex in the boundary
     * @param fids a vector of integers that saves the position indexes of the edges incident in the vertes
     */
    inline void faces_in_vertex_quad(int v_pos, ivect &fids)
    {
        fids.push_back(v_pos);
        fids.push_back((v_pos-1)%this->vertices.size());
    }
    /**
     * @brief A private method that returns the edges position indexes that are NOT incident in a vertex
     * WARNING: The procedure is specific for quads
     *
     * @param v_pos an integer representing the position index of the vertex in the boundary
     * @param fids a vector of integers that saves the position indexes of the edges incident in the vertes
     */
    inline void faces_NOT_in_vertex_quad(int v_pos, ivect &fids)
    {
        fids.push_back(v_pos+1);
        fids.push_back((v_pos+2)%this->vertices.size());
    }

    //////////////// HEXAHEDRA SPECIFIC PROCEDURES ////////////////

    /**
     * @brief A private method that computes the edge at position pos.     
     * The output is a sorted vector
     * WARNING: The procedure is specific for hexahedra
     *
     * @param e a vector of int representing the edge
     * @param pos the position of the edge in the boundary of the hexahedron
     */
    void HE(ivect &e, int pos);
    /**
     * @brief A private method that returns the quad face at position pos in the boundary of a hexahedron
     * The output is a sorted vector
     * WARNING: the procedure is defined only for hexahedra.
     *
     * @param f a vector of integers that saves the vertices in the boundary vertoces forming the quad
     * @param pos an integer representing the position index of the quad in the boundary of the hexahedron
     */
    inline void HF(ivect &f, int pos)
    {
        HF_unsorted(f,pos);
        sort(f.begin(),f.end());
    }
    /**
     * @brief A private method that returns the quad face at position pos in the boundary of a hexahedron
     * The output is a unsorted vector
     * WARNING: the procedure is defined only for hexahedra.
     *
     * @param f a vector of integers that saves the vertices in the boundary vertoces forming the quad
     * @param pos an integer representing the position index of the quad in the boundary of the hexahedron
     *
     * NOTA: the ordering for the d-1 faces is fixed.
     * F0 -> 0 1 4 5
     * F1 -> 0 1 3 6
     * F2 -> 0 2 3 5
     * F3 -> 1 4 6 7
     * F4 -> 2 3 6 7
     * F5 -> 2 4 5 7
     */
    void HF_unsorted(ivect &v, int pos);
    /**
     * @brief A private method that returns the quads position indexes incident in a vertex
     * WARNING: The procedure is specific for hexahedra.
     *
     * @param v_pos an integer representing the position index of the vertex in the boundary
     * @param fids a vector of integers that saves the position indexes of the quads incident in the vertes
     */
    void faces_in_vertex_hex(int v_pos, ivect &fids);
    /**
     * @brief A private method that returns the quads position indexes that are NOT incident in a vertex
     * WARNING: The procedure is specific for hexahedra.
     *
     * @param v_pos an integer representing the position index of the vertex in the boundary
     * @param fids a vector of integers that saves the position indexes of the quads incident in the vertes
     */
    void faces_NOT_in_vertex_hex(int v_pos, ivect &fids);
};

#endif // TOP_CP_CELL_H
