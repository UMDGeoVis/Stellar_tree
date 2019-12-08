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

#ifndef SIMPLEXES_SOUP_H
#define SIMPLEXES_SOUP_H

#include "box.h"
#include "stellar_decomposition/top_simplex.h"
#include "mesh.h"
#include "explicit_cell.h"

using namespace std;

/**
 * @brief A class representing a generic mesh parametrized by the type of top simplexes
 *
 */
class Cells_Soup {
public:
    /**
     * @brief A constructor method
     *
     */
    Cells_Soup()
    {
        domain = Box<COORDBASETYPE>();
    }
    /**
     * @brief A copy-constructor method
     *
     * @param orig
     */
    Cells_Soup(const Cells_Soup& orig)
    {
        this->domain = orig.domain;
        this->cells = orig.cells;
    }
    /**
     * @brief A destructor method
     *
     */
    virtual ~Cells_Soup()
    {
        cells.clear();
    }
    /**
     * @brief A public method that clears the cells array
     *
     */
    inline void clear() { cells.clear(); }
    /**
     * @brief A public method that returns the soup domain
     *
     * @return Box&, the mesh domain
     */
    inline Box<COORDBASETYPE>& get_domain() { return this->domain; }
    /**
     * @brief A public method that sets the soup domain
     *
     * @param d a Box& argument, representing the domain to set
     */
    inline void set_domain(Box<COORDBASETYPE>& d) { this->domain = d; }
    /**
     * @brief A public method that initializes the top cells arrays
     *
     * @param dim a integer argument, representing the top cells types in the soup
     * @param num_top a ivect argument, representing the number of all top d-cells
     */
    inline void init_top_cells_array(int dim, ivect num_top)
    {
        cells.assign(dim,vector<Explicit_Cell<COORDBASETYPE> >());
        for(int d=0; d<dim; d++)
            cells[dim].reserve(num_top[dim]);
    }
    /**
     * @brief A public method that adds a top p-cell to the corresponding array
     *
     * @param pos an integer argument representing the dimension of the top p-cell
     * @param t an ExplicitCell& argument containing the top p-cell to insert
     */
    inline void add_cell(int pos, Explicit_Cell<COORDBASETYPE>& t) { this->cells[pos].push_back(t); }
    /**
     * @brief A public method that returns a top p-cell at a given position
     *
     * @param pos an integer argument representing the dimension of the top p-cell
     * @param id an integer argument representing the position index of the top p-cell in the corresponding array
     * @return ExplicitCell& representing top p-cell at the given position
     */
    inline Explicit_Cell<COORDBASETYPE>& get_cell(int pos, int id) { return this->cells[pos][id-1]; }
    /**
     * @brief A public method that returns the number of top cells types encoded by the soup
     *
     * @return an integer representing the number of top cells types
     */
    inline int get_cells_types_num() { return this->cells.size(); }    
    /**
     * @brief A public methot that returns the number of top p-cells encoded by the soup
     *
     * @param pos an integer argument representing the dimension of the top p-cells
     * @return int an integer containing the number of top p-cells
     */
    inline int get_cells_num(int pos) { return this->cells[pos].size(); }
    /**
     * @brief A public method that adds of type of top p-cells
     *
     * @param type a CellType argument representing the type
     */
    inline void add_type(CellType type) { this->top_type.push_back(type); }
    /**
     * @brief A public method that returns a top cell type at a given dimension
     *
     * @param d an integer argument representing the top cell dimension
     * @return CellType containing the type
     */
    inline CellType get_type(int d) const { return this->top_type[d]; }    
    /**
     * @brief A public method that saves in a string the corresponding top d-cell type
     *
     * @param d an integer argument representing the top cell dimension
     * @return string containing the name of the type
     */
    inline string print_type(int d) { return (this->top_type[d] != UNSET? cell_names[this->top_type[d]] : "UNSET"); }
    /**
     * @brief A public method that checks if the the top cells in a given array are simplices
     *
     * @param d, the position index of the array
     * @return bool, true if the top cells in the d-array are simplices, false otherwise
     */
    inline bool is_simplex(int d) { return (this->top_type[d] == LINE || this->top_type[d] == TRIANGLE || this->top_type[d] == TETRA || this->top_type[d] == nDIM_SIMPL); }

protected:
    ///A private varible representing the soup domain
    Box<COORDBASETYPE> domain;
    ///A private varible representing the top simplexes list of the soup
    vector<vector<Explicit_Cell<COORDBASETYPE> > > cells;
    ///A private variable representing the d-top types
    vector<CellType> top_type;
};

#endif // SIMPLEXES_SOUP_H
