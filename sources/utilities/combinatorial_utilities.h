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

#ifndef COMBINATORIAL_UTILITIES_H
#define COMBINATORIAL_UTILITIES_H

#include "gmp.h"
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <algorithm>

#include "utilities/basic_wrappers.h"

namespace combinatorial_utilities
{
/**
 * @brief A procedure that computes the number of all sub-simplices
 * this function compute the following function to get all the possible d-sub-simplices
 * n! / (d! * (n-d)!)
 * where n is the actual dimension of the top simplex
 * and d is the sub-simplices that we want to extract
 *
 * NOTA: the function requires the GMP library, a library for arbitrary precision arithmetic.
 *
 * @param dim the dimension of the simplex
 * @param numV an integer representing the number of vertices in a simplex
 * @return int the number of sub-simplices
 */
extern int compute_sub_type_number(int dim, int numV);
/**
 * @brief A procedure that computes the factorial of a number
 *
 * @param n the number on which we compute the factorial
 * @param p the final resulting factorial
 */
extern void fact(int n, mpz_t p);
/**
 * @brief A procedure that computes the binomial coefficient of two numbers
 *
 * @param n
 * @param k
 * @return std::size_t
 */
extern std::size_t choose(double n, double k) ;
/**
 * @brief A procedure that returns the largest n such that Choose(n, k) <= pos.
 *
 * @param k
 * @param pos
 * @return int
 */
extern int combination_element(int k, int pos);
/**
 * @brief A procedure that returns an k-combination at position pos
 *
 * @param combination an integer vector containing the combination
 * @param k an integer referring to the size of the combination
 * @param pos the starting position for the combination algorithm
 */
extern void combination(ivect &combination, int k, int pos);
}

#endif // COMBINATORIAL_UTILITIES_H
