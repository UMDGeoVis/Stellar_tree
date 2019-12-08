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

#include "combinatorial_utilities.h"

namespace combinatorial_utilities
{

int compute_sub_type_number(int dim, int numV)
{
    // this function compute the following function to get all the possible d-sub-simplices
    // n! / (d! * (n-d)!)
    // where n is the actual dimension of the top simplex
    // and d is the sub-simplices that we want to extract
    mpz_t res; fact(numV,res);
    mpz_t res2; fact(dim,res2);
    mpz_t res3; fact(numV-dim,res3);
    mpz_t res4; mpz_init(res4); mpz_mul(res4,res2,res3);
    mpz_t res5; mpz_init(res5);
    mpz_tdiv_q(res5,res,res4);

    /// all the mpz_t variables are intrinsically pointers.. must be clared manually!
    mpz_clear(res);
    mpz_clear(res2);
    mpz_clear(res3);
    mpz_clear(res4);

    unsigned long long result = 0;
    mpz_export(&result, 0, -1, sizeof result, 0, 0, res5);
    mpz_clear(res5);
    return result;
}

void fact(int n, mpz_t p)
{
    int i;
    mpz_init_set_ui(p,1); /* p = 1 */
    for (i=1; i <= n ; ++i){
        mpz_mul_ui(p,p,i); /* p = p * i */
    }
}

std::size_t choose(double n, double k)
{
    using boost::math::binomial_coefficient;
    if (n < k) return 0;
    return static_cast<std::size_t>(binomial_coefficient<double>(n, k));
}

// Returns the largest n such that Choose(n, k) <= pos.
int combination_element(int k, int pos)
{
    int n = k;
    int coeff = 1;

    while (coeff <= pos)
    {
        coeff = choose(++n, k);
    }

    return n - 1;
}

// Returns an k-combination at position pos.
void combination(ivect &combination, int k, int pos)
{
    for (; k > 0; --k)
    {
        int n = combination_element(k, pos);
        combination.push_back(n);
        pos -= choose(n, k);
    }
}

}
