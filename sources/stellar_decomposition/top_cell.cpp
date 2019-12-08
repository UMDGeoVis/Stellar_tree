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

#include "top_cell.h"

int Top_Cell::min_v_index()
{
    int min=abs(vertices[0]);
    for(int i=1; i<get_vertices_num(); i++)
    {
        if(min > abs(vertices[i]))
            min=abs(vertices[i]);
    }
    return min;
}

int Top_Cell::max_v_index()
{
    int max=abs(vertices[0]);
    for(int i=1; i<get_vertices_num(); i++)
    {
        if(max < abs(vertices[i]))
            max=abs(vertices[i]);
    }
    return max;
}

bool Top_Cell::has_cell(const ivect &c)
{
    boost::dynamic_bitset<> b(c.size());

    for(unsigned i=0;i<c.size();i++)
    {
        for(int j=0;j<get_vertices_num();j++)
        {
            if(c[i]==abs(vertices[j]))
            {
                b[i] = 1;
                break;
            }
        }
    }
    return b.count() == b.size();
}

bool Top_Cell::are_equal(const ivect &c)
{
    if((unsigned)this->get_vertices_num() != c.size())
        return false;

    boost::dynamic_bitset<> b(this->get_vertices_num());

    for(int i=0;i<this->get_vertices_num();i++)
    {
        for(unsigned j=0;j<c.size();j++)
        {
            if(!b[j] && abs(this->vertices[i])==c[j])
            {
                b[j] = 1;
                break;
            }
        }
    }

    return b.count() == b.size();
}
