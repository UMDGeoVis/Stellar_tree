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

#include "top_simplex.h"

int Top_Simplex::get_sub_types_num(int d)
{
    switch(d)
    {
    case 0:
        return this->get_vertices_num();
    case 1:
        return this->get_edges_num(); //num edges
    case 2:
        if(this->get_vertices_num() == 3) // TRIANGLE
            return 1;
        else if(this->get_vertices_num() == 4) // TETRAHEDRON
            return this->get_dfaces_num(); //num faces
        else
            return this->compute_sub_types_num(d);
    case 3:
        if(this->get_vertices_num() == 4) // TETRAHEDRON
            return 1;
        else
            return this->compute_sub_types_num(d);
    default:
        if(this->get_vertices_num()-1 == d) // generic d-simplex
            return 1;
        else if(this->get_vertices_num()-1 > d) // we are extracting something valid
            return this->compute_sub_types_num(d);
        else
            return -1; // wrong d as argument
    }
}

void Top_Simplex::get_d_cell(ivect &s, int d, int pos)
{
    switch(d)
    {
    case 0:
        TV_vector(s,pos);
        break;
    case 1:
        this->TE(s,pos);
        break;
    case 2:
        if(this->get_vertices_num() == 4) // TETRAHEDRON
            this->TF(s,pos);
        else
            this->compute_d_simplex(s,d,pos); // this catches the TRI.. and returns the sorted vertices list forming the quad
        break;
    default:
        if(this->get_vertices_num()-2 == d)
            this->TF(s,pos);
        else
            this->compute_d_simplex(s,d,pos);
    }
}

int Top_Simplex::compute_sub_types_num(int d)
{
    int dim = d + 1; /// +1 because we need the real number of vertices that form the sub-simplex

    return combinatorial_utilities::compute_sub_type_number(dim,vertices.size());
}


void Top_Simplex::compute_d_simplex(ivect &s, int d, int pos)
{    
    int dim = d + 1; /// +1 because we need the real number of vertices that form the sub-simplex

    if(dim == this->get_vertices_num()) /// we need the simplex itself
        get_positive_sorted_vertices(s);
    else  if(this->get_vertices_num()-2 >= d)
    {
        compute_d_simplex_unsorted(s,dim,pos);
        sort(s.begin(),s.end());
    }
    else
    {
        cerr<<"[ERROR] [compute_d_cell] d is greater than the actual dimension of the simplex. num_v: "<<this->get_vertices_num()<<" d: "<<d<<endl;
    }
}
void Top_Simplex::compute_d_simplex_unsorted(ivect &s, int dim, int pos)
{
    s.clear();
    //modified version of combinatorial_utilities::Combination algorithm
    for (; dim > 0; --dim)
    {
        int n = combinatorial_utilities::combination_element(dim, pos);
        s.push_back(vertices[n]);
        pos -= combinatorial_utilities::choose(n, dim);
    }
}

void Top_Simplex::TetE(ivect &e, int pos)
{
    e.clear();
    switch(pos)
    {
    case(0):
        e.push_back(abs(vertices[0])); e.push_back(abs(vertices[1]));
        break;
    case(1):
        e.push_back(abs(vertices[0])); e.push_back(abs(vertices[2]));
        break;
    case(2):
        e.push_back(abs(vertices[0])); e.push_back(abs(vertices[3]));
        break;
    case(3):
        e.push_back(abs(vertices[1])); e.push_back(abs(vertices[2]));
        break;
    case(4):
        e.push_back(abs(vertices[1])); e.push_back(abs(vertices[3]));
        break;
    case(5):
        e.push_back(abs(vertices[2])); e.push_back(abs(vertices[3]));
        break;
    default:
        cerr<<"[TE] wrong edge position. must be from 0 to 5."<<endl;
        int a; cin>>a;
        break;
    }
    sort(e.begin(),e.end());
}
