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

#include "top_cp_cell.h"

int Top_CP_Cell::get_edges_num()
{
    if(this->get_vertices_num() == 4) // QUAD
        return this->get_vertices_num();
    else if(this->get_vertices_num() == 8) // HEXAHEDRON
        return 12;
    else /// GENERIC CP_CELL
    {
        cerr<<"[Top_CP_Cell](get_edges_num) undefined behavior"<<endl;
        return -1;
    }
}

int Top_CP_Cell::get_dfaces_num()
{
    if(this->get_vertices_num() == 4) // QUAD
        return this->get_vertices_num();
    else if(this->get_vertices_num() == 8) // HEXAHEDRON
        return 6;
    else /// GENERIC CP_CELL
    {
        cerr<<"[Top_CP_Cell](get_dfaces_num) undefined behavior"<<endl;
        return -1;
    }
}

int Top_CP_Cell::get_sub_types_num()
{
    if(this->get_vertices_num() == 4) // QUAD
        return 2; /// edges + vertices
    else if(this->get_vertices_num() == 8) // HEXAHEDRON
        return 3; /// faces + edges + vertices
    else /// GENERIC CP_CELL
    {
        cerr<<"[Top_CP_Cell](get_sub_simplices_types_num) undefined behavior"<<endl;
        return -1;
    }
}

int Top_CP_Cell::get_sub_types_num(int d)
{
    switch(d)
    {
    case 0:
        return this->get_vertices_num();
    case 1:
        return get_edges_num(); ///num edges
    case 2:
        if(this->get_vertices_num() == 4) // QUAD
            return 1;
        else if(this->get_vertices_num() == 8) // HEXAHEDRON
            return 6; ///num faces
    case 3:
        if(this->get_vertices_num() == 8) // HEXAHEDRON
            return 1;
    default:/// GENERIC CP_CELL
        cerr<<"[Top_CP_Cell](get_sub_types_num) undefined behavior"<<endl;
        return -1;
    }
}

void Top_CP_Cell::get_d_cell(ivect &s, int d, int pos)
{
    switch(d)
    {
    case 0:
        TV_vector(s,pos);
        break;
    case 1:
        TE(s,pos);
        break;
    case 2:
        if(this->get_vertices_num() == 8) // HEXAHEDRON
            this->HF(s,pos);
        else if(this->get_vertices_num() == 4) // QUAD
            this->get_positive_sorted_vertices(s); /// this returns the sorted vertices list forming the quad
        break;
    case 3:
        this->get_positive_sorted_vertices(s); // this returns the sorted vertices list forming the hexahedron
    }
}

void Top_CP_Cell::edges_in_face(int f_pos, ivect &eids)
{
    if(this->get_vertices_num()==8)
    {
        eids.assign(4,0);
        switch (f_pos) {
        case 0:
            eids[0] = 0; eids[1] = 2; eids[2] = 3; eids[3] = 9;
            break;
        case 1:
            eids[0] = 0; eids[1] = 1; eids[2] = 4; eids[3] = 8;
            break;
        case 2:
            eids[0] = 1; eids[1] = 2; eids[2] = 5; eids[3] = 6;
            break;
        case 3:
            eids[0] = 3; eids[1] = 4; eids[2] = 10; eids[3] = 11;
            break;
        case 4:
            eids[0] = 5; eids[1] = 7; eids[2] = 7; eids[3] = 11;
            break;
        case 5:
            eids[0] = 6; eids[1] = 7; eids[2] = 9; eids[3] = 10;
            break;
        default:
            break;
        }
    }
}

void Top_CP_Cell::HE(ivect &e, int pos) // hexahedron only
{
    e.assign(2,-1);
    switch (pos)
    {
    case 0:
        e[0] = abs(vertices[0]); e[1] = abs(vertices[1]);
        break;
    case 1:
        e[0] = abs(vertices[0]); e[1] = abs(vertices[3]);
        break;
    case 2:
        e[0] = abs(vertices[0]); e[1] = abs(vertices[5]);
        break;
    case 3:
        e[0] = abs(vertices[1]); e[1] = abs(vertices[4]);
        break;
    case 4:
        e[0] = abs(vertices[1]); e[1] = abs(vertices[6]);
        break;
    case 5:
        e[0] = abs(vertices[2]); e[1] = abs(vertices[3]);
        break;
    case 6:
        e[0] = abs(vertices[2]); e[1] = abs(vertices[5]);
        break;
    case 7:
        e[0] = abs(vertices[2]); e[1] = abs(vertices[7]);
        break;
    case 8:
        e[0] = abs(vertices[3]); e[1] = abs(vertices[6]);
        break;
    case 9:
        e[0] = abs(vertices[4]); e[1] = abs(vertices[5]);
        break;
    case 10:
        e[0] = abs(vertices[4]); e[1] = abs(vertices[7]);
        break;
    case 11:
        e[0] = abs(vertices[6]); e[1] = abs(vertices[7]);
        break;
    default:
        break;
    }
    sort(e.begin(),e.end());
}

void Top_CP_Cell::HF_unsorted(ivect &v, int pos)
{
    v.assign(4,-1);
    switch (pos)
    {
    case 0:
        v[0] = abs(vertices[0]); v[1] = abs(vertices[1]); v[2] = abs(vertices[4]); v[3] = abs(vertices[5]);
        break;
    case 1:
        v[0] = abs(vertices[0]); v[1] = abs(vertices[1]); v[2] = abs(vertices[3]); v[3] = abs(vertices[6]);
        break;
    case 2:
        v[0] = abs(vertices[0]); v[1] = abs(vertices[2]); v[2] = abs(vertices[3]); v[3] = abs(vertices[5]);
        break;
    case 3:
        v[0] = abs(vertices[1]); v[1] = abs(vertices[4]); v[2] = abs(vertices[6]); v[3] = abs(vertices[7]);
        break;
    case 4:
        v[0] = abs(vertices[2]); v[1] = abs(vertices[3]); v[2] = abs(vertices[6]); v[3] = abs(vertices[7]);
        break;
    case 5:
        v[0] = abs(vertices[2]); v[1] = abs(vertices[4]); v[2] = abs(vertices[5]); v[3] = abs(vertices[7]);
        break;
    default:
        break;
    }
}

void Top_CP_Cell::faces_in_vertex_hex(int v_pos, ivect &fids)
{
    fids.assign(3,0);
    switch (v_pos) {
    case 0:
        fids[0] = 0; fids[1] = 1; fids[2] = 2;
        break;
    case 1:
        fids[0] = 0; fids[1] = 1; fids[2] = 3;
        break;
    case 2:
        fids[0] = 2; fids[1] = 4; fids[2] = 5;
        break;
    case 3:
        fids[0] = 1; fids[1] = 2; fids[2] = 4;
        break;
    case 4:
        fids[0] = 0; fids[1] = 3; fids[2] = 5;
        break;
    case 5:
        fids[0] = 0; fids[1] = 2; fids[2] = 5;
        break;
    case 6:
        fids[0] = 1; fids[1] = 3; fids[2] = 4;
        break;
    case 7:
        fids[0] = 3; fids[1] = 4; fids[2] = 5;
        break;
    default:
        break;
    }
}

void Top_CP_Cell::faces_NOT_in_vertex_hex(int v_pos, ivect &fids)
{
    fids.assign(3,0);
    switch (v_pos) {
    case 0:
        fids[0] = 3; fids[1] = 4; fids[2] = 5;
        break;
    case 1:
        fids[0] = 2; fids[1] = 4; fids[2] = 5;
        break;
    case 2:
        fids[0] = 0; fids[1] = 1; fids[2] = 3;
        break;
    case 3:
        fids[0] = 0; fids[1] = 3; fids[2] = 5;
        break;
    case 4:
        fids[0] = 1; fids[1] = 2; fids[2] = 4;
        break;
    case 5:
        fids[0] = 1; fids[1] = 3; fids[2] = 4;
        break;
    case 6:
        fids[0] = 0; fids[1] = 2; fids[2] = 5;
        break;
    case 7:
        fids[0] = 0; fids[1] = 1; fids[2] = 2;
        break;
    default:
        break;
    }
}
