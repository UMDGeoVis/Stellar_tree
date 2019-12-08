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

#ifndef HALFEDGE_H
#define HALFEDGE_H

class HalfEdge
{
public:
    HalfEdge() { v = f = n = o = -1; }

    inline void set_pointed_vertex(int v) { this->v = v; }
    inline void set_reference_face(int f) { this->f = f; }
    inline void set_next_half_edge(int n) { this->n = n; }
    inline void set_opposite_half_edge(int o) { this->o = o; }

    inline int get_pointed_vertex() { return v; }
    inline int get_reference_face() { return f; }
    inline int get_next_half_edge() { return n; }
    inline int get_opposite_half_edge() { return o; }

    inline void print()
    {
        cerr<<"vertex: "<<v<<" face: "<<f<<" next: "<<n<<" opposite: "<<o<<endl;
    }

private:
    int v; /// the vertex pointed by the half edge
    int f; /// the face
    int n; /// the next half edge
    int o; /// the opposite half edge
};

#endif // HALFEDGE_H
