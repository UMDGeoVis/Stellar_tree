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

#include "iastar_generator.h"

void IAstar_Generator::visit_adj_queue_G(int top_id, int d, int real_v_index, bm::bvector<> &bv, Simplicial_Mesh &mesh, IAstar &gia)
{
    iqueue q;
    q.push(top_id);

    while(!q.empty())
    {
        int current = q.front();
        Top_Simplex &t = mesh.get_top_cell(d,current);
        int k = t.vertex_index(real_v_index);


        for(int f=1; f<t.get_dfaces_num(); f++)
        {
            int adj = gia.get_TT(d,current,(k+f)%t.get_dfaces_num()); /// we consider only the incident faces in current
            gia.check_adj(adj,d,q,bv);
        }
        q.pop();
    }
}

void IAstar_Generator::visit_adj_queue_G(int top_id, int d, int real_v_index, bm::bvector<> &bv, CP_Mesh &mesh, IAstar &gia)
{
    iqueue q;
    q.push(top_id);

    while(!q.empty())
    {
        int current = q.front();
        Top_CP_Cell &t = mesh.get_top_cell(d,current);
        int k = t.vertex_index(real_v_index);
        ivect fids;
        t.faces_in_vertex(k,fids);
        for(unsigned f=0; f<fids.size(); f++)
        {
            int adj = gia.get_TT(d,current,fids[f]); /// we consider only the incident faces in current
            gia.check_adj(adj,d,q,bv);
        }
        q.pop();
    }
}

void IAstar_Generator::visit_adj_queue_L(int top_id, int d, int real_v_index, bm::bvector<> &bv, Simplicial_Mesh &mesh, IAstar_local &gia)
{
    iqueue q;
    q.push(top_id);

    while(!q.empty())
    {
        int current = q.front();
        Top_Simplex &t = mesh.get_top_cell(d,gia.get_index(d,current));
        int k = t.vertex_index(real_v_index);

        for(int f=1; f<t.get_dfaces_num(); f++)
        {
            int adj = gia.get_TT(d,current,(k+f)%t.get_dfaces_num()); /// we consider only the incident faces in current
            gia.check_adj(adj,d,q,bv);
        }
        q.pop();
    }
}

void IAstar_Generator::visit_adj_queue_L(int top_id, int d, int real_v_index, bm::bvector<> &bv, CP_Mesh &mesh, IAstar_local &gia)
{
    iqueue q;
    q.push(top_id);

    while(!q.empty())
    {
        int current = q.front();
        Top_CP_Cell &t = mesh.get_top_cell(d,gia.get_index(d,current));
        int k = t.vertex_index(real_v_index);
        ivect fids;
        t.faces_in_vertex(k,fids);
        for(unsigned f=0; f<fids.size(); f++)
        {
            int adj = gia.get_TT(d,current,fids[f]); /// we consider only the incident faces in current
            gia.check_adj(adj,d,q,bv);
        }
        q.pop();
    }
}
