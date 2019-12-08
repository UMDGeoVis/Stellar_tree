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

#include "iastar.h"

void IAstar::save_IAstar(string filename, ivect &original_v_pos, vector<ivect > &original_t_pos)
{
    ofstream output(filename.c_str());
    output << "CONNECTIVITY" << endl;

    /// writing the header
    output << this->VTstar_rel.size() << " ";
    output << this->TT_rel.size() << " ";
    for(unsigned i=0; i<this->TT_rel.size(); i++)
        output << this->TT_rel[i].size() << " ";
    output << this->face_coboundary.size() << " ";
    for(unsigned i=0; i<this->face_coboundary.size();i++)
        output << this->face_coboundary[i].size() << " ";
    output << endl;

    /// init local structure
    vector<Vertex_Coboundary> VTstar_rel_orig;
    VTstar_rel_orig = this->VTstar_rel;
    for(unsigned i=0; i<this->VTstar_rel.size(); i++)
        VTstar_rel_orig[original_v_pos[i]-1] = this->VTstar_rel[i];

    ///writing the partial incidence relations
    for(unsigned i=0; i<VTstar_rel_orig.size(); i++)
    {
        vector<ivect > &pVT = VTstar_rel_orig[i].get_partial_VT();
        for(unsigned d=0; d<pVT.size();d++)
        {
            for(unsigned i=0; i<pVT[d].size(); i++)
            {
                output<<original_t_pos[d][pVT[d][i]-1];
                if(i+1 < pVT[d].size())
                    output<<" ";
            }
            output<<" # ";
        }
        output<<endl;
    }
    VTstar_rel_orig.clear();

    /// init local structure
    vector<vector<Top_Adj> > TT_rel_orig;
    TT_rel_orig = this->TT_rel;
    for(unsigned i=0; i<this->TT_rel.size(); i++)
        for(unsigned j=0; j<this->TT_rel[i].size(); j++)
            TT_rel_orig[i][original_t_pos[i][j]-1] = this->TT_rel[i][j];

    /// writing the adjacency relations
    for(unsigned i=0; i<TT_rel_orig.size(); i++)
    {
        for(unsigned j=0; j<TT_rel_orig[i].size(); j++)
        {
            for(int w=0; w<TT_rel_orig[i][j].size(); w++)
            {
                if(TT_rel_orig[i][j].get_adj(w) >= 0)
                    output<<original_t_pos[i][TT_rel_orig[i][j].get_adj(w)-1]<<" ";
                else /// nm adj
                    output<<TT_rel_orig[i][j].get_adj(w)<<" ";
            }
            output<<endl;
        }
    }

    /// writing the non-manifold adjacency relations (i.e. face coboundary)
    for(unsigned i=0; i<this->face_coboundary.size();i++)
    {
        map<int,ivect > &m = this->face_coboundary[i];
        for(map<int,ivect >::iterator it=m.begin(); it!=m.end(); ++it)
        {
            output << it->first << " ";
            for(ivect_iter it2=it->second.begin(); it2!=it->second.end(); ++it2)
                output << original_t_pos[i][*it2-1] << " ";
            output<< "# " << endl;
        }
    }

    output.close();
}

void IAstar::print_stats()
{
    cerr<<"max VTstar size: "<<this->get_VTstar_Num()<<endl;
    int nm_vertices = 0;
    for(int i=1; i<=this->get_VTstar_Num(); i++)
    {
        if(this->get_VTstar(i).non_manifold_vertex())
            nm_vertices++;
    }
    for(int i=0; i<this->get_TT_Num(); i++)
    {
        cerr<<"max TT"<<i<<" size: "<<this->get_TT_Num(i)<<endl;
    }

    cerr<<"--non-manifold stats--"<<endl;
    cerr<<"   non-manifold vertices: "<<nm_vertices<<endl;
    for(int i=0; i<this->get_nmTT_Num(); i++)
    {
        cerr<<"   non-manifold "<<i+1<<"-simplices: "<<this->get_nmTT_Num(i)<<endl;
    }
}

void IAstar::print_non_manifold_adjacencies()
{
    for(unsigned d=0; d<this->face_coboundary.size(); d++)
    {
        map<int,ivect > &m = this->face_coboundary[d];
        for(map<int,ivect >::iterator it=m.begin(); it!=m.end(); ++it)
        {
            cout<<it->first<<"] ";
            for(ivect_iter it2=it->second.begin(); it2!=it->second.end(); ++it2)
                cout<<*it2<<" ";
            cout<<endl;
        }
    }
}

void IAstar::get_res_and_next_VTop(int v_id, int current, Top_Simplex &top, int d, iqueue &q, bm::bvector<> &bv, VT &res)
{
    int pos = top.vertex_index(v_id);

    for(int i=1; i<top.get_vertices_num(); i++)
        check_adj(this->get_TT(d,current,(int)(pos+i)%top.get_vertices_num()),d,q,bv,res);
}

void IAstar::get_res_and_next_VTop(int v_id, int current, Top_CP_Cell &top, int d, iqueue &q, bm::bvector<> &bv, VT &res)
{
    int pos = top.vertex_index(v_id);
    ivect fids;
    top.faces_in_vertex(pos,fids);
    for(unsigned f=0; f<fids.size(); f++)
    {
        check_adj(this->get_TT(d,current,fids[f]),d,q,bv,res);
    }
}

void IAstar::get_res_and_next_VV(int pos, int current, Top_Simplex &top, int d, iqueue &q, bm::bvector<> &bv, iset &res)
{
    for(int i=1; i<top.get_vertices_num(); i++)
    {
        res.insert(top.TV((int)(pos+i)%top.get_vertices_num()));

        check_adj(this->get_TT(d,current,(int)(pos+i)%top.get_vertices_num()),d,q,bv);
    }
}

void IAstar::get_res_and_next_VV(int pos, int current, Top_CP_Cell &top, int d, iqueue &q, bm::bvector<> &bv, iset &res)
{
    for(int i=0; i<top.get_vertices_num(); i++)
    {
        if(i != pos)
            res.insert(top.TV(i));
    }

    ivect fids;
    top.faces_in_vertex(pos,fids);
    for(unsigned f=0; f<fids.size(); f++)
    {
        check_adj(this->get_TT(d,current,fids[f]),d,q,bv);
    }
}

void IAstar::check_adj(int adj, int d, iqueue &q, bm::bvector<> &bv, vector<ivect > &res)
{
    if(adj > UNSETADJ && !bv[adj]) /// manifold adjacency
    {
        bv.set(adj);
        q.push(adj);
        res[d].push_back(adj);
    }
    else if(adj < UNSETADJ)/// non-manifold adjacency
    {
        for(int ad=0; ad<this->get_nmTT_Num(d,adj); ad++)
        {
            if(!bv[this->get_nmTT(d,adj,ad)])
            {
                bv.set(this->get_nmTT(d,adj,ad));
                q.push(this->get_nmTT(d,adj,ad));
                res[d].push_back(this->get_nmTT(d,adj,ad));
            }
        }
    }
}

