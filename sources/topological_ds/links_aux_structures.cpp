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

#include "links_aux_structures.h"

void s_link::extract_subs(ivect &v_in_v_link, Top_Simplex &)
{
    /// (3) if the p-simplex is at least a triangle, then we add all the sub-simplices in its boundary
    ///     NOTA: we do not consider the 0-simplices and the p-simplex as they have been already processed
    queue<ivect> q;
    for(ivect_iter it=v_in_v_link.begin(); it!=v_in_v_link.end(); ++it)
    {
        int v_id = *it;
        ivect sub_s(v_in_v_link.size()-1);
        copy_if(v_in_v_link.begin(),v_in_v_link.end(),sub_s.begin(), [v_id](int i){return (i!=v_id);} );
        pair<set<ivect>::iterator,bool> ret = insert_simplex(sub_s);
        if(ret.second && sub_s.size() > 2) /// the vertices are already into the link
            q.push(sub_s);
    }

    while(!q.empty())
    {
        ivect current = q.front();
        q.pop();
        for(ivect_iter it=current.begin(); it!=current.end(); ++it)
        {
            int v_id = *it;
            ivect sub_s(current.size()-1);
            copy_if(current.begin(),current.end(),sub_s.begin(), [v_id](int i){return (i!=v_id);} );
            pair<set<ivect>::iterator,bool> ret = insert_simplex(sub_s);
            if(ret.second && sub_s.size() > 2) /// the vertices are already into the link
                q.push(sub_s);
        }
    }
}

void s_link::extract_subs(ivect &v_in_v_link, Top_CP_Cell &)
{
    /// (3) if the p-cell is at least a quad, then we add all the sub-cells in its boundary
    ///     NOTA: we do not consider the 0-cells and the p-cell as they have been already processed
    Top_CP_Cell sub = Top_CP_Cell(v_in_v_link);
    for(int s=1; s<sub.get_sub_types_num(); s++)
    {
        int sub_simpl_num = sub.get_sub_types_num(s);
        for(int i=0; i<sub_simpl_num; i++)
        {
            ivect sub_s;
            sub.get_d_cell(sub_s,s,i);
            insert_cell(sub_s);
        }
    }
}

void s_link::print_link()
{
    cout<<"v in link: "<<this->v_in_link.size()<<endl;
    for(iset_iter it=this->v_in_link.begin(); it!=v_in_link.end(); ++it)
        cout<<*it<<" ";
    cout<<endl;
    int d=0;
    for(vector< set<ivect > >::iterator it=this->s_in_link.begin(); it!=s_in_link.end(); ++it)
    {
        if(it->size() > 0)
        {
            cout<<d<<"-simplices num: "<<it->size()<<endl;
            for(set<ivect >::iterator it2=it->begin(); it2!=it->end(); ++it2)
            {
                for(ivect_const_iter it3=it2->begin(); it3!=it2->end(); ++it3)
                    cout<<*it3<<" ";
                cout<<endl;
            }
        }
    }
}
