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

#include "connectedness_validator_adjs.h"

void connVal_d_1_adj::set_adjacencies()
{
    for(unsigned d=0; d<all_faces.size(); d++)
    {
        face_top_vector &d_faces = all_faces[d];

        if(d_faces.size() > 0)
        {
            sort_container(d_faces);

            for(unsigned j=0; j<d_faces.size(); j++)
            {
                if(j+1<d_faces.size() && d_faces[j] == d_faces[j+1])
                {
                    /// se ho un'adiacenza manifold adj_counter e' 1, altrimenti contiene il numero di tetra
                    int adj_counter = 1;
                    while((j+1+adj_counter)<d_faces.size() && d_faces[j] == d_faces[j+1+adj_counter] )
                        adj_counter++;

                    if(adj_counter==1)
                    {
                        this->get_adj(d_faces[j].get_t_id()).set_adj(d_faces[j].get_f_pos(),d_faces[j+1].get_t_id());
                        this->get_adj(d_faces[j+1].get_t_id()).set_adj(d_faces[j+1].get_f_pos(),d_faces[j].get_t_id());
                    }
                    else
                    {
                        for(unsigned ad=j; ad<=(j+adj_counter); ad++)
                        {
                            for(unsigned bd=ad+1; bd<=(j+adj_counter); bd++)
                            {
                                this->get_adj(d_faces[ad].get_t_id()).set_adj(d_faces[ad].get_f_pos(),d_faces[bd].get_t_id());
                                this->get_adj(d_faces[bd].get_t_id()).set_adj(d_faces[bd].get_f_pos(),d_faces[ad].get_t_id());
                            }
                        }
                    }

                    j+=adj_counter;
                }

                if(j>=d_faces.size())
                    break;
            }
        }
    }

    /// erase the temporary structure
    all_faces.clear();
}

