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

#ifndef CONNECTEDNESS_VALIDATOR_STATS_H
#define CONNECTEDNESS_VALIDATOR_STATS_H

#include <vector>
#include <map>
#include <deque>
#include <iostream>

#include "utilities/basic_wrappers.h"

using namespace std;

class Connectedness_Validator_stats
{
public:
    Connectedness_Validator_stats() { is_pseudo_manifold = true; edges_num = 0; ve_size = 0; }

    ///
    int edges_num, ve_size;
    /// debug
    double time_extraction, time_init, time_visit, time_simpl_gath, time_set_flags/*, time_adding_to_local_struct*/;
    /// array of dimension |V| at runtime that stores at which component of a 1-skeleton each vertex belong (0-connectedness test)
    /// or array of dimension |T| during the d-connectedness test
    ivect cc;
    /// contains the connected components identified and the list of simplices in that conn.comps
    map<int,deque<int> > connected_components;
    ///
    bool is_pseudo_manifold;

    inline void reset_time_variables() { time_extraction = time_init = time_visit = time_simpl_gath = time_set_flags = /*time_adding_to_local_struct =*/ 0; }
    inline void print_time_variables()
    {
        if(time_init > 0)
            cerr << " initialization: " << time_init << endl;
        if(time_extraction > 0)
            cerr << " relations extraction: " << time_extraction << endl;
        if(time_simpl_gath > 0)
            cerr << " structure set-up: " << time_simpl_gath << endl;
        if(time_visit > 0)
            cerr << " visit structure: " << time_visit << endl;
        if(time_set_flags > 0)
            cerr << " setting flags: " << time_set_flags << endl;
    }
    ///
    ///
    inline void print_validation_result()
    {
        if(connected_components.size() == 1)
        {
            cerr<<"    The model is CONNECTED"<<endl;
        }
        else
        {
            cerr<<"    The model has "<<connected_components.size()<<" connected components "<<endl;

            /// check if we have some simplices that are not been visited
            ivect unvisited;
            for(unsigned i=0; i<cc.size(); i++)
            {
                if(cc[i]==0)
                {
                    unvisited.push_back(i+1);
                }
            }
            if(unvisited.size() > 0)
            {
                cout<<"[validation report] "<<unvisited.size()<<" entities that are not touched during validation.."<<endl;
            }
        }
    }

    inline void print_pseudo_manifold_result()
    {
        if(is_pseudo_manifold && connected_components.size() == 1)
            cerr<<"The model is PSEUDO-MANIFOLD"<<endl;
        else
            cerr<<"The model is NOT PSEUDO-MANIFOLD"<<endl;
    }

    inline void print_skeleton_stats(bool debug)
    {
        if(debug)
        {
            cerr<<"MAX edges number: "<<this->edges_num<<" ("<<((4*2*this->edges_num) / (1024.0*1024.0))<<" MBs)"<<endl;
            cerr<<"VE entries: "<<this->ve_size<<" ("<<((4*this->ve_size) / (1024.0*1024.0))<<" MBs)"<<endl;
            cerr<<"==> TOT local skeleton storage: "<<((4*2*this->edges_num+4*this->ve_size) / (1024.0*1024.0))<<" MBs"<<endl;
        }
    }

    inline void reset()
    {
        cc.clear();
        //        s_to_cc.clear();
        connected_components.clear();
        is_pseudo_manifold = true;
    }

};

#endif // CONNECTEDNESS_VALIDATOR_STATS_H
