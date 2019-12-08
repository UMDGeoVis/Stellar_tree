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

#ifndef STATISTICS_IASTAR_H
#define STATISTICS_IASTAR_H

#include <vector>
#include "utilities/basic_wrappers.h"

class IAstar_stats
{
public:
    IAstar_stats(int num_simplices, bool debug) /// local approach
    {
        if(debug)
        {
            this->max_queue_size = 0;
            max_local_vt = 0;
            max_local_tt.assign(num_simplices,0);
            max_local_nm.assign(num_simplices,0);
            ///
            time_tt = 0;
            time_vt = 0;
            time_init = 0;
            time_gath = 0;
        }
    }
    IAstar_stats() /// global approach
    {
        this->max_queue_size = 0;
        ///
        max_local_vt = 0;
        ///
        time_tt = 0;
        time_vt = 0;
        time_init = 0;
        time_gath = 0;
    }

    /// for stats gathering (global and local approaches)
    int max_queue_size;

    /// for stats gathering (local approach only)
    int max_local_vt;
    ivect max_local_tt;
    ivect max_local_nm;

    /// for debug only
    /// we consider the timings for the different part of the application
    double time_init, time_gath, time_tt, time_vt;


    inline void reset_stats(int num_simplices)
    {
        this->max_queue_size = 0;
        max_local_vt = 0;
        max_local_tt.assign(num_simplices,0);
        max_local_nm.assign(num_simplices,0);

        ///
        time_tt = 0;
        time_vt = 0;
        time_init = 0;
        time_gath = 0;
    }

    inline void print_global_stats(IAstar &gia)
    {
        cerr<<"max_queue_size: "<<max_queue_size<<endl;
        gia.print_stats();
    }

    inline void print_local_stats()
    {
        cerr<<"max_queue_size: "<<max_queue_size<<endl;

        cerr<<"max VTstar size: "<<max_local_vt<<endl;

        for(unsigned i=0; i<max_local_tt.size(); i++)
        {
            cerr<<"max TT"<<i<<" size: "<<max_local_tt[i]<<endl;
        }
        for(unsigned i=0; i<max_local_nm.size(); i++)
        {
            cerr<<"max nm faces"<<i<<" size: "<<max_local_nm[i]<<endl;
        }
    }

    inline void print_partial_timings()
    {
        cerr<<"time_init: "<<time_init<<" time_gath; "<<time_gath<<endl;
        cerr<<"time_tt: "<<time_tt<<" time_vt; "<<time_vt<<endl;
    }
};

#endif // STATISTICS_IASTAR_H
