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

#ifndef USAGE_H
#define USAGE_H

#include <sys/types.h>
#include <sys/resource.h>

#include <stdlib.h>
#include <stdio.h>
#include <string>

using namespace std;

/**
 * @brief The MemoryUsage class keeps track of the memory (peak) usage at run-time
 * The class is specified to work on both linux and apple systems
 */
class MemoryUsage
{

private:
    int who = RUSAGE_SELF;
    struct rusage usage;
    int ret;

public:
    /**
     * @brief A public method that return the current memory peak in kilobytes
     * @param print a boolean, if true prints on the standard output a string containing the value of the memory peak
     * @return a float containing the memory peak
     */
    inline float getValue_in_KB(bool print)
    {
         ret=getrusage(who,&usage);
#ifdef __APPLE__
         if(print) cout << "Memory Usage: " << usage.ru_maxrss/(1024.0) << " KB" << endl;
         return usage.ru_maxrss/(1024.0);
#else
         if(print) cout << "Memory Usage: " << usage.ru_maxrss << " KB" << endl;
         return usage.ru_maxrss;
#endif
    }
    /**
     * @brief A public method that return the current memory peak in megabytes
     * @param print a boolean, if true prints on the standard output a string containing the value of the memory peak
     * @return a float containing the memory peak
     */
    inline float getValue_in_MB(bool print)
    {
         ret=getrusage(who,&usage);
#ifdef __APPLE__
         if(print) cout << "Memory Usage: " << usage.ru_maxrss/(1024.0*1024.0) << " MB" << endl;
         return usage.ru_maxrss/(1024.0*1024.0);
#else
         if(print) cout << "Memory Usage: " << usage.ru_maxrss/(1024.0) << " MB" << endl;
         return usage.ru_maxrss/(1024.0);
#endif
    }
    /**
     * @brief A public method that return the current memory peak in gigabytes
     * @param print a boolean, if true prints on the standard output a string containing the value of the memory peak
     * @return a float containing the memory peak
     */
    inline float getValue_in_GB(bool print)
    {
         ret=getrusage(who,&usage);
#ifdef __APPLE__
         if(print) cout << "Memory Usage: " << usage.ru_maxrss/(1024.0*1024.0*1024.0) << " GB" << endl;
         return usage.ru_maxrss/(1024.0*1024.0*1024.0);
#else
         if(print) cout << "Memory Usage: " << usage.ru_maxrss/(1024.0*1024.0) << " GB" << endl;
         return usage.ru_maxrss/(1024.0*1024.0);
#endif
    }
};



#endif // USAGE_H
