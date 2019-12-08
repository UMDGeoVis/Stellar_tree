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

#ifndef _POINT_H
#define	_POINT_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

/**
 * @brief A class representing a generic point in the d-dimensional space
 *
 */
template<class CoordType> class Point {
public:    
    /**
     * @brief A constructor method
     *
     */
    Point() { }    
    /**
     * @brief A copy-constructor method
     *
     * @param orig
     */
    Point(const Point& orig) { this->coords.assign(orig.coords.begin(),orig.coords.end()); }    
    /**
     * @brief A constructor method
     *
     * @param coords, an array containing the coordinates of the point
     */
    Point(vector<CoordType> &coords) { this->coords.assign(coords.begin(),coords.end()); }
    /**
     * @brief A constructor method
     *
     * @param dimension, the dimension of the point
     *
     * NOTA: the point coordinate are set al to 0.0
     */
    Point(int dimension) { this->coords.assign(dimension,0); }
    /**
     * @brief A public method that returns the dimension of the point
     *
     * @return int
     */
    inline int get_dimension() const { return this->coords.size(); }    
    /**
     * @brief A destructor method
     *
     */
    virtual ~Point() {}
    /**
     * @brief
     *
     * @param p
     * @param q
     * @return bool operator
     */
    inline friend bool operator== (const Point& p, const Point &q)
    {
        if(p.get_dimension() != q.get_dimension())
            return false;
        for(int i=0; i<p.get_dimension(); i++)
        {
            if(p.getC(i)!=q.getC(i))
                return false;
        }
        return true;
    }
    /**
     * @brief
     *
     * @param p
     * @param q
     * @return bool operator
     */
    inline friend bool operator!= (const Point& p, const Point &q) { return !(p == q); }
    ///
    /**
     * @brief
     *
     * @param out
     * @param p
     * @return std::ostream &operator
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Point& p)
    {
        out <<"[";
        for(int i=0; i<p.get_dimension(); i++)
        {
            out << p.getC(i) << " ";
        }
        out << "]";
        return out;
    }

    /**
     * @brief A public method that returns the coordinate at the 'pos' position
     *
     * @param pos, the position of the coordinate in the array
     * @return CoordType, the coordinate value at position pos
     */
    inline CoordType getC(int pos) const
    {
        // WARNING: no range checking
        return this->coords[pos];
    }
    ///
    /**
     * @brief A public method that set a coordinate
     *
     * @param pos, the position of the coordinate in the array
     * @param c, the value of the coordinate
     */
    inline void setC(int pos, CoordType c)
    {
        // WARNING: no range checking
        this->coords[pos] = c;
    }
    /**
     * @brief A public method that returns the norm of v - this
     * NOTA: return undefined behavior if the two points have different dimensions
     *
     * @param v, a Point&
     * @return CoordType, the value of the norm
     */
    inline CoordType norm(Point& v)
    {
        if(this->get_dimension() != v.get_dimension())
        {
            cerr<<"[norm] the two points have different dimensions"<<endl;
            return -1;
        }
        CoordType factor = 0;
        for(int i=0; i<this->get_dimension(); i++)
        {
            factor += (v.getC(i)-coords[i])*(v.getC(i)-coords[i]);
        }
        return sqrt(factor);
    }    
    /**
     * @brief A public method that returns the norm of this point
     *
     * @return CoordType, the norm value
     */
    inline CoordType norm()
    {
        CoordType factor = 0;
        for(unsigned i=0; i<this->get_dimension(); i++)
        {
            factor += coords[i]*coords[i];
        }
        return sqrt(factor);
    }
    /**
     * @brief A public method that computes the scalar product between two vectors v1-this and v2-this
     *
     * @param v1, a Point&
     * @param v2, a Point&
     * @return CoordType, the scalar product value
     */
    inline CoordType scalar_product(Point& v1,Point& v2)
    {
        if(v1.get_dimension() != v2.get_dimension() || this->get_dimension() != v1.get_dimension())
        {
            cerr<<"[scalar_product] two points have different dimensions"<<endl;
            return -1;
        }
        CoordType factor = 0;
        for(int i=0; i<v1.get_dimension(); i++)
        {
            factor += (v1.getC(i)-coords[i])*(v2.getC(i)-coords[i]);
        }
        return factor;
    }    
    /**
     * @brief A public method that computes the scalar product of vector v1-this
     *
     * @param v1, a Point&
     * @return CoordType, the scalar product value
     */
    inline CoordType scalar_product(Point& v1)
    {
        if(this->get_dimension() != v1.get_dimension())
        {
            cerr<<"[scalar_product] two points have different dimensions"<<endl;
            return -1;
        }
        CoordType factor = 0;
        for(unsigned i=0; i<v1.get_dimension(); i++)
        {
            factor += v1.getC(i)*coords[i];
        }
        return factor;
    }

    /**
     * @brief A public method that computes the distance between two points
     *
     * @param v, a Point&
     * @return CoordType, the distance value
     */
    inline CoordType distance(Point& v)
    {
        if(this->get_dimension() != v.get_dimension())
        {
            cerr<<"[distance] two points have different dimensions"<<endl;
            return -1;
        }
        CoordType distance_factor = 0;
        for(int i=0; i<v.get_dimension(); i++)
        {
            distance_factor += (coords[i]-v.getC(i))*(coords[i]-v.getC(i));
        }
        return sqrt(distance_factor);
    }

    /**
     * @brief A public method that normalize the coordinates of the current vertex respect to another vertex (usually the domain maximum)
     *
     * @param max, a Point&
     */
    inline void normalize(Point& max)
    {
        if(this->get_dimension() != max.get_dimension())
        {
            cerr<<"[normalize] two points have different dimensions"<<endl;
            return;
        }

        for(unsigned i=0; i<max.get_dimension(); i++)
        {
            coords[i] = coords[i]/max.getC(i);
        }
    }

protected:
    vector<CoordType> coords; /**< the array of point coordinates */
};

#endif	/* _POINT_H */

