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

#ifndef _BOX_H
#define	_BOX_H

#include "point.h"

#include <iostream>
#include <cmath>

/**
 * @brief A template class representing a box in d-dimensional space.
 * The template parameter is CoordType, that define the precision of the lower and upper points of the box.
 */
template<class CoordType> class Box {
public:
    /**
     * @brief A constructor method
     *
     */
    Box();
    ///A copy-constructor method
    /**
     * @brief
     *
     * @param orig is a Box<CoordType> variable
     *
     */
    Box(const Box<CoordType>& orig);
    ///A constructor method
    /**
     * @brief
     *
     * @param min is a Point<CoordType> representing the lower corner of the box
     * @param max is a Point<CoordType> representing the upper corner of the box
     */
    Box(Point<CoordType> &min, Point<CoordType> &max);
    ///A constructor method
    /**
     * @brief
     *
     * @param dimension an integer that set the dimension of the two point forming the box, i.e, the Euclidean space in which they are embedded
     */
    Box(int dimension)
    {
        min = Point<CoordType>(dimension);
        max = Point<CoordType>(dimension);
    }
    ///A destructor method
    /**
     * @brief
     *
     */
    virtual ~Box() {}
    ///Public method that returns the minimum point of the box.    
    /**
     * @brief
     *
     * @return the reference of the lower corner of the box
     */
    inline Point<CoordType>& getMinPoint() { return this->min; }
    ///Public method that returns the maximum point of the box.
    /**
     * @brief
     *
     * @return the reference of the upper corner of the box
     */
    inline Point<CoordType>& getMaxPoint() { return this->max; }
    ///Public method that sets the lower corner of the box
    /**
     * @brief
     *
     * @param min a vector of CoordType containings the coordinates of the lower corner
     */
    inline void setMinPoint(vector<CoordType> &min)
    {
        this->min = min;
    }
    ///Public method that sets the lower corner of the box
    /**
     * @brief
     *
     * @param max a vector of CoordType containings the coordinates of the upper corner
     */
    inline void setMaxPoint(vector<CoordType> &max)
    {
        this->max = max;
    }
    ///Public method that computes the centroid of the box
    /**
     * @brief
     *
     * @param p a Point<CoordType> variable that at the end of processing contains the coordinates of the box centroid
     */
    inline void getCenter(Point<CoordType>& p)
    {        
        for(unsigned i=0; i<min.get_dimension(); i++)
            p.setC(i,((min.getC(i) + max.getC(i)) / 2.0));
    }
    ///Public method that checks if a box intersects the current box.
    /**
     * @brief
     *
     * @param other a Box& argument representing the box we need to test
     * @return a boolean value, true if the box intersects the current box, false otherwise
     */
    bool intersects(Box<CoordType>& other);
    ///Public method that checks if a point is inside the current box, viewed as domain.
    /*!
     * This method is not used to check if a point is inside the current box,
     * as it considers all the faces of the box as closed.
     * @param p a Point<CoordType>& argument, represents the point to check
     * @return a boolean value, true if the point is inside domain, false otherwise
     */
    bool contains_loose(Point<CoordType>& p);
    ///Public method that checks if a point is inside the current box.
    /*!
     * This method considers as closed only three faces, the ones incidentes in the minimum point.
     * It check if the box is at the border of the mesh domain,
     * and in these cases considers as close the corresponding face.
     * @param p a Point& argument, represents the point to check
     * @param max a Point& argument, represent the maximum point of the mesh domain
     * @return a boolean value, true if the point is inside the box, false otherwise
     */
    bool contains(Point<CoordType>& p, Point<CoordType>& max);
    /**
     * @brief Public method that checks if a variable is in the upper half range of a specific coordinate of the box
     *
     * @param c a CoordType variable representing the value to test
     * @param i an integer variable representing the coordinate position in the coordinate system
     * @param domain_max a CoordType variable containing the maximum value of the domain
     * @param new_min a CoordType variable representing the new lower bound of the coordinate
     * @param new_max a CoordType variable representing the new upper bound of the coordinate
     * @return a boolean value, true if c is inside the range, false otherwise
     */
    bool contains_coord_upper_half(CoordType c, int i, CoordType domain_max, CoordType &new_min, CoordType &new_max);
    /**
     * @brief Public method that checks if a variable is in the lower half range of a specific coordinate of the box
     *
     * @param c a CoordType variable representing the value to test
     * @param i an integer variable representing the coordinate position in the cooordinate system
     * @param domain_max a CoordType variable containing the maximum value of the domain
     * @param new_min a CoordType variable representing the new lower bound of the coordinate
     * @param new_max a CoordType variable representing the new upper bound of the coordinate
     * @return a boolean value, true if c is inside the range, false otherwise
     */
    bool contains_coord_lower_half(CoordType c, int i, CoordType domain_max, CoordType &new_min, CoordType &new_max);
    /**
     * @brief Public method that checks if a variable is in range
     *
     *
     * @param c a CoordType variable representing the value to test
     * @param domain_max a CoordType variable containing the maximum value of the domain
     * @param new_min a CoordType variable representing the lower bound of the range
     * @param new_max a CoordType variable representing the upper bound of the range
     * @return a boolean value, true if c is inside the range, false otherwise
     */
    bool coord_in_range(CoordType c, CoordType domain_max, CoordType new_min, CoordType new_max);
    /**
     * @brief Public method that checks if a point is inside the mesh domain, and if the box domain doesn't contain the point it enlarges it accordingly
     *
     * @param p a Point<CoordType>& argument, represents the point to check
     */
    void resize_to_contain(Point<CoordType> &p);
    /**
     * @brief Public method that computes within two boxes
     *
     * @param b a Box<CoordType>& variable, representing the other box
     * @return a CoordType containing the minimum distance between the two boxes
     */
    CoordType min_distance(Box<CoordType>& b);
    /**
     * @brief Public method that computes within the box and a point
     *
     * @param p a Point<CoordType>& variable, representing the point
     * @return a CoordType containing the minimum distance between the box and the point
     */
    CoordType min_distance(Point<CoordType> &p);
    /**
     * @brief
     *
     * @param out
     * @param p
     * @return std::ostream &operator
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Box& p)
    {
        out << p.min << " " << p.max;
        return out;
    }

    inline CoordType get_side(int i) { return fabs(min.getC(i)-max.getC(i)); }
    inline CoordType get_diagonal()
    {
        CoordType diagonal = 0.0;
        for(int i=0; i<min.get_dimension(); i++)
        {
            CoordType s = this->get_side(i);
            diagonal += s*s;
        }
        return sqrt(diagonal);
    }

protected:
    ///A protected variable representing the lower corner of the box
    Point<CoordType> min;
    ///A protected variable representing the upper corner of the box
    Point<CoordType> max;
};

template<class CoordType> Box<CoordType>::Box()
{
    min = Point<CoordType>();
    max = Point<CoordType>();
}

template<class CoordType> Box<CoordType>::Box(const Box<CoordType>& orig)
{
    min = orig.min;
    max = orig.max;
}

template<class CoordType> Box<CoordType>::Box(Point<CoordType>& min, Point<CoordType>& max){
    this->min = min;
    this->max = max;
}

template<class CoordType> bool Box<CoordType>::intersects(Box<CoordType> &other)
{
    if(this->min.get_dimension() != other.min.get_dimension())
    {
        cerr<<"[intersects] the two points have different dimensions"<<endl;
        return false;
    }

    for(int i=0; i<this->min.get_dimension(); i++)
    {
        if (this->max.getC(i) < other.min.getC(i)) return false;
        if (this->min.getC(i) > other.max.getC(i)) return false;
    }
    return true;
}

// WARNING: it is used to check if a point is contained by the box, but this procedure consider all the faces of the domain as closed.
// It is used in a preprocessing phase
template<class CoordType> bool Box<CoordType>::contains_loose(Point<CoordType> &p)
{    
    for(int i=0; i<this->min.get_dimension(); i++)
    {
        if (this->max.getC(i) < p.getC(i)) return false;
        if (this->min.getC(i) > p.getC(i)) return false;
    }
    return true;
}

// WARNING: it is used to check if a point is contained by the box, but this procedure consider as closed the faces that are incident in the lower corner of the domain.
// Moreover it consider the border casistics, and if one (or more) coordinates are on the border, the procedure consider the incident faces as closed.
template<class CoordType> bool Box<CoordType>::contains(Point<CoordType> &p, Point<CoordType> &max)
{
    for(int i=0; i<this->min.get_dimension(); i++)
    {
        if(!coord_in_range(p.getC(i),max.getC(i),this->min.getC(i),this->max.getC(i)))
            return false;
    }

    return true;
}

template<class CoordType> bool Box<CoordType>::contains_coord_upper_half(CoordType c, int i, CoordType domain_max, CoordType &new_min, CoordType &new_max)
{
    new_min = min.getC(i)+(max.getC(i)-min.getC(i))/2.0;
    new_max = max.getC(i);
    return coord_in_range(c,domain_max,new_min,new_max);
}

template<class CoordType> bool Box<CoordType>::contains_coord_lower_half(CoordType c, int i, CoordType domain_max, CoordType &new_min, CoordType &new_max)
{
    new_min = min.getC(i);
    new_max = min.getC(i)+(max.getC(i)-min.getC(i))/2.0;
    return coord_in_range(c,domain_max,new_min,new_max);
}

template<class CoordType> bool Box<CoordType>::coord_in_range(CoordType c, CoordType domain_max, CoordType new_min, CoordType new_max)
{
    if(new_max==domain_max){
        if(new_max < c){
            return false;
        }
    }
    else if(new_max <= c){
        return false;
    }

    if (new_min > c){
        return false;
    }
    return true;
}

//this method is used to resize the mesh domain
template<class CoordType> void Box<CoordType>::resize_to_contain(Point<CoordType>& p)
{
    // If we are big enough to contain the point, do nothing;
    if (contains_loose(p)) return;

    for(int i=0; i<this->min.get_dimension(); i++)
    {
        if (p.getC(i) < this->min.getC(i)) this->min.setC(i,p.getC(i));
        if (p.getC(i) > this->max.getC(i)) this->max.setC(i,p.getC(i));
    }
}

template<class CoordType> CoordType Box<CoordType>::min_distance(Box<CoordType> &b)
{
    //we assume the axis aligned.
    //something will be added to sqrDist only when the projections onto the current axis do not overlap.
    //This means that when the min distance is the distance between two edges, nothing will be added to
    //sqrDist for the axis parallel to it, and so it will in effect be the distance between two corners of a 2D projection of the situation.
    //similarly, if the min distance is between two faces, sqrDist is only added to for one of the axes.
    CoordType sqrDist = 0;

    for(int i=0; i<this->min.get_dimension(); i++)
    {
        CoordType dist;
        if(b.getMaxPoint().getC(i) < this->min.getC(i))
        {
            dist = b.getMaxPoint().getC(i) - this->min.getC(i);
            sqrDist += dist * dist;
        }
        else if(b.getMinPoint().getC(i) > this->max.getC(i))
        {
            dist = b.getMinPoint().getC(i) - this->max.getC(i);
            sqrDist += dist * dist;
        }
    }

    return sqrt(sqrDist);
}

template<class CoordType> CoordType Box<CoordType>::min_distance(Point<CoordType> &p)
{
    CoordType sqrDist = 0;

    for(int i=0; i<this->min.get_dimension(); i++)
    {
        CoordType dist;
        if(p.getC(i) < this->min.getC(i))
        {
            dist = p.getC(i) - this->min.getC(i);
            sqrDist += dist * dist;
        }
        else if(p.getC(i) > this->max.getC(i))
        {
            dist = p.getC(i) - this->max.getC(i);
            sqrDist += dist * dist;
        }
    }

    return sqrt(sqrDist);
}

#endif	/* _BOX_H */

