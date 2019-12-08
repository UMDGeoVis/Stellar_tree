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

#ifndef _VERTEX_H
#define	_VERTEX_H

#include "point.h"
#include <iostream>
#include <cstdio>

/**
 * @brief An inner-class, extending Point, representing a vertex in a mesh
 *
 */
template<class CoordType> class Vertex : public Point<CoordType>
{
public:    
    /**
     * @brief A constructor
     */
    Vertex(): Point<CoordType>() { /*this->field_value = -1;*/ }
    /**
     * @brief A copy-constructor
     * @param orig
     */
    Vertex(const Vertex& orig): Point<CoordType>(orig) { this->field_values = orig.field_values; }
    /**
     * @brief A constructor
     * @param coords, a vector containing the coordinates of the vertex
     *
//     * NOTA: the field value of the vertex is set to the default value -1
     */
    Vertex(vector<CoordType> &coords): Point<CoordType>(coords) { /*this->field_value = -1;*/ }
    /**
     * @brief A constructor
     * @param coords, a vector containing the coordinates of the vertex
     * @param field, the field value associated to the vertex
     */
    Vertex(vector<CoordType> &coords, CoordType field): Point<CoordType>(coords) { this->field_values = { field }; }
    /**
     * @brief A constructor
     * @param coords, a vector containing the coordinates of the vertex
     * @param fields, a vector containing the field values of the vertex
     */
    Vertex(vector<CoordType> &coords, vector<CoordType> fields): Point<CoordType>(coords) { this->field_values = fields; }
    /**
     * @brief A constructor
     * @param dimension, the dimension of the point
     */
    Vertex(int dimension) : Point<CoordType>(dimension) { /*this->field_values = -1;*/ }
    /**
     * @brief A destructor
     */
    virtual ~Vertex() {}

    /**
     * @brief
     *
     * @param p
     * @param q
     * @return bool operator
     */
    inline friend bool operator== (const Vertex& p, const Vertex &q)
    {
        if(p.get_dimension() != q.get_dimension())
            return false;
        for(int i=0; i<=p.get_dimension(); i++)
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
    inline friend bool operator!= (const Vertex& p, const Vertex &q) { return !(p == q); }
    /**
     * @brief
     *
     * @param out
     * @param p
     * @return std::ostream &operator
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Vertex& p)
    {
        out << "v_dim["<<p.get_dimension()<<"] c[";
        for(int i=0; i<p.get_dimension(); i++)
        {
            out << p.getC(i) << " ";
        }
        out << p.field_values << "]";
        return out;
    }
    /**
     * @brief A public method that saves the vertex coordinates in a FILE* variable
     *
     * @param output
     */
    inline void save_to_file_coords(FILE* output)
    {
        for(int i=0; i<this->get_dimension(); i++)
        {
            fprintf(output,"%f ",this->getC(i));
        }
    }
    /**
     * @brief A public method that saves the vertex coordinates and field values in a FILE* variable
     *
     * @param output
     */
    inline void save_to_file(FILE* output)
    {
        for(int i=0; i<this->get_dimension(); i++)
        {
            fprintf(output,"%f ",this->getC(i));
        }
        for(int i=0; i<this->get_fields_num(); i++)
        {
            fprintf(output,"%f ",this->getF(i));
        }
    }
    /**
     * @brief A public method that returns the field value
     *
     * @return CoordType
     */
    inline CoordType setF(int pos, CoordType f) { return this->field_values[pos] = f; }
    /**
     * @brief A public method that returns the field value
     *
     * @return CoordType
     */
    inline CoordType getF(int pos) { return this->field_values[pos]; }
    /**
     * @brief A public method that returns the field value
     *
     * @return CoordType
     */
    inline CoordType get_fields_num() { return this->field_values.size(); }
//    /**
//     * @brief A public method that returns one of the coordinates or the field value of the vertex
//     *
//     * @param pos, the position of the coordinate to return
//     * @return CoordType
//     *
//     * NOTA: if pos is equal to the dimension of the point then, the procedure returns the field value
//     */
//    inline CoordType getC(int pos) const
//    {
//        if(pos < 0 && pos > this->get_dimension()+this->get_fields_num())
//        {
//            std::cerr<<"[getC] Invalid coordinate array position!";
//            int a; cin>>a;
//            return -1;
//        }
//        else if(pos == this->get_dimension())
//            return this->field_values;
//        else
//            return this->coords[pos];
//    }
//    /**
//     * @brief A public method that set the value of a coordinate or the field value.
//     *
//     * @param pos, the position to set
//     * @param c, the new value of the coordinate / field value
//     */
//    inline void setC(int pos, CoordType c)
//    {
//        if(pos < 0 && pos > this->get_dimension()+1)
//            std::cerr<<"Invalid coordinate array position!";
//        else if(pos == this->get_dimension())
//            this->field_values = c;
//        else
//            this->coords[pos] = c;
//    }
    /**
     * @brief A public method that compute the distance between two vertices
     * considering also the field values
     *
     * @param v, the other vertex
     * @return CoordType, the distance value
     */
    inline CoordType distance_withF(Vertex& v)
    {
        if(this->get_dimension() != v.get_dimension())
        {
            cerr<<"[distance] two points have different dimensions"<<endl;
            return -1;
        }
        CoordType distance_factor = 0;
        for(unsigned i=0; i<v.get_dimension(); i++)
        {
            distance_factor += (this->coords[i]-v.getC(i))*(this->coords[i]-v.getC(i));
        }
        for(unsigned i=0; i<v.get_fields_num(); i++)
            distance_factor += (this->getF(i)-v.getF(i))*(this->getF(i)-v.getF(i));
        return sqrt(distance_factor);
    }
//    /**
//     * @brief  A public method that returns the norm of v-this considering in the computation also the field value
//     *
//     * NOTA: return undefined behavior if the two points have different dimensions
//     *
//     * @param v, a Vertex&
//     * @return CoordType, the value of the norm
//     */
//    inline CoordType norm_withF(Vertex& v)
//    {
//        if(this->get_dimension() != v.get_dimension())
//        {
//            cerr<<"[norm] the two points have different dimensions"<<endl;
//            return -1;
//        }
//        CoordType factor = 0;
//        for(int i=0; i<this->get_dimension(); i++)
//        {
//            factor += (v.getC(i)-this->coords[i])*(v.getC(i)-this->coords[i]);
//        }
//        factor += (v.get_field_value()-this->field_values)*(v.get_field_value()-this->field_values);
//        return sqrt(factor);
//    }
//    /**
//     * @brief A public method that returns the norme of this vertex considering in the computation also the field value
//     *
//     * @return CoordType, the norm value
//     */
//    inline CoordType norm_withF()
//    {
//        CoordType factor = 0;
//        for(unsigned i=0; i<this->get_dimension(); i++)
//        {
//            factor += this->coords[i]*this->coords[i];
//        }
//        factor += field_values*field_values;
//        return sqrt(factor);
//    }
//    /**
//     * @brief A public method that computes the scalar product between two vectors v1-this and v2-this
//     *  considering in the computation also the field value
//     *
//     * @param v1, a Point&
//     * @param v2, a Point&
//     * @return CoordType, the scalar product value
//     */
//    inline CoordType scalar_product_withF(Vertex& v1,Vertex& v2) { return(this->scalar_product(v1,v2)+((v1.get_field_value()-this->field_values)*(v2.get_field_value()-this->field_values))); }
//    /**
//     * @brief A public method that computes the scalar product of vector v1-this
//     *  considering in the computation also the field value
//     *
//     * @param v1, a Point&
//     * @return CoordType, the scalar product value
//     */
//    inline CoordType scalar_product_withF(Vertex& v1) {  return(this->scalar_product(v1) + ((v1.get_field_value())*(this->field_values))); }

private:
    ///
    vector<CoordType> field_values; /**< A private variable encoding the field value of a vertex */
};

#endif	/* _VERTEX_H */

