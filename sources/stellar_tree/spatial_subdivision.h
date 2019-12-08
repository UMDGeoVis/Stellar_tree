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

#ifndef _SPATIALSUBDIVISION_H
#define	_SPATIALSUBDIVISION_H

#include <boost/dynamic_bitset.hpp>
#include "box.h"
#include "vertex.h"

///
/**
 * @brief A class representing a division criterion that is applied to the tree.
 *
 * It executes a n-quadtree subdivision up to 6D (included), while from 7D it scales to kD-tree subdivion.
 */
class Spatial_Subdivision
{
public:
    /**
     * @brief A constructor method
     *
     */
    Spatial_Subdivision() {}
    /**
     * @brief A constructor method
     *
     * @param dimension an integer representing Euclidean space dimension
     */
    Spatial_Subdivision(int dimension)
    {
        this->dimension = dimension;
         if(this->dimension < 7)
             this->son_number = pow(2,dimension); /// generalized quad tree subdivision
         else
             this->son_number = 2; /// kd-tree subdivision
    }
    /**
     * @brief A copy-constructor method
     *
     * @param orig
     */
    Spatial_Subdivision(const Spatial_Subdivision& orig) {this->son_number = orig.son_number; this->dimension = orig.dimension; }    
    /**
     * @brief A destructor method
     *
     */
    virtual ~Spatial_Subdivision() {}
    /**
     * @brief A public method that sets the Euclidean space dimension and the number of sons to be created at each split operation
     * Up to 6D it is maintained a nd-quadtree like subdivision, while from 7D the subdivision switch to kD-tree
     *
     * @param dimension
     */
    inline void set(int dimension)
    {
        this->dimension = dimension;
         if(this->dimension < 7)
             this->son_number = pow(2,dimension); /// generalized quad tree subdivision
         else
             this->son_number = 2; /// kd-tree subdivision
    }
    /**
     * @brief Public method that returns the number of son nodes
     *
     * @return an integer value representing the number of son nodes
     */
    inline int get_son_number() { return son_number; }    
    /**
     * @brief Public method that computes the box domain of a son node
     *
     * @param parent_dom a Box& argument, representing the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * @param child_ind a integer argument, representing the index position of the son in the sub-tree
     * @return Box<C> representing the domain of the son
     */
    template<class C> Box<C> compute_domain(Box<C>& parent_dom, int level, int child_ind);
    /**
     * @brief A public method that identifies the son node that indexes the vertex v and return the son domain and its position in the sons array
     *
     * @param parent_dom a Box& argument, representing the node domain
     * @param level an integer argument representing the level of n in the hierarchy
     * @param v a Vertex& argument, representing the target vertex
     * @param max_point a Point& argument, representing the upper corner of the mesh domain
     * @param child_ind an integer argument that saves the child position index in the sons array
     * @param sonDom a Box& argument that saves the domain of the son
     * @return bool true if the son is identified, false otherwise
     */
    template<class C> bool get_son_index_and_domain(Box<C>& parent_dom, int level, Vertex<C> &v, Point<C> &max_point, int &child_ind, Box<C> &sonDom);

private:
    int son_number; /**< A private variable encoding the number of son nodes that each internal node of the tree has */
    int dimension; /**< A private variable encoding the Euclidean space dimension that embeds the indexed mesh */

    /**
     * @brief A private method that computes the domain of the son at position child_ind following the n-quadtree subdivision
     *
     * @param parent_dom a Box& argument representing the domain of the parent node
     * @param child_ind an integer argument representing the position index of the child in the array of sons of the parent
     * @return Box representing the domain of the son
     */
    template<class C> Box<C> compute_domain_nquad(Box<C>& parent_dom, int child_ind);
    /**
     * @brief A private method that computes the domain of the son node at position child_ind following the kD-tree subdivision
     *
     * @param parent_dom a Box& argument representing the domain of the parent node
     * @param level an integer argument representing the level of n in the hierarchy
     * @param child_ind an integer argument representing the position index of the child in the array of sons of the parent
     * @return Box representing the domain of the son
     */
    template<class C> Box<C> compute_domain_kd(Box<C>& parent_dom, int level, int child_ind);
    /**
     * @brief A private method that compute the position index and the domain of son node starting from the coordinates of a vertex -- n-quadtree version
     *
     * @param parent_dom a Box& argument representing the domain of the parent node
     * @param v a Vertex& argument representing the target vertex used in the computations
     * @param max_point a Point& argument representing the upper corner of the mesh domain
     * @param child_ind an integer argument representing the position index of the child in the array of sons of the parent
     * @param sonDom representing the domain of the son
     * @return bool true if the computation ends successfully, false otherwise
     */
    template<class C> bool get_son_index_and_domain_nquad(Box<C>& parent_dom, Vertex<C> &v, Point<C> &max_point, int &child_ind, Box<C> &sonDom);
    /**
     * @brief A private method that compute the position index and the domain of son node starting from the coordinates of a vertex -- kD-tree version
     *
     * @param parent_dom a Box& argument representing the domain of the parent node
     * @param level an integer argument representing the level of n in the hierarchy
     * @param v a Vertex& argument representing the target vertex used in the computations
     * @param max_point a Point& argument representing the upper corner of the mesh domain
     * @param child_ind an integer argument representing the position index of the child in the array of sons of the parent
     * @param sonDom representing the domain of the son
     * @return bool true if the computation ends successfully, false otherwise
     */
    template<class C> bool get_son_index_and_domain_kd(Box<C>& parent_dom, int level, Vertex<C> &v, Point<C> &max_point, int &child_ind, Box<C> &sonDom);
};

template<class C> Box<C> Spatial_Subdivision::compute_domain(Box<C> &parent_dom, int level, int child_ind)
{
    /// we change the subdivision type depending on the dimension in which the points lies
    if(dimension < 7)
        return compute_domain_nquad(parent_dom,child_ind);
    else
        return compute_domain_kd(parent_dom,level,child_ind);
}

template<class C> Box<C> Spatial_Subdivision::compute_domain_nquad(Box<C> &parent_dom, int child_ind)
{
    boost::dynamic_bitset<> son_id_bits(dimension,child_ind);
    Box<C> sonDom = Box<C>(dimension);

    Point<C> &min = sonDom.getMinPoint();
    Point<C> &max = sonDom.getMaxPoint();

    for(int i=0; i<dimension; i++)
    {
        if(son_id_bits[i])
        {
            min.setC(i,parent_dom.getMinPoint().getC(i)+(parent_dom.getMaxPoint().getC(i)-parent_dom.getMinPoint().getC(i))/2.0);
            max.setC(i,parent_dom.getMaxPoint().getC(i));
        }
        else
        {
            min.setC(i,parent_dom.getMinPoint().getC(i));
            max.setC(i,parent_dom.getMinPoint().getC(i)+(parent_dom.getMaxPoint().getC(i)-parent_dom.getMinPoint().getC(i))/2.0);
        }
    }

    return sonDom;
}

template<class C> Box<C> Spatial_Subdivision::compute_domain_kd(Box<C>& parent_dom, int level, int child_ind)
{
    int coord_to_change = level % parent_dom.getMinPoint().get_dimension();

    Box<C> sonDom = parent_dom;

    if(child_ind == 1)
    {
        sonDom.getMinPoint().setC(coord_to_change,parent_dom.getMinPoint().getC(coord_to_change)+(parent_dom.getMaxPoint().getC(coord_to_change)-parent_dom.getMinPoint().getC(coord_to_change))/2.0);
    }
    else if(child_ind == 0)
    {
        sonDom.getMaxPoint().setC(coord_to_change,parent_dom.getMinPoint().getC(coord_to_change)+(parent_dom.getMaxPoint().getC(coord_to_change)-parent_dom.getMinPoint().getC(coord_to_change))/2.0);
    }

    return sonDom;
}

template<class C> bool Spatial_Subdivision::get_son_index_and_domain(Box<C>& parent_dom, int level, Vertex<C> &v, Point<C> &max_point, int &child_ind, Box<C> &sonDom)
{
    /// we change the subdivision type depending on the dimension in which the points lies
    if(dimension < 7)
        return get_son_index_and_domain_nquad(parent_dom,v,max_point,child_ind,sonDom);
    else
        return get_son_index_and_domain_kd(parent_dom,level,v,max_point,child_ind,sonDom);
}

template<class C> bool Spatial_Subdivision::get_son_index_and_domain_nquad(Box<C>& parent_dom, Vertex<C> &v, Point<C> &max_point, int &child_ind, Box<C> &sonDom)
{
    boost::dynamic_bitset<> son_id_bits(dimension,0);
    sonDom = Box<C>(dimension);

    Point<C> &min = sonDom.getMinPoint();
    Point<C> &max = sonDom.getMaxPoint();
    double new_min, new_max;

    for(int i=0; i<dimension; i++)
    {
        if(parent_dom.contains_coord_upper_half(v.getC(i),i,max_point.getC(i),new_min,new_max))
        {
            son_id_bits.set(i);
        }
        else if(!parent_dom.contains_coord_lower_half(v.getC(i),i,max_point.getC(i),new_min,new_max))
        {
            return false;
        }
        min.setC(i,new_min);
        max.setC(i,new_max);
    }

    child_ind = son_id_bits.to_ulong();
    return true;
}

template<class C> bool Spatial_Subdivision::get_son_index_and_domain_kd(Box<C>& parent_dom, int level, Vertex<C> &v, Point<C> &max_point, int &child_ind, Box<C> &sonDom)
{
    int coord_to_change = level % parent_dom.getMinPoint().get_dimension();

    sonDom = parent_dom;

    C new_min = 0, new_max = 0;

    if(parent_dom.contains_coord_upper_half(v.getC(coord_to_change),coord_to_change,max_point.getC(coord_to_change),new_min,new_max))
    {
       child_ind = 1;
    }
    else if(parent_dom.contains_coord_lower_half(v.getC(coord_to_change),coord_to_change,max_point.getC(coord_to_change),new_min,new_max))
    {
        child_ind = 0;
    }
    else
        return false;

    sonDom.getMinPoint().setC(coord_to_change,new_min);
    sonDom.getMaxPoint().setC(coord_to_change,new_max);

    return true;
}

#endif	/* _SPATIALSUBDIVISION_H */


