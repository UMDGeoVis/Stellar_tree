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

#ifndef _NODE_H
#define	_NODE_H

#include <vector>
using namespace std;

/**
 * @brief A super-class, not instantiable, that represents a generic node of the tree
 *
 */
template<class N> class Node
{
public:
    /**
     * @brief A typename alias for accessing the iterator to the children nodes vector
     *
     */
    typedef typename std::vector<N*>::iterator child_iterator;
    /**
     * @brief A destructor method
     *
     */
    virtual ~Node() {}    
    /**
     * @brief A public method that sets a node son
     *
     * @param n a N* argument, represents the son to set
     * @param pos an integer, represents the son position
     */
    inline void set_son(N* n, int pos) { (*this->sons)[pos] = n; }
    /**
     * @brief A public method that checks if the node is a leaf node
     *
     * @return a boolean, true if the node is a leaf, false otherwise
     */
    inline bool is_leaf() const { return (this->sons==NULL); }
    /**
     * @brief A public method that returns the son at the i-th position
     *
     * @param i an integer argument, represents the son position into the list
     * @return N* representing the son at the i-th position
     */
    inline N* get_son(int i) { return (*this->sons)[i]; }    
    /**
     * @brief A public method that initialize the sons array
     *
     * @param num_sons an integer argument representing the number of sons
     */
    inline void init_sons(int num_sons) { this->sons = new vector<N*>(); this->sons->assign(num_sons,NULL); }    
    /**
     * @brief A public method that clear the son array
     *
     */
    inline void delete_sons() { sons->clear(); delete sons; sons=NULL; }

    /**
     * @brief A public method that returns the begin iterator to the son array
     *
     * @return child_iterator an iterator to the first position of the array
     */
    inline child_iterator begin() { return sons->begin(); }
    /**
     * @brief A public method that returns the end iterator to the son array
     *
     * @return child_iterator an iterator to the last position of the array
     */
    inline child_iterator end() { return sons->end(); }

protected:
    /**
     * @brief A constructor method
     *
     */
    Node()
    {
        this->sons = NULL;
    }
    /**
     * @brief A copy-constructor method
     *
     * @param orig
     */
    Node(const Node& orig)
    {
        this->sons = orig.sons;
    }

    ///A protected variable representing the list of node sons
    vector<N*> *sons;

private:

};

#endif	/* _NODE_H */

