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

#include "stellar_tree.h"

template<typename F, class M, class P> void Stellar_Tree::visit(F function_name, Node_Stellar &n, M &mesh, P &parameters)
{
    if (n.is_leaf())
    {
        function_name(n,mesh,parameters);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->visit(function_name,**it,mesh,parameters);
        }
    }
}

template<typename F, class P> void Stellar_Tree::visit(F function_name, Node_Stellar &n, P &parameters)
{
    if (n.is_leaf())
    {
        function_name(n,parameters);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->visit(function_name,**it,parameters);
        }
    }
}

template<typename F, class M, class P> void Stellar_Tree::visit(F function_name, Node_Stellar &n, Box<COORDBASETYPE> &dom, int level, M &mesh, P &parameters)
{
    if (n.is_leaf())
    {
        function_name(n,dom,mesh,parameters);
    }
    else
    {
        Spatial_Subdivision &subdiv = this->get_subdivision();
        int son_level = level + 1;
        for (int i = 0; i < subdiv.get_son_number(); i++)
        {
            if(n.get_son(i) != NULL)
            {
                Box<COORDBASETYPE> son_dom = subdiv.compute_domain(dom,level,i);

                this->visit(function_name,*n.get_son(i),son_dom,son_level,mesh,parameters);
            }
        }
    }
}

template<typename F, class M, class P> void Stellar_Tree::parallel_visit(F function_name, M &mesh, P &parameters)
{
#pragma omp parallel for
    for(unsigned i=0; i<this->get_leaves_number(); i++)
    {
        Node_Stellar* leaf = this->get_leaf(i);

        function_name(*leaf,mesh,parameters);
    }
}

template<typename F, class P> void Stellar_Tree::parallel_visit(F function_name, P &parameters)
{
#pragma omp parallel for
    for(unsigned i=0; i<this->get_leaves_number(); i++)
    {
        Node_Stellar* leaf = this->get_leaf(i);

        function_name(*leaf,parameters);
    }
}

template<typename F, class M, class C, class P> void Stellar_Tree::visit_with_cache(F function_name, Node_Stellar &n, M &mesh, C &cache, P &parameters)
{
    if (n.is_leaf())
    {
        function_name(n,mesh,cache,parameters);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->visit_with_cache(function_name,**it,mesh,cache,parameters);
        }
    }
}

template<typename F, class C, class P> void Stellar_Tree::visit_with_cache(F function_name, Node_Stellar &n, C &cache, P &parameters)
{
    if (n.is_leaf())
    {
        function_name(n,cache,parameters);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->visit_with_cache(function_name,**it,cache,parameters);
        }
    }
}

template<typename F, class M, class C, class P> void Stellar_Tree::parallel_visit_with_cache(F function_name, M &mesh, C &cache, P &parameters)
{
#pragma omp parallel for
    for(unsigned i=0; i<this->get_leaves_number(); i++)
    {
        Node_Stellar* leaf = this->get_leaf(i);

        function_name(*leaf,mesh,cache,parameters);
    }
}

template<typename F, class C, class P> void Stellar_Tree::parallel_visit_with_cache(F function_name, C &cache, P &parameters)
{
#pragma omp parallel for
    for(unsigned i=0; i<this->get_leaves_number(); i++)
    {
        Node_Stellar* leaf = this->get_leaf(i);

        function_name(*leaf,cache,parameters);
    }
}
