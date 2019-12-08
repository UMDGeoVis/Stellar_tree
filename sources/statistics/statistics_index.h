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

#include <deque>

#include "statistics.h"

template<class C,class T> void Statistics::init_vector(Mesh<C,T>& mesh)
{
    for(int i=0; i<mesh.get_top_cells_types(); i++)
    {
        this->indexStats.numLeafForTop[i].assign(mesh.get_top_cells_num(i),0);

        if(i<mesh.get_top_cells_types()-1)
        {
            this->indexStats.numLeafForTop.push_back(ivect());

            this->indexStats.minLeafForTop.push_back(INT_MAX);
            this->indexStats.maxLeafForTop.push_back(0);

            this->indexStats.min_Top_PerLeaf.push_back(INT_MAX);
            this->indexStats.avg_Top_PerLeaf.push_back(0);
            this->indexStats.max_Top_PerLeaf.push_back(0);

            this->indexStats.min_Reindexed_TopPerLeaf.push_back(INT_MAX);
            this->indexStats.avg_Reindexed_TopPerLeaf.push_back(0);
            this->indexStats.max_Reindexed_TopPerLeaf.push_back(0);

            this->indexStats.avgWeightedLeafForTop.push_back(0);

            this->indexStats.t_list_length.push_back(0);
            this->indexStats.real_t_list_length.push_back(0);

            this->indexStats.min_run_length.push_back(INT_MAX);
            this->indexStats.max_run_length.push_back(0);
            this->indexStats.avg_run_length.push_back(0);
            this->indexStats.tot_number_of_run.push_back(0);
        }
    }
}

template<class C, class T> void Statistics::get_index_statistics(Stellar_Tree &tree, Mesh<C,T> &mesh, bool compute_run_histogram)
{    
    init_vector(mesh);
    visit_tree(tree.get_root(),0,mesh,compute_run_histogram); // the root is at level 0
    calc_remaining_index_stats();
    check_inconsistencies();
    Writer::write_index_stats(this->indexStats,mesh,compute_run_histogram);
    this->indexStats = IndexStatistics();
    return;
}

template<class C,class T> void Statistics::visit_tree(Node_Stellar &n, int n_level, Mesh<C,T> &mesh, bool compute_run_histogram)
{
    this->indexStats.numNode++;

    if(n.is_leaf())
    {
        if(this->indexStats.minTreeDepth==INT_MAX || this->indexStats.minTreeDepth > n_level)
            this->indexStats.minTreeDepth = n_level;
        if(this->indexStats.maxTreeDepth < n_level)
            this->indexStats.maxTreeDepth = n_level;
        this->indexStats.avgTreeDepth += n_level;

        this->get_leaf_stats(n,mesh,compute_run_histogram);
    }
    else
    {
        int sons_level = n_level + 1;

        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
                this->visit_tree(**it,sons_level,mesh,compute_run_histogram);
        }
    }
    return;
}

template<class C, class T> void Statistics::get_leaf_stats(Node_Stellar &n, Mesh<C, T> &mesh, bool compute_run_histogram)
{
    int num_t_completely;

    /// all the leaf nodes are full
    this->indexStats.numFullLeaf++;

    if(compute_run_histogram)
    {
        /// compute the statistics concerning the run efficiency
        n.compute_run_histogram(this->indexStats);
    }

    for(int i=0; i<n.get_num_top_cells_encoded();i++)
    {
        num_t_completely = 0;

        this->indexStats.t_list_length[i] += n.get_t_array_size(i);
        this->indexStats.real_t_list_length[i] += n.get_real_t_array_size(i);
        this->indexStats.avg_run_length[i] += n.get_t_in_run_size(i,this->indexStats.tot_number_of_run[i],this->indexStats.min_run_length[i],this->indexStats.max_run_length[i]);

        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& runIt = itPair.first;
            T& top = mesh.get_top_cell(i,*runIt);
            if((n.completely_indexes_top(top)) )
                num_t_completely++;
            this->indexStats.numLeafForTop[i][*runIt-1]++; ///MODIFICA 27 FEBBRAIO 2013 (README MAIL 8 GIUGNO 2011)
        }

        /// we count the number of leaf nodes that have an internal run of top k-simplexes
        if(num_t_completely > 0)
            this->indexStats.numLeaf_InternalRun++;

        if((num_t_completely) > 0)
        {
            int num_vertex = n.get_real_v_array_size();
            if(this->indexStats.minVertexInFullLeaf==INT_MAX || this->indexStats.minVertexInFullLeaf > num_vertex)
                this->indexStats.minVertexInFullLeaf = num_vertex;
            if(this->indexStats.maxVertexInFullLeaf < num_vertex)
                this->indexStats.maxVertexInFullLeaf = num_vertex;
            this->indexStats.avgVertexInFullLeaf += num_vertex;

            if(this->indexStats.min_Reindexed_TopPerLeaf[i]==INT_MAX || this->indexStats.min_Reindexed_TopPerLeaf[i] > num_t_completely)
                this->indexStats.min_Reindexed_TopPerLeaf[i] =  num_t_completely;
            if(this->indexStats.max_Reindexed_TopPerLeaf[i] < num_t_completely)
                this->indexStats.max_Reindexed_TopPerLeaf[i] = num_t_completely;
            this->indexStats.avg_Reindexed_TopPerLeaf[i] += num_t_completely;

            /// we consider the real statistic (i.e. no compression, no reindexing)
            if(this->indexStats.min_Top_PerLeaf[i]==INT_MAX || this->indexStats.min_Top_PerLeaf[i] > n.get_real_t_array_size(i))
                this->indexStats.min_Top_PerLeaf[i] = n.get_real_t_array_size(i);
            if(this->indexStats.max_Top_PerLeaf[i] < n.get_real_t_array_size(i))
                this->indexStats.max_Top_PerLeaf[i] = n.get_real_t_array_size(i);
            this->indexStats.avg_Top_PerLeaf[i] += n.get_real_t_array_size(i);
        }
    }
}

template<class C, class T> void Statistics::get_unique_cp_cells_counters(Stellar_Tree &tree, Mesh<C, T> &mesh)
{
    long long tot_simplices = 0;

    ///we exclude the vertices and the top simplices
    ///then we visit all the implicitly encoded simplexes
    cerr<<"unique number of d-simplexes"<<endl;

    cerr<<"0-simplices num: "<<mesh.get_vertices_num()<<endl;
    tot_simplices += mesh.get_vertices_num();

    cerr<<mesh.get_implicitly_encoded_cells_num()<<"-simplices num: "<<mesh.get_top_cells_num(mesh.get_top_cells_types()-1)<<endl;
    tot_simplices += mesh.get_top_cells_num(mesh.get_top_cells_types()-1);

    // first the simplest simplices that can be extracted, i.e., edges and d-1 simplices
    pair<int,long long> p;

//    p.first = 1; p.second = 0;
//    this->get_unique_cells_counter(tree,mesh,p);
//    tot_simplices += p.second;

//    if(!mesh.get_implicitly_encoded_cells_num()-1 != 1) // handling triangle meshes
//    p.first = mesh.get_implicitly_encoded_cells_num()-1; p.second = 0;
//    this->get_unique_cells_counter(tree,mesh,p);
//    tot_simplices += p.second;


    for(int d=1; d<mesh.get_implicitly_encoded_cells_num(); d++) // from triangles to d-2 simplices
    {        
        p.first = d; p.second = 0;
        this->get_unique_cells_counter(tree,mesh,p);
        tot_simplices += p.second;
    }

    cerr<<"[TOT] "<<tot_simplices<<endl;
}

template<class C, class T> void Statistics::get_unique_dcells_counter(Node_Stellar &n, Mesh<C,T> &mesh, int d, long long &num)
{
    deque<leaf_p_faces> unique_simplexes;
    unique_simplexes.assign(n.get_v_end()-n.get_v_start(),set<ivect>());

    ivect dsimplex;

    /// we cannot place d instead of 0 because we do not know here if we are encoding the tops using a verbose encoding or a compact one
    for(int i=0; i<n.get_num_top_cells_encoded();i++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& runIt = itPair.first;
            T& top = mesh.get_top_cell(i,*runIt);

            int num_s = top.get_sub_types_num(d);

            if(num_s == -1)
            {
                break; /// we stop this for as we are analyzing a top simplex that is lower dimensional than the one that we want to extract
            }

            for(int s=0; s<num_s; s++)
            {
                top.get_d_cell(dsimplex,d,s);
                if(n.indexes_vertex(dsimplex[0]))
                {
                    int v = dsimplex[0];
                    dsimplex.erase(dsimplex.begin());
                    unique_simplexes[v-n.get_v_start()].insert(dsimplex);
                }
                dsimplex.clear();
            }
        }
    }

    int local_counter = 0;
    for(unsigned i=0; i<unique_simplexes.size(); i++)
        local_counter += unique_simplexes[i].size();

#pragma omp critical
    {
        num += local_counter;
    }

    unique_simplexes.clear();
}

