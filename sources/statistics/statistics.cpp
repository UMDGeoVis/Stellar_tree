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

#include "statistics.h"

void Statistics::calc_remaining_index_stats()
{
    /// versione generalizzata al tipo di simplessi indicizzato
    for(unsigned j=0;j<this->indexStats.numLeafForTop.size();j++) /// ciclo sui tipi di simplessi
    {
        int numTinMoreLeaf = 0;
        for(unsigned i=0;i<this->indexStats.numLeafForTop[j].size();i++)
        {
            if(this->indexStats.minLeafForTop[j] == INT_MAX || this->indexStats.minLeafForTop[j] > this->indexStats.numLeafForTop[j][i])
                this->indexStats.minLeafForTop[j] = this->indexStats.numLeafForTop[j][i];
            if(this->indexStats.maxLeafForTop[j] < this->indexStats.numLeafForTop[j][i])
                this->indexStats.maxLeafForTop[j] = this->indexStats.numLeafForTop[j][i];

            //effective chi value for partial
            if(this->indexStats.numLeafForTop[j][i]!=1)
                this->indexStats.avgWeightedLeafForTop[j] += this->indexStats.numLeafForTop[j][i];

            if(this->indexStats.numLeafForTop[j][i]>1)
                numTinMoreLeaf++;
        }

        if(this->indexStats.numFullLeaf > 0)
        {
            this->indexStats.avg_Top_PerLeaf[j] /= this->indexStats.numFullLeaf;
            if(this->indexStats.avg_Reindexed_TopPerLeaf[j] != 0)
                this->indexStats.avg_Reindexed_TopPerLeaf[j] /= this->indexStats.numFullLeaf;
        }

        if(this->indexStats.numLeafForTop[j].size() > 0)
        {
            this->indexStats.avgWeightedLeafForTop[j] /= numTinMoreLeaf;
        }

        this->indexStats.avg_run_length[j] /= this->indexStats.tot_number_of_run[j];
    }

    if(this->indexStats.numNode > 0)
    {
        this->indexStats.avgTreeDepth /= this->indexStats.numFullLeaf;
    }

    if(this->indexStats.numFullLeaf > 0)
    {
        this->indexStats.avgVertexInFullLeaf /= this->indexStats.numFullLeaf;
    }
}

void Statistics::check_inconsistencies()
{
    if(this->indexStats.numFullLeaf == 0){
        this->indexStats.minVertexInFullLeaf = 0;
    }
    if(this->indexStats.numNode == 0)
        this->indexStats.minTreeDepth = 0;

    if(this->indexStats.minVertexInFullLeaf==INT_MAX)
        this->indexStats.minVertexInFullLeaf=0;
    return;
}

void Statistics::get_unique_cells_counter(Stellar_Tree &tree, CP_Mesh &mesh, pair<int, long long> &p)
{
    cerr<<p.first<<"-simplices num: ";
    tree.visit(Statistics::get_unique_dcells_counter_CP_wrapper,tree.get_root(),mesh,p);
    cerr<<p.second<<endl;
}

void Statistics::get_unique_cells_counter(Stellar_Tree &tree, Simplicial_Mesh &mesh, pair<int,long long> &p)
{
    cerr<<p.first<<"-simplices num: ";
    tree.visit(Statistics::get_unique_dcells_counter_Simplicial_wrapper,tree.get_root(),mesh,p);
    cerr<<p.second<<endl;
}
