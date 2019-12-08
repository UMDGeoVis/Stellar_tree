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


#ifndef _INDEXSTATISTICS_H
#define	_INDEXSTATISTICS_H

#include <vector>
#include <set>
#include <map>
#include <climits>
#include <iostream>

#include "utilities/basic_wrappers.h"

using namespace std;

/**
 * @brief A container class used to store the statistics of a spatial index
 *
 */
class IndexStatistics{
public:    
    /**
     * @brief A constructor method
     *
     */
    IndexStatistics()
    {
        numNode=0;
        numFullLeaf=0;
        numLeaf_InternalRun = 0;

        minTreeDepth=INT_MAX;
        avgTreeDepth=maxTreeDepth=0;
        minVertexInFullLeaf=INT_MAX;
        avgVertexInFullLeaf=maxVertexInFullLeaf=0;

        min_Top_PerLeaf.push_back(INT_MAX);
        avg_Top_PerLeaf.push_back(0);
        max_Top_PerLeaf.push_back(0);

        minLeafForTop.push_back(INT_MAX);
        maxLeafForTop.push_back(0);

        min_Reindexed_TopPerLeaf.push_back(INT_MAX);
        avg_Reindexed_TopPerLeaf.push_back(0);
        max_Reindexed_TopPerLeaf.push_back(0);

        avgWeightedLeafForTop.push_back(0);

        /// default only one entry..
        /// for variable entry number work outside the constructor
        t_list_length.push_back(0);
        real_t_list_length.push_back(0);

        numLeafForTop.push_back(ivect());

        min_run_length.push_back(INT_MAX);
        max_run_length.push_back(0);
        avg_run_length.push_back(0);
        tot_number_of_run.push_back(0);

        run_histogram.assign(13,0); ///
        run_strings.push_back("no_run");
        run_strings.push_back("3");
        run_strings.push_back("4");
        run_strings.push_back("5");
        run_strings.push_back("6->10");
        run_strings.push_back("11->15");
        run_strings.push_back("16->20");
        run_strings.push_back("21->30");
        run_strings.push_back("31->40");
        run_strings.push_back("41->50");
        run_strings.push_back("51->75");
        run_strings.push_back("76->100");
        run_strings.push_back(">100");
    }

    /**
     * @brief A public method that increment the counter associated to a run range
     *
     * @param run_size a integer representing a run
     */
    inline void set_run_histogram(int run_size)
    {
        if(run_size <= 2)
            run_histogram[0]++;
        else if(run_size == 3)
            run_histogram[1]++;
        else if(run_size == 4)
            run_histogram[2]++;
        else if(run_size == 5)
            run_histogram[3]++;
        else if(run_size > 5 && run_size <= 10)
            run_histogram[4]++;
        else if(run_size > 10 && run_size <= 15)
            run_histogram[5]++;
        else if(run_size > 15 && run_size <= 20)
            run_histogram[6]++;
        else if(run_size > 20 && run_size <= 30)
            run_histogram[7]++;
        else if(run_size > 30 && run_size <= 40)
            run_histogram[8]++;
        else if(run_size > 40 && run_size <= 50)
            run_histogram[9]++;
        else if(run_size > 50 && run_size <= 75)
            run_histogram[10]++;
        else if(run_size > 75 && run_size <= 100)
            run_histogram[11]++;
        else if(run_size > 100)
            run_histogram[12]++;
    }

    /**
     * @brief A public method that prints the histograms of run ranges
     *
     */
    inline void print_run_histogram()
    {
        cerr<<"run histogram"<<endl;
        for(unsigned i=0; i<run_histogram.size(); i++)
            cerr<<"["<<run_strings[i]<<"] "<<run_histogram[i]<<endl;
    }

    ///A public variable representing the number of tree nodes
    int numNode;
    ///A public variable representing the number of leaf containing significant information
    int numFullLeaf;

    /// a variable that keep track of the number of leaf that have a run of internal top k-simplexes
    /// (i.e. those indexed in a single leaf node)
    int numLeaf_InternalRun;

    ///A public variable representing the minimum tree depth
    int minTreeDepth;
    ///A public variable representing the maximum tree depth
    int maxTreeDepth;
    ///A public variable representing the average tree depth
    double avgTreeDepth;
    ///A public variable representing the minimum number of vertex per full leaf node
    int minVertexInFullLeaf;
    ///A public variable representing the maximum number of vertex per full leaf node
    int maxVertexInFullLeaf;
    ///A public variable representing the average number of vertex per full leaf node
    double avgVertexInFullLeaf;

    ///A public variable representing the minimum number of top-cells per full leaf node
    ivect min_Top_PerLeaf;
    ///A public variable representing the maximum number of top-cells per full leaf node
    ivect max_Top_PerLeaf;
    ///A public variable representing the average number of top-cells per full leaf node
    dvect avg_Top_PerLeaf;

    ///A public list representing, for each tetrahedron, the number of leaf that contain the specified tetrahedron
    vector<ivect > numLeafForTop;
    ///A public variable representing the minimum number of leaf that contains a single tetrahedron
    ivect minLeafForTop;
    ///A public variable representing the maximum number of leaf that contains a single tetrahedron
    ivect maxLeafForTop;

    /// A public associative array that contains for each leaf the number of leaves (without duplicates) which share at least a tetrahedron
    map<int,set<pair<int,int> > > adjacent_leaf_histogram;

    ///A public variable representing the minimum number of top cells that are indexed in a leaf block
    ivect min_Reindexed_TopPerLeaf;
    ///A public variable representing the maximum number of top cells that are indexed in a leaf block
    ivect max_Reindexed_TopPerLeaf;
    ///A public variable representing the average number of top cells that are in indexed in a leaf block
    dvect avg_Reindexed_TopPerLeaf;
    ///A public variable representing the average number of top cells that are in a run but are not completely indexed by a single leaf
    dvect avgWeightedLeafForTop;

    /// A public variable that contains the summation of the top cells lists per top cells type
    ivect t_list_length;
    /// A public variable that contains the summation of the top cells lists (real size obtained expending the runs) per top cells type
    ivect real_t_list_length;

    /// A public variable containing the minimum run size per top cells type
    ivect min_run_length;
    /// A public variable containing the maximum run size per top cells type
    ivect max_run_length;
    /// A public variable containing the average run size per top cells type
    dvect avg_run_length;
    /// A public variable containing the total number of runs per top cells type
    ivect tot_number_of_run;

    /// A public variable containing the run histogram describing the distribution of the run in a given tree
    ivect run_histogram;
    /// A public variable containing the run histogram names (for outputting the results only)
    vector<string> run_strings;
};

#endif	/* _INDEXSTATISTICS_H */
