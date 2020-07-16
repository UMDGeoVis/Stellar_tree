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

#ifndef _WRITER_H
#define	_WRITER_H

#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <queue>
#include <iostream>

#include "statistics/index_statistics.h"
#include "stellar_tree/box.h"
#include "stellar_tree/mesh.h"
#include "stellar_decomposition/run_iterator.h"
#include "stellar_decomposition/node_stellar.h"

using namespace std;

///A class that writes to file or standard output some library structure or statistics executed over the spatial index
/**
 * @brief
 *
 */
class Writer {
public:
    /**
     * @brief A public method that writes to file the tree structure
     *
     * @param fileName a string argument, represents the output path where the tree structure will be written
     * @param root a Node_Stellar& argument, represents the root of the tree to visit
     * @param original_vertices_ordering
     */
    static void write_tree(string fileName, Node_Stellar& root, ivect & original_vertices_ordering);
    /**
     * @brief A public method that writes to standard output the spatial index statistics
     *
     * @param indexStats an IndexStatistics& argument, represents the statistics to output
     * @param mesh representing the mesh indexed by the tree
     */
    template<class C, class T> static void write_index_stats(IndexStatistics& indexStats, Mesh<C,T> &mesh, bool verbose);
    /**
     * @brief A public method that writes to file a simplicial mesh in OFF format
     *
     * @param mesh_name a string argument represents the output path of the mesh
     * @param operation_type a string containing the name of the application from which we have obtained the new mesh file
     * @param mesh representing the simplicial mesh to save
     */
    template<class C, class T> static void write_mesh(string mesh_name, string operation_type, Mesh<C,T> &mesh); /// OFF format
    /**
     * @brief A public method that writes to file a simplicial mesh in VTK format
     *
     * @param mesh_name a string argument represents the output path of the mesh
     * @param operation_type a string containing the name of the application from which we have obtained the new mesh file
     * @param vertices_per_leaf an integer representing the leaf block threshold used during the tree generation phase
     * @param mesh representing the simplicial mesh to save
     */
    static void write_mesh_VTK(string mesh_name, Simplicial_Mesh& mesh);/// VTK FORMAT FOR (MANIFOLD and NON-MANIFOLD) TRIANGLE AND TETRAHEDRAL MESHES
//    static void write_2D_3D_mesh_VTK(string mesh_name, Simplicial_Mesh& mesh, dvect &field); /// VTK FORMAT FOR TRIANGLE AND TETRAHEDRAL MESHES
    static void write_mesh_VTK(string mesh_name, string operation_type, int vertices_per_leaf, Simplicial_Mesh& mesh);/// VTK FORMAT FOR (MANIFOLD and NON-MANIFOLD) TRIANGLE AND TETRAHEDRAL MESHES


private:
    /**
     * @brief A constructor method
     *
     */
    Writer() {}
    /**
     * @brief A constructor method
     *
     * @param
     */
    Writer(const Writer&) {}
    /**
     * @brief A destructor method
     *
     */
    virtual ~Writer() {}
    /**
     * @brief A private method that writes to output stream a tree node information (Generic Version)
     *
     * @param output an ofstream& argument, represents the stream where the tree node structure will be saved
     * @param n a Node_Stellar* argument, represents the node to save
     * @param original_vertices_ordering an integer vector containing the original vertices position indexes (prior the exploiting of spatial coherence)
     */
    static void write_node(ofstream& output, Node_Stellar* n, ivect & original_vertices_ordering);

};

template<class C, class T> void Writer::write_mesh(string mesh_name, string operation_type, Mesh<C,T> &mesh)
{
    stringstream stream;
    stream<<mesh_name<<"_"<<operation_type<<".off";
    FILE *output;
    output = fopen(stream.str().c_str(),"w");

    fprintf(output,"OFF\n");
    fprintf(output,"%d %d 0\n",mesh.get_vertices_num(),mesh.count_top_cells_num());
    for(int v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex<C>& vert = mesh.get_vertex(v);
        vert.save_to_file_coords(output);
        fprintf(output,"\n");

    }
    for(int i=0; i<mesh.get_top_cells_types(); i++)
    {
        for(int t=1; t<=mesh.get_top_cells_num(i); t++)
        {
            T &top = mesh.get_top_cell(i,t);
            fprintf(output,"%d ",top.get_vertices_num());
            top.save_to_file(output);
            fprintf(output,"\n");
        }
    }
    fclose(output);
}

template<class C, class T> void Writer::write_index_stats(IndexStatistics& indexStats, Mesh<C,T> &mesh, bool verbose)
{
    cerr << indexStats.numNode << " ";
    cerr << indexStats.numFullLeaf << " ";
    cerr << indexStats.numLeaf_InternalRun << " ";

    if(verbose)
    {
        cerr << indexStats.minTreeDepth << " ";
        cerr << indexStats.avgTreeDepth << " ";
        cerr << indexStats.maxTreeDepth << " ";
    }
    if(indexStats.avg_Reindexed_TopPerLeaf.size()==1)
    {
        double chi = indexStats.real_t_list_length[0] / (double)indexStats.numLeafForTop[0].size();
        double mu = indexStats.t_list_length[0] / (double)indexStats.numLeafForTop[0].size();

        if(verbose)
        {
            cerr << chi << " ";
            cerr << indexStats.avgWeightedLeafForTop[0] << " ";
            cerr << indexStats.minVertexInFullLeaf << " ";
            cerr << indexStats.avgVertexInFullLeaf << " ";
            cerr << indexStats.maxVertexInFullLeaf << " ";

            cerr << indexStats.min_Top_PerLeaf[0] << " " << indexStats.avg_Top_PerLeaf[0] << " " << indexStats.max_Top_PerLeaf[0] << " ";
        }

        cerr << indexStats.t_list_length[0] << " ";
        cerr << indexStats.real_t_list_length[0] << " ";
        cerr << indexStats.avg_run_length[0] << " ";
        cerr << indexStats.tot_number_of_run[0] << " ";
        cerr << endl;

        if(verbose)
        {
            ivect numTinXLeaf; numTinXLeaf.assign(indexStats.maxLeafForTop[0]+1,0);
            for(unsigned j=0; j<indexStats.numLeafForTop[0].size(); j++)
            {
                numTinXLeaf[indexStats.numLeafForTop[0][j]]++;
            }
            cerr<<"[T_in_L] ";
            for(unsigned j=1; j<numTinXLeaf.size(); j++)
            {
                if(numTinXLeaf[j] > 0)
                    cerr<<"["<<j<<"] "<<numTinXLeaf[j]<<" ";
            }
            cerr<<endl;
        }

        cerr << "  explicit/vertex_encoding: "<< chi << " " <<indexStats.real_t_list_length[0] << endl;
        cerr << "  compressed_encoding: " << mu << " " << indexStats.t_list_length[0] << endl;

    }
    else
    {
        cerr << endl;
        for(unsigned i=0; i<indexStats.avg_Reindexed_TopPerLeaf.size(); i++)
        {
            if(mesh.get_top_cells_num(i) > 0)
            {
                double chi = indexStats.real_t_list_length[i] / (double)indexStats.numLeafForTop[i].size();
                double mu = indexStats.t_list_length[i] / (double)indexStats.numLeafForTop[i].size();

                cerr << "[TOP_POS_"<<i<<"] ";

                if(verbose)
                {
                    cerr << chi << " ";
                    cerr << indexStats.avgWeightedLeafForTop[i] << " ";

                    cerr << indexStats.min_Top_PerLeaf[i] << " " << indexStats.avg_Top_PerLeaf[i] << " " << indexStats.max_Top_PerLeaf[i] << " ";
                }

                cerr << indexStats.t_list_length[i] << " ";
                cerr << indexStats.real_t_list_length[i] << " ";
                cerr << indexStats.avg_run_length[i] << " ";
                cerr << indexStats.tot_number_of_run[i] << " ";

                if(verbose)
                {
                    ivect numTinXLeaf; numTinXLeaf.assign(indexStats.maxLeafForTop[i]+1,0);
                    for(unsigned j=0; j<indexStats.numLeafForTop[i].size(); j++)
                    {
                        numTinXLeaf[indexStats.numLeafForTop[i][j]]++;
                    }
                    cerr<<"[T_in_L] ";
                    for(unsigned j=1; j<numTinXLeaf.size(); j++)
                    {
                        if(numTinXLeaf[j] > 0)
                            cerr<<"["<<j<<"] "<<numTinXLeaf[j]<<" ";
                    }
                    cerr<<endl;
                }

                cerr << "  explicit_encoding: "<< chi << " " <<indexStats.real_t_list_length[i]<< " ";
                cerr << "  compressed_encoding: " << mu << " " << indexStats.t_list_length[i] << endl;
            }
        }

    }

    if(verbose)
        indexStats.print_run_histogram();

    cerr << "=== Stellar tree storage costs (in MBs) ===" << endl;
    long indexed_mesh_cost = 0;
    long hierarchy_cost = 48 * indexStats.numNode;
    long compr_refs_cost = 0, real_refs_cost = 0;

    for(int i=0; i<mesh.get_top_cells_types(); i++)
    {
        if(mesh.get_top_cells_num(i) > 0)
            indexed_mesh_cost += 4 * mesh.get_top_cell(i,1).get_vertices_num() * mesh.get_top_cells_num(i); /// int_cost * tv_size * top_i_cells_num
    }
    for(unsigned i=0; i<indexStats.real_t_list_length.size(); i++)
        real_refs_cost += 4 * indexStats.real_t_list_length[i];
    for(unsigned i=0; i<indexStats.t_list_length.size(); i++)
        compr_refs_cost += 4 * indexStats.t_list_length[i];


    std::cerr<<"indexed mesh cost: "<<indexed_mesh_cost / (1024.0*1024.0) <<" MBs"<<std::endl;
    std::cerr<<"hierarchy cost: "<<hierarchy_cost / (1024.0*1024.0) <<" MBs"<<std::endl;
    std::cerr<<"explicit references cost: "<<real_refs_cost / (1024.0*1024.0) <<" MBs"<<std::endl;
    std::cerr<<"compressed references cost: "<<compr_refs_cost / (1024.0*1024.0) <<" MBs"<<std::endl;
    std::cerr<<"TOTAL EXPLICIT Stellar treer cost: "<<(indexed_mesh_cost + hierarchy_cost + real_refs_cost) / (1024.0*1024.0) <<" MBs"<<std::endl;
    std::cerr<<"TOTAL COMPRESSED Stellar treer cost: "<<(indexed_mesh_cost + hierarchy_cost + compr_refs_cost) / (1024.0*1024.0) <<" MBs"<<std::endl;
    std::cerr<<"=== === === === === ==="<<std::endl;

    return;
}

#endif	/* _WRITER_H */
