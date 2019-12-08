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

#ifndef _READER_H
#define	_READER_H

#include <vector>
#include <string>
#include <set>
#include <limits>
#include <queue>
#include <cstdlib>
#include <iterator>
#include <fstream>
#include <iostream>

#include "stellar_tree/point.h"
#include "stellar_tree/box.h"
#include "stellar_decomposition/top_simplex.h"
#include "stellar_tree/mesh.h"
#include "stellar_tree/cells_soup.h"
#include "utilities/string_management.h"
#include "stellar_tree/stellar_tree.h"
#include "topological_ds/skeleton.h"

using namespace std;
///
/**
 * @brief A static class that reads file and initializes the main structures of the library
 *
 */
class Reader {
public:
    /**
     * @brief A public method that reads a file containing a tetrahedral mesh
     *
     * @param mesh represents the tetrahedral mesh to initialize
     * @param path a string argument, representing the path to the mesh file
     * @param verbose_encoding a boolean that activates the verbose encoding of the simplices of the mesh. see NOTA2 below.
     * @return bool value, true if the file is correctly readed, false otherwise
     *
     * NOTA: the tetrahedral mesh file must be in TS format
     * NOTA2: the verbose encoding initializes a vector for each p-simplex, with p>0, i.e. even those that are not explicitly encoded
     * (in this case edges and triangular faces). The verbose encoding can be a requirement for a simplification procedure
     */
    static bool read_tetra_mesh(Simplicial_Mesh& mesh, string path, bool verbose_encoding);
    /**
     * @brief A public method that reads a file containing a triangle mesh
     *
     * @param mesh represents the triangle mesh to initialize
     * @param path a string argument, representing the path to the mesh file
     * @param verbose_encoding a boolean that activates the verbose encoding of the simplices of the mesh. see NOTA2 below.
     * @return bool value, true if the file is correctly readed, false otherwise
     *
     * NOTA: the triangle mesh file must be in TRI format
     * NOTA2: the verbose encoding initializes a vector for each p-simplex, with p>0, i.e. even those that are not explicitly encoded
     * (in this case edges). The verbose encoding can be a requirement for a simplification procedure
     */
    static bool read_triangle_mesh(Simplicial_Mesh& mesh, string path, bool verbose_encoding);
    /**
     * @brief A public method that reads a file containing a n-dimensional mesh
     *
     * @param mesh represents the n-dimensional mesh to initialize
     * @param path a string argument, representing the path to the mesh file
     * @return bool value, true if the file is correctly readed, false otherwise
     *
     * NOTA: the n-dimensional mesh file must be in SIERP format
     */
    static bool read_sierp_mesh(Simplicial_Mesh& mesh, string path);
    /**
     * @brief A public method that reads a file containing a non-manifold mesh
     *
     * @param mesh represents the non-manifold mesh to initialize
     * @param path a string argument, representing the path to the mesh file
     * @return bool value, true if the file is correctly readed, false otherwise
     *
     * NOTA: the non-manifold mesh file must be in GMV format
     */
    static bool read_gmv_mesh(Simplicial_Mesh& mesh, string path);
    /**
     * @brief A public method that reads a file containing a terrain
     *
     * @param mesh represents the terrain to initialize
     * @param path a string argument, representing the path to the mesh file
     * @param verbose_encoding a boolean that activates the verbose encoding of the simplices of the mesh. see NOTA2 below.
     * @return bool value, true if the file is correctly readed, false otherwise
     *
     * NOTA: the terrain file must be in TRI format
     * NOTA2: the verbose encoding initializes a vector for each p-simplex, with p>0, i.e. even those that are not explicitly encoded
     * (in this case edges). The verbose encoding can be a requirement for a simplification procedure
     */
    static bool read_terrain_mesh(Simplicial_Mesh& mesh, string path, bool verbose_encoding);
    /**
     * @brief A public method that reads a file containing a CP complex
     *
     * @param mesh represents the CP complex to initialize
     * @param path a string argument, representing the path to the mesh file
     * @return bool value, true if the file is correctly readed, false otherwise
     *
     * NOTA: the quad mesh file must be in OFF format
     */
    static bool read_cell_mesh(CP_Mesh& mesh, string path);
    /**
     * @brief A public method that reads a file containing a hexahedral mesh
     *
     * @param mesh represents the hexahedral mesh to initialize
     * @param path a string argument, representing the path to the mesh file
     * @return bool value, true if the file is correctly readed, false otherwise
     *
     * NOTA: the hexahedral mesh file must be in HEX format
     */
    static bool read_hexa_mesh(CP_Mesh& mesh, string path);
    /**
     * @brief A public method that reads a file containing a soup of top cells
     *
     * @param mesh represents the soup to initialize
     * @param path a string argument, representing the path to the mesh file
     * @return bool value, true if the file is correctly readed, false otherwise
     *
     * NOTA: the soup file must be in SOUP format
     */
    static bool read_soup(Cells_Soup& soup, string path);
    /**
     * @brief A public method that reads a file containing a simplicial mesh in arbitrary dimensions
     *
     * @param mesh represents the simplicial mesh to initialize
     * @param path a string argument, representing the path to the mesh file
     * @param verbose_encoding a boolean that activates the verbose encoding of the simplices of the mesh. see NOTA3 below.
     * @return bool value, true if the file is correctly readed, false otherwise
     *
     * NOTA: the simplicial mesh file must be in OFF format
     * NOTA2: the top simplices in the input file must be ordered considering their dimension
     *        thus, first 1-simplices (if any), then 2-simplices (if any) and so on up to p-simplices.
     * NOTA3: the verbose encoding initializes a vector for each p-simplex, with p>0, i.e. even those that are not explicitly encoded
     * (in this case edges). The verbose encoding can be a requirement for a simplification procedure
     */
    static bool read_simplicial_mesh(Simplicial_Mesh& mesh, string path, bool verbose_encoding, bool force3D=false);
    /**
     * @brief A public method that read a list of vertices
     *    
     * @param mesh represents the soup to initialize
     * @param path a string argument, representing the path to the mesh file
     * @return bool value, true if the file is correctly readed, false otherwise
     *
     * NOTA: the input file can be in TS, TRI or OFF format
     *
     */
    template<class C, class T> static bool read_vertices(Mesh<C,T>& mesh, string path/*, bool force3D=false*/);

    /**
     * @brief A public method that reads a file containing a tree
     *
     * @param root represents the root node of the tree
     * @param division represents the spatial subdivision used by the tree
     * @param fileName a string argument, represents the path to the mesh file
     * @return a boolean value, true if the file is correctly readed, false otherwise
     */
    static bool read_tree(Node_Stellar& root, Spatial_Subdivision &division, string fileName);

    template<class C, class T> static bool read_skeleton(skeleton_graph &graph, Mesh<C,T>& mesh, string path);

private:    
    /**
     * @brief A constructor method
     *
     * NOTA: the class is not instantiable
     *
     */
    Reader() {}    
    /**
     * @brief A constructor method
     *
     * NOTA: the class is not instantiable
     * @param
     */
    Reader(const Reader& ) {}
    /**
     * @brief A destructor method
     *
     */
    virtual ~Reader() {}

    /**
     * @brief A private method that reads a list of vertices with field value
     *
     * @param mesh represents the mesh to initialize
     * @param input represents the input stream
     * @param num_vertices represents the number of vertices of the mesh
     * @param v_dimension represents the dimension of the vertices
     *
     * NOTA: the procedure also set the mesh domain
     */
    template<class C, class T> static void read_vertices_list_with_field(Mesh<C,T>& mesh, ifstream &input, int num_vertices, int v_dimension);
    /**
     * @brief A private method that reads a list of vertices with just the coordinates
     *
     * @param mesh represents the mesh to initialize
     * @param input represents the input stream
     * @param num_vertices represents the number of vertices of the mesh
     * @param v_dimension represents the dimension of the vertices
     *
     * NOTA: the procedure also set the mesh domain
     */
    template<class C, class T> static void read_vertices_list_no_field(Mesh<C,T>& mesh, ifstream &input, int num_vertices);
    /**
     * @brief A private method that reads a list of vertices with a list of fields of which only the first 3 are considered as coordinates.
     * The others are encoded in the field values array
     *
     * @param mesh represents the mesh to initialize
     * @param input represents the input stream
     * @param num_vertices represents the number of vertices of the mesh
     * @param v_dimension represents the dimension of the vertices
     *
     * NOTA: the procedure also set the mesh domain
     */
    template<class C, class T> static void read_vertices_list_force3D(Mesh<C,T>& mesh, ifstream &input, int num_vertices);
    /**
     * @brief A private method that reads a list of top cells from file
     *
     * @param mesh represents the mesh to initialize
     * @param input represents the input stream
     * @param num_top_cells represents the number of top cells of the mesh
     * @param numV represents the number of vertices of a top cell
     * @param top_cells_pos represents the position of the top cells in the mesh array
     */
    static void read_top_cells(Simplicial_Mesh &mesh, ifstream &input, int num_top_cells, int numV, int top_cells_pos);
    /**
     * @brief A private method that reads a node into the file and saves the information readed (Generic Version)
     *
     * @param n a Node_Stellar* argument, representing the empty node to be set
     * @param input represents the input stream
     * @param division represents the spatial subdivision used by the tree
     */
    static void read_node(Node_Stellar* n, ifstream& input, Spatial_Subdivision& division);
    /**
     * @brief A protected method that reads a leaf into the file and saves the information readed
     *
     * @param n a Node_Stellar* argument, representing the empty node to be set
     * @param input represents the input stream
     * @param tokens a vector<string>& argument, representing the a row of the file
     */
    static void read_leaf_node(Node_Stellar* n, ifstream& input, vector<string>& tokens);
};

template<class C, class T> bool Reader::read_vertices(Mesh<C,T> &mesh, string path/*, bool force3D*/)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string extension = string_management::get_file_extension(path);

    string line, trash;
    int num_vertices;

    if(extension == "ts")
    {
        getline(input, line);
        int tpos = line.find_first_of(' ');

        int num_vertices = atoi(line.substr(0, tpos).c_str());

        if (num_vertices == 0) {
            cerr << "This is not a valid .ts file: " << path << endl;
            return false;
        }

        Reader::read_vertices_list_with_field(mesh,input,num_vertices,3);
    }
    else if(extension == "tri")
    {
        input >> num_vertices;

        if (num_vertices == 0)
        {
            cerr << "This is not a valid .tri file: " << path << endl;
            return false;
        }

        Reader::read_vertices_list_no_field(mesh,input,num_vertices/*,3*/);
    }
    else if(extension == "off")
    {
        /// INIT PHASE
        input >> trash;
        input >> num_vertices;
        input >> trash;
        input >> trash;

        if (num_vertices == 0)
        {
            cerr << "This is not a valid .off file: " << path << endl;
            return false;
        }

//        if(force3D)
//            Reader::read_vertices_list_force3D(mesh,input,num_vertices);
//        else
            Reader::read_vertices_list_no_field(mesh,input,num_vertices/*,3*/);
    }

    return true;
}

template<class C, class T> void Reader::read_vertices_list_with_field(Mesh<C,T> &mesh, ifstream &input, int num_vertices, int v_dimension)
{
    mesh.reserve_vertices_space(num_vertices);

    vector<C> coords;
    coords.assign(v_dimension,-1);
    C field;

    //legge i vertici aggiustando il dominio..
    for (int i = 0; i < num_vertices; i++)
    {
        for(int j=0; j<v_dimension; j++)
            input >> coords[j];
        input >> field;

        if (input.eof())
            break;

        Vertex<C> v = Vertex<C>(coords,field);
        mesh.add_vertex(v);
        if (i == 0) {
            Box<C> b = Box<C>(v, v);
            mesh.set_domain(b);
        } else {
            mesh.get_domain().resize_to_contain(v);
        }
    }
}

template<class C, class T> void Reader::read_vertices_list_no_field(Mesh<C,T>& mesh, ifstream &input, int num_vertices)
{
    mesh.reserve_vertices_space(num_vertices);

    string line;
    getline(input,line); /// to clear up any '\n' symbol

    for (int i = 0; i < num_vertices; i++)
    {
        getline(input,line);
        istringstream iss(line);
        string coord;
        vector<C> coords;
        while( getline(iss, coord, ' '))
        {
            if(coord != "")
                coords.push_back(atof(coord.c_str()));
        }

        if (input.eof())
            break;

        //we set a fictious field value
        Vertex<C> v = Vertex<C>(coords,coords[coords.size()-1]);
        mesh.add_vertex(v);
        if (i == 0)
        {
            Box<C> b = Box<C>(v, v);
            mesh.set_domain(b);
        }
        else
        {
            mesh.get_domain().resize_to_contain(v);
        }
    }
}

template<class C, class T> void Reader::read_vertices_list_force3D(Mesh<C,T>& mesh, ifstream &input, int num_vertices)
{
//    cout<<"ci entro"<<endl;
    mesh.reserve_vertices_space(num_vertices);

    string line;
    getline(input,line); /// to clear up any '\n' symbol

    for (int i = 0; i < num_vertices; i++)
    {
        getline(input,line);
        istringstream iss(line);
        string coord;
        vector<C> all_values;
        while( getline(iss, coord, ' '))
        {
            if(coord != "")
                all_values.push_back(atof(coord.c_str()));
        }

        if (input.eof())
            break;

        /// THIS IS NOT WORKING ON 2D DATA
        vector<C> coords(all_values.begin(),all_values.begin()+3); /// first three positions are coordinates used for the decomposition
        vector<C> fields(all_values.begin()+3,all_values.end());
//        cout<<"coords_size: "<<coords.size()<<" -- fields: "<<fields.size()<<endl;
//        int a; cin>>a;
        if(fields.empty()) /// if we have a 3D dataset we add a dummy field value
            fields = { coords[2] };

        Vertex<C> v = Vertex<C>(coords,fields);
        mesh.add_vertex(v);
        if (i == 0)
        {
            Box<C> b = Box<C>(v, v);
            mesh.set_domain(b);
        }
        else
        {
            mesh.get_domain().resize_to_contain(v);
        }
    }
}

template<class C, class T> bool Reader::read_skeleton(skeleton_graph &graph, Mesh<C,T>& mesh, string path)
{
    ifstream input(path.c_str());
    graph = skeleton_graph(mesh.get_vertices_num());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    int v1, v2;
    while(input.good())
    {
        input >> v1;
        input >> v2;
//        cout<<v1<<" "<<v2<<endl;
//        cout<<mesh.get_vertices_num()<<endl;

        if(input.eof())
            break;

        graph.insert(v1,v2+1); // NOTA: the vertices range of the input file is different from
        graph.insert(v2,v1+1); //       the one used internally by the library
//        cout<<"inserted"<<endl;
    }

    cout<<"end of reading"<<endl;
    graph.check_skeleton_graph();

    return true;
}

#endif	/* _READER_H */

