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

#include "reader.h"
#include <sstream>
#include <algorithm>

bool Reader::read_tetra_mesh(Simplicial_Mesh &mesh, string path, bool verbose_encoding)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string line;
    getline(input, line);
    int tpos = line.find_first_of(' ');

    int num_vertices = atoi(line.substr(0, tpos).c_str());
    int num_top_simplices = atoi(line.substr(tpos).c_str());

    if (num_vertices == 0 || num_top_simplices == 0) {
        cerr << "This is not a valid .ts file: " << path << endl;
        return false;
    }

    Reader::read_vertices_list_with_field(mesh,input,num_vertices,3);

    int t_pos;
    if(verbose_encoding)
    {
        mesh.init_top_cells_array(3);
        mesh.add_type(LINE);
        mesh.add_type(TRIANGLE);
        mesh.add_type(TETRA);
        mesh.reserve_top_cells_space(0,num_top_simplices); ///
        mesh.reserve_top_cells_space(1,num_top_simplices); ///
        mesh.reserve_top_cells_space(2,num_top_simplices); ///
        t_pos = 2;
    }
    else
    {
        mesh.init_top_cells_array(1);
        mesh.add_type(TETRA);
        mesh.reserve_top_cells_space(0,num_top_simplices);
        t_pos = 0;
    }

    //legge i top simplessi
    Reader::read_top_cells(mesh,input,num_top_simplices,4,t_pos);

    return true;
}

bool Reader::read_triangle_mesh(Simplicial_Mesh &mesh, string path, bool verbose_encoding)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    int num_vertices;
    input >> num_vertices;

    if (num_vertices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    Reader::read_vertices_list_no_field(mesh,input,num_vertices/*,3*/);

    int num_top_simplices;
    input >> num_top_simplices;

    if(num_top_simplices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    int t_pos;
    if(verbose_encoding)
    {
        mesh.init_top_cells_array(2);
        mesh.add_type(LINE);
        mesh.add_type(TRIANGLE);
        mesh.reserve_top_cells_space(0,num_top_simplices); ///
        mesh.reserve_top_cells_space(1,num_top_simplices); ///
        t_pos = 1;
    }
    else
    {
        mesh.init_top_cells_array(1);
        mesh.add_type(TRIANGLE);
        mesh.reserve_top_cells_space(0,num_top_simplices);
        t_pos = 0;
    }

    //legge i top simplessi
    Reader::read_top_cells(mesh,input,num_top_simplices,3,t_pos);

    return true;
}

bool Reader::read_terrain_mesh(Simplicial_Mesh& mesh, string path, bool verbose_encoding)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    int num_vertices;
    input >> num_vertices;

    if (num_vertices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    Reader::read_vertices_list_with_field(mesh,input,num_vertices,2);

    int num_top_simplices;
    input >> num_top_simplices;

    if(num_top_simplices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    int t_pos;
    if(verbose_encoding)
    {
        mesh.init_top_cells_array(2);
        mesh.add_type(LINE);
        mesh.add_type(TRIANGLE);
        mesh.reserve_top_cells_space(0,num_top_simplices); ///
        mesh.reserve_top_cells_space(1,num_top_simplices); ///
        t_pos = 1;
    }
    else
    {
        mesh.init_top_cells_array(1);
        mesh.add_type(TRIANGLE);
        mesh.reserve_top_cells_space(0,num_top_simplices);
        t_pos = 0;
    }

    //legge i top simplessi
    Reader::read_top_cells(mesh,input,num_top_simplices,3,t_pos);

    return true;
}

void Reader::read_top_cells(Simplicial_Mesh &mesh, ifstream &input, int num_top_cells, int numV, int top_cells_pos)
{
    for (int i = 0; i < num_top_cells; i++) {
        ivect v;
        v.assign(numV,0);
        for (int j = 0; j < numV; j++) {
            int index;
            input >> index;
            v[j] = index+1;// the vertices position indexes starts from 1 (thus we have to increment the actual index of the file)
        }
        Top_Simplex t = Top_Simplex(v);
        mesh.add_top_cell(top_cells_pos,t);
    }
}

/// OFF FORMAT!
bool Reader::read_cell_mesh(CP_Mesh& mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    /// INIT PHASE
    int num_vertices, num_cells;
    string trash;
    input >> trash;
    input >> num_vertices;
    input >> num_cells;
    input >> trash;


    if (num_vertices == 0 || num_cells == 0)
    {
        cerr << "This is not a valid .off file: " << path << endl;
        return false;
    }

    Reader::read_vertices_list_no_field(mesh,input,num_vertices/*,3*/);

    mesh.init_top_cells_array(1);

    mesh.reserve_top_cells_space(0,num_cells);

    int cell_num_v;

    for (int i = 0; i < num_cells; i++) {
        ivect v;
        input >> cell_num_v; //recover the dimension
        v.assign(cell_num_v,0);
        for (int j = 0; j < cell_num_v; j++) {
            int index;
            input >> index;
            v[j] = index+1;// the vertices position indexes starts from 1 (thus we have to increment the actual index of the file)
        }
        Top_CP_Cell t = Top_CP_Cell(v);
        mesh.add_top_cell(0,t);
    }

    if(cell_num_v==4) /// we have a quad mesh
        mesh.add_type(QUAD);
    else if(cell_num_v==8) /// we have a quad mesh
        mesh.add_type(HEX);

    return true;
}

bool Reader::read_sierp_mesh(Simplicial_Mesh& mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    int dimension;
    input >> dimension;

    int num_vertices;
    input >> num_vertices;
    int num_topSimplexes;
    input >> num_topSimplexes;

    if (num_vertices == 0 || num_topSimplexes == 0) {
        cerr << "This is not a valid .sierp file: " << path << endl;
        return false;
    }

    Reader::read_vertices_list_no_field(mesh,input,num_vertices/*,dimension*/);

    mesh.init_top_cells_array(1);
    mesh.add_type(nDIM_SIMPL);
    mesh.reserve_top_cells_space(0,num_topSimplexes);

    int index;
    //legge i top simplessi
    for (int i = 0; i < num_topSimplexes; i++)
    {
        ivect v;
        v.assign(dimension+1,-1);

        input >> index; //trash

        for (unsigned j = 0; j < v.size(); j++) {
            input >> index;
            v[j] = index+1;// the vertices position indexes starts from 1 (thus we have to increment the actual index of the file)
        }
        Top_Simplex t = Top_Simplex(v/*,nDIM_SIMPL*/);

        mesh.add_top_cell(0,t); /// we have just a single top simplex which is stored in the first position of the global array
    }

    return true;
}

bool Reader::read_gmv_mesh(Simplicial_Mesh &mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string trash = "";
    while( trash != "nodes" )
        input >> trash;
    int num_vertices;
    input>>num_vertices;
    vector<COORDBASETYPE> xs, ys, zs;
    xs.assign(num_vertices,0.0);
    ys.assign(num_vertices,0.0);
    zs.assign(num_vertices,0.0);
    //read the coordinates from file
    for(int i=0; i<num_vertices; i++)
        input>>xs[i];
    for(int i=0; i<num_vertices; i++)
        input>>ys[i];
    for(int i=0; i<num_vertices; i++)
        input>>zs[i];
    //and save them in the structure
    vector<COORDBASETYPE> coords; coords.assign(3,0.0);
    for(int i=0; i<num_vertices; i++)
    {
        coords[0]=xs[i];
        coords[1]=ys[i];
        coords[2]=zs[i];

        Vertex<COORDBASETYPE> v = Vertex<COORDBASETYPE>(coords,coords[2]);
        mesh.add_vertex(v);

        if (i == 0) {
            Box<COORDBASETYPE> b = Box<COORDBASETYPE>(v, v);
            mesh.set_domain(b);
        } else {
            mesh.get_domain().resize_to_contain(v);
        }
    }

    input>>trash;
    int num_topSimplexes=0;
    if(trash == "cells")
        input>>num_topSimplexes;
    else
        cerr<<"something wrong here..."<<endl;

    vector<Top_Simplex> top_simplex;
    top_simplex.reserve(num_topSimplexes);

    vector<CellType> temp_top_types;
    temp_top_types.assign(num_topSimplexes,UNSET);

    int index=0;
    bool has_line = false, has_tri = false, has_quad = false, has_tetra = false;

    //read the top simplices from file
    for(int i=0; i<num_topSimplexes; i++)
    {
        input>>trash;

        if(trash == "line")
        {
            has_line = true;
            temp_top_types[i] = LINE;
        }
        else if(trash == "tri")
        {
            has_tri = true;
            temp_top_types[i] = TRIANGLE;
        }
        else if(trash == "quad")
        {
            has_quad = true;
            temp_top_types[i] = QUAD;
        }
        else if(trash == "tet")
        {
            has_tetra = true;
            temp_top_types[i] = TETRA;
        }

        input>>index; //get num of vertex
        ivect v;
        v.assign(index,-1);

        for (unsigned j = 0; j < v.size(); j++) {
            input >> index;
            v[j] = index;// the vertices position indexes starts from 1 (thus we have to increment the actual index of the file)
            /// nota 12 settembre 2013 -> sembra che questo formato faccia partire gli indici dei vertici da 1 quindi non serve mettere il +1 come negli altri casi
        }

        Top_Simplex t = Top_Simplex(v);
        top_simplex.push_back(t);
    }

    //and we save them in the structure
    //we must know the number of non-manifold non-top simplices
    //in order to understand in which array to place them
    int pos_quad = -1, pos_tri = -1, pos_line =-1, pos_tetra = -1; /// the default pos_max value is 0 because if we have just a type of top simplex we have it at position 0

    int pos_counter = 0;
    if(has_line)
    {
        mesh.add_type(LINE);
        pos_line = pos_counter;
        pos_counter++;
    }
    if(has_tri)
    {
        mesh.add_type(TRIANGLE);
        pos_tri = pos_counter;
        pos_counter++;
    }
    if(has_quad)
    {
        mesh.add_type(QUAD);
        pos_quad = pos_counter;
        pos_counter++;
    }
    if(has_tetra)
    {
        mesh.add_type(TETRA);
        pos_tetra = pos_counter;
        pos_counter++;
    }
    mesh.init_top_cells_array(pos_counter);

    //in the end we place the top simplices in the correct list
    for(unsigned i=0; i<top_simplex.size(); i++)
    {
        if(temp_top_types[i] == TETRA)
            mesh.add_top_cell(pos_tetra,top_simplex[i]);
        else if(temp_top_types[i] == QUAD)
        {
            mesh.add_top_cell(pos_quad,top_simplex[i]);
        }
        else if(temp_top_types[i] == TRIANGLE)
        {
            mesh.add_top_cell(pos_tri,top_simplex[i]);
        }
        else if(temp_top_types[i] == LINE)
        {
            mesh.add_top_cell(pos_line,top_simplex[i]);
        }
    }
    return true;
}

bool Reader::read_hexa_mesh(CP_Mesh& mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string trash = "";
    while( trash != "v" )
        input >> trash;

    vector<COORDBASETYPE> coords;
    coords.assign(3,-1);
    bool is_first = true;

    while( trash != "#" )
    {
        for(unsigned j=0; j<coords.size(); j++)
        {
            input >> coords[j];
        }

        input >> trash;

        Vertex<COORDBASETYPE> v = Vertex<COORDBASETYPE>(coords,coords[2]);

        mesh.add_vertex(v);
        if (is_first)
        {
            Box<COORDBASETYPE> b = Box<COORDBASETYPE>(v, v);
            mesh.set_domain(b);
            is_first = false;
        } else {
            mesh.get_domain().resize_to_contain(v);
        }
    }

    mesh.init_top_cells_array(1);
    mesh.add_type(HEX);

    int index;
    while(input.good())
    {
        input >> trash;
        if(trash != "h")
        {
            while( trash != "h" && input.good())
                input >> trash;
        }

        if(input.eof())
            break;

        ivect v;
        v.assign(8,-1);

        for (unsigned j = 0; j < v.size(); j++) {
            input >> index;
            v[j] = index;
        }
        Top_CP_Cell t = Top_CP_Cell(v);
        mesh.add_top_cell(0,t);

        if(input.eof())
            break;
    }

    return true;
}

bool Reader::read_soup(Cells_Soup& soup, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    int v_dim, num_top_types;
    ivect top_num;
    input >> v_dim;
    input >> num_top_types;

    for(int i=0; i<num_top_types; i++)
    {
        int t; input >> t;
        top_num.push_back(t);
    }

    bool is_first = true;

    soup.init_top_cells_array(num_top_types,top_num);

    for(int d=0; d<num_top_types; d++)
    {
        for(int i=0; i<top_num[d]; i++)
        {
            Explicit_Cell<COORDBASETYPE> cell;
            int num_v=0;
            string type;
            input >> type;

            if(type == "line")
            {
                num_v = 2;
                if(i==0)
                    soup.add_type(LINE);
            }
            else if(type == "tri")
            {
                num_v = 3;
                if(i==0)
                    soup.add_type(TRIANGLE);
            }
            else if(type == "quad")
            {
                num_v = 4;
                if(i==0)
                    soup.add_type(QUAD);
            }
            else if(type == "tet")
            {
                num_v = 4;
                if(i==0)
                    soup.add_type(TETRA);
            }
            else if(type == "hex")
            {
                num_v = 8;
                if(i==0)
                    soup.add_type(HEX);
            }
            else if(type == "sierp")
            {
                num_v = 1 + v_dim;
                if(i==0)
                    soup.add_type(nDIM_SIMPL);
            }

            for(int v=0; v<num_v; v++)
            {
                Vertex<COORDBASETYPE> tmp(v_dim-1);
                double coord;
                for(int c=0; c<v_dim; c++)
                {
                    input >> coord;
                    tmp.setC(c,coord);
                }
                cell.add_vertex(tmp);

                if (is_first)
                {
                    Box<COORDBASETYPE> b = Box<COORDBASETYPE>(tmp, tmp);
                    soup.set_domain(b);
                    is_first = false;
                } else {
                    soup.get_domain().resize_to_contain(tmp);
                }
            }
            soup.add_cell(d,cell);
        }
    }

    input.close();

    return true;
}

/// OFF FORMAT!
bool Reader::read_simplicial_mesh(Simplicial_Mesh& mesh, string path, bool verbose_encoding, bool force3D)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    /// INIT PHASE
    int num_vertices, num_simplices;
    string trash;
    input >> trash;
    input >> num_vertices;
    input >> num_simplices;
    input >> trash;

    if (num_vertices == 0 || num_simplices == 0)
    {
        cerr << "This is not a valid .off file: " << path << endl;
        return false;
    }

    if(force3D)
        Reader::read_vertices_list_force3D(mesh,input,num_vertices);
    else
        Reader::read_vertices_list_no_field(mesh,input,num_vertices/*,3*/);

    map<int,int> top_pos_map; // associative array that keep track of the position index of top d-simplices
    map<int,int> dtop_num; // associative array that keep track of the number of top d-simplices
    int num_v = 0, current_num_v = 0;

    if(verbose_encoding)
    {
        //we have to exploit the maximum top-simplex dimension
        int max_top_dim = 0;
        for (int i = 0; i < num_simplices; i++)
        {
            input >> num_v;

            for (int j = 0; j < num_v; j++)
                input >> trash;

            if(current_num_v != num_v)
            {
                top_pos_map.insert(make_pair(num_v,top_pos_map.size()));

                current_num_v = num_v;
            }

            dtop_num[num_v]++;

            if(num_v > max_top_dim)
                max_top_dim = num_v;
        }
        mesh.init_top_cells_array(max_top_dim-1); // init the top simplices lists (-1 as we encode the vertices in a separate array)
        // then init the top simplices types
        for(int i=0; i<mesh.get_top_cells_types(); i++)
        {
            if(i==0)
                mesh.add_type(LINE);
            else if(i==1)
                mesh.add_type(TRIANGLE);
            else if(i==2)
                mesh.add_type(TETRA);
            else
                mesh.add_type(nDIM_SIMPL);
            top_pos_map[i+2]=i;
        }
    }
    else // just encode what it is encoded in the input file
    {
        for (int i = 0; i < num_simplices; i++)
        {
            input >> num_v;

            for (int j = 0; j < num_v; j++)
                input >> trash;

            if(current_num_v != num_v)
            {
                top_pos_map.insert(make_pair(num_v,top_pos_map.size()));
                mesh.add_top_cells_array();
                current_num_v = num_v;

                if(current_num_v==2)
                    mesh.add_type(LINE);
                else if(current_num_v==3)
                    mesh.add_type(TRIANGLE);
                else if(current_num_v==4)
                    mesh.add_type(TETRA);
                else
                    mesh.add_type(nDIM_SIMPL);
            }

            dtop_num[num_v]++;
        }
    }

    for(map<int,int>::iterator it=dtop_num.begin(); it!=dtop_num.end(); ++it)
    {
        mesh.reserve_top_cells_space(top_pos_map[it->first],it->second);
    }

    input.close();
    input.clear();
    input.open(path.c_str());

    string_management::go_to_line(input,num_vertices+2); // we set the iterator back to the beginning of the top simplices list

    for (int i = 0; i < num_simplices; i++) {
        ivect v;
        input >> num_v; //recover the dimension

        v.assign(num_v,0);
        for (int j = 0; j < num_v; j++) {
            int index;
            input >> index;
            v[j] = index+1;// the vertices position indexes starts from 1 (thus we have to increment the actual index of the file)
        }
        Top_Simplex t = Top_Simplex(v);
        mesh.add_top_cell(top_pos_map[num_v],t);
    }

//    mesh.print_mesh_stats(cout);
//    int a; cin>>a;

    return true;
}

bool Reader::read_tree(Node_Stellar &root, Spatial_Subdivision &division, string fileName)
{
    ifstream input(fileName.c_str());

    if (input.is_open() == false)
    {
        cerr << "Error opening tree file. The file does not exist or is unreadable." << endl;
        return false;
    }

    queue<Node_Stellar*> coda;
    Node_Stellar* current;

    coda.push(&root);

    while (input.good())
    {
        current = coda.front();
        coda.pop();

        Reader::read_node(current,input,division);

        if (!current->is_leaf())
        {
            for(Node_Stellar::child_iterator it=current->begin(); it!=current->end(); ++it)
            {
                if(*it != NULL)
                    coda.push(*it);
            }
        }

        if (input.eof() || coda.empty())
            break;
    }

    return true;
}

void Reader::read_node(Node_Stellar *n, ifstream &input, Spatial_Subdivision &division)
{
    string line;
    vector<string> tokens;
    getline(input, line);
    istringstream iss(line);
    copy(istream_iterator<string > (iss), istream_iterator<string > (), back_inserter<vector<string> >(tokens));

    if (tokens[0] == "N")
    {
        n->init_sons(division.get_son_number());
        for(unsigned i=1; i<tokens.size(); i++)
        {
            Node_Stellar* s = new Node_Stellar();

            int child_pos = std::atol(tokens[i].c_str());
            if(child_pos < 0)
            {
                cerr<<"[readNode] the child position cannot be negative.. somehing wrong during reading happen"<<endl;
                cerr<<"original string value: "<<tokens[i]<<endl;
                int a; cin>>a;
            }
            n->set_son(s,child_pos);
        }
    }
    else if (tokens[0] == "L")
    {
        Reader::read_leaf_node(n,input,tokens);
    }
}

void Reader::read_leaf_node(Node_Stellar *n, ifstream& input, vector<string>& tokens)
{
    string line;
    int vertices = atoi(tokens[1].c_str());

    ivect tops;
    for(unsigned i=2;i<tokens.size();i++)
    {
        tops.push_back(atoi(tokens[i].c_str()));
    }

    if(vertices > 0)
    {
        getline(input, line);
        vector<string> tokens2;
        istringstream iss2(line);
        copy(istream_iterator<string > (iss2),
             istream_iterator<string > (),
             back_inserter<vector<string> >(tokens2));
        for (unsigned int i = 1; i < tokens2.size(); i++)
            n->add_vertex(atoi(tokens2[i].c_str()));
    }

    n->init_sub_cells_vectors(tops.size());

    for(unsigned j=0;j<tops.size();j++)
    {
        if(tops[j] > 0)
        {
            getline(input, line);
            vector<string> tokens2;
            istringstream iss2(line);
            copy(istream_iterator<string > (iss2),
                 istream_iterator<string > (),
                 back_inserter<vector<string> >(tokens2));
            for (unsigned i = 1; i < tokens2.size(); i++)
                n->add_top_cell(j,atoi(tokens2[i].c_str()));
        }
    }
}
