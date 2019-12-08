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

#include "writer.h"

void Writer::write_tree(string fileName, Node_Stellar& root, ivect & original_vertices_ordering)
{
    ofstream output(fileName.c_str());

    queue<Node_Stellar*> coda;
    Node_Stellar* visited;

    bool is_root = true;

    coda.push(&root);

    while (!coda.empty()) {
        visited = coda.front();
        coda.pop();

        string begin;

        if (visited->is_leaf())
            begin = "L";
        else
            begin = "N";

        if(!is_root)
            output << endl;
        output << begin << " ";

        is_root = false;

        Writer::write_node(output, visited, original_vertices_ordering);

        if (visited->is_leaf() == false)
        {
            int counter = 0;
            for(Node_Stellar::child_iterator it=visited->begin(); it!=visited->end(); ++it)
            {
                if(*it != NULL)
                {
                    output << counter << " ";
                    coda.push(*it);
                }
                counter++;
            }
        }
    }
    output.close();
}

void Writer::write_node(ofstream& output, Node_Stellar* n, ivect &original_vertices_ordering)
{
    if (n->is_leaf())
    {
        output << n->get_real_v_array_size() << " ";

        for(int i=0; i<n->get_num_top_cells_encoded(); i++)
        {
            output << n->get_real_t_array_size(i)<< " ";
        }

        if (n->get_real_v_array_size() > 0)
        {
            output << endl << "  V ";
            for(RunIteratorPair itPair = n->make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
            {
                RunIterator const& v_id = itPair.first;
                output << original_vertices_ordering[*v_id-1] << " ";
            }
        }

        for(int i=0; i<n->get_num_top_cells_encoded(); i++)
        {
            if (n->get_real_t_array_size(i) > 0)
            {
                output << endl << "  T"<<i<<" ";
                for(RunIteratorPair itPair = n->make_t_array_iterator_pair(i); itPair.first != itPair.second; ++itPair.first)
                {
                    RunIterator const& t_id = itPair.first;
                    output << *t_id << " ";
                }
            }
        }
    }
}

void Writer::write_mesh_VTK(string mesh_name, Simplicial_Mesh &mesh)
{
    // temporary check... (to-do) extend to quad and hexahedral meshes as well)
    for(int t=0; t<mesh.get_top_cells_types(); t++)
    {
        if(mesh.get_top_cells_num(t) > 0)
        {
            if(mesh.get_type(t)==QUAD || mesh.get_type(t)==HEX || mesh.get_type(t)==nDIM_SIMPL || mesh.get_type(t)==UNSET)
            {
                cerr<<"[ERROR] write_mesh_VTK is compatible only for simplicial complexes in dimension up to 3."<<endl;
                cerr<<"        Not compatible with cell complexes and higher-dimensional CP complexes."<<endl;
                return;
            }
        }
    }

    stringstream stream;
    stream<<mesh_name<<".vtk";   

    FILE *output;
    output = fopen(stream.str().c_str(),"w");

    fprintf(output,"# vtk DataFile Version 2.0 \n \nASCII\nDATASET UNSTRUCTURED_GRID\n\n");
    fprintf(output,"POINTS %d float \n",mesh.get_vertices_num());

    for(int v=1; v<=mesh.get_vertices_num(); v++)
    {
        Vertex<COORDBASETYPE>& vert = mesh.get_vertex(v);
        vert.save_to_file_coords(output);
        fprintf(output,"\n");
    }

//    cerr<<"qui sono dentro"<<endl;
//    mesh.print_mesh_stats(cerr);

//    int t_pos = 0;
    size_t num_entries = 0;
    size_t top_cells_num = 0;
    int value = 0;
    for(int t=0; t<mesh.get_top_cells_types(); t++)
    {
//        cout<<"num_top_cell "<<mesh.get_top_cells_num(t)<<endl;
        if(mesh.get_top_cells_num(t) > 0)
        {
//            t_pos++;
//            break;
            if(mesh.get_type(t)==TETRA)
                value = 5;
            else if(mesh.get_type(t)==TRIANGLE)
                value = 4;
            else if(mesh.get_type(t)==LINE)
                value = 3;
            num_entries += mesh.get_top_cells_num(t)*value;
        }
        top_cells_num += mesh.get_top_cells_num(t);
    }

//    fprintf(output,"\nCELLS %d %d \n",mesh.get_top_cells_num(t_pos),(mesh.get_top_cells_num(t_pos)*value));
    fprintf(output,"\nCELLS %zu %zu \n",top_cells_num,num_entries); // zu for unsigned long type


    for(int t=0; t<mesh.get_top_cells_types(); t++)
    {
        if(mesh.get_top_cells_num(t) > 0)
        {
            for(int v=1; v<=mesh.get_top_cells_num(t); v++)
            {
                Top_Simplex &top = mesh.get_top_cell(t,v);
                fprintf(output,"%d ",top.get_vertices_num());
                top.save_to_file(output);
                fprintf(output,"\n");
            }
        }
    }

//    int value = 0;
//    if(mesh.get_type(t_pos)==TETRA)
//        value = 5;
//    else if(mesh.get_type(t_pos)==TRIANGLE)
//        value = 4;

//    cerr<<"t_pos: "<<t_pos<<" - top_type: " << mesh.get_type(t_pos) << endl;
//    int a; cin>>a;



//    for(int v=1; v<=mesh.get_top_cells_num(t_pos); v++)
//    {
//        Top_Simplex &t = mesh.get_top_cell(t_pos,v);
//        fprintf(output,"%d ",t.get_vertices_num());
//        t.save_to_file(output);
//        fprintf(output,"\n");
//    }

//    if(mesh.get_type(t_pos)==TETRA)
//        value = 10;
//    else if(mesh.get_type(t_pos)==TRIANGLE)
//        value = 6;

//    fprintf(output,"\nCELL_TYPES %d \n",mesh.get_top_cells_num(t_pos));
    fprintf(output,"\nCELL_TYPES %zu \n",top_cells_num);
    for(int t=0; t<mesh.get_top_cells_types(); t++)
    {
        if(mesh.get_top_cells_num(t) > 0)
        {
            if(mesh.get_type(t)==TETRA)
                value = 10;
            else if(mesh.get_type(t)==TRIANGLE)
                value = 5;
            else if(mesh.get_type(t)==LINE)
                value = 3;
            for(int v=1; v<=mesh.get_top_cells_num(t); v++)
            {
                fprintf(output,"%d ",value);
            }
        }
    }
//    for (int i = 0; i < mesh.get_top_cells_num(t_pos); ++i)
//        fprintf(output,"%d ",value);
    fprintf(output,"\n");


    fprintf(output,"POINT_DATA %d \n\n",mesh.get_vertices_num());
    fprintf(output,"FIELD FieldData 1 \n\n");
    fprintf(output,"fieldvalue 1 %d float \n",mesh.get_vertices_num());

    for (int i=1; i <= mesh.get_vertices_num(); ++i)
        fprintf(output,"%f ",mesh.get_vertex(i).getF(0)); // we save in the vtk file the first field value for visualization purposes
    fprintf(output,"\n");

    fclose(output);
}

//void Writer::write_2D_3D_mesh_VTK(string mesh_name, Simplicial_Mesh& mesh, dvect &field)
//{
//    stringstream stream;
//    stream<<mesh_name<<"_custom_field.vtk";

//    int t_pos = 0;
//    for(int t=0; t<mesh.get_top_cells_types(); t++)
//    {
//        if(mesh.get_top_cells_num(t) > 0)
//        {
//            t_pos++;
//            break;
//        }
//    }

//    FILE *output;
//    output = fopen(stream.str().c_str(),"w");

//    fprintf(output,"# vtk DataFile Version 2.0 \n \nASCII\nDATASET UNSTRUCTURED_GRID\n\n");
//    fprintf(output,"POINTS %d float \n",mesh.get_vertices_num());

//    for(int v=1; v<=mesh.get_vertices_num(); v++)
//    {
//        Vertex<COORDBASETYPE>& vert = mesh.get_vertex(v);
//        vert.save_to_file_coords(output);
//        fprintf(output,"\n");
//    }

//    int value = 0;
//    if(mesh.get_type(t_pos)==TETRA)
//        value = 5;
//    else if(mesh.get_type(t_pos)==TRIANGLE)
//        value = 4;

//    fprintf(output,"\nCELLS %d %d \n",mesh.get_top_cells_num(t_pos),(mesh.get_top_cells_num(t_pos)*value));

//    for(int v=1; v<=mesh.get_top_cells_num(t_pos); v++)
//    {
//        Top_Simplex &t = mesh.get_top_cell(t_pos,v);
//        fprintf(output,"%d ",t.get_vertices_num());
//        t.save_to_file(output);
//        fprintf(output,"\n");
//    }

//    if(mesh.get_type(t_pos)==TETRA)
//        value = 10;
//    else if(mesh.get_type(t_pos)==TRIANGLE)
//        value = 6;

//    fprintf(output,"\nCELL_TYPES %d \n",mesh.get_top_cells_num(t_pos));
//    for (int i = 0; i < mesh.get_top_cells_num(t_pos); ++i)
//        fprintf(output,"%d ",value);
//    fprintf(output,"\n");


//    fprintf(output,"POINT_DATA %d \n\n",mesh.get_vertices_num());
//    fprintf(output,"FIELD FieldData 1 \n\n");
//    fprintf(output,"fieldvalue 1 %d float \n",mesh.get_vertices_num());

//    for (unsigned i=0; i <= field.size(); ++i)
//        fprintf(output,"%f ",field[i]);
//    fprintf(output,"\n");

//    fclose(output);
//}

void Writer::write_mesh_VTK(string mesh_name, string operation_type, int vertices_per_leaf, Simplicial_Mesh &mesh)
{
    stringstream stream;
//    stream<<mesh_name<<"_kv_"<<vertices_per_leaf<<"_"<<operation_type<<".vtk";
    stream<<mesh_name<<"_kv_"<<vertices_per_leaf<<"_"<<operation_type;
    Writer::write_mesh_VTK(stream.str(),mesh);

//    int t_pos = 0;
//    for(int t=0; t<mesh.get_top_cells_types(); t++)
//    {
//        if(mesh.get_top_cells_num(t) > 0)
//        {
//            t_pos++;
//            break;
//        }
//    }

//    FILE *output;
//    output = fopen(stream.str().c_str(),"w");

//    fprintf(output,"# vtk DataFile Version 2.0 \n \nASCII\nDATASET UNSTRUCTURED_GRID\n\n");
//    fprintf(output,"POINTS %d float \n",mesh.get_vertices_num());

//    for(int v=1; v<=mesh.get_vertices_num(); v++)
//    {
//        Vertex<COORDBASETYPE>& vert = mesh.get_vertex(v);
//        vert.save_to_file_coords(output);
//        fprintf(output,"\n");
//    }

//    int value = 0;
//    if(mesh.get_type(t_pos)==TETRA)
//        value = 5;
//    else if(mesh.get_type(t_pos)==TRIANGLE)
//        value = 4;

//    fprintf(output,"\nCELLS %d %d \n",mesh.get_top_cells_num(t_pos),(mesh.get_top_cells_num(t_pos)*value));

//    for(int v=1; v<=mesh.get_top_cells_num(t_pos); v++)
//    {
//        Top_Simplex &t = mesh.get_top_cell(t_pos,v);
//        fprintf(output,"%d ",t.get_vertices_num());
//        t.save_to_file(output);
//        fprintf(output,"\n");
//    }

//    if(mesh.get_type(t_pos)==TETRA)
//        value = 10;
//    else if(mesh.get_type(t_pos)==TRIANGLE)
//        value = 6;

//    fprintf(output,"\nCELL_TYPES %d \n",mesh.get_top_cells_num(t_pos));
//    for (int i = 0; i < mesh.get_top_cells_num(t_pos); ++i)
//        fprintf(output,"%d ",value);
//    fprintf(output,"\n");


//    fprintf(output,"POINT_DATA %d \n\n",mesh.get_vertices_num());
//    fprintf(output,"FIELD FieldData 1 \n\n");
//    fprintf(output,"fieldvalue 1 %d float \n",mesh.get_vertices_num());

//    for (int i=1; i <= mesh.get_vertices_num(); ++i)
//        fprintf(output,"%f ",mesh.get_vertex(i).get_field_value());
//    fprintf(output,"\n");

//    fclose(output);
}
