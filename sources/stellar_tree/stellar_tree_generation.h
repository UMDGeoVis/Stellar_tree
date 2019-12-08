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

template<class C, class T, class R> void Stellar_Tree::build(Mesh<C,T> &mesh, string filename, R &reindexer, bool save_original_positions)
{
    Timer t;
    double tot = 0;

    cerr<<"--Generation Timings--"<<endl;
    t.start();
    add_vertices(mesh);
    t.stop();
    t.print_elapsed_time("   Add Vertices: ");
    tot += t.get_elapsed_time();

//    cout<<"entering reorder vertices"<<endl;

    t.start();
    reindexer.reorganize_mesh_vertices(*this,mesh,true); // always true because we need to save the tree structure considering the original vertices position indexes
    t.stop();
    t.print_elapsed_time("   Exploit Vertices Spatial Locality: ");
//    tot += t.get_elapsed_time();

    t.start();
    add_top_cells(mesh);
    t.print_elapsed_time("   Add Top-cells: ");
    t.stop();
    tot += t.get_elapsed_time();

    // writing the tree structure with the ORIGINAL vertices ordering
    Writer::write_tree(filename,root,reindexer.get_original_vertices_ordering());
    if(!save_original_positions)
        reindexer.reset_vertices_variables();

    t.start();
    reindexer.reorganize_mesh_top_simplices(*this,mesh,save_original_positions);
    reindexer.reset();
    t.print_elapsed_time("   Exploit Top-cells Spatial Locality: ");
    t.stop();
    tot += t.get_elapsed_time();

    cerr<< kv << " Stellar tree [Build and Reordering] "<<tot<<endl;
}

template<class C, class T, class R> void Stellar_Tree::build(Cells_Soup &soup, Mesh<C,T> & mesh, R &reindexer, bool save_original_positions)
{
    Timer t;
    double tot = 0;
    t.start();
    cerr<<"--Generation from SOUP Timings--"<<endl;
    t.start();
    add_vertices(soup,mesh);
    t.stop();
    t.print_elapsed_time("   Add Vertices: ");
    tot += t.get_elapsed_time();
    soup.clear();

    t.start();
    reindexer.reorganize_mesh_vertices(*this,mesh,save_original_positions);
    t.stop();
    t.print_elapsed_time("   Exploit Vertices Spatial Locality: ");
    tot += t.get_elapsed_time();

    t.start();
    add_top_cells(mesh);
    t.print_elapsed_time("   Add Top-cells: ");
    t.stop();
    tot += t.get_elapsed_time();

    t.start();
    reindexer.reorganize_mesh_top_simplices(*this,mesh,save_original_positions);
    t.stop();
    t.print_elapsed_time("   Exploit Top-cells Spatial Locality: ");
    tot += t.get_elapsed_time();
    cerr<< kv << " Stellar tree [Build and Reordering] "<<tot<<endl;
}

template<class C, class T, class R> void Stellar_Tree::generate_spatial_subdivision(Mesh<C,T> &mesh, R &reindexer)
{
    Timer t;
    double tot = 0;

    cerr<<"--Generation Timings--"<<endl;
    t.start();
    add_vertices(mesh);
    t.stop();
    t.print_elapsed_time("   Add Vertices: ");
    tot += t.get_elapsed_time();

    t.start();
    reindexer.reorganize_mesh_vertices(*this,mesh,true);
    t.stop();
    t.print_elapsed_time("   Exploit Vertices Spatial Locality: ");
    tot += t.get_elapsed_time();

    cerr<< kv << " Stellar tree [Build and Reordering] "<<tot<<endl;
}

template<class C, class T, class R> void Stellar_Tree::add_and_resort_top_cells(Mesh<C, T> &mesh, R &reindexer)
{
    Timer t;
    double tot = 0;

    // we have to init these lists, when we generate the index in two distinct phases,
    // after generating the Vietori Rips complex
    this->init_tops_lists(this->root,mesh.get_top_cells_types());

    cerr<<"--Generation Timings--"<<endl;
    t.start();
    add_top_cells(mesh);
    t.print_elapsed_time("   Add Top-cells: ");
    t.stop();
    tot += t.get_elapsed_time();

    t.start();
    reindexer.reorganize_mesh_top_simplices(*this,mesh,true);
    t.print_elapsed_time("   Exploit Top-cells Spatial Locality: ");
    t.stop();
    tot += t.get_elapsed_time();

    cerr<< kv << " Stellar tree [Build and Reordering] "<<tot<<endl;
}

template<class C,class T> void Stellar_Tree::add_vertices(Mesh<C, T> &mesh)
{
    for(int i=1;i<=mesh.get_vertices_num();i++)
    {
        if (mesh.get_domain().contains(mesh.get_vertex(i),mesh.get_domain().getMaxPoint()))
        {
            this->insert_vertex(this->root,mesh.get_domain(),0,mesh.get_vertex(i),i,mesh);
        }
    }
}

template<class C, class T> void Stellar_Tree::add_vertices(Cells_Soup &soup, Mesh<C,T> &mesh)
{
    // init mesh domain
    mesh.set_domain(soup.get_domain());
    int counter = 1;
    // init top array
    mesh.init_top_cells_array(soup.get_cells_types_num());

    // we loop over the top simplexes, and we get the pr-star structure by inserting the vertices
    // during the index generation we have to avoid the insertion of a vertex several time
    // inline we also generate the mesh indexed representation
    for(int d=0; d<soup.get_cells_types_num(); d++)
    {
        for(int i=1; i<=soup.get_cells_num(d); i++)
        {
            Explicit_Cell<C> &cell = soup.get_cell(d,i);
            ivect indexed_simplex;
            for(int v=0; v<cell.get_vertices_num(); v++)
            {
                if (mesh.get_domain().contains(cell.get_vertex(v),mesh.get_domain().getMaxPoint()))
                {
                    bool first_time = false;
                    this->insert_vertex_from_soup(this->root,mesh.get_domain(),0,cell.get_vertex(v),counter,indexed_simplex,first_time,mesh);
                    if(first_time)
                        counter++;
                }
            }
            mesh.add_top_cell(d,indexed_simplex);
        }
    }
}

template<class C,class T> void Stellar_Tree::add_top_cells(Mesh<C,T> &mesh)
{
    // init phase for the root node only in case it is a leaf node (i.e. index with just one node)
    if(this->root.is_leaf())
        this->root.init_sub_cells_vectors(mesh.get_top_cells_types());

    for(int w=0; w<mesh.get_top_cells_types(); w++)
    {
        for(int i=1;i<=mesh.get_top_cells_num(w);i++)
        {
            for(int j=0; j<mesh.get_top_cell(w,i).get_vertices_num(); j++)
            {
                if(mesh.get_domain().contains(mesh.get_vertex(abs(mesh.get_top_cell(w,i).TV(j))),mesh.get_domain().getMaxPoint()))
                {
                    this->insert_top_cell(this->root,mesh.get_top_cell(w,i),i,w);
                    break;
                }
            }
        }
    }
}

template<class C,class T> void Stellar_Tree::insert_vertex(Node_Stellar &n, Box<C> &domain, int level, Vertex<C> &v, int vertex_index, Mesh<C, T> &mesh)
{
    if (n.is_leaf())
    {
        n.add_vertex(vertex_index);
        if (n.is_full(this->kv))
            this->split(n,domain,level,mesh);
    }
    else
    {
        // instead of visiting all the sons until we get the one that contains the vertex
        // now we get in one pass the node that contain the vertex, executing range comparison for each coordinate axis

        int ind = 0;

        Box<C> sd;

        if(this->subdivision.get_son_index_and_domain(domain,level,v,mesh.get_domain().getMaxPoint(),ind,sd))
        {
            int son_level = level+1;
            /// if the son is not set.. then init it
            if(n.get_son(ind)==NULL)
            {
                Node_Stellar* s = new Node_Stellar();
                s->init_sub_cells_vectors(mesh.get_top_cells_types());
                n.set_son(s,ind);
            }
            this->insert_vertex(*n.get_son(ind),sd,son_level,v,vertex_index,mesh);
        }
        else
        {
            cerr<<"[addvertex] a vertex is not in any branch"<<endl;
            int a; cin>>a;
        }
    }
}

template<class C,class T> void Stellar_Tree::insert_vertex_from_soup(Node_Stellar &n, Box<C> &domain, int level, Vertex<C> &v, int vertex_index,
                                                                     ivect &indexed_cell, bool &first_time, Mesh<C, T> &mesh)
{
    if (n.is_leaf())
    {
        int ind = this->is_already_inserted(n,v,mesh);
        if(ind == -1)
        {
            ind = vertex_index;
            mesh.add_vertex(v);
            first_time = true;
            n.add_vertex(ind);
            /// in case of split we call the "usual" addVertex because we have all we need into the mesh variable
            if (n.is_full(this->kv))
            {
                this->split(n,domain,level,mesh);
            }
        }
        /// we have to add the index of the vertex to the indexed representation of the index
        indexed_cell.push_back(ind);
    }
    else
    {
        /// instead of visiting all the son until we get the one that contains the vertex
        /// now we get in one pass the node that contain the vertex, executing range comparison for each coordinate axis

        int ind = 0;
        Box<C> sd;
        if(this->subdivision.get_son_index_and_domain(domain,level,v,mesh.get_domain().getMaxPoint(),ind,sd))
        {
            int son_level = level+1;
            /// if the son is not set.. then init it
            if(n.get_son(ind)==NULL)
            {
                Node_Stellar* s = new Node_Stellar();
                s->init_sub_cells_vectors(mesh.get_top_cells_types());
                n.set_son(s,ind);
            }
            return this->insert_vertex_from_soup(*n.get_son(ind),sd,son_level,v,vertex_index,indexed_cell,first_time,mesh);
        }
        else
        {
            cerr<<"[addvertex] a vertex is not in any branch"<<endl;
            int a; cin>>a;
        }
    }
}

template<class T> void Stellar_Tree::insert_top_cell(Node_Stellar& n, T& t, int topsimplex_index, int vector_position)
{
    if(n.is_leaf())
    {
        // we do not need to execute this test. when we arrive in a leaf node, this encodes surely the top simplex
        n.add_top_cell(vector_position,topsimplex_index);
    }
    else
    {
        for(Node_Stellar::child_iterator it=n.begin(); it!=n.end(); ++it)
        {
            if(*it != NULL)
            {
                for(int j=0; j<t.get_vertices_num(); j++)
                {
                    if((*it)->indexes_vertex(abs(t.TV(j))))
                    {
                        this->insert_top_cell(**it,t,topsimplex_index,vector_position);
                        break;
                    }
                }
            }
        }
    }
}

template<class C,class T> void Stellar_Tree::split(Node_Stellar &n, Box<C> &domain, int level, Mesh<C,T> &mesh)
{
    n.init_sons(this->subdivision.get_son_number());

    for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        RunIterator const& v_id = itPair.first;
        this->insert_vertex(n,domain,level,mesh.get_vertex(*v_id),*v_id,mesh);
    }

    n.clear_lists();
}

template<class C,class T> int Stellar_Tree::is_already_inserted(Node_Stellar& n, Vertex<C> &v, Mesh<C, T> &mesh)
{
    int ind = -1;
    for(RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first)
    {
        Vertex<C> &vert = mesh.get_vertex(*itPair.first);
        if(v == vert)
        {
            ind = *itPair.first;
            break;
        }
    }
    return ind;
}
