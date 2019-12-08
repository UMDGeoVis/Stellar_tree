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

#include "vietoris_rips_generator.h"

void VietorisRips_Generator::extract_VietorisRips_complex_GLOBAL(Stellar_Tree &tree, Simplicial_Mesh &mesh, double eps, bool debug)
{
    cerr<<"[VietorisRips Generation] GLOBAL approach"<<endl;

    Timer time;
    COORDBASETYPE tot_time;

    pair<skeleton_graph,Stellar_Tree*> p = make_pair(skeleton_graph(eps,mesh.get_vertices_num()),&tree);

    /// (1) extract skeleton ///
    time.start();
    tree.visit(VietorisRips_Generator::extract_neighborhood_graph,tree.get_root(),mesh.get_domain(),0,mesh,p);
    time.stop();
    if(!debug)
        time.print_elapsed_time("VR-Gen --> extracting skeleton: ");
    tot_time = time.get_elapsed_time();

//    if(debug)
    p.first.compute_skeleton_stats();

    /// (2) extract top simplices
    time.start();
    this->get_maximal_cliques(p.first,mesh,debug);
    time.stop();
    if(!debug)
        time.print_elapsed_time("VR-Gen --> extract top simplices: ");
    tot_time += time.get_elapsed_time();

    if(!debug)
        cerr<<"VR-Gen --> TOT: "<<tot_time<<endl;
}

void VietorisRips_Generator::extract_VietorisRips_complex_GLOBAL(Simplicial_Mesh &mesh, skeleton_graph &graph, bool debug)
{
    cerr<<"[VietorisRips Generation] GLOBAL approach"<<endl;

    Timer time;
    if(debug)
        graph.compute_skeleton_stats();

    /// (2) extract top simplices
    time.start();
    this->get_maximal_cliques(graph,mesh,debug);
    time.stop();
    if(!debug)
        time.print_elapsed_time("VR-Gen --> extract top simplices: ");
}

void VietorisRips_Generator::extract_VietorisRips_complex_LOCAL(Stellar_Tree &tree, Simplicial_Mesh &mesh, double eps, bool debug)
{
    cerr<<"[VietorisRips Generation] LOCAL approach"<<endl;

    Timer time;

    local_generation_parameters p = local_generation_parameters(eps,&tree,debug);

    time.start();
    tree.visit(VietorisRips_Generator::extract_VR_complex,tree.get_root(),mesh.get_domain(),0,mesh,p);
    time.stop();
    if(debug)
        p.print_local_stats();
    else
        time.print_elapsed_time("VR-Gen --> TOT: ");
}

void VietorisRips_Generator::extract_VietorisRips_complex_HYBRID(Stellar_Tree &tree, Simplicial_Mesh &mesh, double eps, bool debug, string mesh_name)
{
    cerr<<"[VietorisRips Generation] HYBRID approach (global skeleton and local cliques computation)"<<endl;

    Timer time;
    COORDBASETYPE tot_time;

    pair<skeleton_graph,Stellar_Tree*> p = make_pair(skeleton_graph(eps,mesh.get_vertices_num()),&tree);

    /// (1) extract skeleton ///
    time.start();
    tree.visit(VietorisRips_Generator::extract_neighborhood_graph,tree.get_root(),mesh.get_domain(),0,mesh,p);
    time.stop();
    time.print_elapsed_time("VR-Gen --> extracting skeleton: ");
    tot_time = time.get_elapsed_time();

    p.first.compute_skeleton_stats();
    if(debug)
    {
        stringstream ss;
        ss << mesh_name << "_VRgen_eps_" << eps << "_proximity_graph";
        p.first.save_skeleton_graph(ss.str());
    }

    /// (2) extract top simplices
    time.start();
    tree.visit(VietorisRips_Generator::get_local_maximal_cliques_v2,tree.get_root(),mesh,p.first);
    time.stop();
    time.print_elapsed_time("VR-Gen --> extract top simplices: ");
    tot_time += time.get_elapsed_time();

    cerr<<"VR-Gen --> TOT: "<<tot_time<<endl;
}

void VietorisRips_Generator::extract_VietorisRips_complex(Stellar_Tree &tree, Simplicial_Mesh &mesh, skeleton_graph &graph, bool debug)
{
    cerr<<"[VietorisRips Generation] HYBRID approach (global skeleton and local cliques computation)"<<endl;
    Timer time;
    if(debug)
        graph.compute_skeleton_stats();
    time.start();
    tree.visit(VietorisRips_Generator::get_local_maximal_cliques_v2,tree.get_root(),mesh,graph);
    time.stop();
    if(!debug)
        time.print_elapsed_time("VR-Gen --> extract top simplices: ");
}

void VietorisRips_Generator::extract_VietorisRips_complex_PARALLEL(Stellar_Tree &tree, Simplicial_Mesh &mesh, double eps)
{
    cerr<<"[VietorisRips Generation] PARALLEL approach (global skeleton and parallel cliques computation)"<<endl;

    Timer time;
    COORDBASETYPE tot_time;

    /// (0) extract the leaf blocks domain ///
    map<int,Box<COORDBASETYPE> > domains;
    time.start();
    tree.visit(VietorisRips_Generator::get_leaf_block_domain,tree.get_root(),mesh.get_domain(),0,mesh,domains);
    time.stop();
    time.print_elapsed_time("VR-Gen --> extracting leaves domains: ");

    tuple<skeleton_graph,Stellar_Tree*,map<int,Box<COORDBASETYPE> > > p = make_tuple(skeleton_graph(eps,mesh.get_vertices_num()),&tree,domains);

    /// (1) extract skeleton ///
    time.start();
    tree.parallel_visit(VietorisRips_Generator::extract_neighborhood_graph_PARALLEL,mesh,p);
    time.stop();
    time.print_elapsed_time("VR-Gen --> extracting skeleton: ");
    tot_time = time.get_elapsed_time();

    get<0>(p).compute_skeleton_stats();
    get<0>(p).check_skeleton_graph();

    /// (2) parallel extraction of top simplices
    time.start();
    tree.parallel_visit(VietorisRips_Generator::get_local_maximal_cliques_v2,mesh,get<0>(p));
    time.stop();
    time.print_elapsed_time("VR-Gen --> parallel extract top simplices: ");
    tot_time += time.get_elapsed_time();

    cerr<<"VR-Gen --> TOT: "<<tot_time<<endl;
}

void VietorisRips_Generator::extract_VietorisRips_complex_PARALLEL(Stellar_Tree &tree, Simplicial_Mesh &mesh, skeleton_graph &graph)
{
    cerr<<"[VietorisRips Generation] PARALLEL approach (global skeleton and parallel cliques computation)"<<endl;
    Timer time;
    time.start();
    tree.parallel_visit(VietorisRips_Generator::get_local_maximal_cliques_v2,mesh,graph);
    time.stop();
    time.print_elapsed_time("VR-Gen --> parallel extract top simplices: ");
}

void VietorisRips_Generator::extract_neighborhood_graph(Node_Stellar &n, Box<COORDBASETYPE> &n_domain, Simplicial_Mesh &mesh, pair<skeleton_graph,Stellar_Tree*> &p)
{
    for(int v_id=n.get_v_start(); v_id<n.get_v_end(); v_id++)
    {
        Vertex<COORDBASETYPE> &v1 = mesh.get_vertex(v_id);

        /// (1) search in the current leaf block
        VietorisRips_Generator::get_inner_nearest_vertices(v_id,v1,p.first,n,mesh,1);
    }

    /// (2) search outside the current leaf block
    VietorisRips_Generator::get_outer_nearest_vertices(n,n_domain,p.first,p.second->get_root(),mesh.get_domain(),0,mesh,p.second->get_subdivision(),1);
}

void VietorisRips_Generator::extract_neighborhood_graph_PARALLEL(Node_Stellar &n, Simplicial_Mesh &mesh, tuple<skeleton_graph, Stellar_Tree *, map<int, Box<COORDBASETYPE> > > &p)
{
    for(int v_id=n.get_v_start(); v_id<n.get_v_end(); v_id++)
    {
        Vertex<COORDBASETYPE> &v1 = mesh.get_vertex(v_id);

        /// (1) search in the current leaf block
        VietorisRips_Generator::get_inner_nearest_vertices(v_id,v1,get<0>(p),n,mesh,1);
    }

    /// (2) search outside the current leaf block
    VietorisRips_Generator::get_outer_nearest_vertices(n,get<2>(p)[n.get_v_start()],get<0>(p),get<1>(p)->get_root(),mesh.get_domain(),0,mesh,get<1>(p)->get_subdivision(),1);
}

void VietorisRips_Generator::get_leaf_block_domain(Node_Stellar &n, Box<COORDBASETYPE> &n_domain, Simplicial_Mesh &mesh, map<int,Box<COORDBASETYPE> > &domains)
{
    domains.insert(make_pair(n.get_v_start(),n_domain));
}

void VietorisRips_Generator::extract_VR_complex(Node_Stellar &n, Box<COORDBASETYPE> &n_domain, Simplicial_Mesh &mesh, local_generation_parameters &p)
{
    skeleton_graph n_graph = skeleton_graph(p.eps,n.get_v_end()-n.get_v_start());
    unordered_map<int,iset> outer_vertices;

    /// STEP 1: build the local neighborhood graph
    for(int v_id=n.get_v_start(); v_id<n.get_v_end(); v_id++)
    {
        Vertex<COORDBASETYPE> &v1 = mesh.get_vertex(v_id);

        /// (a) search in the current leaf block
        VietorisRips_Generator::get_inner_nearest_vertices(v_id,v1,n_graph,n,mesh,n.get_v_start());
    }

    /// (b) search outside the current leaf block
    VietorisRips_Generator::get_all_outer_nearest_vertices(n,n_domain,n_graph,p.tree->get_root(),mesh.get_domain(),0,mesh,p.tree->get_subdivision(),n.get_v_start(),outer_vertices);

    unordered_map<int,iset>::iterator next;
    /// (c) add the edges of the neighborhood graph from vertices completely outside
    for(unordered_map<int,iset>::iterator it=outer_vertices.begin(); it!=outer_vertices.end(); ++it)
    {
        Vertex<COORDBASETYPE> &v1 = mesh.get_vertex(it->first);

        next = it;
        ++next;

        for(unordered_map<int,iset>::iterator it2=next; it2!=outer_vertices.end(); ++it2)
        {
            Vertex<COORDBASETYPE> &v2 = mesh.get_vertex(it2->first);
            if(v1.distance(v2) < n_graph.get_epsilon())
            {
                outer_vertices[it->first].insert(it2->first);
                outer_vertices[it2->first].insert(it->first);
            }
        }
    }

    if(p.debug)
    {
        n_graph.set_skeleton_stats(p.maxV,p.maxE);
        if(p.maxOuterV < outer_vertices.size())
            p.maxOuterV = outer_vertices.size();
        long long outer_entities = 0;
        for(unordered_map<int,iset>::iterator it=outer_vertices.begin(); it!=outer_vertices.end(); ++it)
            outer_entities += it->second.size();
        if(p.maxOuterE < outer_entities)
            p.maxOuterE = outer_entities;
    }

    /// STEP 2: extract the local maximal cliques
    VietorisRips_Generator::get_local_maximal_cliques(n_graph,outer_vertices,mesh,n,p);
}

void VietorisRips_Generator::get_inner_nearest_vertices(int v_id, Vertex<COORDBASETYPE> &v, skeleton_graph &graph, Node_Stellar &n, Simplicial_Mesh &mesh, int offset)
{
    for(int v_id2=v_id+1; v_id2<n.get_v_end(); v_id2++)
    {
        Vertex<COORDBASETYPE> &v2 = mesh.get_vertex(v_id2);
        if(v.distance(v2) < graph.get_epsilon())
        {
            #pragma omp critical
            {
                graph.insert(v_id-offset,v_id2); /// -offset for the position index in the n_graph
                graph.insert(v_id2-offset,v_id);
            }
        }
    }
}

void VietorisRips_Generator::get_outer_nearest_vertices(Node_Stellar &pivot_n, Box<COORDBASETYPE> &pivot_domain, skeleton_graph &graph,
                                                        Node_Stellar &n, Box<COORDBASETYPE> &dom, int level,
                                                        Simplicial_Mesh &mesh, Spatial_Subdivision &subdiv, int offset)
{
    if (n.is_leaf())
    {
        if(n.is_after_vertex(pivot_n.get_v_start()) && pivot_domain.min_distance(dom) < graph.get_epsilon())
        {
            for(int v_id=pivot_n.get_v_start(); v_id<pivot_n.get_v_end(); v_id++)
            {
                Vertex<COORDBASETYPE> &v = mesh.get_vertex(v_id);

                for(int v_id2=n.get_v_start(); v_id2<n.get_v_end(); v_id2++)
                {
                    Vertex<COORDBASETYPE> &v2 = mesh.get_vertex(v_id2);
                    if(v.distance(v2) < graph.get_epsilon())
                    {
                        #pragma omp critical
                        {
                            graph.insert(v_id-offset,v_id2);
                            graph.insert(v_id2-offset,v_id);
                        }
                    }
                }
            }
        }
    }
    else
    {
        /// we only look "forward" to leaf blocksat distance < epsilon with respect to the pivot block
        if(n.indexes_vertex(pivot_n.get_v_start()) || (n.is_after_vertex(pivot_n.get_v_start()) && pivot_domain.min_distance(dom) < graph.get_epsilon()))
        {

            int sl = level + 1;
            for (int i = 0; i < subdiv.get_son_number(); i++)
            {
                if(n.get_son(i) != NULL)
                {
                    Box<COORDBASETYPE> son_dom = subdiv.compute_domain(dom,level,i);
                    VietorisRips_Generator::get_outer_nearest_vertices(pivot_n,pivot_domain,graph,*n.get_son(i),son_dom,sl,mesh,subdiv,offset);
                }
            }
        }
    }
}

void VietorisRips_Generator::get_all_outer_nearest_vertices(Node_Stellar &pivot_n, Box<COORDBASETYPE> &pivot_domain, skeleton_graph &graph, Node_Stellar &n, Box<COORDBASETYPE> &dom, int level,
                                                            Simplicial_Mesh &mesh, Spatial_Subdivision &subdiv, int offset, unordered_map<int,iset> &outer_vertices)
{
    if (n.is_leaf())
    {
        /// we avoid the pivot leaf
        if(!n.indexes_vertex(pivot_n.get_v_start()) && pivot_domain.min_distance(dom) < graph.get_epsilon())
        {
            for(int v_id=pivot_n.get_v_start(); v_id<pivot_n.get_v_end(); v_id++)
            {
                Vertex<COORDBASETYPE> &v = mesh.get_vertex(v_id);

                for(int v_id2=n.get_v_start(); v_id2<n.get_v_end(); v_id2++)
                {
                    Vertex<COORDBASETYPE> &v2 = mesh.get_vertex(v_id2);
                    if(v.distance(v2) < graph.get_epsilon())
                    {
                        #pragma omp critical
                        {
                            graph.insert(v_id-offset,v_id2);
                            outer_vertices[v_id2].insert(v_id);
                        }
                    }
                }
            }
        }
    }
    else
    {
        if(n.indexes_vertex(pivot_n.get_v_start()) || pivot_domain.min_distance(dom) < graph.get_epsilon())
        {
            int sl = level + 1;
            for (int i = 0; i < subdiv.get_son_number(); i++)
            {
                if(n.get_son(i) != NULL)
                {
                    Box<COORDBASETYPE> son_dom = subdiv.compute_domain(dom,level,i);
                    VietorisRips_Generator::get_all_outer_nearest_vertices(pivot_n,pivot_domain,graph,*n.get_son(i),son_dom,sl,mesh,subdiv,offset,outer_vertices);
                }
            }
        }
    }
}

void VietorisRips_Generator::get_maximal_cliques(skeleton_graph &graph, Simplicial_Mesh &mesh, bool debug)
{
    iset points;
    for(int i=1; i<=mesh.get_vertices_num(); i++)
    {
        points.insert(i);
    }

    map<int, list<ivect > > top_simplices_lists;
    iset setR, setX;
    unsigned maxX = 0;
    VietorisRips_Generator::get_maximal_cliques(graph,1,setR, points, setX, top_simplices_lists,debug,maxX);

    /// analyze the top simplexes
    int max_dim = top_simplices_lists.rbegin()->first;

    if(debug)
        cerr<<"max-setR: "<<max_dim<<" -- max-setP: "<<mesh.get_vertices_num()<<" -- max-setX: "<<maxX<<endl;

    mesh.init_top_cells_array(max_dim);
    mesh.add_types(max_dim,nDIM_SIMPL);

    for(top_simplices_map::iterator it=top_simplices_lists.begin(); it!=top_simplices_lists.end(); ++it)
    {
        mesh.add_top_cells(it->first-1,it->second); /// -1 to get the effective position
    }
}

void VietorisRips_Generator::get_maximal_cliques(skeleton_graph &graph, int offset, iset &setR, iset setP, iset &setX,
                                                 map<int, list<ivect > > &top_simplexes_local, bool debug, unsigned &maxX)
{
    if(setP.size() == 0 && setX.size() == 0)
    {
        //set setR as maximal
        if(setR.size() > 1)
        {
            ivect new_top(setR.begin(), setR.end());
            if(top_simplexes_local.find(new_top.size()-1) == top_simplexes_local.end()){
                top_simplexes_local[new_top.size()-1]= list<ivect >();
            }
            top_simplexes_local[new_top.size()-1].push_back(new_top);
        }
    }
    else if(setP.size() > 0)
    {
        //choose here the pivot vertex
        int pivot=*(setP.begin());
        for(iset_iter it=setP.begin(); it!=setP.end(); ++it)
        {
            if(graph.adjacent_vertices_num(*it-offset) > graph.adjacent_vertices_num(pivot-offset))
                pivot=*it;
        }

        iset diff;
        difference_of_sets(setP,graph.get_graph(pivot-offset),diff); // diff = P \ N(pivot)

        for(iset_iter it=diff.begin(); it!=diff.end(); ++it)
        {
            unsigned ui=*it;

            // setP = setP - {ui}
            setP.erase(ui);
            iset adjui = graph.get_graph(ui-offset);
            iset newR = setR;
            newR.insert(ui);
            // newP = N(ui) [intersection] setP
            iset newP;
            intersect_sets(adjui,setP,newP);
            // newX = N(ui) [intersection] setX
            iset newX;
            intersect_sets(adjui,setX,newX);

            get_maximal_cliques(graph,offset,newR, newP, newX, top_simplexes_local,debug,maxX);

            // setX = setX + {ui}
            setX.insert(ui);

            if(debug)
                if(maxX < setX.size())
                    maxX = setX.size();
        }
    }
}

void VietorisRips_Generator::get_local_maximal_cliques(skeleton_graph &graph, unordered_map<int, iset> &outer_vertices,
                                                       Simplicial_Mesh &mesh, Node_Stellar &n, local_generation_parameters &p)
{
    iset points;
    for(int v_id=n.get_v_start(); v_id<n.get_v_end(); v_id++)
        points.insert(v_id);
    for(unordered_map<int, iset>::iterator it=outer_vertices.begin(); it!=outer_vertices.end(); ++it)
        points.insert(it->first);

    if(p.debug)
    {
        if(p.maxP < (int)points.size())
            p.maxP = points.size();
    }

    top_simplices_map tops;
    iset setR, setX;
    VietorisRips_Generator::get_local_maximal_cliques_rec(graph,outer_vertices,setR, points, setX, tops, n, n.get_v_start(),p);

    /// analyze the top simplexes
    if(tops.size() > 0)
    {
        int max_dim = tops.rbegin()->first;

        if(p.debug)
        {
            if(p.maxR < max_dim)
                p.maxR = max_dim;
        }

        if(mesh.get_top_cells_types() == 0)
        {
            mesh.init_top_cells_array(max_dim);
            mesh.add_types(max_dim,nDIM_SIMPL);
        }
        else if(max_dim > mesh.get_top_cells_types())
        {
            int diff = max_dim - mesh.get_top_cells_types();
            mesh.add_top_cells_arrays(diff);
            mesh.add_types(diff,nDIM_SIMPL);
        }

        for(top_simplices_map::iterator it=tops.begin(); it!=tops.end(); ++it)
        {
            mesh.add_top_cells(it->first-1,it->second); /// -1 to get the effective position
        }
    }
}

void VietorisRips_Generator::get_local_maximal_cliques_rec(skeleton_graph &graph, unordered_map<int,iset> &outer_vertices, iset &setR, iset setP, iset &setX,
                                                           map<int, list<ivect > > &top_simplexes_local, Node_Stellar &n, int offset, local_generation_parameters &p)
{
    if(setP.size() == 0 && setX.size() == 0)
    {
        //set setR as maximal only if it is at least an edge
        if(setR.size() > 1)
        {
            /// NOTA: encoding the graph globally lead to the extraction of the same maximal clique
            /// several times, thus we add to the top simplices list the first time we find it
            /// (i.e. when we process the leaf block indexing its minimum vertex)
            ///
            if(n.indexes_vertex(*setR.begin()))
            {
                ivect new_top(setR.begin(), setR.end());
                if(top_simplexes_local.find(new_top.size()-1) == top_simplexes_local.end()){
                    top_simplexes_local[new_top.size()-1]= list<ivect >();
                }
                top_simplexes_local[new_top.size()-1].push_back(new_top);
            }
        }
    }
    else if(setP.size() > 0)
    {

        // chose here the pivot vertex for P
        int pivot=*(setP.begin());
        unsigned pivot_size;
        if(n.indexes_vertex(pivot))
            pivot_size = graph.adjacent_vertices_num(pivot-offset);
        else
            pivot_size = outer_vertices[pivot].size();

        for(iset_iter it=setP.begin(); it!=setP.end(); ++it)
        {
            if(n.indexes_vertex(*it))
            {
                if(graph.adjacent_vertices_num(*it-offset) > pivot_size)
                {
                    pivot = *it;
                    pivot_size = graph.adjacent_vertices_num(*it-offset);
                }
            }
            else
            {
                if(outer_vertices[*it].size() > pivot_size)
                {
                    pivot=*it;
                    pivot_size = outer_vertices[*it].size();
                }
            }
        }

        if(pivot_size == 0) // nothing to do
            return;

        iset diff;
         // diff = P \ N(pivot)
        if(n.indexes_vertex(pivot))
            difference_of_sets(setP,graph.get_graph(pivot-offset),diff);
        else
            difference_of_sets(setP,outer_vertices[pivot],diff);

        for(iset_iter it=diff.begin(); it!=diff.end(); ++it)
        {
            unsigned ui=*it;

            // setP = setP - {ui}
            setP.erase(ui);

            iset newR, newP, newX, adjui;
            // newR = setR + {ui}
            newR = setR;
            newR.insert(ui);

            /// NOTA: we have to discretize between the local indexed vertices
            /// and those (external) that only belong to the neighborhood of a local vertex
            if(n.indexes_vertex(ui))
                adjui = graph.get_graph(ui-offset);
            else
                adjui = outer_vertices[ui];

            // newP = N(ui) [intersection] setP
            intersect_sets(adjui,setP,newP);
            // newX = N(ui) [intersection] setX
            intersect_sets(adjui,setX,newX);

            VietorisRips_Generator::get_local_maximal_cliques_rec(graph,outer_vertices,newR, newP, newX, top_simplexes_local,n,offset,p);

            // setX = setX + {ui}
            setX.insert(ui);

            if(p.debug)
            {
                if(p.maxX < (int)setX.size())
                    p.maxX = setX.size();
            }
        }
    }
}

void VietorisRips_Generator::get_local_maximal_cliques_v2(Node_Stellar &n, Simplicial_Mesh &mesh, skeleton_graph &graph)
{
    int offset = 1; /// the graph is stored globally
    iset points;
    for(int v_id=n.get_v_start(); v_id<n.get_v_end(); v_id++)
    {
        points.insert(v_id);
        points.insert(graph.begin(v_id-offset),graph.end(v_id-offset));
    }

    top_simplices_map tops;
    iset setR, setX;
    VietorisRips_Generator::get_local_maximal_cliques_v2_rec(graph, setR, points, setX, tops, n, offset);

    /// analyze the top simplexes
    if(tops.size() > 0)
    {
        int max_dim = tops.rbegin()->first;

        if(mesh.get_top_cells_types() == 0)
        {
#pragma omp critical
            {
                mesh.init_top_cells_array(max_dim);
                mesh.add_types(max_dim,nDIM_SIMPL);
            }
        }
        else if(max_dim > mesh.get_top_cells_types())
        {
            int diff = max_dim - mesh.get_top_cells_types();
#pragma omp critical
            {
                mesh.add_top_cells_arrays(diff);
                mesh.add_types(diff,nDIM_SIMPL);
            }
        }

#pragma omp critical
        {
            for(top_simplices_map::iterator it=tops.begin(); it!=tops.end(); ++it)
            {
                mesh.add_top_cells(it->first-1,it->second); /// -1 to get the effective position
            }
        }
    }
}

void VietorisRips_Generator::get_local_maximal_cliques_v2_rec(skeleton_graph &graph, iset &setR, iset setP, iset &setX,
                                                           map<int, list<ivect > > &top_simplexes_local, Node_Stellar &n, int offset)
{
    if(setP.size() == 0 && setX.size() == 0)
    {
        //set setR as maximal only if it is at least an edge
        if(setR.size() > 1)
        {
            /// NOTA: encoding the graph globally lead to the extraction of the same maximal clique
            /// several times, thus we add to the top simplices list the first time we find it
            /// (i.e. when we process the leaf block indexing its minimum vertex)
            ///
            if(n.indexes_vertex(*setR.begin()))
            {
                ivect new_top(setR.begin(), setR.end());
                if(top_simplexes_local.find(new_top.size()-1) == top_simplexes_local.end()){
                    top_simplexes_local[new_top.size()-1]= list<ivect >();
                }
                top_simplexes_local[new_top.size()-1].push_back(new_top);
            }
        }
    }
    else if(setP.size() > 0)
    {
        // chose here the pivot vertex for P
        int pivot=*(setP.begin());
        unsigned pivot_size = graph.adjacent_vertices_num(pivot-offset);

        for(iset_iter it=setP.begin(); it!=setP.end(); ++it)
        {
            if(graph.adjacent_vertices_num(*it-offset) > pivot_size)
            {
                pivot = *it;
                pivot_size = graph.adjacent_vertices_num(*it-offset);
            }
        }

        if(pivot_size == 0) // nothing to do
            return;

        iset diff;
         // diff = P \ N(pivot)
        difference_of_sets(setP,graph.get_graph(pivot-offset),diff);

        for(iset_iter it=diff.begin(); it!=diff.end(); ++it)
        {
            unsigned ui=*it;

            // setP = setP - {ui}
            setP.erase(ui);

            iset newR, newP, newX, adjui;
            // newR = setR + {ui}
            newR = setR;
            newR.insert(ui);

            adjui = graph.get_graph(ui-offset);

            // newP = N(ui) [intersection] setP
            intersect_sets(adjui,setP,newP);
            // newX = N(ui) [intersection] setX
            intersect_sets(adjui,setX,newX);

            VietorisRips_Generator::get_local_maximal_cliques_v2_rec(graph, newR, newP, newX, top_simplexes_local, n, offset);

            // setX = setX + {ui}
            setX.insert(ui);
        }
    }
}

void VietorisRips_Generator::validity_check(Stellar_Tree &tree, Simplicial_Mesh &mesh)
{
    /// here we check if the VRips generator algorithm is implemented correctly
    /// starting from a full simplicial mesh and a stellar tree built on top of it
    /// we first extract the 1-skeleton
    /// and then we extract the maximal cliques on it
    /// if the number of maximal cliques is equal to the number of top simplices then the algorithm is correct

    cout<<"[validity_check]"<<endl;

    skeleton_graph n_graph = skeleton_graph(mesh.get_vertices_num());

    /// (1) extract skeleton ///
    tree.visit(VietorisRips_Generator::extract_skeleton,tree.get_root(),mesh,n_graph);

    /// (2) extract top simplices
    this->check_maximal_cliques(n_graph,mesh);
}

void VietorisRips_Generator::extract_skeleton(Node_Stellar &n, Simplicial_Mesh &mesh, skeleton_graph &graph)
{
    for(int d=0; d<n.get_num_top_cells_encoded();d++)
    {
        for(RunIteratorPair itPair = n.make_t_array_iterator_pair(d); itPair.first != itPair.second; ++itPair.first)
        {
            RunIterator const& t_id = itPair.first;

            if(!mesh.is_top_cell_removed(d,*t_id))
            {
                Top_Simplex& t = mesh.get_top_cell(d,*t_id);
                ivect f;
                for(int v=0; v<t.get_vertices_num(); v++)
                {
                    int v_id = abs(t.TV(v));

                    if(n.indexes_vertex(v_id))
                    {
                        t.TF(f,v);
                        graph.insert(v_id-1,f);
                    }
                }
            }
        }
    }
}

void VietorisRips_Generator::check_maximal_cliques(skeleton_graph &graph, Simplicial_Mesh &mesh)
{
    iset points;
    for(int i=1; i<=mesh.get_vertices_num(); i++)
    {
        points.insert(i);
    }

    top_simplices_map tops;
    iset setR, setX;
    unsigned trash = 0;
    get_maximal_cliques(graph,1,setR, points, setX, tops,false,trash);

    for(top_simplices_map::iterator it=tops.begin(); it!=tops.end(); ++it)
    {
        cout<<it->first<<"] --> "<<it->second.size()<<endl;
    }

    mesh.print_mesh_stats(cout);
}

