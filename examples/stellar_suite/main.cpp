#include "main_utility_functions.h"
#include "main_simplicial_functions.h"
#include "main_cp_functions.h"

template<class C, class T> void load_index(Stellar_Tree& tree, Mesh<C,T> &mesh, Statistics& stats, Reindexer &reindexes, cli_parameters &cli);
template<class M> void exec_topo_queries(Stellar_Tree& tree, M &mesh, cli_parameters &cli);

template<class M> void exec_alpha_shape_generation(Stellar_Tree &tree, M &mesh, cli_parameters &cli);

// NOTA: old approach, the visit function is implemented into the Generator class
template<class C, class T> int extract_half_edges(Stellar_Tree& tree, Mesh<C,T>& mesh, cli_parameters &cli);
// NOTA: old approach, the visit function is implemented into the Generator class
template<class C, class T> int extract_incidence_graph(Stellar_Tree& tree, Mesh<C,T>& mesh, cli_parameters &cli);
// NOTA: old approach, the visit function is implemented into the Generator class
int extract_vietoris_rips_complex(Stellar_Tree& tree, Statistics& stats, Reindexer &reindexer, cli_parameters &cli);

//EXPERIMENTAL FEATURE..
// NOTA: old approach, the visit function is implemented into the Repair class
template<class T,class C, class Top> int repair_mesh(T& tree, Mesh<C,Top>& mesh, cli_parameters &cli);

int main(int argc, char** argv)
{
    if(argc == 1)
    {
        print_help();
        return 0;
    }

    cli_parameters cli = cli_parameters();

    if (read_arguments(argc, argv,cli) == -1)
    {
        print_usage();
        return (EXIT_FAILURE);
    }

    /// if kv=INT_MAX we do not want to have a hierarchical decomposition
    int kv = (cli.vertices_per_leaf == 0) ? INT_MAX : cli.vertices_per_leaf;

    Stellar_Tree tree = Stellar_Tree(/*cli.vertices_per_leaf*/kv);
    Reindexer reindexer = Reindexer(/*variables.reindex_type*/);
    Statistics stats = Statistics();

    if(cli.kind_of_operation == "vrG" || cli.kind_of_operation == "vrL"
            || cli.kind_of_operation == "vrH" || cli.kind_of_operation == "vrP")
        // if we want to build a Vietoris-Rips complex, then we have a dedicated main
    {
        extract_vietoris_rips_complex(tree,stats,reindexer,cli);
    }
    // SIMPLICIAL MESHES -> TRI - TETRA - SIERP - GMV - nDim
    else if(is_simplicial_mesh(cli))
    {
        cerr << "==Simplicial_Mesh==" <<endl;

        Simplicial_Mesh mesh = Simplicial_Mesh();

        load_mesh(mesh,cli);
        cerr << "RAM peak for Indexing the mesh: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
//        mesh.print_mesh_stats(cerr);
//        cerr << mesh.get_domain() << endl;
        load_index(tree,mesh,stats,reindexer,cli);
        cerr << "RAM peak for encoding the Stellar tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

        /// get the adjacency graph -- DISABLE AFTER DEBUGGING ///////////////////////////////
//        if
//        {
//            stringstream ss;
//            ss << get_file_name(cli.mesh_path) << "_kv_";
//            if(cli.vertices_per_leaf == 0)
//                ss << "inf";
//            else
//                ss << cli.vertices_per_leaf;
//    //            mesh.print_mesh_adj_graph(ss2.str());
////            /// prints the 1-skeleton indexed by each block
////            string f = ss.str();
////            tree.visit(Adjacency_Graph_Extractor::extract_local_1skeleton_graph_Simplicial_wrapper,tree.get_root(),mesh,f);

//            /// prints the 1-skeleton of the complex
//            cerr<<"[OUTPUT] compute the GLOBAL 1-skeleton within each leaf block of the INITIAL mesh"<<endl;
//            skeleton_graph n_graph = skeleton_graph(mesh.get_vertices_num());
//            tree.visit(skeleton_generator::extract_skeleton_Simplicial_wrapper,tree.get_root(),mesh,n_graph);
//            n_graph.write_gv_file(get_file_name(cli.mesh_path));

////            /// prints an adjacency graph for each leaf block
////            Local_Adjacency_Graph g(ss.str(),mesh);
////            tree.visit(Adjacency_Graph_Extractor::extract_local_adjacency_graph_Simplicial_wrapper,tree.get_root(),mesh,g);

//            /// prints an adjacency graph for the complete complex
//            cerr<<"[OUTPUT] computing the GLOBAL adjacency matrix of the INITIAL mesh"<<endl;
//            Adjacency_Graph gg(mesh);
//            tree.visit(Adjacency_Graph_Extractor::extract_lite_adjacency_graph_Simplicial_wrapper,tree.get_root(),mesh,gg);
//            gg.write_graph_file(get_file_name(cli.mesh_path),mesh);
//            gg.write_fixed_graph_file(get_file_name(cli.mesh_path),mesh);

//            /// TODO compute graphs statistics
//            cerr<<"===[TODO] COMPUTE GRAPHS STATISTICS==="<<endl;
//        }
        ///////////////////////////////////////////////////////////////////////////////////

        exec_alpha_shape_generation(tree, mesh, cli);

        // TOPOLOGICAL RELATION EXTRACTION
        exec_topo_queries(tree,mesh,cli);
        // TOPOLOGICAL DATA STRUCTURE GENERATION
        extract_adjacencies(tree,mesh,reindexer,cli);
        if(get_file_extension(cli.mesh_path) == "tri")
            extract_half_edges(tree,mesh,cli);
        extract_incidence_graph(tree,mesh,cli);
        // VALIDATOR
        validate_connectivity(tree,mesh,cli);
        if(get_file_extension(cli.mesh_path) == "ts" || get_file_extension(cli.mesh_path) == "tri")
        {
            validate_link_conditions(tree,mesh,cli); // the "classical" link conditions make sense only on low dimensional simplices
            /// to-do: add gmv support
        }
        //SIMPLIFICATION
        if(get_file_extension(cli.mesh_path) == "ts" ||
                get_file_extension(cli.mesh_path) == "tri" ||
                get_file_extension(cli.mesh_path) == "off")
        {
            simplify(tree,reindexer,mesh,cli,stats);
        }        
        // FOR DEBUGGING THE VRIPS GENERATOR
        if(cli.kind_of_operation == "vrDebug")
        {
            VietorisRips_Generator vr;
            vr.validity_check(tree,mesh);
        }
        //EXPERIMENTAL FEATURE
        repair_mesh(tree,mesh,cli);
        //
        extract_global_skeleton(tree,mesh,cli);        
    }
    // CP MESHES -> QUAD - HEXA
    else if(is_cp_mesh(cli))
    {
        cerr << "==CP_Mesh==" <<endl;

        CP_Mesh mesh = CP_Mesh();

        load_mesh(mesh,cli);
        cerr << "RAM peak for Indexing the mesh: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        load_index(tree,mesh,stats,reindexer,cli);
        cerr << "RAM peak for encoding the Stellar tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;
        // TOPOLOGICAL RELATION EXTRACTION
        exec_topo_queries(tree,mesh,cli);
        // TOPOLOGICAL DATA STRUCTURE GENERATION
        extract_adjacencies(tree,mesh,reindexer,cli);
        if(get_file_extension(cli.mesh_path) == "off")
            extract_half_edges(tree,mesh,cli);
        extract_incidence_graph(tree,mesh,cli);
        // VALIDATOR
        validate_connectivity(tree,mesh,cli);
        if(get_file_extension(cli.mesh_path) == "off") // JUST FOR QUAD-MESHES
        {
            validate_link_conditions(tree,mesh,cli);
            //EXPERIMENTAL FEATURE
            repair_mesh(tree,mesh,cli);
        }
    }
    // GET INDEXED REPRESENTATION FROM SOUP + APPS
    else if(get_file_extension(cli.mesh_path) == "soup")
    {
        Cells_Soup soup;
        Simplicial_Mesh mesh;

        if (!Reader::read_soup(soup, cli.mesh_path))
        {
            cerr << "Error Loading soup file. Execution Stopped." << endl;
            return (EXIT_FAILURE);
        }

        tree.set_subdivision(soup.get_cell(0,1).get_vertex(0).get_dimension());
        tree.build(soup,mesh,reindexer,(cli.kind_of_operation == "adjG"));

        if (cli.get_index_stats)
            stats.get_index_statistics(tree,mesh,true);
        if (cli.count_sub_simplices)
            stats.get_unique_cp_cells_counters(tree,mesh);

        exec_topo_queries(tree,mesh,cli);

        extract_adjacencies(tree,mesh,reindexer,cli);
        extract_incidence_graph(tree,mesh,cli);

        validate_connectivity(tree,mesh,cli);
    }

    return (EXIT_SUCCESS);
}

template<class C, class T> void load_index(Stellar_Tree &tree, Mesh<C,T> &mesh, Statistics &stats, Reindexer &reindexer, cli_parameters &cli)
{
    // init the division variable in tree
    tree.set_subdivision(mesh);

    stringstream out;
    string kv = (cli.vertices_per_leaf == 0) ? "inf" : to_string(cli.vertices_per_leaf);
    // NOTA: the verbose encoding is needed/used by the contraction and reduction simplifiers
    out << get_path_without_file_extension(cli.mesh_path) << "_v_" << kv << "_" << get_file_extension(cli.mesh_path);
    if(cli.verbose_encoding)
        out << "_verbose";
    if(cli.force3D)
        out << "_force3D";


    out << "_.stellar";
    cout << out.str() << endl;

    // first we check if exists a basic Stellar tree already built using the given kv
    if (!Reader::read_tree(tree.get_root(), tree.get_subdivision(), out.str()))
    {
        // if not we have to build the index with the given threshold
        tree.build(mesh,out.str(),reindexer,cli.save_original_indices);
    }
    else
    {
        // we need to explicitly exploit the spatial coherence
        stringstream base_info;
        base_info << kv << " Stellar-tree [Reordering Index] ";
        Timer t;
        t.start();
        // if we want to build the global IA/IA* data structure
        // we need the original position indexes of vertices and of top cells in order to save them consistently in the output file
        reindexer.reorganize_index_and_mesh(tree,mesh,cli.save_original_indices);
        t.stop();
        t.print_elapsed_time(base_info.str());
    }

    tree.init_leaves_list(tree.get_root()); // NOTA: enables the parallel OpenMP-based processing of the leaf block

    if (cli.get_index_stats)
        stats.get_index_statistics(tree,mesh,cli.debug_mode);
    if (cli.count_sub_simplices)
        stats.get_unique_cp_cells_counters(tree,mesh);
}

template<class M> void exec_topo_queries(Stellar_Tree &tree, M &mesh, cli_parameters &cli)
{
    if(cli.kind_of_operation == "vt")
        extract_VTop_relations(tree,mesh);
    if(cli.kind_of_operation == "vv")
        extract_VV_relations(tree,mesh);
    if(cli.kind_of_operation == "link")
        extract_links_relations(tree,mesh);
    /// pick one of these two strategies for computing p-faces ///
    if(cli.kind_of_operation == "pfaces")
    {
        extract_p_faces(tree,mesh);
    }
    if(cli.kind_of_operation == "pfaces2")
        extract_p_faces_v2(tree,mesh);
}

COORDBASETYPE calculate_top_simplex_radius(Top_Simplex &t, Simplicial_Mesh &mesh)
{
    int dim = t.get_vertices_num();
    if (dim < 2)
    {
        return -1;
    }

    if (dim == 2) // LINE
    {
        ivect &vertices = t.get_vertices_array();
        Vertex<COORDBASETYPE> &v1 = mesh.get_vertex(vertices[0]);
        Vertex<COORDBASETYPE> &v2 = mesh.get_vertex(vertices[1]);
        return v1.distance(v2) / 2.0;
    }
    else if (dim == 3) // TRIANGLE
    {
        ivect &vertices = t.get_vertices_array();
        Vertex<COORDBASETYPE> &v1 = mesh.get_vertex(vertices[0]);
        Vertex<COORDBASETYPE> &v2 = mesh.get_vertex(vertices[1]);
        Vertex<COORDBASETYPE> &v3 = mesh.get_vertex(vertices[2]);
        COORDBASETYPE d1 = v1.distance(v2);
        COORDBASETYPE d2 = v1.distance(v3);
        COORDBASETYPE d3 = v2.distance(v3);

        COORDBASETYPE circum_radius_sq = pow(d1 * d2 * d3, 2) / ((d1 + d2 + d3) * (d2 + d3 - d1) * (d3 + d1 - d2) * (d1 + d2 - d3));
        return sqrt(circum_radius_sq);
    }
    else if (dim == 4) // TETRAHEDRON
    {
        ivect &vertices = t.get_vertices_array();
        Vertex<COORDBASETYPE> &v1 = mesh.get_vertex(vertices[0]);
        Vertex<COORDBASETYPE> &v2 = mesh.get_vertex(vertices[1]);
        Vertex<COORDBASETYPE> &v3 = mesh.get_vertex(vertices[2]);
        Vertex<COORDBASETYPE> &v4 = mesh.get_vertex(vertices[3]);

        {
            // http://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter
            // https://github.com/czlowiekimadlo/Mag/blob/master/src/tetrahedron.cpp

            // double circumcenter[3];
            double *xi;
            double *eta;
            double *zeta;

            double denominator;
            // Use coordinates relative to point `a' of the tetrahedron.

            // ba = b - a
            double ba_x = v2.getC(0) - v1.getC(0);
            double ba_y = v2.getC(1) - v1.getC(1);
            double ba_z = v2.getC(2) - v1.getC(2);

            // ca = c - a
            double ca_x = v3.getC(0) - v1.getC(0);
            double ca_y = v3.getC(1) - v1.getC(1);
            double ca_z = v3.getC(2) - v1.getC(2);

            // da = d - a
            double da_x = v4.getC(0) - v1.getC(0);
            double da_y = v4.getC(1) - v1.getC(1);
            double da_z = v4.getC(2) - v1.getC(2);

            // Squares of lengths of the edges incident to `a'.
            double len_ba = ba_x * ba_x + ba_y * ba_y + ba_z * ba_z;
            double len_ca = ca_x * ca_x + ca_y * ca_y + ca_z * ca_z;
            double len_da = da_x * da_x + da_y * da_y + da_z * da_z;

            // Cross products of these edges.

            // c cross d
            double cross_cd_x = ca_y * da_z - da_y * ca_z;
            double cross_cd_y = ca_z * da_x - da_z * ca_x;
            double cross_cd_z = ca_x * da_y - da_x * ca_y;

            // d cross b
            double cross_db_x = da_y * ba_z - ba_y * da_z;
            double cross_db_y = da_z * ba_x - ba_z * da_x;
            double cross_db_z = da_x * ba_y - ba_x * da_y;

            // b cross c
            double cross_bc_x = ba_y * ca_z - ca_y * ba_z;
            double cross_bc_y = ba_z * ca_x - ca_z * ba_x;
            double cross_bc_z = ba_x * ca_y - ca_x * ba_y;

            // Calculate the denominator of the formula.
            denominator = 0.5 / (ba_x * cross_cd_x + ba_y * cross_cd_y + ba_z * cross_cd_z);

            // Calculate offset (from `a') of circumcenter.
            double circ_x = (len_ba * cross_cd_x + len_ca * cross_db_x + len_da * cross_bc_x) * denominator;
            double circ_y = (len_ba * cross_cd_y + len_ca * cross_db_y + len_da * cross_bc_y) * denominator;
            double circ_z = (len_ba * cross_cd_z + len_ca * cross_db_z + len_da * cross_bc_z) * denominator;

            // circumcenter[0] = circ_x + v1.getC(0);
            // circumcenter[1] = circ_y + v1.getC(1);
            // circumcenter[2] = circ_z + v1.getC(2);

            double circum_radius_sq = pow(circ_x, 2) + pow(circ_y, 2) + pow(circ_z, 2);
            if (circum_radius_sq <= 0) {
                cerr << "Wrong circumradius calculation for tetrahedron" << endl;
                exit(0);
            }
            return sqrt(circum_radius_sq);
        }
    }
    else
    {
        cerr << "Unsupported dimension for circumradius calculation: " << dim << std::endl;
        exit(0);
    }
}

double alpha_min = 1e-8;
// double alpha_value = 0.1; // 100; for terrain data
set<Top_Simplex> final_tetra, final_tri, final_edges;

void alpha_test(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool, int> &p)
{
    // determine node density
    COORDBASETYPE x_min = std::numeric_limits<COORDBASETYPE>::infinity(),
                  y_min = std::numeric_limits<COORDBASETYPE>::infinity(),
                  z_min = std::numeric_limits<COORDBASETYPE>::infinity(),
                  x_max = -std::numeric_limits<COORDBASETYPE>::infinity(),
                  y_max = -std::numeric_limits<COORDBASETYPE>::infinity(),
                  z_max = -std::numeric_limits<COORDBASETYPE>::infinity();
    int num_of_vertices = 0;
    for (RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first, ++num_of_vertices)
    {
        RunIterator const &v_id = itPair.first;
        Vertex<COORDBASETYPE> &v = mesh.get_vertex(*v_id);

        if (v.getC(0) < x_min)
        {
            x_min = v.getC(0);
        }
        if (v.getC(1) < y_min)
        {
            y_min = v.getC(1);
        }
        if (v.getC(2) < z_min)
        {
            z_min = v.getC(2);
        }

        if (v.getC(0) > x_max)
        {
            x_max = v.getC(0);
        }
        if (v.getC(1) > y_max)
        {
            y_max = v.getC(1);
        }
        if (v.getC(2) > z_max)
        {
            z_max = v.getC(2);
        }
    }
    COORDBASETYPE volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);
    double alpha_value = volume / num_of_vertices / 1.2 * 0.4 + 0.1;

    // alpha test
    set<Top_Simplex> decomposed_tri, decomposed_edge;
    if (n.get_num_top_cells_encoded() == 0)
    {
        return;
    }
    RunIterator t_id = n.t_array_begin_iterator(0);
    Top_Simplex &t = mesh.get_top_cell(0, *t_id);
    int dim = t.get_vertices_num();

    // 4: tetrahedron
    // 3: triangle
    for (int d = dim; d > 1; d--)
    {
        if (dim - d < n.get_num_top_cells_encoded())
        {
            for (RunIteratorPair itPair = n.make_t_array_iterator_pair(dim - d); itPair.first != itPair.second; ++itPair.first)
            {
                RunIterator const &t_id = itPair.first;
                Top_Simplex &t = mesh.get_top_cell(dim - d, *t_id);
                COORDBASETYPE circumradius = calculate_top_simplex_radius(t, mesh);

                if (circumradius < alpha_min)
                {
                    continue;
                }

                if (circumradius > alpha_value)
                {
                    for (int s = 0; s < t.get_sub_types_num(d - 2); s++)
                    {
                        ivect sub;
                        t.get_d_cell(sub, d - 2, s);

                        // for (auto x : sub) {
                        //     cerr << x << " ";
                        // }
                        // Top_Simplex tmp = Top_Simplex(sub);
                        // cerr <<circumradius << " "<< calculate_top_simplex_radius(tmp, mesh) << std::endl;
                        // exit(0);

                        if (d == 4)
                        {
                            decomposed_tri.insert(Top_Simplex(sub));
                        }
                        else if (d == 3)
                        {
                            decomposed_edge.insert(Top_Simplex(sub));
                        }
                    }
                }
                else
                {
                    if (d == 4)
                    {
                        final_tetra.insert(t);
                    }
                    else if (d == 3)
                    {
                        final_tri.insert(t);
                    }
                    else if (d == 2)
                    {
                        final_edges.insert(t);
                    }
                }
            }
        }

        if (d == 3)
        {
            for (Top_Simplex tri : decomposed_tri)
            {
                COORDBASETYPE circumradius = calculate_top_simplex_radius(tri, mesh);

                if (circumradius < alpha_min)
                {
                    continue;
                }
                if (circumradius > alpha_value)
                {
                    for (int s = 0; s < 3; s++)
                    {
                        ivect sub;
                        tri.get_d_cell(sub, 1, s);
                        decomposed_edge.insert(Top_Simplex(sub));
                    }
                }
                else
                {
                    final_tri.insert(tri);
                }
            }
        }
        else if (d == 2)
        {
            for (Top_Simplex edge : decomposed_edge)
            {
                COORDBASETYPE circumradius = calculate_top_simplex_radius(edge, mesh);

                if (circumradius < alpha_min)
                {
                    continue;
                }
                if (circumradius < alpha_value)
                {
                    final_edges.insert(edge);
                }
            }
        }
    }
}

// build connectivity graph for alpha shape test result
//    every top simplex is considered as a collection of edges carrying connectivity info for each vertex
std::unordered_map<int, int> connected_component_labels;
void dfs_on_connectivity_graph(int vid, std::unordered_map<int, iset> &connectivity, std::unordered_map<int, bool> &is_visited, int label)
{
    is_visited[vid] = true;
    connected_component_labels[vid] = label;
    for (int adj_vid : connectivity[vid])
    {
        if (!is_visited[adj_vid])
        {
            dfs_on_connectivity_graph(adj_vid, connectivity, is_visited, label);
        }
    }
}

void build_connectivity_graph()
{
    std::unordered_map<int, iset> connectivity;
    std::unordered_map<int, bool> is_visited;
    iset all_vertices;

    std::vector<Top_Simplex> alpha_shape;
    std::copy(final_tetra.begin(), final_tetra.end(), std::back_inserter(alpha_shape));
    std::copy(final_tri.begin(), final_tri.end(), std::back_inserter(alpha_shape));
    std::copy(final_edges.begin(), final_edges.end(), std::back_inserter(alpha_shape));

    for (Top_Simplex &simplex : alpha_shape)
    {
        for (int s = 0; s < simplex.get_sub_types_num(1); s++)
        {
            ivect sub;
            simplex.get_d_cell(sub, 1, s);
            connectivity[sub[0]].insert(sub[1]);
            connectivity[sub[1]].insert(sub[0]);
            is_visited[sub[0]] = false;
            is_visited[sub[1]] = false;
        }

        ivect &vertices = simplex.get_vertices_array();
        std::copy(vertices.begin(), vertices.end(), std::inserter(all_vertices, all_vertices.end()));
    }


    // for (auto x : connectivity) {
    //     std::cout << x.first << ": \t";
    //     for (auto y : x.second) {
    //         std::cout << y << " ";
    //     }
    //     std::cout << std::endl;
    // }

    int connected_component_label = 0;
    for (int vid : all_vertices)
    {
        if (!is_visited[vid])
        {
            dfs_on_connectivity_graph(vid, connectivity, is_visited, connected_component_label++);
        }
    }
    std::cout << "Number of connected component: " << connected_component_label << std::endl;
}

void save_alpha_shape(const char *ofname, Simplicial_Mesh &mesh)
{
    int pts_num = mesh.get_vertices_num();
    int cells_num = final_tetra.size() + final_tri.size() + final_edges.size();

    std::cout << "Save alpha shape: number of cells " << final_tetra.size() << " " << final_tri.size() << " " << final_edges.size() << std::endl;

    ofstream ofs(ofname);
    ofs << std::fixed << std::setprecision(8);
    ofs << "# vtk DataFile Version 2.0\n\n";
    ofs << "ASCII \n";
    ofs << "DATASET UNSTRUCTURED_GRID\n\n";
    ofs << "POINTS " << pts_num << " float\n";
    //  points
    for (int i = 1; i <= pts_num; ++i)
    {
        Vertex<COORDBASETYPE> &v = mesh.get_vertex(i);
        for (int j = 0; j < v.get_dimension(); j++)
        {
            ofs << v.getC(j) << " ";
        }
        ofs << std::endl;
    }
    //  cells
    ofs << "CELLS " << cells_num << " " << 5 * final_tetra.size() + 4 * final_tri.size() + 3 * final_edges.size() << std::endl;
    // tetras
    for (auto tt : final_tetra)
    {
        // vector<int> pids = tt.getVertices();
        ofs << tt.get_vertices_num() << " ";
        for (auto &pid : tt.get_vertices_array())
        {
            ofs << pid - 1 << " ";
        }
        ofs << std::endl;
    }
    // triangles
    for (auto tt : final_tri)
    {
        // vector<int> pids = tt.getVertices();
        ofs << tt.get_vertices_num() << " ";
        for (auto &pid : tt.get_vertices_array())
        {
            ofs << pid - 1 << " ";
        }
        ofs << std::endl;
    }
    // edges
    for (auto tt : final_edges)
    {
        // vector<int> pids = tt.getVertices();
        ofs << tt.get_vertices_num() << " ";
        for (auto &pid : tt.get_vertices_array())
        {
            ofs << pid - 1 << " ";
        }
        ofs << std::endl;
    }

    //  cell types
    ofs << "CELL_TYPES " << cells_num << endl;
    //  VTK_TETRA (=10)
    //  tetras
    for (const auto &tt : final_tetra)
    {
        ofs << "10 ";
    }
    //  VTK_TRIANGLE(=5)
    for (const auto &tt : final_tri)
    {
        ofs << "5 ";
    }
    //  VTK_LINE (=3)
    for (const auto &tt : final_edges)
    {
        ofs << "3 ";
    }

    // alpha_sq for each simplex on alpha shape
    // ofs << "CELL_DATA " << cells_num << endl;
    // ofs << "FIELD FieldData 1\n";
    // ofs << "a_sq 1 " << cells_num << " float\n";
    // //  tetras
    // for (Top_Simplex tt : final_tetra)
    // {
    //     ofs << calculate_top_simplex_radius(tt, mesh) << "\n";
    // }
    // //  triangles
    // for (Top_Simplex tt : final_tri)
    // {
    //     ofs << calculate_top_simplex_radius(tt, mesh) << "\n";
    // }
    // //  edges
    // for (Top_Simplex tt : final_edges)
    // {
    //     ofs << calculate_top_simplex_radius(tt, mesh) << "\n";
    // }

    // connected component labeling
    ofs << "POINT_DATA " << pts_num << std::endl;
    ofs << "FIELD FieldData 1\n"
        << "label 1 " << pts_num << " int\n";
    for (int vid = 1; vid <= pts_num; vid++)
    {
        if (connected_component_labels.count(vid) == 1)
        {
            ofs << connected_component_labels[vid] << " ";
        }
        else
        {
            ofs << "-1 ";
        }
    }
    ofs << std::endl;

    ofs.close();
}

/* Test functions start */
std::vector<std::vector<COORDBASETYPE>> bboxes;
void nodes_vis(Node_Stellar &n, Simplicial_Mesh &mesh, pair<bool, int> &p)
{
    int node_id = rand() % 1000;
    COORDBASETYPE x_min = std::numeric_limits<COORDBASETYPE>::infinity(),
                  y_min = std::numeric_limits<COORDBASETYPE>::infinity(),
                  z_min = std::numeric_limits<COORDBASETYPE>::infinity(),
                  x_max = -std::numeric_limits<COORDBASETYPE>::infinity(),
                  y_max = -std::numeric_limits<COORDBASETYPE>::infinity(),
                  z_max = -std::numeric_limits<COORDBASETYPE>::infinity();
    int num_of_vertices = 0;
    for (RunIteratorPair itPair = n.make_v_array_iterator_pair(); itPair.first != itPair.second; ++itPair.first, ++num_of_vertices)
    {
        RunIterator const &v_id = itPair.first;
        Vertex<COORDBASETYPE> &v = mesh.get_vertex(*v_id);
        v.setF(1, node_id);

        if (v.getC(0) < x_min)
        {
            x_min = v.getC(0);
        }
        if (v.getC(1) < y_min)
        {
            y_min = v.getC(1);
        }
        if (v.getC(2) < z_min)
        {
            z_min = v.getC(2);
        }

        if (v.getC(0) > x_max)
        {
            x_max = v.getC(0);
        }
        if (v.getC(1) > y_max)
        {
            y_max = v.getC(1);
        }
        if (v.getC(2) > z_max)
        {
            z_max = v.getC(2);
        }
    }
    COORDBASETYPE volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);
    bboxes.push_back({x_min, y_min, z_min, x_max, y_max, z_max, volume / num_of_vertices});
}

void save_node_vis(string mesh_path, Simplicial_Mesh &mesh)
{
    int pts_num = mesh.get_vertices_num();
    ofstream ofs((get_path_without_file_extension(mesh_path) + "_test.vtk").c_str());
    ofs << std::fixed << std::setprecision(8);
    ofs << "# vtk DataFile Version 2.0\n\n";
    ofs << "ASCII \n";
    ofs << "DATASET UNSTRUCTURED_GRID\n\n";
    ofs << "POINTS " << pts_num << " float\n";

    //  points
    for (int i = 1; i <= pts_num; ++i)
    {
        Vertex<COORDBASETYPE> &v = mesh.get_vertex(i);
        for (int j = 0; j < v.get_dimension(); j++)
        {
            ofs << v.getC(j) << " ";
        }
        ofs << std::endl;
    }

    ofs << "POINT_DATA " << pts_num << std::endl;
    ofs << "FIELD FieldData 1\n";
    ofs << "node_id 1 " << pts_num << " int\n";
    for (int i = 1; i <= pts_num; ++i)
    {
        Vertex<COORDBASETYPE> &v = mesh.get_vertex(i);
        ofs << int(v.getF(1)) << std::endl;
    }
    ofs.close();

    // save bboxes to another vtk file
    size_t num_of_bboxes = bboxes.size();
    ofstream ofs2((get_path_without_file_extension(mesh_path) + "_test_bbox.vtk").c_str());
    ofs2 << std::fixed << std::setprecision(8);
    ofs2 << "# vtk DataFile Version 2.0\n\n";
    ofs2 << "ASCII \n";
    ofs2 << "DATASET POLYDATA\n\n";
    ofs2 << "POINTS " << num_of_bboxes * 8 << " float\n";
    for (size_t i = 0; i < num_of_bboxes; i++)
    {
        std::vector<COORDBASETYPE> &bbox = bboxes[i];
        ofs2 << bbox[0] << " " << bbox[1] << " " << bbox[2] << "\n"
             << bbox[3] << " " << bbox[1] << " " << bbox[2] << "\n"
             << bbox[3] << " " << bbox[4] << " " << bbox[2] << "\n"
             << bbox[0] << " " << bbox[4] << " " << bbox[2] << "\n"
             << bbox[0] << " " << bbox[1] << " " << bbox[5] << "\n"
             << bbox[3] << " " << bbox[1] << " " << bbox[5] << "\n"
             << bbox[3] << " " << bbox[4] << " " << bbox[5] << "\n"
             << bbox[0] << " " << bbox[4] << " " << bbox[5] << "\n";
    }

    ofs2 << "POLYGONS " << num_of_bboxes * 6 << " " << num_of_bboxes * 30 << std::endl;
    for (size_t i = 0; i < num_of_bboxes; i++)
    {
        size_t offset = i * 8;
        ofs2 << "4 " << offset + 0 << " " << offset + 1 << " " << offset + 2 << " " << offset + 3 << "\n"
             << "4 " << offset + 4 << " " << offset + 5 << " " << offset + 6 << " " << offset + 7 << "\n"
             << "4 " << offset + 0 << " " << offset + 1 << " " << offset + 5 << " " << offset + 4 << "\n"
             << "4 " << offset + 2 << " " << offset + 3 << " " << offset + 7 << " " << offset + 6 << "\n"
             << "4 " << offset + 0 << " " << offset + 4 << " " << offset + 7 << " " << offset + 3 << "\n"
             << "4 " << offset + 1 << " " << offset + 2 << " " << offset + 6 << " " << offset + 5 << "\n";
    }
    ofs2 << std::endl;

    ofs2 << "CELL_DATA " << num_of_bboxes * 6 << "\n"
        << "FIELD FieldData 1\n"
        << "density 1 " << num_of_bboxes * 6 << " float" << std::endl;
    for (size_t i = 0; i < num_of_bboxes; i++)
    {
        ofs2 << bboxes[i][6] << " " << bboxes[i][6] << " " << bboxes[i][6] << " " << bboxes[i][6] << " " << bboxes[i][6] << " " << bboxes[i][6] << " ";
    }
    ofs2 << std::endl;
    ofs2.close();
}
/* Test functions end */

template <class M>
void exec_alpha_shape_generation(Stellar_Tree &tree, M &mesh, cli_parameters &cli)
{
    pair<bool, int> p;
    p.first = false;
    p.second = 0;
    Timer time;

    time.start();
    tree.visit(alpha_test, tree.get_root(), mesh, p);
    time.stop();
    time.print_elapsed_time("Alpha Shape Generation ");

    // time.start();
    // build_connectivity_graph();
    // time.stop();
    // time.print_elapsed_time("Connected component labeling ");

    save_alpha_shape((get_path_without_file_extension(cli.mesh_path) + "_as.vtk").c_str(), mesh);

    // time.start();
    // tree.visit(nodes_vis, tree.get_root(), mesh, p);
    // time.stop();
    // time.print_elapsed_time("nodes_vis function ");
    // save_node_vis(cli.mesh_path, mesh);
}

// NOTA: old approach, the visit function is implemented into the HalfEdge_Generator class
template<class C, class T> int extract_half_edges(Stellar_Tree& tree, Mesh<C,T>& mesh, cli_parameters &cli)
{
    HalfEdge_Generator he_gen;
    Timer time;

    if(cli.kind_of_operation == "heG")
    {
        he_gen = HalfEdge_Generator(mesh);
        time.start();
        he_gen.build_global_HalfEdge(tree.get_root(),mesh,cli.debug_mode);//,debug);
        time.stop();
        time.print_elapsed_time("Build global HalfEdge ");
        he_gen.print_stats();

    }
    else if(cli.kind_of_operation == "heL")
    {
        he_gen = HalfEdge_Generator();
        time.start();
        he_gen.build_HalfEdge(tree.get_root(),mesh);
        time.stop();
        time.print_elapsed_time("Build local HalfEdge ");
        he_gen.print_stats();

    }

    return (EXIT_SUCCESS);
}

// NOTA: old approach, the visit function is implemented into the IG_Generator class
template<class C, class T> int extract_incidence_graph(Stellar_Tree& tree, Mesh<C,T>& mesh, cli_parameters &variables)
{
    IG_Generator ig_generator(mesh.get_implicitly_and_explicitly_encoded_cells_num());
    Timer time;

    if(variables.kind_of_operation == "igG")
    {
        IG ig = IG(mesh.get_implicitly_and_explicitly_encoded_cells_num());
        time.start();
        ig_generator.build_IG(tree.get_root(),mesh,ig);
        time.stop();
        time.print_elapsed_time("Build global IG ");
        ig_generator.print_global_stats(ig);
        ig.clear_IG();
    }
    else if(variables.kind_of_operation == "igL")
    {
        time.start();
        ig_generator.build_IG(tree.get_root(),mesh,variables.debug_mode);
        time.stop();
        time.print_elapsed_time("Build local IG ");
        if(variables.debug_mode)
            ig_generator.print_local_stats();
    }

    return (EXIT_SUCCESS);
}

int extract_vietoris_rips_complex(Stellar_Tree& tree, Statistics& stats, Reindexer &reindexer, cli_parameters &cli)
{
    Simplicial_Mesh mesh = Simplicial_Mesh();
    VietorisRips_Generator vr_generator;

    Reader::read_vertices(mesh,cli.mesh_path);
    /// init the division variable in tree
    tree.set_subdivision(mesh);

    double diagonal = mesh.get_domain().get_diagonal();
//    for(int i=0; i< mesh.get_domain().getMinPoint().get_dimension(); i++)
//        cerr << mesh.get_domain().getMinPoint().getC(i) << " " << mesh.get_domain().getMaxPoint().getC(i) << " "
//             << i << "_side: " << mesh.get_domain().get_side(i) << endl;
    cerr << "[STATS] espilon: " << cli.eps << " diagonal: " << diagonal << " -> ratio: " << cli.eps / diagonal << endl;

    cerr << "RAM peak for loading the vertices: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    tree.generate_spatial_subdivision(mesh,reindexer);
    tree.init_leaves_list(tree.get_root());

    cerr << "RAM peak for generating the spatial subdivision + leaves list: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    ///for debug only
    if(cli.debug_mode)
    {
        mesh.get_vertices_distance_stats();
    }

    if(cli.eps > 0) /// we generate the Vietoris-Rips complex from a point cloud
    {
        if(cli.kind_of_operation == "vrG")
            // GLOBAL graph + GLOBAL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_GLOBAL(tree,mesh,cli.eps,cli.debug_mode);
        else if(cli.kind_of_operation == "vrL")
            // LOCAL graph + LOCAL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_LOCAL(tree,mesh,cli.eps,cli.debug_mode);
        else if(cli.kind_of_operation == "vrH")
            // GLOBAL graph + LOCAL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_HYBRID(tree,mesh,cli.eps,cli.debug_mode,get_path_without_file_extension(cli.mesh_path));
        else if(cli.kind_of_operation == "vrP")
            // GLOBAL graph + PARALLEL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_PARALLEL(tree,mesh,cli.eps);
    }
    else /// otherwise we extract the maximal cliques from a 1-skeleton
    {
        skeleton_graph graph;
        stringstream g_path;
        g_path << get_path_without_file_extension(cli.mesh_path) << ".edges";
//        cout<<"A"<<endl;
        Reader::read_skeleton(graph,mesh,g_path.str());
//        cout<<"A"<<endl;

        if(cli.kind_of_operation == "vrG")
            // GLOBAL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_GLOBAL(mesh,graph,cli.debug_mode);
        else if(cli.kind_of_operation == "vrL" || cli.kind_of_operation == "vrH")
            // LOCAL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex(tree,mesh,graph,cli.debug_mode);
        else if(cli.kind_of_operation == "vrP")
            // PARALLEL top-simplices extraction procedure
            vr_generator.extract_VietorisRips_complex_PARALLEL(tree,mesh,graph);
    }

    cerr << "RAM peak for generating a V-Rips complex: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    stringstream ss;
    ss<<"VRgen_eps_"<<cli.eps;

    if(!cli.debug_mode) //NOT IN DEBUG MODE
    {
        stringstream ss2;
        ss2 << get_path_without_file_extension(cli.mesh_path) << "_" << ss.str() << "_mesh_stats.txt";
        ofstream output(ss2.str());
        mesh.print_mesh_stats(output);
        Writer::write_mesh(get_path_without_file_extension(cli.mesh_path),ss.str(),mesh);
    }

//    if(cli.debug_mode)
//    {
//        cerr<<"computing the adjacency matrix of the GENERATED mesh"<<endl;
//        stringstream ss2;
//        ss2 << get_file_name(cli.mesh_path) << "_" << ss.str();
//        Adjacency_Graph g(mesh);
//        tree.visit(adj_graph_extractor::extract_adjacency_graph_Simplicial_wrapper,tree.get_root(),mesh,g);
//        g.write_graph_file(ss2.str(),mesh);
////        mesh.print_mesh_adj_graph(ss2.str());
//    }

    tree.add_and_resort_top_cells(mesh,reindexer);

    cerr << "RAM peak for add/resort the top cells in the tree: " << to_string(MemoryUsage().getValue_in_MB(false)) << " Mbs" << std::endl;

    if (cli.get_index_stats)
        stats.get_index_statistics(tree,mesh,false);
    if (cli.count_sub_simplices)
        stats.get_unique_cp_cells_counters(tree,mesh);

    if(cli.debug_mode) //IN DEBUG MODE
    {
        tuple<IAstar,bool,IAstar_stats> params;
        get<0>(params) = IAstar(mesh);
        get<1>(params) = cli.debug_mode;
        get<2>(params) = IAstar_stats();

        Timer time;
        time.start();
        tree.visit(IAstar_Generator::global_generation_Simplicial_wrapper,tree.get_root(),mesh,params);
        time.stop();
        time.print_elapsed_time("Build global IA-star ");

        stringstream stream2;
        stream2<<get_path_without_file_extension(cli.mesh_path)<<"_"<<ss.str()<<".iastar";
        get<0>(params).save_IAstar(stream2.str(),reindexer.get_original_vertices_ordering(),reindexer.get_original_tops_ordering());
        /// *** ///
    }

    return (EXIT_SUCCESS);
}

// NOTA: old approach, the visit function is implemented into the Mesh_Repairer class
template<class T,class C, class Top> int repair_mesh(T& tree, Mesh<C,Top>& mesh, cli_parameters &cli)
{
    if(cli.kind_of_operation == "fix")
    {
        cerr<<"[WARNING] The mesh repair is a limited and experimental feature.."<<endl;
        Mesh_Repairer<C,Top> repairer;
        repairer.repair_mesh(tree.get_root(),mesh);
        repairer.print_repair_output(get_path_without_file_extension(cli.mesh_path),mesh);
    }
    return (EXIT_SUCCESS);
}
