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

#include "main_utility_functions.h"

bool checkParameters(cli_parameters &cli)
{

    if (cli.vertices_per_leaf == DEFAULT_X_PER_LEAF)
    {
        cout << "Error initializing vertices_per_leaf. Execution Stopped." << endl;
        print_usage();
        return false;
    }
    return true;
}

int read_arguments(int argc, char** argv, cli_parameters &cli)
{
    vector<string> tokens;
    string tok;
    cli.get_index_stats = false;
    cli.exe_name = argv[0];
    cli.exe_name = strip_path(cli.exe_name);

    for(int i=1; i<(argc); i++)
    {
        char* tag = argv[i];
        if(strcmp(tag, "-i") == 0)
        {
            cli.mesh_path = argv[i+1];
            i++;
        }
        else if(strcmp(tag, "-v") == 0)
        {
            cli.vertices_per_leaf = atoi(argv[i+1]);
            if (cli.vertices_per_leaf < 0) { /// NOTA we accept also 0 -> that we use as a key for the hierarchies with kv=infinite
                cerr << "Error: the limit of vertices per leaf must a positive integer "<< endl;
                return -1;
            }
            i++;
        }
        else if(strcmp(tag, "-s") == 0)
        {
            cli.get_index_stats = true;
        }
        else if(strcmp(tag, "-c") == 0)
        {
            cli.count_sub_simplices = true;
        }
        else if(strcmp(tag, "-q") == 0)
        {
            cli.kind_of_operation = argv[i+1];
            tokenize(cli.kind_of_operation,tokens,"_");

            if(tokens.size() == 0)
            {
                cerr << "[error] invalid number of parameters. 0 tokens error."<<endl;
                return -1;
            }
            else if((tokens[0] == "vt" || tokens[0] == "vv" || tokens[0] == "link" || tokens[0] == "pfaces" || tokens[0] == "pfaces2"
                     || tokens[0] == "adjL" || tokens[0] == "adjG" || tokens[0] == "adjLP" || tokens[0] == "adjGP"
                     || tokens[0] == "heL" || tokens[0] == "heG"
                     || tokens[0] == "igL" || tokens[0] == "igG"
                     || tokens[0] == "validate"
                     || tokens[0] == "fix"
                     || tokens[0] == "contractW" || tokens[0] == "contractT"
                     || tokens[0] == "skel"
                     )
                    && tokens.size() > 1)
            {
                cerr << "[error] too many parameters."<<endl;
                return -1;
            }
            else if((tokens[0] == "vrG" || tokens[0] == "vrL" || tokens[0] == "vrH" || tokens[0] == "vrP") && tokens.size() == 2)
            {
                cli.kind_of_operation = tokens[0];
                cli.eps = atof(tokens[1].c_str());
            }
            i++;

            /// if we reduce/contract our model, then we have to enable the verbose encoding in the mesh and in the leaf blocks
            if(cli.kind_of_operation == "contractW" || cli.kind_of_operation == "contractT")
            {
                cerr << "[NOTICE] Verbose Encoding Enabled." << endl;
                cli.verbose_encoding = true;
            }
        }
        else if(strcmp(tag, "-debug") == 0)
        {
            cerr<<"[WARNING] --> DEBUG MODE ACTIVATED <--"<<endl;
            cerr<<"          computing partial timings and the storage requirements of local structures"<<endl;
            cli.debug_mode = true;
        }
        else if(strcmp(tag, "-debugprint") == 0)
        {
            cerr<<"[WARNING] --> DEBUG PRINTS MODE ACTIVATED <--"<<endl;
            cli.debug_prints_mode = true;
        }
        else if(strcmp(tag, "-cp") == 0)
        {
            cli.is_cp_mesh = true;
        }
        else if(strcmp(tag, "-cache") == 0)
        {
            tok = argv[i+1];
            cli.cache_size = atoi(tok.c_str());
            i++;
        }
        else if(strcmp(tag, "-w") == 0)
        {
            cli.weights_file = argv[i+1];
            i++;
        }
        else if(strcmp(tag, "-blockersM") == 0)
        {
            cli.extract_blockers = "blockersM";
        }
        else if(strcmp(tag, "-blockersC") == 0)
        {
            cli.extract_blockers = "blockersC";
        }
//        else if(strcmp(tag, "-qem") == 0)
//        {
//            cli.qem = true;
//            cli.eps = atof(argv[i+1]);
//            i++;
//        }
        else if(strcmp(tag, "-len") == 0)
        {
            cli.len = true;
            cli.eps = atof(argv[i+1]);
            i++;
        }
        else if(strcmp(tag, "-force3D") == 0)
        {
            cli.force3D = true;
            cerr << "[NOTICE] Forcing a 3D Spatial Decomposition" << endl;
        }
    }    

    if(cli.kind_of_operation == "adjG" || cli.kind_of_operation == "adjGP" ||
            (cli.kind_of_operation == "contractT" && cli.debug_mode) ||
            (cli.kind_of_operation == "contractT" && cli.weights_file!=DEFAULT) ||
            (cli.kind_of_operation == "contractW" && cli.weights_file!=DEFAULT) )
        cli.save_original_indices = true;

    return 0;
}

bool is_simplicial_mesh(cli_parameters &cli)
{
    return (get_file_extension(cli.mesh_path) == "ts" || get_file_extension(cli.mesh_path) == "tri" ||
            get_file_extension(cli.mesh_path) == "sierp" || get_file_extension(cli.mesh_path)=="gmv" ||
            ((get_file_extension(cli.mesh_path) == "off" && !cli.is_cp_mesh)));
}

bool is_cp_mesh(cli_parameters &cli)
{
    return ((get_file_extension(cli.mesh_path) == "off" && cli.is_cp_mesh) || get_file_extension(cli.mesh_path) == "hex");
}

void print_usage()
{
    cerr<<"Wrong Usage. Run ./stellar_tree for detailed instructions."<<endl;
}

void print_help()
{
    FILE* fp;
    char path[1035];

    int cols;
    fp = popen("tput cols", "r");
    if(fp != NULL){
        fgets(path, sizeof(path)-1, fp);
        cols = atoi(path);
    }
    else{
        cols = 80;
    }

    printf(BOLD "\n  NAME:\n" RESET);
    printf("\t Stellar tree library.\n\n" RESET);

    printf(BOLD "  USAGE: \n" RESET);
    printf(BOLD "    ./stellar_tree -v [kv] -q [application] -s -c -cp -cache [cache_value] -i [mesh_file] -debug -debugtime\n\n" RESET);

    printf(BOLD "    -v [kv]\n" RESET);
    print_paragraph("'kv' is the vertex bucketing threshold per leaf block. If a stellar hierarchy with that kv, and that indexes the input mesh exists, then we load the hierarchy in memory, "
                    "otherwise, we build a brand new index with that threshold kv.", cols);
    printf(BOLD "    -q [application]\n" RESET);
    print_paragraph("'application' can be: vt - vv - link - adjL/adjG/adjGP/adjLP - heL/heG - igL/igG - vrG/vrL/vrH/vrP - validate - contractW/contractT - skel - fix", cols);
    printf(BOLD "        vt \n" RESET);
    print_paragraph(" computes the VTop relations, keeping track of the net timings and/or of the auxiliary data structure requirements (if debug mode is enabled).", cols);
    printf(BOLD "        vv \n" RESET);
    print_paragraph(" computes the VV relations, keeping track of the net timings and/or of the auxiliary data structure requirements (if debug mode is enabled).", cols);
    printf(BOLD "        link \n" RESET);
    print_paragraph(" computes the links of vertices and edges in the mesh, keeping track of the net timings and/or of the auxiliary data structure requirements (if debug mode is enabled).", cols);
    printf(BOLD "        adjL - adjG - adjGP - adjLP\n" RESET);
    print_paragraph(" generate the IA/IA* data structure. The available approaches are a global (G key), that explicitly encodes the topological data structure, "
                    " and a local approach (L key), that just generate a locally coherent data structure, that will be discarded after processing a leaf block."
                    " Both approaches can be executed in parallel (P key)", cols);
    printf(BOLD "        heL - heG\n" RESET);
    print_paragraph(" computes the half edge data structure, using a local (L key) or a global (G key) approach. (implemented only for triangle and quad meshes).", cols);
    printf(BOLD "        igL - igG\n" RESET);
    print_paragraph(" computes the Incidence Graph (IG), using a local (L key) or a global (G key) approach.", cols);
    printf(BOLD "        vrG_eps - vrL_eps - vrH_eps - vrP_eps\n" RESET);
    print_paragraph(" generate a Vietoris-Rips complex, with maximum distance of two points 'eps', from the points of the initial complex. The procedure is splitted in two main operation: "
                    "(1) the extraction of the neighborhood graph and (2) the maximal cliques computation. We provide four strategies, that generate the same output: "
                    "a fully local strategy (L key), that computes the neighborhood graph and the maximal cliques locally within each leaf block, "
                    "a global strategy (G key), that computes, using the hierarchy, and stores globally the neighborhood graph and then, computes the maximal cliques using a global procedure, "
                    "an hibrid strategy (H key), that computes, using the hierarchy, and stores globally the neighborhood graph and then, computes the maximal cliques locally within the hierarchy, and"
                    "a parallel strategy (P key), that computes, using the hierarchy, the neighborhood graph and stores globally the neighborhood graph and, then, computes in parallel the maximal cliques.", cols);
    printf(BOLD "        validate\n" RESET);
    print_paragraph(" validate the input model. It checks the 0-connectness, the d-connecteness, "
                    "the pseudo-manifoldness and the link conditions. The link conditions is available only on triangle,"
                    " quad and tetrahedral meshes. In debug mode in the program prints a verbose output if singularities occur.", cols);
    printf(BOLD "        contractW - contractT \n" RESET);
    print_paragraph(" executes the edge-contraction operator allover the indexed complex. We provide two strategies: "
                    "(1) one based on the weak-link-condition operator (contractW key) and (2) based on the top-based link condition operator (contractT key). "
                    "If debug mode is enabled the procedure outputs also extra details about times and memory peak.", cols);
    printf(BOLD "        skel \n" RESET);
    print_paragraph(" computes the 1-skeleton, keeping track of the net timings and/or of the auxiliary data structure requirements (if debug mode is enabled).", cols);
    printf(BOLD "        fix \n" RESET);
    print_paragraph(" execute an experimental repair operation, that checks the presence of duplicated vertices, duplicated/degenerate top simplices/cells."
                    " As this feature is still experimental is available only on triangle, tetrahedral and quad meshes.", cols);

    printf(BOLD "    -s\n" RESET);
    print_paragraph("computes the index statistics.", cols);

    printf(BOLD "    -c\n" RESET);
    print_paragraph("computes the number of sub-simplices in the imput mesh/soup.", cols);

    printf(BOLD "    -cp\n" RESET);
    print_paragraph("specify that the current input mesh is a constant polytope CP mesh. This is mandatory if the input file is an .off file representing a quad mesh.", cols);

    printf(BOLD "    -cache [cache_value]\n" RESET);
    print_paragraph("updates the maximum cache size. 'cache_value' is an integer containing this value. NOTA: the cache is not required by all applications. ", cols);

    printf(BOLD "    -i [mesh_file]\n" RESET);
    print_paragraph("reads the indexed mesh representation in 'mesh_file'. The supported formats are .ts (for tetrahedral meshes), .tri (for triangle meshes), "
                    ".gmv (for non-pure non-manifold meshes), .sierp (for sierpinski-based synthetic complexes), .off (for quad-meshes and arbitrary simplicial complexes), "
                    ".hex (for hexahedral meshes) and .soup for arbitrary soup of complexes.", cols);

    printf(BOLD "    -debug\n" RESET);
    print_paragraph("enables the debugging mode.", cols);

    printf(BOLD "    -debugtime\n" RESET);
    print_paragraph("enables the debugging timing mode.", cols);


    printf(BOLD "  EXAMPLE[1]: \n" RESET);
    printf("          ./stellar_tree -v 20 -q adjGP -s -i mesh.ts\n\n");
    print_paragraph("read the input file 'mesh.ts, containing a tetrahedral mesh. Then, it builds/load a stellar hierarchy (based on the octree subdivision) with kv = 20. "
                    "After it computes the index statistics (-s) and, finally, it extract the IA/IA* data structure following a global-parallel approach.", cols);

    printf(BOLD "  EXAMPLE[2]: \n" RESET);
    printf("          ./stellar_tree -v 50 -q vt -i mesh.tri\n\n");
    print_paragraph("read the input file 'mesh.tri, containing a triangle mesh. Then, it builds/load a stellar hierarchy (based on the octree subdivision) with kv = 50. "
                    "Finally, it extracts the local VTop relations for all the vertices, outputting the net timings.", cols);

    printf(BOLD "  IMPLEMENTATION:\n" RESET);
    printf("          Author: Riccardo Fellegara\n");
    printf("          Group: G3 Geometry and Graphics Group\n");
    printf("          Man-page Last Update: May 2016\n\n");

    printf(BOLD "  DESCRIPTION: \n" RESET);
    print_paragraph("The efficient representation and management of simplicial and cell complexes "
                    "is an active research topic in several fields, including geometric modeling, "
                    "computer graphics, scientific visualization, and geographic data processing. "
                    "Managing complexes in three dimensions and higher is not a simple task, since "
                    "the topological data structures proposed in the literature do not scale when "
                    "the size and the dimension of a complex become high. "
                    "We propose the Stellar library as a topological C++ framework "
                    "for performing efficient topological queries on simplicial and non-simplicial meshes. "
                    "The Stellar library provides a scalable and compact representation "
                    "that encodes the minimal information to locally reconstruct the topological "
                    "connectivity of its indexed elements. "
                    "This provides the flexibility to efficiently construct the optimal data "
                    "structures to solve the task at hand using a fraction of the memory "
                    "required for a corresponding topological data structure on the global mesh. "
                    "The efficiency of the Stellar library increases with the execution of "
                    "successive queries, as the construction costs of these runtime data structures "
                    "are amortized over multiple accesses while processing each node. "
                    "This project has been supported by the Italian Ministry of Education and "
                    "Research under the PRIN 2009 program, and by the National Science Foundation "
                    "under grant number IIS-1116747.", cols);
}

void print_paragraph(string stringa, int cols)
{
    if((int)stringa.size() < cols-20){
        printf("          %s\n\n",stringa.c_str());
    }
    else{
        float dim = (float)(stringa.size()/((float)cols-20));
        int dim_int = dim;
        if(dim > dim_int) dim_int++;
        for(int i=0; i<dim_int; i++){
            if(stringa.at((cols-20)*i) == ' ')
                printf("         %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
            else printf("          %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
        }
        printf("\n");
    }
}
