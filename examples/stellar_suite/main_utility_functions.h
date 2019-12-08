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

#ifndef MAIN_UTILITY_FUNCTIONS_H
#define MAIN_UTILITY_FUNCTIONS_H

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <sstream>
#include <boost/tuple/tuple.hpp>

/// STELLAR CORE LIBRARY ///
#include "io/reader.h"
#include "io/writer.h"
#include "stellar_decomposition/node_stellar.h"
#include "stellar_decomposition/reindexer.h"
#include "stellar_tree/spatial_subdivision.h"
#include "stellar_tree/stellar_tree.h"
#include "topological_queries/topological_query_extractor.h"
#include "statistics/statistics.h"
#include "utilities/cli_parameters.h"
#include "utilities/string_management.h"
#include "utilities/timer.h"
#include "utilities/usage.h"

/// STELLAR APPS LIBRARY ///
#include "topological_ds_generator/iastar_generator.h"
#include "topological_ds_generator/halfedge_generator.h"
#include "topological_ds_generator/ig_generator.h"
#include "topological_ds_generator/skeleton_generator.h"
#include "topological_ds_generator/adjacency_graph_extractor.h"

#include "vietoris_rips/vietoris_rips_generator.h"
#include "vietoris_rips/blockers_extractor.h"

#include "validator/connectedness_validator.h"
#include "validator/links_validator.h"
#include "validator/mesh_repairer.h"

#include "simplification/contraction_simplifier_weak.h"
#include "simplification/contraction_simplifier_top.h"

using namespace std;
using namespace string_management;


#define BOLD  "\033[1m\033[33m" //for dark background shell
#define RESET   "\033[0m"

int read_arguments(int argc, char** argv, cli_parameters &cli);
void print_usage();
bool checkParameters(cli_parameters &cli);
void print_help();
void print_paragraph(string stringa, int cols);
bool is_simplicial_mesh(cli_parameters &cli);
bool is_cp_mesh(cli_parameters &cli);

#endif // MAIN_UTILITY_FUNCTIONS_H
