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

#ifndef MAIN_CP_FUNCTIONS_H
#define MAIN_CP_FUNCTIONS_H

#include "main_utility_functions.h"

void load_mesh(CP_Mesh& mesh, cli_parameters &cli);

int extract_VTop_relations(Stellar_Tree& tree, CP_Mesh &mesh);
int extract_VV_relations(Stellar_Tree& tree, CP_Mesh &mesh);
int extract_links_relations(Stellar_Tree& tree, CP_Mesh &mesh);
int extract_p_faces(Stellar_Tree& tree, CP_Mesh &mesh);
int extract_p_faces_v2(Stellar_Tree& tree, CP_Mesh &mesh);

int extract_adjacencies(Stellar_Tree &tree, CP_Mesh &mesh, Reindexer &reindexer, cli_parameters &cli);

int validate_connectivity(Stellar_Tree& tree, CP_Mesh& mesh, cli_parameters &cli);
int validate_link_conditions(Stellar_Tree& tree, CP_Mesh& mesh, cli_parameters &cli);

#endif // MAIN_CP_FUNCTIONS_H
