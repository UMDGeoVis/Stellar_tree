#include "skeleton_generator.h"

namespace skeleton_generator
{

void extract_global_directed_1skeleton_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, skeleton_graph &graph)
{
    extract_global_directed_1skeleton(n,mesh,graph);
}

void extract_global_directed_1skeleton_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, skeleton_graph &graph)
{
    extract_global_directed_1skeleton(n,mesh,graph);
}

void extract_global_undirected_1skeleton_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, skeleton_graph &graph)
{
    extract_global_undirected_1skeleton(n,mesh,graph);
}

void extract_global_undirected_1skeleton_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, skeleton_graph &graph)
{
    extract_global_undirected_1skeleton(n,mesh,graph);
}


}
