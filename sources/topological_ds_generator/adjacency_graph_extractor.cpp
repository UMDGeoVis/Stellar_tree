#include "adjacency_graph_extractor.h"

void Adjacency_Graph_Extractor::extract_lite_adjacency_graph_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, Adjacency_Graph &graph)
{
//    Adjacency_Graph_Extractor::extract_adjacency_graph_lite(n,mesh,graph);
    Adjacency_Graph_Extractor::extract_adjacency_graph_blockers_lite(n,mesh,graph);
}

void Adjacency_Graph_Extractor::extract_lite_adjacency_graph_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, Adjacency_Graph &graph)
{
    Adjacency_Graph_Extractor::extract_adjacency_graph_lite(n,mesh,graph);
//    Adjacency_Graph_Extractor::extract_adjacency_graph_blockers_lite(n,mesh,graph);
}

void Adjacency_Graph_Extractor::extract_full_adjacency_graph_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, Adjacency_Graph &graph)
{
    Adjacency_Graph_Extractor::extract_adjacency_graph(n,mesh,graph);
}

void Adjacency_Graph_Extractor::extract_full_adjacency_graph_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, Adjacency_Graph &graph)
{
    Adjacency_Graph_Extractor::extract_adjacency_graph(n,mesh,graph);
}

void Adjacency_Graph_Extractor::extract_local_lite_adjacency_graph_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, Local_Adjacency_Graph &graph)
{
    Adjacency_Graph_Extractor::extract_local_adjacency_graph_lite(n,mesh,graph);
}

void Adjacency_Graph_Extractor::extract_local_lite_adjacency_graph_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, Local_Adjacency_Graph &graph)
{
    Adjacency_Graph_Extractor::extract_local_adjacency_graph_lite(n,mesh,graph);
}

void Adjacency_Graph_Extractor::extract_local_full_adjacency_graph_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, Local_Adjacency_Graph &graph)
{
    Adjacency_Graph_Extractor::extract_local_adjacency_graph(n,mesh,graph);
}

void Adjacency_Graph_Extractor::extract_local_full_adjacency_graph_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, Local_Adjacency_Graph &graph)
{
    Adjacency_Graph_Extractor::extract_local_adjacency_graph(n,mesh,graph);
}

void Adjacency_Graph_Extractor::extract_local_1skeleton_graph_Simplicial_wrapper(Node_Stellar &n, Simplicial_Mesh &mesh, string mesh_name)
{
    Adjacency_Graph_Extractor::extract_local_1skeleton_graph(n,mesh,mesh_name);
}

void Adjacency_Graph_Extractor::extract_local_1skeleton_graph_CP_wrapper(Node_Stellar &n, CP_Mesh &mesh, string mesh_name)
{
    Adjacency_Graph_Extractor::extract_local_1skeleton_graph(n,mesh,mesh_name);
}
