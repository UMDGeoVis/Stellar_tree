set(INCLUDE_DIRS		
		io
		simplification				
		statistics
		stellar_decomposition
		stellar_tree
		topological_ds
		topological_ds_generator		
		topological_queries
		utilities
		validator
		vietoris_rips
	)
	
INCLUDE_DIRECTORIES( ${INCLUDE_DIRS} )

set(IO_SRC
		io/reader.h
		io/reader.cpp
		io/writer.h
		io/writer.cpp
	)

set(SIMPLIFICATION_SRC
		simplification/contraction_simplifier_top.h
		simplification/contraction_simplifier_top.cpp
        simplification/contraction_simplifier_top_length.cpp
		simplification/contraction_simplifier_weak.h
		simplification/contraction_simplifier_weak.cpp
        simplification/contraction_simplifier_weak_length.cpp
		simplification/edges_aux_structures.h
		simplification/simplification_aux_structures.h
		simplification/simplification_aux_structures.cpp
	)
    
set(STATISTICS_SRC		
		statistics/index_statistics.h
		statistics/statistics.h
		statistics/statistics_index.h
		statistics/statistics.cpp		
	)
	
set(STELLAR_DECOMPOSITION_SRC
		stellar_decomposition/node.h
		stellar_decomposition/node_stellar.h
		stellar_decomposition/node_stellar.cpp
		stellar_decomposition/node_v.h
		stellar_decomposition/reindexer.h
		stellar_decomposition/reindexer.cpp
		stellar_decomposition/run_iterator.h	
		stellar_decomposition/top_cell.h
		stellar_decomposition/top_cell.cpp
		stellar_decomposition/top_cp_cell.h
		stellar_decomposition/top_cp_cell.cpp
		stellar_decomposition/top_simplex.h
		stellar_decomposition/top_simplex.cpp		
	)
	
set(STELLAR_TREE_SRC
		stellar_tree/box.h
		stellar_tree/cells_soup.h
		stellar_tree/explicit_cell.h
		stellar_tree/mesh.h
		stellar_tree/mesh_updater.h
		stellar_tree/point.h
		stellar_tree/vertex.h	
		stellar_tree/spatial_subdivision.h
		stellar_tree/stellar_tree.h
		stellar_tree/stellar_tree_generation.h
		stellar_tree/stellar_tree_visit.h
		stellar_tree/stellar_tree_update.h
		stellar_tree/stellar_tree.cpp
	)
	
set(TOPOLOGICAL_DS_SRC
		topological_ds/adjacency_aux_structure.h
		topological_ds/adjacency_graph.h
		topological_ds/adjacency_graph.cpp
		topological_ds/halfedge.h
		topological_ds/iastar.h
		topological_ds/iastar.cpp
		topological_ds/ig.h
		topological_ds/ig_node.h
		topological_ds/links_aux_structures.h
		topological_ds/links_aux_structures.cpp
		topological_ds/skeleton.h
		topological_ds/skeleton.cpp
	)
	
set(TOPOLOGICAL_DS_GEN_SRC	
		topological_ds_generator/adjacency_graph_extractor.h
		topological_ds_generator/adjacency_graph_extractor.cpp
		topological_ds_generator/halfedge_generator.h
		topological_ds_generator/halfedge_generator.cpp		
		topological_ds_generator/iastar_generator.h
		topological_ds_generator/iastar_generator.cpp
		topological_ds_generator/iastar_statistics.h		
		topological_ds_generator/ig_generator.h
		topological_ds_generator/ig_generator.cpp		
		topological_ds_generator/ig_statistics.h		
		topological_ds_generator/skeleton_generator.h
		topological_ds_generator/skeleton_generator.cpp
    )

set(TOPOLOGICAL_QUERIES_SRC		
		topological_queries/topological_query_extractor.h
		topological_queries/topological_query_extractor.cpp
	)

set(UTILITIES_SRC 
		utilities/basic_wrappers.h
		utilities/cli_parameters.h
		utilities/color.h
		utilities/color.cpp
		utilities/combinatorial_utilities.h
		utilities/combinatorial_utilities.cpp
		utilities/container_utilities.h
		utilities/lru_cache.h
		utilities/string_management.h
		utilities/string_management.cpp
		utilities/timer.h
		utilities/timer.cpp
		utilities/usage.h
	)
	
set(VALIDATOR_SRC
		validator/connectedness_validator.h
		validator/connectedness_validator.cpp
		validator/connectedness_validator_adjs.h
		validator/connectedness_validator_adjs.cpp
		validator/connectedness_validator_skeleton.h
		validator/connectedness_validator_skeleton.cpp
		validator/connectedness_validator_stats.h
		validator/links_validator.h
		validator/links_validator.cpp
		validator/mesh_repairer.h
	)
	
set(VRIPS_SRC
		vietoris_rips/blockers_extractor.h
		vietoris_rips/blockers_extractor.cpp
		vietoris_rips/vietoris_rips_aux_structures.h
		vietoris_rips/vietoris_rips_generator.h
		vietoris_rips/vietoris_rips_generator.cpp
	)

set(STELLAR_ALL_SRC
		${BASIC_TYPES_SRC}
		${IO_SRC}
		${SIMPLIFICATION_SRC}
		${STATISTICS_SRC}
		${STELLAR_DECOMPOSITION_SRC}
		${STELLAR_TREE_SRC}
		${TOPOLOGICAL_DS_SRC}
		${TOPOLOGICAL_DS_GEN_SRC}
		${TOPOLOGICAL_QUERIES_SRC}
		${UTILITIES_SRC}
		${VALIDATOR_SRC}
		${VRIPS_SRC}
	)
	
add_library(stellar ${STELLAR_ALL_SRC})


