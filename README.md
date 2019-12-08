# Stellar Library #

The efficient representation and management of simplicial and cell complexes is an active research
topic in several fields, including geometric modeling, computer graphics, scientific visualization, and geographic data processing. Managing complexes in three dimensions and higher is not a simple task, since the topological data structures proposed in the literature do not scale when the size and the dimension of a complex become high. We propose the Stellar library as a topological C++
framework for performing efficient topological queries on simplicial and non-simplicial meshes. 
The Stellar library provides a scalable and compact representation that encodes the minimal information to locally reconstruct the topological connectivity of its indexed elements. This provides the flexibility to efficiently construct the optimal data structures to solve the task at hand using a fraction of the memory required for a corresponding topological data structure on the global mesh. The efficiency of the Stellar library increases with the execution of successive queries, as the construction costs of these runtime data structures are amortized over multiple accesses while processing each node. This project has been supported by the Italian Ministry of Education and Research under the PRIN 2009 program, and by the National Science Foundation under grant number IIS-1116747.

### Features ###

+ One spatial decomposition based on
    * point threshold (PR tree)
+ Two spatial decomposition
    * nD-quadtree
    * kD-tree
+ Topological queries
    * batched co-boundary queries
    * batched adjacency queries
    * batched link extraction for vertices
+ Topological data structures generation
    * halfedge data structure
    * IA* data structure
    * Incidence Graph
    * 1-skeleton
+ Complex Validation
    * 0-connectedness
    * d-connectedness
    * pseudo-manifoldness
    * link conditions (only on 2D and 3D complexes only)
+ Vietoris-Rips complexes generation
+ Homology preserving simplification based on two operators
    * weak link condition
    * top-based link condition
+ Soup to indexed mesh conversion

### How to compile ###

The library requires the [boost library](http://www.boost.org/) (for dynamic_bitset class), the [GMP library](https://gmplib.org/) and [cmake](https://cmake.org/) installed in your system.

Once in the root of the repository type from the command line
```
#!

cmake CMakeList.txt
```
and once configured
```
#!

make
```
This command generates a portable library file, located into *lib* folder, as well as some executables in `dist` folder.

The compilation has been test on linux and mac systems.

### Execute unit tests ###

Once compiled, it is possible to test the main functionalities running a script located into the `dist` folder.

From command line executing 
```
#!

sh run_tests.sh
```
checks the functionalities of the main implemented features.
The output files are saved into the `data` folder (where the input datasets are located).

It is possible to clean the output files running from command line (from `data` folder) the following command
```
#!

sh clean_up.sh
```

### Use the main library ###

In the `dist` folder there is another executable file named `stellar_suite` that contains the main library interface. For a complete list of the command line options refer the [wiki](https://bitbucket.org/riccardo_fellegara/stellar-library/wiki/Command%20line%20parameters) page.