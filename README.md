# Stellar Library #

The efficient representation and management of simplicial and cell complexes is an active research
topic in several fields, including geometric modeling, computer graphics, scientific visualization, and geographic data processing. Managing complexes in three dimensions and higher is not a simple task, since the topological data structures proposed in the literature do not scale when the size and the dimension of a complex become high. We propose the Stellar library as a topological C++ framework for performing efficient topological queries on simplicial and non-simplicial meshes. 

The Stellar library provides a scalable and compact representation that encodes the minimal information to locally reconstruct the topological connectivity of its indexed elements. This provides the flexibility to efficiently construct the optimal data structures to solve the task at hand using a fraction of the memory required for a corresponding topological data structure on the global mesh. The efficiency of the Stellar library increases with the execution of successive queries, as the construction costs of these runtime data structures are amortized over multiple accesses while processing each node. 

### Acknowledgments ###

This work has been partially supported by the US National Science Foundation under grant number IIS-1910766 and by the University of Maryland under the 2017-2018 BSOS Dean Research Initiative Program. It has also been performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.

### Publications ###

**Main paper (describing the Stellar decomposition idea and Stellar tree library)**

- **The Stellar tree: a Compact Representation for Simplicial Complexes and Beyond**  
R. Fellegara, K. Weiss, and L. De Floriani  
*Arxiv e-prints, 2019 (v2)* - [doi](https://arxiv.org/abs/1707.02211)

**Papers using Stellar tree library (or some older version of it)**

- **Efficient Homology‐Preserving Simplification of High‐Dimensional Simplicial Shapes**  
R. Fellegara, F. Iuricich, L. De Floriani, and U. Fugacci  
*Computer Graphics Forum, 39: 244-259, 2019* - [doi](http://dx.doi.org/10.1111/cgf.13764)
- **An efficient approach for verifying manifold properties of simplicial complexes**  
R. Fellegara, K. Weiss, and L. De Floriani  
*Proceedings of the 25th International Meshing Roundtable, 2016* - [doi](http://imr.sandia.gov/papers/abstracts/Fe830.html)
- **Analysis of Geolocalized Social Networks Based on Simplicial Complexes**  
R. Fellegara, U. Fugacci, F. Iuricich, L. De Floriani  
*Proceedings of the 9th ACM SIGSPATIAL Workshop on Location-based Social Networks, 5:1--5:8, 2016* - [doi](https://dl.acm.org/doi/10.1145/3021304.3021309)
- **Efficient Computation and Simplification of Discrete Morse Decompositions on Triangulated Terrains**  
R. Fellegara, F. luricich, L. De Floriani, and K. Weiss  
*Proceedings of the 22Nd ACM SIGSPATIAL International Conference on Advances in Geographic Information Systems, 223--232, 2014* - [doi](https://doi.org/10.1145/2666310.2666412)
- **A primal/dual representation for discrete Morse complexes on tetrahedral meshes**  
K. Weiss, F. Iuricich, R. Fellegara, and L. De Floriani  
*Computer Graphics Forum (Proceedings of Eurovis 2013), 32(3), 361--370, 2013* - [doi](https://doi.org/10.1111/cgf.12123)
- **A Spatial Approach to Morphological Feature Extraction from Irregularly Sampled Scalar Fields**  
L. De Floriani, R. Fellegara, F. Iuricich, and K. Weiss  
*Proceedings of the Third ACM SIGSPATIAL International Workshop on GeoStreaming, 40--47, 2012* - [doi](https://doi.org/10.1145/2442968.2442974)
- **The PR-star Octree: A Spatio-topological Data Structure for Tetrahedral Meshes**  
K. Weiss, R. Fellegara, L. De Floriani, and M. Velloso  
*Proceedings of the 19th ACM SIGSPATIAL International Conference on Advances in Geographic Information Systems, 92--101, 2011* - [doi](https://doi.org/10.1145/2093973.2093987)

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

The library requires the [boost library](http://www.boost.org/) (for dynamic_bitset class), the [GMP library](https://gmplib.org/), the [BitMagic library](http://bitmagic.io/) and [cmake](https://cmake.org/) installed in your system.

In Debian-based OSes (including Ubuntu and Linux Mint) these requirements can be satisfied by running the following command in the terminal:
```
#!

sudo apt install cmake libboost-all-dev libgmp3-dev bmagic
```

Once installed, execute from the command line in the root folder of the project the following command:
```
#!

cmake CMakeList.txt
```
and once configured execute:
```
#!

make
```
This latter command generates a portable library file, located into `lib` folder, as well as some executables in `dist` folder.

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
