## To compile and run ##
```
g++ -c factorial.cpp PrimePowers.cpp tenJ.cpp test.cpp tet.cpp theta.cpp 
g++ -o 10j *.o
./10j
```
## Genesis of Spin Networks ##

In the late 60's, Sir Roger Penrose created the theory of spin networks as an approach to quantum
geometry. Musing on the nature of the continuum concept, which had served as the bona fide representation of spacetime, Penrose sought to vanquish this continuum as purely a mathematical utility devoid of any grounding in physical reality, in favor of combinatorial principles. The barebones mathematical framework date back to the work of _Eudoxus_ in 4th century B.C.E. He argues against the spacetime constructions:
- lattice space-time of Schild
- discrete causal space 
- structure based on Ahmavaara's large finite field 

> If points are not to be the basic elements of discrete spacetime, then how are we to decide what these basic elements should in fact be?



 the total angular momentum(j-values) can be represented in terms of cominbatorial terms. A spin-network evaluates to a nonnegative integer, the norm. 

> If a spin-network has a number of large end units such that the angle between any two of them is well defined in the above sense, then these angles can be consistently interpreted as angles between directions in a Euclidean 3-d space.


His seminal work has bloomed into the theory of spin foams and loop quantum gravity, bequeathed by Rovelli and Smolin. 

## preliminaries ##
The _antisymmetrizer_ is a diagram sum described by a bundle of $N$ lines. This diagrammatic sum is the key ingredient for making spin netowrk calculations of Clebsch-Gordon coefficients, __3j and 6j symbols__, and other representations of angular momentum.
It is defined by:

$ \sum_{\sigma \in S_N} sgn(\sigma)$ where $\sigma$ iterates over all permutations in the group $S_N$ 




def 2-skeleton

def 4-simplex

_ A 10j symbol is a Spin(4) spin network with 5 vertices and 10 edges, which are labeled by spins._

A spin foam is the 2-d analogue of a Feynman diagram; A 2-d cell complex with polygonal faces labelled by representations and edges labelled by intertwinin goperators. The time evolution of the system is given by a linear combination of quantum histories, with weights given by amplitudes. The 1-d spin networks give rise to the 2-d nature of spin foams, just as point particles (0-d) give rise to 1-d Feynman diagrams.

The computation of an amplitude for any spin foam is the product of __face__, __edge__, and __vertex__ amplitudes. The partition function is then computed as a sum of these amplitudes.

Lorentzian quantum gravity: any QFT whose partition function is:
$\int e^{iS}$ where $S$ is the Einstein-Hilbert action for a lOrentizan metric on spacetime.

Riemannian QG seems to have limited relevance to real-world physics. rotation group for Riemannian models are compact thus the irreducible unitary representations are indexed by discrete paramters, allowing one to show convergence of a single spin foam amplitude.

1. Ponzano-Regge model of 3-d Riemannian QG. 
    triangulate a given 3-manifold
    express $\int e^{iS}$ as a sum over spin foams in the dual 2-skeleton of resulting triangulation
    gauge group: Spin(3) = SU(2) of 3d rotation group

2. Turaev and Viro q-deformed model: regularized 1. by replacing SU(2) with corresponding quantum group $SU_q(2)$
    - triangulation independent
    - deformation parameter $q$ related to the cosmological constant by :

3. Barrett and Crane spin foam model of 4-d riemannian QG   
    - partition function computed as sum over spin foams in the dual 2-skeleton of a triangulated 4-manifold
    - faces labelled by representations $Spin(4) = SU(2) \times SU(2)$, $j \otimes j$
    - edges: barrett-crane intertwiner
    - first motivation for using 10j symbols for vertex amplitudes 

    was impossible to determine if partition function converges or not since did not give formulas for edge and face amplitudes. Can q-deform this model based on the quantum group $SU_q(2) \times SU_{\bar(q)2}$

4. DFKR proposal: BC model arises from QFT on product of 4 copies of 3-sphere 
    - feynman diagrams correspond to the spin foams in BC model, vertex amplitudes same
    - group field theory gives formulas for edge and face amplitudes 
    - computes partition function by summing over spin foams lying in dual 2-skeleton of ALL triangulations of ALL compact 4-manifolds 
    - ___extends BC model to incorporate sum over triangulations and sum over topologies__


    for any simplicial complex:
        finite set of 4-simplices
        routine:
            attach distinct ones pairwise along tetrahedral faces
            stop: all faces are paired 

5. Perez and Rovelli modification:
    - modified edge and face amplitudes 
    - goal: eliminate divergences from model 

6. New version..

__need to compute expectation value of observables, not partition functions__

this then leads to the question of whether a spin foam model reduces to GR in large scale limit.

the amplitude of a single spin foam is given by:

$\prod_{f\in\Delta_2}A(f) \prod_{e\in\Delta_3} A(e) \prod_{v\in\Delta_4} A(v)$
with the amplitudes of face, edge, vertex computed only using the spin $j(f)$ labelling that face. 

where the partition function is $Z(M) = \sum_F Z(F)$,

which, when imposing an ``infared" cutoff which rules out larg areas, $\abs(F) \leq J$ the above partition function becomes a __finite sum__

$Z_J(M) = \sum_{\abs(F)\leq J} Z(F)$,

In algorithmic treatments of low-level topology, manifolds are typically represented by triangulations(see Regina). A d-manifold triangulation consists of a set of d-simplices along with instructions on how some or all of their $(d-1)$-dimensional facets should be glued together in pairs. Let $M$ be a triangulated compact 4-manifold and let $Delta_n$ be the set of n-simplices in the triangulation.

For example, take $M$ to be a 4-sphere triangulated as the boundary of a 5-simplex.


Thus, what is a spin network? A graph whose vertices are labeled by tensors, connected by edges which contain the tensor contraction "rules". This graph evaluates to
a complex number with input as ten spins. 

Algorithm:

Objects: vertices of type tensor; The barrett-crane intertwiners

morphisms: edges with tensor contraction capabilities

Space: $Spin(4) \equiv SU(2) \times SU(2)$

Class Triangulation of a 4-manifold
    Subclass dual 2-skeleton




evaluation of spin network: v1, v2

v1: works directly with the spin(4) network 
1. choose bases for the spin representations labelling the edges 
2. compute components of the tensors representing the five vertices of the network, 
aka the Barrett-Crane intertwiners 
3. compute contraction of of the tensors in same way 
    3a. direct contraction
    3b. staged contraction
    iteratively contract the tensors in the decagonal su2 network. 
    3c. 3cut
    3d. 2cut 
    trace of resulting operator from ray from center, cut the two edges

v2: converts spin(4) netowrk into a 5-fold sum over SU(2) networks via expansion
of each Barrett-Crane intertwiner:
Note: the two decagonal networks obtained by expanding each BC intertwiner are the same, only need to evaluate one of them then square the result. 

## References ##

https://math.ucr.edu/home/baez/penrose/
doi:10.1142/4256 part II, sec.12



