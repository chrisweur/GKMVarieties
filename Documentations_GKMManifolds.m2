


-------------------------------------------------------------------------------------------
--------------------------------------< Documentation >------------------------------------
-------------------------------------------------------------------------------------------


beginDocumentation()


doc ///
	Key
		GKMManifolds
	Headline
		a package for computations with GKM manifolds and moment graphs
	Description
		Text
			A GKM manifold is a variety X, often assumed smooth and complete, with an 
			action of an algebraic torus T satisfying the following conditions:
			(i) X is equivariantly formal with respect to the the action of T,
			(ii) X has finitely many T-fixed points, and (iii) X has finitely
			many one-dimensional T-orbits.  The data of the zero and one dimensional
			T-orbits of X defines a moment graph, with which one can carry out
			T-equivariant cohomology and T-equivariant K-theory computations by
			combinatorial means.  This package provides methods for these computations
			in Macaulay2.
			
		Text
			For mathematical background see:
			
			@UL{
			{"T. Braden and R. MacPherson", EM "From moment graphs to intersection cohomology", "Math. Ann. 321 (2001), 533-551."},
			{"E. Bolker, V. Guillemin, and T. Holm", EM "How is a graph like a manifold?", 	"arXiv:math/0206103"},
			{"M. Goresky, R. Kottwitz, and R. MacPherson", EM "Equivariant cohomology, Koszul duality, and the localization theorem", " Invent. Math. 131 (1998), no. 1, 25-83."},
			{"I. Rosu", EM "Equivariant K-theory and equivariant cohomology", "with an Appendix by I. Rosu and A. Knutson",  "Math. Z. 243 (2003), 423-448."},
			{"J. Tymoczko", EM "An introduction to equivariant cohomology and homology, following Goresky, Kottwitz, and MacPherson", "Contemp. Math. 388 (2005), 169-188."},
			{"G. Vezzosi and A. Vistoli", EM "Higher algebraic K-theory for actions of diagonalizable groups", "Invent. Math. 153 (2003), no. 1, 1–44."}
			}@ 
		    
		SUBSECTION "Examples"
		Text
			"smooth toric varieties"
			"generalized flag varieties"
			--"Example: generalized Schubert varieties"

		SUBSECTION "Contributors"
		Text
			The following people have contributed code, improved existing code, or enhanced the documentation:
			@UL{
			{HREF("https://www.mis.mpg.de/combag/members/tim-seynnaeve.html","Tim Seynnaeve")}
			}@


///


doc ///
	Key
		TVariety
	Headline
		the class of all T-varieties
	Description
		Text
			To see how to specify a T-variety, see @TO tVariety@.

			Describe basic functionality of the package...		
		Example
			Here....

///


doc ///
	Key
		tVariety
		(tVariety, List, List, Ring)
	Headline
		constructs a T-variety
	Usage
		X = tVariety(P,T,L,R)
	Inputs
		P:List
			of Torus fixed points
		T:HashTable
			encoding the one-dimensional Torus fixed orbits.
			The values are the orbits and the corresponding key is the pair
			of Torus fixed points in the orbit
		L:List
			of lists; each entry consists of a Torus fixed point along with
			the characters of the contracting affine chart around it
		R:Ring
			representing the characteristic ring the torus
	Outputs
		X:TVariety
	Description
		Text
			Here..
			
		Example
			Here...

	Caveat
		This function does not check if X defines a T-variety
	
	SeeAlso
		(symbol **, TVariety, TVariety)
		tFlagVariety
		tMap
///


doc ///
	Key
		(symbol **, TVariety, TVariety)
	Headline
		product of T-varieties
	Usage
		X ** Y
	Inputs
		X:TVariety
		Y:TVariety
	Outputs
		X:TVariety
		    product of X and Y
	Description
		Text
			Here..
			
		Example
			Here...

	Caveat
		This function does not check if X defines a T-variety
	
///


doc ///
	Key
		TKClass
	Headline
		a type of HashTable
	Description
		Text
			Describe it...
			
	SeeAlso
	    

///


doc ///
	Key
		tKClass
		(tKClass, TVariety, List)
	Headline
		constructs a equivariant K-class
	Usage
		C = tKClass(X,L)
	Inputs
		X:TVariety
		L:List
			of Laurent polynomials corresponding to each T-invariant point
	Outputs
		C:TKClass
	Description
		Text
			Here..
			
		Example
			Here...

	Caveat
		This function does not check if X defines a T-variety - see 
		@TO2{(isWellDefined, TKClass), "isWellDefined"}@.
	
	SeeAlso
		(isWellDefined, TKClass)
		(symbol *, TKClass, TKClass)
		(symbol +, TKClass, TKClass)
		pullback
		pushforward
///

