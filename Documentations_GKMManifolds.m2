-*----
Things to ask Justin:

Why does UL not work?
The main page will not work
What's up with the @something@
How to send version 0.1 to be on the M2 website
Printing on isWellDefined



----*-


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
			A GKM manifold is a variety $X$, often assumed to be smooth and complete, with an 
			action of an algebraic torus $T$ satisfying the following conditions:
			(i) X is equivariantly formal with respect to the the action of $T$,
			(ii) X has finitely many $T$-fixed points, and (iii) X has finitely
			many one-dimensional $T$-orbits.  The data of the zero and one dimensional
			$T$-orbits of $X$ defines a moment graph, with which one can carry out
			$T$-equivariant cohomology and $T$-equivariant $K$-theory computations by
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
			"Example: generalized flag varieties"
			"Example: smooth toric varieties"

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
		the class of all GKM manifolds
	Description
		Text
			A @TO "TVariety" @TT "X"@ is a @TO "MutableHashTable"@ representing a GKM manifold $X$ with an action of a torus $T$.
			Its keys include:
			@UL{
			{@TT "points"@, whose value is a list representing the $T$-fixed points of $X$},
			{@TT "charRing"@, whose value is a ring representing the character ring of $T$},
			{@TT "momentGraph"@, whose value is the @TO "MomentGraph"@ of $X$},
			{@TT "charts"@, whose value is a @TO "HashTable"@ representing the (negatives of) characters of the action of $T$
			on each $T$-invariant affine chart around a $T$-fixed point.  The keys of @TT "X.charts"@ are @TT "X.points"@, and the values are lists consisting of lists of integers}.
			}@

	SeeAlso
		tVariety
		"Example: generalized flag varieties"
		"Example: smooth toric varieties"




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
			
		--Example
			Here...

	Caveat
		This function does not check if X defines a T-variety
	
	SeeAlso
		(symbol **, TVariety, TVariety)
		tGeneralizedFlagVariety
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
			
		--Example
			Here...

	Caveat
		This function does not check if X defines a T-variety
	
///


doc ///
	Key
		TKClass
	Headline
		the class of all $T$-equivariant $K$-classes
	Description
		Text
			For $X$ a GKM manifold with an action of a torus $T$ whose character ring is $R$,
			a $T$-equivariant $K$-class $C \in K_T^0(X)$ of is encoded by its image in $K_T^0(X^T) \simeq \prod_{x\in X^T} R$
			under the injective restriction map $K_T^0(X) \to K_T^0(X^T)$.  See REFERENCE HERE.

		Text
			A @TO "TKClass"@ @TT "C"@ is a @TO "HashTable"@
			consisting of two keys:
			@UL{
			{@TT "tvar"@, whose value is a @TO "TVariety"@ of which @TT "C"$ is a $K$-class of},
			{@TT "hilb"@, whose value is a @TO "HashTable"@; its keys are @TT "X.points"@ and the values are
			Laurent polynomials in the character ring $R$ representing the values of the $C$ under the restriction map}.
			}@

	SeeAlso
		tKClass
		pushforward
		pullback
		tChi	    

///



doc ///
	Key
		tKClass
		(tKClass, TVariety, List)
	Headline
		constructs a $T$-equivariant $K$-class
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
			This method creates a @TO TKClass@ given a @TO TVariety@ @TT "X"@ and a list @TT "L"@ of Laurent polynomials in its
			character ring.  The order of Laurent polynomials in the list must correspond to the order of the list
			of $T$-fixed points @TT "X.points"@.

		Text
			The following example is the class of $O(1)$ on the projective space $\mathbb P^3$.
			
		Example
			PP3 = tProjectiveSpace 3;
			R = PP3.charRing;
			L = gens R
			C = tKClass(PP3,L) --the class of O(1) on PP3
			assert (C === ampleTKClass P33)
			assert (isWellDefined C)




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

doc ///
	Key
		(isWellDefined, TKClass)
	Headline
		whether the input is a well-defined T-equivariant K-class
	Usage
		isWellDefined C
	Inputs
		C:TKClass
	Outputs
		:Boolean
			whether or not a list of Laurent polynomials satisfies edge compatibility condition
			--prints out the edges of the moment graph for which @TT "C"@ fails the compatibility condition
	Description
		Text
			If $\{f_x \mid x\in X^T\}$ is a collection of Laurent polynomials in the
			character ring $\mathbb Z[T_0, \ldots, T_n]$ of a @TO TVariety@ $X$, one per each $T$-fixed point, representing an element $C$ of $K_T^0(X^T)$,
			then $C$ is in the image of $K_T^0(X)$ under the injective restriction map $K_T^0(X)\to K_T^0(X^T)$ if and only if
			it satisfies the following "edge compatibility condition":

			For each one-dimensional $T$-orbit-closure in $X$ with boundary points $x$ and $x'$, one has
			$f_x \equiv f_{x'} \mod 1 - T^\lambda_{x,x'}$ where $\lambda_{x,x'}$ is the character of the action of $T$ on the
			one-dimensional orbit.

		Example
			PP3 = tProjectiveSpace 3
			isWellDefined ampleTKClass PP3 --the O(1) class on PP3 is well-defined
			badC = tKClass(PP3, reverse gens PP3.charRing) --reverse the order of Laurent polynomials defining the O(1) class
			isWellDefined badC --no longer well-defined

	Caveat
		A @TO MomentGraph@ must be defined on the @TO TVariety@ on which the @TO TKClass@ is a $K$-class of.

	SeeAlso
		TKClass
		tKClass		


///

doc ///
	Key
		(symbol *, TKClass, TKClass)
	Headline
		computes the product of two T-equivariant K-classes
	Usage
		C1 * C2
	Inputs
		C1:TKClass
		C2:TKClass
	Output
		:TKClass
			the product of C1 and C2
	Description
		Text
			This method computes the product of two $T$-equivariant $K$-classes.

		Example
			Gr24 = tGeneralizedFlagVariety(3,{2}); --the Grassmannian of projective lines in projective 3-space
			O1 = ampleTKClass Gr24 -- the O(1) bundle on Gr24 as a T-equivariant K-class
			O2 = O1 * O1
			peek O2

	SeeAlso
		tKClass
		(symbol ^, TKClass, ZZ)
		(symbol +, TKClass, TKClass)
///

doc ///
	Key
		(symbol ^, TKClass, ZZ)
	Headline
		computes powers of a T-equivariant K-classes
	Usage
		C^n
	Inputs
		C:TKClass
		n:ZZ
	Output
		:TKClass
			the n-th power of C
	Description
		Text
			This method computes the $n$-th power of a $T$-equivariant $K$-class $C$.

		Example
			Gr24 = tGeneralizedFlagVariety(3,{2}); --the Grassmannian of projective lines in projective 3-space
			O1 = ampleTKClass Gr24 -- the O(1) bundle on Gr24 as a T-equivariant K-class
			O2 = O1 * O1
			peek O

	Caveat
		$n$ is allowed to be negative only when $C$ is a line bundle, or a direct sum of copies of a line bundle.

	SeeAlso
		tKClass
		(symbol *, TKClass, TKClass)
		(symbol +, TKClass, TKClass)
///


doc ///
	Key
		(symbol +, TKClass, TKClass)
	Headline
		computes the sum of two T-equivariant K-classes
	Usage
		C1 + C2
	Inputs
		C1:TKClass
		C2:TKClass
	Output
		:TKClass
			the sum of C1 and C2
	Description
		Text
			This method computes the sum of two $T$-equivariant $K$-classes.

		Example
			Gr24 = tGeneralizedFlagVariety(3,{2}); --the Grassmannian of projective lines in projective 3-space
			O1 = ampleTKClass Gr24 -- the O(1) bundle on Gr24 as a T-equivariant K-class
			E = O1 + (O1*O1)
			peek E

	SeeAlso
		tKClass
		(symbol *, TKClass, TKClass)
///





