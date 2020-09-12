-*----
Things to ask Justin:


How to send version 0.1 to be on the M2 website

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
			(i) $X$ is equivariantly formal with respect to the the action of $T$,
			(ii) $X$ has finitely many $T$-fixed points, and (iii) $X$ has finitely
			many one-dimensional $T$-orbits.  The data of the zero and one dimensional
			$T$-orbits of $X$ defines the moment graph of $X$, with which one can carry out
			$T$-equivariant cohomology and $T$-equivariant $K$-theory computations by
			combinatorial means.  This package provides methods for these computations
			in Macaulay2.
			
		Text
			For mathematical background see:
			
			@UL{
			{"T. Braden and R. MacPherson.  From moment graphs to intersection cohomology.  Math. Ann. 321 (2001), 533-551."},
			{"E. Bolker, V. Guillemin, and T. Holm.  How is a graph like a manifold?  arXiv:math/0206103."},
			{"R. Dinu, C. Eur, and T. Seynnaeve.  K-theoretic Tutte polynomials of morphisms of matroids.  arXiv:math/2004.00112."},
			{"A. Fink and S. Speyer.  K-classes for matroids and equivariant localization.  Duke Math. J. 161 (2012), no. 14, 2699-2723."},
			{"M. Goresky, R. Kottwitz, and R. MacPherson.  Equivariant cohomology, Koszul duality, and the localization theorem. Invent. Math. 131 (1998), no. 1, 25-83."},
			{"I. Rosu.  Equivariant K-theory and equivariant cohomology.  With an Appendix by I. Rosu and A. Knutson.  Math. Z. 243 (2003), 423-448."},
			{"J. Tymoczko.  An introduction to equivariant cohomology and homology, following Goresky, Kottwitz, and MacPherson.  Contemp. Math. 388 (2005), 169-188."},
			{"G. Vezzosi and A. Vistoli.  Higher algebraic K-theory for actions of diagonalizable groups. Invent. Math. 153 (2003), no. 1, 1–44."}
			}@

		Text
			@SUBSECTION "Contributors"@
		Text
			The following people have contributed code, improved existing code, or enhanced the documentation:
			@HREF("https://www.mis.mpg.de/combag/members/tim-seynnaeve.html","Tim Seynnaeve")@.
		    
	SeeAlso
		"Example: generalized flag varieties"
		"Example: smooth toric varieties"

		

///


doc ///
	Key
		"Example: generalized flag varieties"
	Description
		Text
			Generalized flag vareities are GKM manifolds.  If $G$ be a reductive complex Lie group and $P$ a 
			parabolic subgroup containing a maximal torus $T$, then the generalized flag variety $G/P$ is a GKM manifold 
			with the action of $T$.  This package allows users to create a generalized flag variety for classical Lie types 
			($A$, $B$, $C$, and $D$) as a @TO "TVariety"@ with conventions explicitly laid out as follows.

		Text
			For type $A_{n-1}$, the group $G$ is $GL_{n}$, and the torus $T$ is $diag(t_1, \ldots, t_n)$, the group of invertible diagonal matrices.

		Text
			For type $B_n$, the group $G$ is $SO_{2n+1}$, where the symmetric bilinear form on $\mathbb C^{2n+1}$ is given by the matrix $\begin{pmatrix} 0 & I_n & 0 \\ I_n & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}$, and the torus $T$ is $diag(t_1, \ldots,t_n, t_1^{-1}, \ldots, t_n^{-1}, 1)$.

		Text
			For type $C_n$, the group $G$ is $Sp_{2n}$, where the alternating bilinear form on $\mathbb C^{2n}$ is given by the matrix $\begin{pmatrix} 0 & -I_n \\ I_n & 0 \end{pmatrix},$ and the torus $T$ is $diag(t_1, \ldots,t_n, t_1^{-1}, \ldots, t_n^{-1})$.

		Text
			For type $D_n$, the group $G$ is $SO_{2n}$, where the symmetric bilinear form on $\mathbb C^{2n}$ is given by the matrix $\begin{pmatrix} 0 & I_n \\ I_n & 0 \end{pmatrix}$, and the torus $T$ is $diag(t_1, \ldots,t_n, t_1^{-1}, \ldots, t_n^{-1})$.

		Text
			In all the cases, the standard action of $(\mathbb C^*)^n$ on $\mathbb C^n$ is defined by $(t_1, \ldots, t_n) \cdot (x_1, \ldots, x_n) = (t_1^{-1}x_1, \ldots, t_n^{-1}x_n)$.  

		Text
			(With a Borel subgroup of $G$ fixed), let $\{w_1, \ldots, w_n\}$ be the fundamental weights, 
			ordered for each type in the same way as in Example 3.7 of [ACEP20]. 
			For a sequence $(a_1, \ldots, a_n)\in \mathbb N^n$ of nonnegative integers, 
			let $I = \{i \mid a_i \neq 0\}$ and $P_I$ the corresponding parabolic subgroup of $G$.
			Then the generalized flag variety $G/P_I$ is embedded in the irreducible representation of $G$ 
			with the highest weight $a_1w_1 + \cdots a_nw_n$.    

			These generalized flag varieties can be created as a @TO "TVariety"@ using the method 
			@TO (tGeneralizedFlagVariety, String, ZZ, List)@.  For instance, the Grassmannian $Gr(2,4)$ of
			2-dimensional subspaces in $\mathbb C^4$, embedded in $\mathbb P^5$ by the usual Plucker embedding, 
			can be created as follows.

		Example
			Gr24 = tGeneralizedFlagVariety("A",3,{2})
			peek Gr24

		Text
			The moment graph of $Gr(2,4)$ is the 1-skeleton of the hypersimplex $\Delta(2,4)$, a.k.a. the octahedron.

		Example
			G = momentGraph Gr24
			underlyingGraph G

		Text
			The $O(1)$ line bundle on $Gr(2,4)$ via its Plucker embedding can be accessed as follows, and 
			its Lefschetz trace (a.k.a. T-equivariant Euler characteristic) is a Laurent polynomial in 
			the character ring of the torus $T$, whose terms correspond to be weights of the second
			exterior power of the standard representation of $GL_4$.

		Example
			O1 = ampleTKClass Gr24 --the O(1) bundle on Gr24 via its Plucker embedding
			tChi O1

		Text
			If $Gr(2,4)$ is embedded differently, say by $O(2)$ line bundle instead, the Lefschetz trace changes 
			accordingly, and its coefficients record the multiplicities of the associated weight spaces in the second
			symmetric power of the second exterior power of the standard representation of $GL_4$.

		Example
			tChi (O1^2)

		Text
			The following example features the isotropic Grassmannian $SpGr(2,4)$ consisting of 
			2-dimensional subspaces in $\mathbb C^4$ that are isotropic with respect to the standard alternating form.
			Its moment graph is the a complete graph on 4 vertices.

		Example
			SpGr24 = tGeneralizedFlagVariety("C",2,{2})
			peek SpGr24
			underlyingGraph momentGraph SpGr24
			tChi ampleTKClass SpGr24


		
			

///


doc ///
	Key
		TVariety
	Headline
		the class of all GKM manifolds
	Description
		Text
			A @TO "TVariety"@ TT "X" is a @TO "MutableHashTable"@ representing a GKM manifold $X$ with an action of a torus $T$.
			Its keys include:
			UL{
			{TT "points", "whose value is a list representing the $T$-fixed points of $X$"},
			{TT "charRing", "whose value is a ring representing the character ring of $T$"},
			{TT "momentGraph", "whose value is the @TO "MomentGraph"@ of $X$"},
			{TT "charts", "whose value is a @TO "HashTable"@ representing the (negatives of) characters of the action of $T$
			on each $T$-invariant affine chart around a $T$-fixed point.  The keys of TT "X.charts" are TT "X.points", and the values are lists consisting of lists of integers."}
			}

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
		the class of all T-equivariant K-classes
	Description
		Text
			For $X$ a GKM manifold with an action of a torus $T$ whose character ring is $R$,
			a $T$-equivariant $K$-class $C \in K_T^0(X)$ of is encoded by its image in $K_T^0(X^T)$, which is isomorphic
			to $\prod_{x\in X^T} R$,
			under the injective restriction map $K_T^0(X) \to K_T^0(X^T)$.  See REFERENCE HERE.

		Text
			A @TO "TKClass"@ C is a @TO "HashTable"@
			consisting of two keys:

			@UL{
			{TT "tvar", ", whose value is a ", TO "TVariety", " of which C is a K-class of"},
			{TT "hilb", ", whose value is a ", TO "HashTable", "; its keys are ", TT "X.points", " and the values are
			Laurent polynomials in the character ring representing the values of the K-class under the restriction map."}
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
		constructs a T-equivariant K-class
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
			assert (C === ampleTKClass PP3)
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
	Outputs
		:TKClass
			the product of C1 and C2
	Description
		Text
			This method computes the product of two $T$-equivariant $K$-classes.

		Example
			Gr24 = tGeneralizedFlagVariety("A",3,{2}); --the Grassmannian of projective lines in projective 3-space
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
	Outputs
		:TKClass
			the n-th power of C
	Description
		Text
			This method computes the $n$-th power of a $T$-equivariant $K$-class $C$.

		Example
			Gr24 = tGeneralizedFlagVariety("A",3,{2}); --the Grassmannian of projective lines in projective 3-space
			O1 = ampleTKClass Gr24 -- the O(1) bundle on Gr24 as a T-equivariant K-class
			O2 = O1 * O1
			peek O2

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
	Outputs
		:TKClass
			the sum of C1 and C2
	Description
		Text
			This method computes the sum of two $T$-equivariant $K$-classes.

		Example
			Gr24 = tGeneralizedFlagVariety("A",3,{2}); --the Grassmannian of projective lines in projective 3-space
			O1 = ampleTKClass Gr24 -- the O(1) bundle on Gr24 as a T-equivariant K-class
			E = O1 + (O1*O1)
			peek E

	SeeAlso
		tKClass
		(symbol *, TKClass, TKClass)
///

doc ///
	Key
		(tGeneralizedFlagVariety, String, ZZ, List)
	Headline
		makes a generalized flag variety
	Usage
		tGeneralizedFlagVariety(LT,d,L)
	Inputs
		LT:String
			"A", "B", "C", or "D"
		d:ZZ
			the dimension of the root system
		L:List
			of integers strictly between 0 and d+1
	Outputs
		:TVariety
		    the T-variety representing the corresponding generalized flag variety.

	Description
		Text
			Let $w = a_1w_1 + \cdots + a_nw_n$ be a weight in the root system of the type $LT_d$ where $a_i$ is the number 
			of times $i$ appears in the list $L$.  Let $G$ be the Lie group corresponding to $LT_d$.
			This method outputs the T-variety representing the generalized flag variety $G/P$ embedded in the irreducible 
			representation of $G$ with the highest weight $w$.

	SeeAlso
		"Example: generalized flag varieties"

///

-*--
doc ///
	Key
		(pullback, TMap)
	Headline
		computes the pullback map
	Usage
		pullback(f)

///
--*-


TEST ///
PP3 = tProjectiveSpace 3
assert (3 == #gens PP3.charRing)
///
