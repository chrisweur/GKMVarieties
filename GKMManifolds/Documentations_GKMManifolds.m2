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
			For type $B_n$, the group $G$ is $SO_{2n+1}$, where we set the standard symmetric bilinear form on $\mathbb C^{2n+1}$ 
			to be is given by the matrix $\begin{pmatrix} 0 & I_n & 0 \\ I_n & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}$, and the torus $T$ is $diag(t_1, \ldots,t_n, t_1^{-1}, \ldots, t_n^{-1}, 1)$.

		Text
			For type $C_n$, the group $G$ is $Sp_{2n}$, where we set the standard the alternating bilinear form on $\mathbb C^{2n}$ 
			to be given by the matrix $\begin{pmatrix} 0 & -I_n \\ I_n & 0 \end{pmatrix},$ and the torus $T$ is 
			$diag(t_1, \ldots,t_n, t_1^{-1}, \ldots, t_n^{-1})$.

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
			@TO tGeneralizedFlagVariety@.  For instance, the Grassmannian $Gr(2,4)$ of
			2-dimensional subspaces in $\mathbb C^4$, embedded in $\mathbb P^5$ by the usual Plucker embedding, 
			can be created as follows.

		Example
			Gr24 = tGeneralizedFlagVariety("A",3,{2})
			peek Gr24

		Text
			The @TO MomentGraph@ of $Gr(2,4)$ is the 1-skeleton of the hypersimplex $\Delta(2,4)$, a.k.a. the octahedron.

		Example
			G = momentGraph Gr24
			underlyingGraph G

		Text
			The $O(1)$ line bundle on $Gr(2,4)$ via its Plucker embedding can be accessed by @TO ampleTKClass@.
			The method @TO (tChi, TKClass)@ computes its Lefschetz trace (a.k.a. T-equivariant Euler characteristic),
			which in this case is the Laurent polynomial in the character ring of the torus $T$ 
			whose terms correspond to be weights of the second exterior power of the standard representation of $GL_4$.

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
			The Schubert decomposition of $Gr(2,4)$, and more generally the Bruhat decomposition of $G/P$, can be accessed 
			by the method @TO (bruhatOrder, TVariety)@, which outputs the poset of the Bruhat order.  Moreover, the Schubert
			varieties can be created via the method @TO tGeneralizedSchubertVariety@.

		Example
			P1 =  bruhatOrder Gr24
			Sch = tGeneralizedSchubertVariety(Gr24,{set{1,2}})
			P2 = bruhatOrder Sch
			--{P1,P2}/displayPoset --to view the posets in pdf

		Text
			The "forgetful" map from the full flag variety $Fl(4)$ of full flags in $\mathbb C^4$ to $Gr(2,4)$, 
			given by forgetting the subpsaces in the flag except for the 2-dimensional one, and be created as 
			a @TO TMap@ by the method @TO tFlagMap@.

		Example
			Fl4 = tGeneralizedFlagVariety("A",3,{1,2,3},Gr24.charRing) --Fl(4) with the torus having the same character ring as Gr24
			f = tFlagMap(Fl4,Gr24)
			Fl4 === f.source and Gr24 === f.target
			(trivialTKClass Gr24) === (pushforward f)(trivialTKClass Fl4)

		Text
			The (derived) pushforward of the structure sheaf of $Fl(4)$ is the structure sheaf of $Gr(2,4)$ since the higher direct images vanish under the forgetful map.

		Example
			(trivialTKClass Gr24) === (pushforward f)(trivialTKClass Fl4)

		Text
			The following example features the isotropic Grassmannian $SpGr(2,6)$ consisting of 
			2-dimensional subspaces in $\mathbb C^6$ that are isotropic with respect to the standard alternating form.
			The vertices of its moment graph can be considered as the vertices of the cuboctahedron.

		Example
			SpGr26 = tGeneralizedFlagVariety("C",3,{2})
			peek SpGr26
			momentGraph SpGr26

		Text
			The second fundamental representation of $Sp_{6}$ is 14-dimensional with 12 extremal weights.

		Example
			tChi ampleTKClass SpGr26

		Text
			The following example features the isotropic Grassmannian $OGr(2,5)$ consisting of
			3-dimensional subspaces in $\mathbb C^5$ that are isotropic with respect to the standard symmetric form.
			Its moment graph is the a complete graph on 4 vertices.
			Note that Spin groups and their representations are not implemented, so for the type $B_n$ the coefficient
			$a_n$ need be a multiple of 2.

		Example
			OGr25 = tGeneralizedFlagVariety("B",2,{2,2}) --inputing {2} instead of {2,2} results in error: spin groups not implemented yet
			peek OGr25
			tChi ampleTKClass OGr25

	SeeAlso
		tGeneralizedFlagVariety
		tFlagMap





///


doc ///
	Key
		"Example: smooth toric varieties"
	Description
		Text
			G

		Text
			In

		Example
		    	hi	

       		Text
			The following example shows the third Hizerburch surface as a @TO "NormalToricVariety"@ and as a @TO "TVariety"@ 
		Example
			FF3 = hirzebruchSurface 3; X = tVariety FF3;
			peek FF3
			peek X
		Text 
		        We can recover the Toric variety as follows
		Example
		        Y = normalToricVariety(X);
		    	peek Y


	SeeAlso
		tVariety
		tProjectiveSpace

///


doc ///
	Key
		TVariety
	Headline
		the class of all GKM manifolds
	Description
		Text
		        A @TO "TVariety"@ $X$ is  a @TO "MutableHashTable"@ representing a GKM manifold $X$ with an action of a torus $T$.
			Its keys include:
		    	
			@UL{
			{TT "points", ", whose value is a list representing $T$-fixed points of $X$"},
			{TT "charRing", ", whose value is a ring representing the character ring of $T$"},
			{TT "momentGraph", ", whose value is the ", TO "MomentGraph", " of  $X$"},
			{TT "charts", ", whose value is a ", TO "HashTable", " representing the (negatives of) characters of the action of $T$
			on each $T$-invariant affine chart around a $T$-fixed point. The keys of ", TT "X.charts", " are ", TT "X.points", " and the
			values are lists consisting of lists of integers."}
			}@

		   			
	SeeAlso
	    	tVariety
		"Example: generalized flag varieties"
		"Example: smooth toric varieties"




///


doc ///
	Key
		tVariety
		(tVariety, List, Ring)
		(tVariety, List, List, Ring)
		(tVariety, MomentGraph)
		(tVariety, MomentGraph, Ring)
		(tVariety, NormalToricVariety)
		(tVariety, NormalToricVariety,Ring)
		(tVariety, TKClass)
	Headline
		constructs a T-variety
	Usage
		X = tVariety(L,R)
		X = tVariety(L,M,R)
		X = tVariety(G)
		X = tVariety(G,R)
		X = tVariety(R)
		X = tVariety(Y,R)
		X = tVariety(C)
	Inputs
		L:List
			of $T$-fixed points of $X$
		M:List
			of lists; each list consists of the (negatives of) characters of the 
			action of $T$ on each $T$-invariant affine chart around the $T$-fixed point 
			corresponding to L_i
		G:MomentGraph
			representing the one dimensional $T$-orbits of $X$
		R:Ring
			representing the characteristic ring of $T$
		Y:NormalToricVariety
		C:TKClass	
	Outputs
		X:TVariety
	Description
		Text
			The minimum data needed to make a GKM Manifold are the set of $T$-fixed points
			and the character ring. Here is an example with projective space	
		Example
			L = {0,1,2,3}; R = ZZ[T_0..T_3];
			X = tVariety(L,R)
		Text	
			If necessary, we can add the (negatives of) characters of the action of $T$ on each 
			$T$-invariant chart of $X$. Note that the i-th entry of the list below corresponds to
			the i-th entry of L.
		Example
		    	M = {{{-1, 1, 0, 0}, {-1, 0, 1, 0}, {-1, 0, 0, 1}},
			    {{1, -1, 0, 0}, {0, -1, 1, 0}, {0, -1, 0, 1}},
			    {{1, 0, -1, 0}, {0, 1, -1, 0}, {0, 0, -1, 1}},
			    {{1, 0, 0, -1}, {0, 1, 0, -1}, {0, 0, 1, -1}}};
			Y = tVariety(L,M,R);
			peek Y
		Text
			To produce one of the generalized flag varieties we use the method @TO tGeneralizedFlagVariety@
			Here is an example of the Lagrangian Grassmannian $SpGr(2,4)$ consisting of 2-dimensional subspaces
			in $\mathbb C^4$ that are isotropic with respect to the standard alternating form.

		Example
			SpGr24 = tGeneralizedFlagVariety("C",2,{2})
			peek SpGr24
       		Text
		       Here is the complete flag variety for $SP_4$
		Example
			FlGr24 = tGeneralizedFlagVariety("C",2,{1,2})
			peek FlGr24
		
		Text
		    	The following example produces the Orthogonal Grassmaninnian OGr(2,5) from its
			moment graph.
		Example
		    V = {{set {0, 1}}, {set {0, "1*"}}, {set {"0*", 1}}, {set {"0*", "1*"}}};
		    edgs = {{{set {"0*", 1}}, {set {"0*", "1*"}}},
                            {{set {0, "1*"}}, {set {"0*", "1*"}}},
			    {{set {0, "1*"}}, {set {"0*", 1}}},
			    {{set {0, "1*"}}, {set {0, 1}}},
			    {{set {0, 1}}, {set {"0*", "1*"}}},
			    {{set {0, 1}}, {set {"0*", 1}}}};
		    wghts = {{0,-1},{-1,0},{-1,1},{0,1},{-1,-1},{-1,0}}
		    E = hashTable(apply(edgs, v-> (v,wghts)));
		    G = momentGraph(V,E,ZZ[T_0..T_2]);
		    Z = tVariety(G);
		    peek Z			

	Caveat
		This function does not check if $X$ is a well defined $T$-variety
	
	SeeAlso
		(symbol **, TVariety, TVariety)
		tProjectiveSpace
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
			Given two $T$-varieties $X$ and $Y$ with an action of a common torus $T$, the
			product is $X \times Y$ with the structure of a $T$-variety given by the 
			diagonal action of $T$.
			
		Text
		    	The following example exhibits the product of the Orthogonal Grassmannian OGr(2,5)
			with the Lagrangian Grassmannian SpGr(2,4) 
    	    	Example
		        X = tGeneralizedFlagVariety("C",2,{2});
			Y = tGeneralizedFlagVariety("B",2,{2,2});
			try(Z = X ** Y);
			try(peek Z)
	SeeAlso
		tVariety

	
///


doc ///
	Key
		TKClass
	Headline
		the class of all T-equivariant K-classes
	Description
		Text
			For $X$ a GKM manifold with an action of a torus $T$ whose character ring is $R$,
			a $T$-equivariant $K$-class $C \in K_T^0(X)$ of is encoded by its image in $K_T^0(X^T) = \prod_{x\in X^T} R$,
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
		tGeneralizedFlagVariety
		(tGeneralizedFlagVariety, String, ZZ, List)
		(tGeneralizedFlagVariety, String, ZZ, List, Ring)
	Headline
		makes a generalized flag variety as a T-variety
	Usage
		X = tGeneralizedFlagVariety(LT,d,L)
		X = tGeneralizedFlagVariety(LT,d,L,R)
	Inputs
		LT:String
			one of "A", "B", "C", and "D"
		d:ZZ
			the dimension of the root system
		L:List
			of integers strictly between 1 and d (inclusive)
		R:Ring
			the character ring of the torus acting on the generalized flag variety

	Outputs
		X:TVariety
		    representing the corresponding generalized flag variety

	Description
		Text
			Let $w = a_1w_1 + \cdots + a_nw_n$ be a weight in the root system of the type $LT_d$ where $a_i$ is the number 
			of times $i$ appears in the list $L$.  Let $G$ be the Lie group corresponding to $LT_d$.
			This method outputs the T-variety representing the generalized flag variety $G/P$ embedded in the irreducible 
			representation of $G$ with the highest weight $w$.

	SeeAlso
		"Example: generalized flag varieties"

///



doc ///
	Key
		TMap
	Headline
		the class of all T-equivariant morphisms between GKM manifolds
	Description
		Text
			Given two GKM manifolds $X$ and $Y$, a T-equivariant morphism from $X$ to $Y$
			induces a map from the $T$-fixed points of $X$ to the $T$-fixed points of $Y$.
			 			
		Text
		    	A @TO "TMap"@ C is a @TO "HashTable"@
			consisting of three keys:

			@UL{
			{TT "source", ", whose value is a ", TO "TVariety", " corresponding to the domain of f"},
			{TT "target", ", whose value is a ", TO "TVariety", " corresponding to the codomain of f"},
			{TT "ptsMap", ", whose value is a ", TO "HashTable", "; its keys are ", TT "X.points", " and the values are
			points of ", TT "Y.points", " that the key maps to."}
			}@
						
	SeeAlso
		(symbol **, TMap, TMap)
		(compose, TMap, TMap)
	        tMap
		tFlagMap
		pullback
		pushforward
		tChi
///


doc ///
	Key
		(symbol **, TMap, TMap)
	Headline
		computes the product of two T-equvariant morphisms
	Usage
		f ** g
	Inputs
		f:TMap
		g:TMap
	Outputs
		:TMap
			the product of f and g
	Description
		Text
			This method computes the cartesian product of two $T$-equivariant morphisms.

		Example
		    	X1 = tGeneralizedFlagVariety("A",4,{1,3})
			X2 = tGeneralizedFlagVariety("A",4,{2,3})
			Y = tGeneralizedFlagVariety("A",4,{3})
			f = tFlagMap(X1,Y); --the projection of Fl(1,3;5) onto Gr(3,5)
			g = tFlagMap(X2,Y); --the projection of Fl(2,3;5) onto Gr(3,5)
			try(f ** g)

	SeeAlso
		tMap
		tFlagMap
		(compose, TMap, TMap)
///



doc ///
	Key
		(compose, TMap, TMap)
	Headline
		computes the composition of two T-equvariant morphisms
	Usage
		compose(f,g)
	Inputs
		f:TMap
		g:TMap
	Outputs
		:TMap
			the composition of f and g
	Description
		Text
			This method computes the composition of two $T$-equivariant morphisms.

		Example
		    	R = ZZ[T_0..T_3]
		    	X = tGeneralizedFlagVariety("A",3,{1,2,3},R)
			Y = tGeneralizedFlagVariety("A",3,{2,3},R)
			Z = tGeneralizedFlagVariety("A",3,{3},R)
			f = tFlagMap(X,Y); --the projection of Fl(1,2,3;4) onto Gr(2,3;4)
			g = tFlagMap(Y,Z); --the projection of Fl(2,3;3) onto Gr(3;4)
			try compose(g,f)


	SeeAlso
	    	(symbol **, TMap, TMap)
		tMap
		tFlagMap
		
///





doc ///
	Key
		tMap
		(tMap, TVariety, TVariety, List)
	Headline
		creates a TMap
	Usage
		f = tMap(X,Y,L)
	Inputs
		X:TVariety
			the source T-variety of the map
		Y:TVariety
			the target T-variety of the map
		L:List
			of pairs (x,y) where x and y are members of @TT "X.points"@ and @TT "Y.points"@, respectively.
	Outputs
		f:TMap
	Description
		Text
			This method creates a @TO TMap@ given a T-Variety $X$, a T-variety $Y$,
			and a list @TT "L"@ of pairs (x,y) where x and y are members of @TT "X.points"@ and @TT "Y.points"@ (respectively),
			indicating that the T-fixed point x of X is sent to the T-fixed point y of Y under the map.

		Text
			The following describes the projection from the third Hizerbruch surface to the projective 
			line.
			
		Example
		    	R = ZZ[T_0,T_1]
			F3 = tVariety(hirzebruchSurface 3,R);
			PP1 = tProjectiveSpace(1,R);
			L = {({0,1},set {0}), ({0,3}, set{0}), ({1,2}, set{1}), ({2,3}, set{1})};
			f = tMap(F3,PP1,L)
	Caveat
		This does not check that the morphism is well defined. In particular, it does not
		verify that the map on T-fixed points is induced by a morphism of T-varieties.						
	SeeAlso
		diagonalTMap
		tFlagMap
		(pullback, TMap)
		pushforward
		tChi
///

doc ///
	Key
		tFlagMap
		(tFlagMap, TVariety, TVariety)
	Headline
		creates TMaps arising from forgetful maps bewteen generalized flag varieties
	Usage
		f = tFlagMap(X,Y)
	Inputs
		X:TVariety
			the source T-variety which is a generalized flag variety
		Y:TVariety
			the target T-variety which is a generalized flag variety
	Outputs
		f:TMap
	Description
		Text
			Given two generalized flag vareities X = Fl(k_1,..,k_m;n) and Y = Fl(k_r,..,k_m;n) of the same
			lie type, this method produced the canonical projection from $X$ to $Y$.
		Example
		    	FlOG = tGeneralizedFlagVariety("B",3,{1,2})
			OGr17 = tGeneralizedFlagVariety("B",3,{1})
			peek tFlagMap(FlOG,OGr17)
	SeeAlso
		tFlagMap
		diagonalTMap
		tGeneralizedFlagVariety

///




doc ///
	Key
		(pullback, TMap)
	Headline
		computes the pullback map of T-equivariant K-classes of a T-map
	Usage
		pullback(f)
	Inputs
		f:TMap
	Outputs
		:FunctionClosure
			whose input is a TKClass on the target T-variety of f and output is its pullback along f
	Description
		Text
			Given two $T$-varieties $X$ and $Y$, this method computes the pullback of a @TO TKClass@ on $Y$ 
			along a T-equivariant morphism $X \to Y$.

		Example
		    	R = ZZ[T_0..T_3];
			FlGr = tGeneralizedFlagVariety("A",3,{1,2},R);
			Gr24 = tGeneralizedFlagVariety("A",3,{2},R);
			f = tFlagMap(FlGr,Gr24);
			O1 = ampleTKClass Gr24;
			try (pullback f)(O1)
			
	SeeAlso
		tFlagMap
		pushforward
///


doc ///
	Key
		pushforward
		(pushforward, TMap)
	Headline
		computes the pushforward map of T-equivariant K-classes of a T-map
	Usage
		pushforward(f)
	Inputs
		f:TMap
	Outputs
		:FunctionClosure
			whose input is a TKClass on the source T-variety of f and output is its pushforward along f
	Description
		Text
			Given two $T$-varieties $X$ and $Y$, this method computes the pushforward of a @TO TKClass@ on $X$ 
			along a T-equivariant morphism $X \to Y$.

		Example
		    	R = ZZ[T_0..T_3];
			FlGr = tGeneralizedFlagVariety("A",3,{1,2},R);
			Gr24 = tGeneralizedFlagVariety("A",3,{2},R);
			f = tFlagMap(FlGr,Gr24);
			O1 = ampleTKClass FlGr;
			try (pushforward f)(O1)

	SeeAlso
		tFlagMap
		(pullback, TMap)
		pushforward
///

doc ///
	Key
		tChi
		(tChi, TKClass)
	Headline
		computes the T-equivariant Euler characteristic of a T-equivariant K-class
	Usage
		tChi C
	Inputs
		C:TKClass
	Outputs
		:RingElement
			in the character ring of the torus of the T-variety on which C is defined
	Description
		Text
			This method computes the pushforward of a @TO TKClass@ on a @TO TVariety@ $X$ along the structure map 
			$X \to pt$, where $pt$ is a point with trivial $T$-action.
			
		Example
			PP3 = tProjectiveSpace 3
			O1 = ampleTKClass PP3
			tChi O1

///




doc ///
	Key
		tProjectiveSpace
		(tProjectiveSpace, ZZ)
		(tProjectiveSpace, ZZ, Ring) 
	Headline
		constructs Projective space as a GKM manifold
	Usage
		tProjectiveSpace n
		tProjectiveSpace(n,R)
	Inputs
		n:ZZ
		R:Ring
	Outputs
		:TVariety	
	Description
		Text
			Given an integer $n$ this method constructs the n-dimensional Projective space, $\mathbb P^n$, as a $T$-variety. The action
			of $(\mathbb C^*)^{n+1}$ on $\mathbb P^n$ is defined by 
			$(t_0, \ldots, t_n) \cdot (x_0, \ldots, x_n) = (t_0^{-1}x_0, \ldots, t_n^{-1}x_n)$.

		Example
			PP4 = tProjectiveSpace 4;
			peek PP4
	SeeAlso
		tFlagMap
		pushforward
///



doc ///
	Key
		makeCharRing
		(makeCharRing, ZZ) 
	Headline
		constructs the character ring of a torus
	Usage
		makeCharRing n
	Inputs
		n:ZZ
	Outputs
		:Ring	
	Description
		Text
			Given an integer n, this method outputs the character ring of T = $(\mathbb C^*)^n$.

		Example
			makeCharRing 4

///



-----------------------------------------------------------------------------------
------------------------------------< TESTS >--------------------------------------
-----------------------------------------------------------------------------------
TEST ///
PP3 = tProjectiveSpace 3
assert (3 == #gens PP3.charRing)
///
