restart
-- We can use this reference for the statements on T-fixed points and orbits : https://hal.archives-ouvertes.fr/hal-01882469/document 
-- We need to modify it for the case of the flags  
 
-- Creates the list of all 0,1 vectors of a fixed size
vecList = k -> (
    if k <= 0 then return {};
    L = {{0},{1}};
    i = 1;
    while i < k do  (
	L = join(apply(L, v -> append(v,0)), apply(L,v -> append(v,1)));
	i = i+1);
    return L
) 

--- Takes two lists, L1 and L2 of equal size, and outputs the list whose i-th element is (L1#i, L2#i)     
pairwiseList = (L1,L2) -> (
    if not (#L1 == #L2) then return print "Error";
    a = #L1;
    return apply(a, v -> (L1#v, L2#v))
)

--- Computes the character of the torus action for each type of matrix entry.
pmFunction = (L,n) -> (
    apply(L_0, w-> (
      apply(n, i -> (
	      if (i == w_0 and L_1 == 1) then -1 
	      else if (i == w_0 and L_1 == 2) then -1
	      else if (i == w_0 and L_1 == 3) then 1
	      else if (i == w_0 and L_1 == 4) then 1
	      else 0)) - 
      apply(n, j -> (
	      if (j == w_1 and L_1 == 1) then -1 
	      else if (j == w_1 and L_1 == 2) then 1
	      else if (j == w_1 and L_1 == 3) then -1
	      else if (j == w_1 and L_1 == 4) then 1
	      else 0))
      )
  )
)
  

-- Let K = {k_1,...,k_s} and consider a T-fixed point P = (V_1,..,V_s) in flag-SGr(K,2n). Consider the usual
-- affine neighbourhood of P in the standard flag variety. This corresponds to a sequence of matrices of variables
-- (M_1,..,M_s) in row reduced echelon form with compatibility conditions. To obtain the neighbourhood of P in
-- the flag-Symplectic grassmannian we need to cut out by the conditions that make V_s isotropic.
-- This amounts to dropping some of the coordinates of M_s (if we index the columns of M_s by 1,..,n,1*,..,n* we
-- are essentially taking the "lower triangular" portion" of the columns corresponding to 1*,..,n*). 
--The coordinates of the M_i descend to coordiantes in the Symplectic grassmannian.

-- The function partialTchart takes in V_{i-1}, V_{i} and computes M_{i-1} for i < s.
-- The function topTchart computes the chart around V_s in SGr(K,2n)

partialTchart = (A,B,n) -> (
    p = # select(toList A, v -> v <= n-1);
    A' = sort toList A;
    rowCol = toList(set toList(0..#A-1) ** (B-A));
    L1 = apply(select(rowCol, v -> isSubset({v}, toList(set toList(0..p-1) ** set toList(0..n-1)))),   v-> (A'#(v_0),v_1));
    L2 = apply(select(rowCol, v -> isSubset({v}, toList(set toList(0..p-1) ** set toList(n..2*n-1)))), v-> (A'#(v_0),v_1-n));
    L3 = apply(select(rowCol, v -> isSubset({v}, toList(set toList(p..k-1) ** set toList(0..n-1)))),   v-> (A'#(v_0)-n,v_1));
    L4 = apply(select(rowCol, v -> isSubset({v}, toList(set toList(p..k-1) ** set toList(n..2*n-1)))), v-> (A'#(v_0)-n,v_1-n) );
    return join(pmFunction({L1,1},n), pmFunction({L2,2},n), pmFunction({L3,3},n), pmFunction({L4,4},n))    
    )
 
 
topTchart = (A,n) -> (
    I0' = sort select(toList A, v -> v<= n-1);
    I0 = sort apply(select(toList A, v -> v<= n-1), v-> v+n);
    p =  #I0;
    I1' = sort apply(select(toList A, v -> v>= n), v-> v-n);
    I1 = sort select(toList A, v -> v >= n);
    I' = join(I0',I1');
    I = join(I0,I1);
    k = #I;
    A11 = flatten apply(I0, v -> apply(toList(0..position(I, w-> w == v)), z -> (z,v)));
    A12 = flatten apply(I1, v -> apply(toList(0..position(I, w-> w == v)), z -> (z,v-n)));
    A1 = join(A11,A12);
    A2 = toList ((set toList(0..k-1) ) ** (set toList(0..n-1) - set I')); 
    A3 = toList ((set toList(0..k-1)) ** (set toList(n..2*n-1) - set I) );
    L = join(A1,A2,A3);
    L1 = apply(select(L, v -> isSubset({v}, toList(set toList(0..p-1) ** set toList(0..n-1)))),   v-> (I0'#(v_0),v_1));
    L2 = apply(select(L, v -> isSubset({v}, toList(set toList(0..p-1) ** set toList(n..2*n-1)))), v-> (I0'#(v_0),v_1-n));
    L3 = apply(select(L, v -> isSubset({v}, toList(set toList(p..k-1) ** set toList(0..n-1)))),   v-> (I1'#(v_0-p),v_1));
    L4 = apply(select(L, v -> isSubset({v}, toList(set toList(p..k-1) ** set toList(n..2*n-1)))), v-> (I1'#(v_0-p),v_1-n));        
    return join(pmFunction({L1,1},n), pmFunction({L2,2},n), pmFunction({L3,3},n), pmFunction({L4,4},n))    
   )





-- Given a list K = {k_1,...,k_s}, in weakly increasing order, and a number n 
-- defining the symplectic flag variety flag-Sgr(k_1,...,k_s;2n),
-- returns the TVariety where the action of T on the vector space kk^{2n} is given by
-- (t_1,...,t_n).(v_1,...,v_{2n}) = (t_1^(-1)v_1,...t_n^(-1)v_n,t_1v_{n+1},...,t_nv_{2n})
-- Since M2 doesn't like negative powers, we just a polynomial ring in n-variables
tSymplecticFlagVariety = method();
tSymplecticFlagVariety(List,ZZ,Ring) := TVariety => (K,n,R) -> (
    if not #(gens R) == n then << "check character ring" << return error;
    --if not max K < n then << "check rank sequence K" << return error;
    vecL = vecList(K#-1);
    topPtsSigned = flatten apply(subsets(n,K#-1), v -> apply(vecL, w -> pairwiseList(v,w)));
    --- Computes the fixed points of SGr(k_s;n)
    topPts = apply(apply(topPtsSigned, v -> apply(v, w -> (first w) + n*(last w))), v -> {set v});
    --- Computes the fixed points of the desired flag
    i = #K-2;
    pts = topPts;
    while i >=0 do (
	pts = flatten apply(pts, v -> apply(subsets(toList v_0,K_i), w -> flatten join({set w}, {v})));
	i = i - 1);
    --- Computes the characters of the affine chart:
    chrts = {};
    for point in pts do (
    	--- First computet the chart of SGr(k_s;2n);
	L = topTchart(point#-1,n);
	j = #K-1;
	while j >=1 do (
	    L = join(L,partialTchart(point#(j-1),point#j,n));
	    j = j -1);
	chrts = join(chrts,{L}));
    X := tVariety(pts,chrts,R);
    -- The class of O(1)
    L := apply(X.points, p -> (
	    Exps := flagToVec(p,2*n); 
	    product(#Exps, i -> if i < n then R_i^(Exps_i) else ((frac R)_(i-n))^(-Exps_i)))
	);
    -- The following doesn't work. tKClass doesn't like negative powers. We can just 
    -- modify tKClass later.
     X.cache.ampleTKClass = tKClass(X,L);
     X
)



-- Test
tSymplecticFlagVariety({1},2,QQ[T_1,T_2])
tSymplecticFlagVariety({1,1},2,QQ[T_1,T_2])
tSymplecticFlagVariety({2},3,QQ[T_1,T_2,T_3])
tSymplecticFlagVariety({1,2},3,QQ[T_1,T_2,T_3])


-- Comment:
-- I think O(1) in the case of Fl(1,1,n) (the standard flag variety) is wrong or I might just be confused?
-- Fl(1,1,n) is isomorphic to P^n, but the code for O(1) outputs {t_0^2,...,t_n^2}
-- instead of outputting {t_0,..,t_n}. Since my class of O(1) is just the restriction
-- of that, I think it's also incorrect for the symplectic flag. For ex:
((tFlagVariety({1,1},2,QQ[T_1,T_2])).cache.ampleTKClass).hilb
-- versus
((tFlagVariety({1},2,QQ[T_1,T_2])).cache.ampleTKClass).hilb

-- The final class of O(1) is incorrect in the symplecitc case anyway. The class is 
-- fine as a map Gr^T(K,2n) -> Z[char T] (on line 124). But tKClass doesn't work as intended
-- when it converts it to the numerator of the Hilbert series.
((tSymplecticFlagVariety({1},2,QQ[T_1,T_2])).cache.ampleTKClass).hilb


------------------------------------------------------------------------------------ 
--Code for 1-dimensional T-orbits
------------------------------------------------------------------------------------ 

-- Note: In my comments I interchange between using {1,..,n,1*,..,n*} and {0,..,2n-1} to denote [n].

-- 1) Standard flag variety.

-- Base case of Gr(k,2n)
-- A pair of T-fixed points (I,J) lie in the boundary of a 1-dim orbit iff
-- the following is true:
--  (a) I - J  = {i} and J - I = {j}

-- More generally a pair ({I_1,..,I_s},{J_1,..,J_s}) is the boundary of a 1-dim orbit of the
-- flag variety X = Fl({k_1,..,k_s},n) iff there exists i,j such that for all l
-- either I_l = J_l or the pair (I_l,J_l) satisfies condition (a) above

-- The set of edges of the moment graph 
fixedPts = (K,n) -> (unique permutations sum(K/(k -> toList(k:1) | toList(n-k:0))))/vecToFlag



boundaryPts = (K,n) -> (
    flatten apply(fixedPts(K,n), w -> (
	    ptsAtInfinity := apply(toList(w#-1 ** (set toList(0..n-1) - w#-1)), v -> (
		    apply(w, u ->  if isSubset(set{v_0},u) then u - set{v_0} + set{v_1} else u))
		);
	    apply(ptsAtInfinity, v -> (w,v))
	    )
	)
    )


-- The key is a pair (I,J) such that I - J = {i} and J - I = {j}. The value is the 
-- e_i - e_j correspnding to the character of the orbit, t_j^{-1}t_i. Alternatively,
-- we can just create a list of triples (I,J, char)...   
orbits = (K,n) -> hashTable apply(boundaryPts(K,n), u -> (
	d := first select(#(u_0), v -> not u_0_v === u_1_v);
	(u,setIndicator(u_0_d - u_1_d,n) - setIndicator(u_1_d - u_0_d,n)))
    )

--- Test
netList boundaryPts({1},4)
orbits({1},4)
netList boundaryPts({1,2},4)
orbits({1,2},4)

--- 2) Symplectic variety

-- Base case of SGr(k,2n)
-- A pair of T-fixed points (I,J) lie in the boundary of a 1-dim orbit iff
-- the following is true:
-- (a) I - J  = {i} and J - I = {j} 
-- (b) I - J = {i,j} and J - I = {n+i,n+j} and vice-versa (with i,j <= n)


-- More generally a pair ({I_1,..,I_s},{J_1,..,J_s}) is the boundary of a 1-dim orbit of the
-- flag variety flag-SGr({k_1,..,k_s},2n) iff there exists i,j such that for all l
-- either I_l = J_l or the pair (I_l,J_l) satisfies condition (a) or (b) above.
-- Note the following can happen does lie in a 1-dim boundary: ({3}, {0,1,2}, {5}, {1,5,6}) of SGr(3,6)




-- Computing the Torus fixed points of the flag symplectic grassmannian, flag-SGr(K,2n)
-- (This has already been used inside tSymplecticFlagVariety)
fixedPts = (K,n) -> (
    vecL = vecList(K#-1);
    topPtsSigned = flatten apply(subsets(n,K#-1), v -> apply(vecL, w -> pairwiseList(v,w)));
    -- Computes the fixed points of SG(k_s;n) first:
    topPts = apply(apply(topPtsSigned, v -> apply(v, w -> (first w) + n*(last w))), v -> {set v});
    i = #K-2;
    pts = topPts;
    while i >=0 do (
	pts = flatten apply(pts, v -> apply(subsets(toList v_0,K_i), w -> flatten join({set w}, {v})));
	i = i - 1);
    return pts
    )



-- Given a subset A of {1,..,n,1*,..,n*} returns the "bar" of A i.e. apply * to every element A.
bar = (A,n) -> set apply(toList A, v -> if v < n then v + n else v-n)

pairBar = (A,n) -> apply(toList A, v -> if v < n then (v,v + n) else (v,v-n))

-- Set of edges of the moment graph
boundaryPts = (K,n) -> (
    flatten apply(fixedPts(K,n), w -> (
	    -- Satisfying condition (a) (need to modify the flag code to make sure it's still admissible)
	    ptsAtInfinity1 := apply(join(toList(w#-1 ** (set toList(0..2*n-1) - w#-1 - bar(w#-1,n))),pairBar(w#-1,n)), v -> (
		    apply(w, u ->  if isSubset(set{v_0},u) then u - set{v_0} + set{v_1} else u))
		);
	    -- Satisfying condition (b)
	    ptsAtInfinity2 = apply(subsets(w#-1,2), u -> apply(#w, a -> (
			    p = toList(u);
			    if  u * w_a === u then w_a - u + bar(u,n) 
			    else if  u * w_a === set{p_0} then w_a - set{p_0} + bar(set{p_0},n)
			    else if  u * w_a === set{p_1} then w_a - set{p_1} + bar(set{p_1},n)
			    else w_a)
			)
		    );
	    apply(join(ptsAtInfinity1,ptsAtInfinity2), v -> (w,v))
	    )
	)
    )



-- Hash table of characters (this can be computed from the affine chart)
-- Follow easily from the code in tSymplecticFlagVariety



-- Test
netList fixedPts({1},2)
netList boundaryPts({1},2)


netList boundaryPts({1,2},2)
netList boundaryPts({2},3)





-- Prototype code (once we fix the data type of moment graph)
-- Hasn't been tested yet.
-- Assumings X.orbits is a list of triples (I,J, character) where I is 0 and J is infty
isWellDefined(TKClass) := Boolean => C -> (
    X = C.tvar;
    Hilb = C.hilb; 
    all(X.orbits, v -> (
	    ind0 := position(X.pts, w -> w == v_0);
	    ind1 := position(X.pts, w -> w == v_1);
	    use X.charRing;
	    div0 := product apply(select(#v_2, w -> w > 0), u -> R_u);
	    div1 := product apply(select(#v_2, w -> w < 0), u -> R_u);
	    (Hilb#indo0 - Hilb#ind1) % (div0 - div1)  == 0)
	)
    )
	    
	    










--------------------------------------------------------
-- Code for all isotoripic flag varieties at once (Haven't tested it yet)
----------------------------------------------------



-- Generalized version of topTchart defined above 
-- m = 1,2,3,4 corresponds to Type A, B, C,D respectively 
topTchart = (A,n,m) -> (
    I0' = sort select(toList A, v -> v<= n-1);
    I0 = sort apply(select(toList A, v -> v<= n-1), v-> v+n);
    p =  #I0;
    I1' = sort apply(select(toList A, v -> v>= n), v-> v-n);
    I1 = sort select(toList A, v -> v >= n);
    I' = join(I0',I1');
    I = join(I0,I1);
    k = #I;
    if (m ==1 or m ==3) then (
	A11 = flatten apply(I0, v -> apply(toList(0..position(I, w-> w == v)-1), z -> (z,v)));
	A12 = flatten apply(I1, v -> apply(toList(0..position(I, w-> w == v)-1), z -> (z,v-n))) ) else (
	A11 = flatten apply(I0, v -> apply(toList(0..position(I, w-> w == v)), z -> (z,v)));
	A12 = flatten apply(I1, v -> apply(toList(0..position(I, w-> w == v)), z -> (z,v-n))) );	
    coordinates1 = join(A11,A12);
    coordinates2 = toList ((set toList(0..k-1) ) ** (set toList(0..n-1) - set I')); 
    coordinates3 = toList ((set toList(0..k-1)) ** (set toList(n..2*n-1) - set I) );
    L = join(coordinates1,coordinates2,coordinates3);
    L1 = apply(select(L, v -> isSubset({v}, toList(set toList(0..p-1) ** set toList(0..n-1)))),   v-> (I0'#(v_0),v_1));
    L2 = apply(select(L, v -> isSubset({v}, toList(set toList(0..p-1) ** set toList(n..2*n-1)))), v-> (I0'#(v_0),v_1-n));
    L3 = apply(select(L, v -> isSubset({v}, toList(set toList(p..k-1) ** set toList(0..n-1)))),   v-> (I1'#(v_0-p),v_1));
    L4 = apply(select(L, v -> isSubset({v}, toList(set toList(p..k-1) ** set toList(n..2*n-1)))), v-> (I1'#(v_0-p),v_1-n));
    -- There are a few more characters when m = 1 i.e. in the case of OG(k,2n+1).
    L5 = if m == 1 then apply(toList A, v -> 
	(if v <= n-1 then toList(v:0) | {-1} | toList(n-1-v:0) else toList(v-n:0) | {1} | toList(2*n-1-v:0) )) else {};
    return join(pmFunction({L1,1},n), pmFunction({L2,2},n), pmFunction({L3,3},n), pmFunction({L4,4},n),L5)    
   )



tIsotropicFlagVariety = method();
tIsotropicFlagVariety(List,ZZ,Ring,Integer) := TVariety => (K,n,R,m) -> (
    -- Checking presentation of input
    if not #(gens R) == n then << "check character ring" << return error;
    if not all(K, w -> 1<= w and w <=n) then << "check rank sequence K" << return error;
    if not sort(K, DegreeOrder => Ascending) == K then << "rank sequence K not in increasing order" << return error; 
    -- Computing the fixed points
    vecL = vecList(K#-1);
    topPtsSigned = flatten apply(subsets(n,K#-1), v -> apply(vecL, w -> pairwiseList(v,w)));
    topPts = apply(apply(topPtsSigned, v -> apply(v, w -> (first w) + n*(last w))), v -> {set v});
    i = #K-2;
    pts = topPts;
    while i >=0 do (
	pts = flatten apply(pts, v -> apply(subsets(toList v_0,K_i), w -> flatten join({set w}, {v})));
	i = i - 1);
    -- Computing the characters of the affine chart:
    chrts = {};
    for point in pts do (
	L = topTchart(point#-1,n,m);
	j = #K-1;
	while j >=1 do (
	    L = join(L,partialTchart(point#(j-1),point#j,n));
	    j = j -1);
	chrts = join(chrts,{L}));
    X := tVariety(pts,chrts,R);
    -- The class of O(1)
--    L := 
--    assignAmpleTKClass tKClass(X,L);
    X
)








