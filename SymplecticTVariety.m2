restart
-- Symplectic grassmannian    
 
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
  


-- Consider a T-fixed point  with constituents A \subseteq B. The following computes the characters
-- of the affine chart around A \subseteq B corresponding to A (this is explained sloppily)
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
   
   
-- Computes the characters of the affine chart around a T-fixed point, A, in SGr(k;n)   
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
-- defining the symplectic flag variety SFl(k_1,...,k_s;2n),
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
    L := apply(X.points, p -> (Exps := flagToVec(p,2*n); product(#Exps, i -> 
		 if i < n then R_i^(Exps_i) else 1/(R_(i-n)^(Exps_i)))));
    assignAmpleTKClass tKClass(X,L);
    X
)



------------------------------------------------------------------------------------ 
--Code for 1-dimensional T-orbits
------------------------------------------------------------------------------------ 

--- 1) Symplectic variety

-- Base case of SGr(k,2n)
-- A pair of T-fixed points (I,J) lie in the boundary of a 1-dim orbit iff
-- the following is true:
--  i) I - J  = {i} and J - I = {j}
-- ii) I - J = {i,j} and J - I = {k+i,k+j}


-- More generally for a flag FGr({k_1,..,k_s},2n), the projection onto Gr(k_i,2n), denoted Pr_i,
-- is T-equivariant. A pair of T-fixed points (I,J) lie in the boundary of a 1-dim orbit 
-- in the *flag* iff 
-- a) The image, (Pr_i(I), Pr_j(J)) in Gr(k_i,2n) satisfies conditions (i), (ii) above or 
-- the Pr_i(I) = Pr_j(J)
-- b) There is some i for which Pr_i(I) \ne Pr_i(J)



-- Computing the Torus fixed points of the flag symplectic grassmannian, FSGr(K,2n)

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


-- Set of edges of the moment graph
boundaryPts = (K,n) -> (
    flatten apply(fixedPts(K,n), w -> (
	    ptsAtInfinity := select(fixedPts(K,n), v -> 
	    	all(apply(#w, u-> (
			    setDiff := toList(v_u - w_u) - toList(w_u-v_u);
			    (#(v_u * w_u) >= K#u-1) or (#setDiff == 2 and abs(first setDiff) == n))
		    	), b -> b == true));
	    apply(delete(w,ptsAtInfinity), v -> (w,v))
	    )    
    	)
    )

-- Hash table of characters (TBD)


-- Test
netList boundaryPts({1,2},2)




-- 2) Standard flag variety.

-- Base case of Gr(k,2n)
-- A pair of T-fixed points (I,J) lie in the boundary of a 1-dim orbit iff
-- the following is true:
--  i) I - J  = {i} and J - I = {j}

-- The general case is then similar to 1)

-- The set of edges of the moment graph 
fixedPts = (K,n) -> (unique permutations sum(K/(k -> toList(k:1) | toList(n-k:0))))/vecToFlag


boundaryPts = (K,n) -> (
    flatten apply(fixedPts(K,n), w -> (
	    ptsAtInfinity := select(fixedPts(K,n), v -> 
	    	all(apply(#w, u-> #(v_u * w_u) >= K#u-1), b -> b == true));
	    apply(delete(w,ptsAtInfinity), v -> (w,v))
	    )    
    	)
    )



-- The key is a pair (I,J) such that I - J = {i} and J - I = {j}. The value is the 
-- e_i - e_j correspnding to the character of the orbit, t_j^{-1}t_i. This is the character
-- *if* we restrict to the open set P_I \ne 0. 
orbits = hashTable apply(boundaryPts(k,n), u -> (u,{setIndicator(toList(u_0 - u_1),n) - setIndicator(toList(u_1 - u_0),n)}))



--- Test
netList boundaryPts({2,3},4)













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









