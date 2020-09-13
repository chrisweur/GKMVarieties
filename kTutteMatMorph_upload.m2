-*----------------------------------------------------------------------------------------------
M2 codes for "K-theoretic Tutte polynomials of morphisms of matroids"

Authors: Christopher Eur, Tim Seynnaeve
Last update: 3/28/2020

-----------------------------------------------------------------------------------------------*-

needsPackage "Matroids"

--load "Polyhedra" if using the "kTutteEhrhart" command
--needsPackage "Polyhedra"

----------------------------------< some matroid functions >--------------------------------------

--shortcut for uniform matroid typing
U = (r,n) -> uniformMatroid(r,n)


--tests whether N is a matroid quotient of M;
isQuot = method();
isQuot(Matroid,Matroid) := Boolean => (N,M) -> isSubset(flats N, flats M)

--truncation (principal truncation by the ground set) of a matroid
truncation = method();
truncation(Matroid) := Matroid => M -> (
    B := independentSets(M,rank M -1);
    matroid(elements M.groundSet, B/elements)
)


----------------------------------< ordinary flag matroids >--------------------------------------

FlagMatroid = new Type of HashTable

flagMatroid = method();

--Given a list L of concordant matroids (M_1, ... , M_k) returns the flag matroid object
--Convention: M_i is quotient of M_(i+1)
--Does not check concordance of matroids in the list, but checks same cardinality of ground set.
flagMatroid(List) := FlagMatroid => L -> (
    E := (first L).groundSet; n := #E;
    if any(L, m -> not #m.groundSet == n) then << "ground set not all same size" << return error;
    new FlagMatroid from {
	symbol groundSet => E,
	symbol constituents => L,
	cache => new CacheTable
    }
)

--Given l x n matrix A and a list r of ranks, outputs the associated flag matroid
flagMatroid(Matrix,List) := FlagMatroid => (A,r) -> (
    k := #r; l := numcols A;
    if r == {} or any(k-1, i -> r_i > r_(i+1)) or r_(k-1) > l then << "check rank sequence" << return error;
    ML := apply(k, i -> matroid A^(apply(r_i, j -> j)) );
    flagMatroid ML
)

--checks that the constituents of a flag matroid M are concordant
isWellDefined(FlagMatroid) := Boolean => M -> (
    L := M.constituents; k := #L;
    all(k-1, i -> isQuot(L_i,L_(i+1))) 
)

--for a flag matroid F, outputs a list of chains of bases of the constituent matroids
bases(FlagMatroid) := List => M -> (
    if not M.cache.?bases then M.cache.bases = (
	ML := M.constituents; k := #ML;
	BL := flatten apply(k, i -> (bases ML_i)/(b -> (i,b)) );
	rel := (a,b) -> first a < first b and isSubset(last a, last b);
	P := poset(BL, rel , AntisymmetryStrategy => "none");
	(maximalChains P)/(c -> c/last)
    );
    M.cache.bases
)


--computes the lattice points of the base polytope of a flag matroid M
latticePts = method();
latticePts(FlagMatroid) := M -> (
    n := #M.groundSet;
    BL := M.constituents/bases/(B -> B/(b -> setIndicator(b,n)));
    k := #BL;
    unique apply((k:0)..toSequence(BL/(B -> #B-1)), i -> sum(k, j -> (BL_j)_(i_j)))
)


--direct sum, deletion, restriction, and contraction of a flag matroid is done by doing
--respective operation on the operation on each constituents
FlagMatroid ++ FlagMatroid := FlagMatroid => (M,N) -> (
    ML := M.constituents; NL := N.constituents;
    if not #ML == #NL then << "two flag matroids have different number of constituents" << return error;
    flagMatroid apply(#ML, i -> ML_i ++ NL_i)
)

restriction(FlagMatroid,Set) := FlagMatroid => (M,S) -> (
    flagMatroid(M.constituents/(m -> restriction(m,S)))
)

restriction(FlagMatroid,List) := FlagMatroid => (M,L) -> (
    flagMatroid(M.constituents/(m -> restriction(m,L)))
)

contraction(FlagMatroid,Set) := FlagMatroids => (M,S) -> (
    flagMatroid(M.constituents/(m -> contraction(m,S)))
)

contraction(FlagMatroid,List) := FlagMatroids => (M,L) -> (
    flagMatroid(M.constituents/(m -> contraction(m,L)))
)

deletion(FlagMatroid,Set) := FlagMatroid => (M,S) -> (
    flagMatroid(M.constituents/(m -> deletion(m,S)))
)

deletion(FlagMatroid,List) := FlagMatroid => (M,L) -> (
    flagMatroid(M.constituents/(m -> deletion(m,L)))
)

FlagMatroid | Set := (M,S) -> restriction(M,S)
FlagMatroid \ Set := (M,S) -> deletion(M,S)
FlagMatroid / Set := (M,S) -> contraction(M,S)
FlagMatroid | List := (M,L) -> restriction(M,L)
FlagMatroid \ List := (M,L) -> deletion(M,L)
FlagMatroid / List := (M,L) -> contraction(M,L)


dual(FlagMatroid) := FlagMatroid => {} >> opts -> M -> (
    ML := M.constituents;
    flagMatroid apply(reverse ML, m -> dual m)
)

truncation(FlagMatroid) := FlagMatroid => M -> (
    flagMatroid(M.constituents/(m -> truncation m))
    )


face(Matroid,Set) := (M,S) -> (M | S) ++ (M / S)

face(FlagMatroid,Set) := (M,S) -> (M | S) ++ (M / S)


--rank (of a subset) in a flag matroid
rank(FlagMatroid) := ZZ => M -> sum(M.constituents/rank)

rank(FlagMatroid,Set) := ZZ => (M,A) -> sum(M.constituents, m -> rank(m,A))


--given matroid (M,N) where M is quotient of N, returns the Higgs lift matroid of M to N
higgsLift = method();
higgsLift(Matroid,Matroid) := Matroid => (M,N) -> (
    E := M.groundSet; r := rank M;
    if not #E == #N.groundSet then << "ground sets are different" << return error;
    if not r < rank N then << "rank of first need be smaller than first" << return error;
    I := independentSets N;
    B := select(subsets(E,r+1), s -> member(s,I) and rank(M,s) == r);
    matroid(toList E,B/toList)
)

--outputs the full flag matroid obtained by Higgs lifts of a matroid M
higgsLift(Matroid) := FlagMatroid => M -> (
    n := #M.groundSet; r := rank M;
    ML1 := {uniformMatroid(0,n)};
    apply(r-1, i -> ML1 = ML1 | {higgsLift(last ML1, M)} );
    ML2 := {M}; U := uniformMatroid(n,n);
    apply(n-r-1, i -> ML2 = ML2 | {higgsLift(last ML2, U)} );
    flagMatroid drop(ML1 | ML2,1)
)

--given a flag matroid F with two constituents (M,N),
--outputs the flag matroid of the higgs factorization
higgsLift(FlagMatroid) := FlagMatroid => F -> (
    if not #F.constituents == 2 then << "need be 2-step flag matroid" << return error;
    M := first F.constituents;
    N := last F.constituents;
    n := #F.groundSet;
    r1 := rank M;
    r2 := rank N;
    ML := {M};
    apply(r2-r1-1, i -> ML = ML | {higgsLift(last ML, N)} );
    flagMatroid (ML | {N})
)


--the Las Vergnas Tutte polynomial of a matroid quotient M1 <<-- M2
LVTutte = method();
LVTutte(Matroid,Matroid) := RingElement => (M1,M2) -> (
    x := symbol x; y := symbol y; z := symbol z;
    E := M1.groundSet; r1 := rank M1; r2 := rank M2;
    R := QQ[x,y,z];
    sum(subsets E, s -> (R_0-1)^(r1-rank(M1,s))*(R_1-1)^(#s-rank(M2,s))*R_2^(r2-r1-rank(M2,s)+rank(M1,s)))
)

LVTutte(FlagMatroid) := RingElement => F -> (
    if #F.constituents != 2 then << "LVTutte only for two-step flag matroids" << return error;
    LVTutte(first F.constituents, last F.constituents)
    )



-----------------------------< T-equivariant K-class computations >-----------------------------

--Given m x n matrix A with integers entries, where each column represents a character,
--outputs the ideal in QQ[x_1..x_n] of the associated AFFINE toric variety
matToIdeal = method();
matToIdeal(Matrix) := Ideal => M -> (
    n := numcols M;
    K := gens ker M; l:= numcols K;
    x := symbol x; R := QQ[x_1..x_n];
    if l == 0 then return ideal(0_R);
    ColtoBinom := C -> product(positions(C, c -> c>0), i -> R_i^(C_i)) - product(positions(C, c -> c<0), i -> R_i^(-C_i));
    I := saturate(ideal(apply(l, j -> ColtoBinom(entries K_j))), product(gens R))
)

--I just didn't want to call MonomialAlgebras package
--given a list L of list of degrees, outputs the associated multigraded ring
monomialRing = method();
monomialRing(List) := Ring => L -> (
    d := #L; x := symbol x;
    R := QQ[x_1..x_d, Degrees => L]
)



--a TVariety X has data of:
--X.points : torus-invariant points
--X.charts : a hash table whose keys are X.points and values are the characters of chart at point
--X.charRing : the character ring of the torus
TVariety = new Type of HashTable


--tVariety created from list L of points and corresponding list M of their charts' characters,
--and R the character ring of the torus
tVariety = method();
tVariety(List,List,Ring) := TVariety => (L,M,R) -> (
    new TVariety from {
	symbol points => L,
	symbol charts => hashTable apply(#L, i -> (L_i,M_i)),
	symbol charRing => R,
	cache => new CacheTable 
    }
)


--product of two TVarieties X,Y
TVariety ** TVariety := TVariety => (X,Y) -> (
    R := X.charRing;
    if not R === Y.charRing then << "character rings need be same" << return error;
    L := toList((set X.points) ** (set Y.points));
    M := apply(L, l -> (X.charts)#(first l) | (Y.charts)#(last l));
    XxY := tVariety(L,M,R)
)



--for setting up a character lattice ring: given n outputs the character ring of T^n
makeCharRing = method();
makeCharRing(ZZ) := Ring => n -> (
    T := symbol T;
    ZZ[T_0..T_(n-1), MonomialOrder => GLex, Inverses=>true]
)



--auxiliary function: makes a Laurent polynomial f into an element of X.charRing
toCharRing = (X,f) -> (
    if f == 1 then return 1_(X.charRing);
    R := ring f; n := #(gens R);
    if not #(gens X.charRing) == n then
        << "different number of variables" << return error;
    sub(f,apply(n, i -> R_i=>X.charRing_i))    
)

--Given a list L of characters of torus of X, outputs the denominator in the Hilbert series
--of the associated monomial subalgebra
tHilbDenom = method();
tHilbDenom(TVariety,List) := RingElement => (X,L) -> (
    R := X.charRing;
    product(L, l -> 1 - product(#l, i -> R_i^(l_i)))
)

--Given a list L of characters of torus of X, outputs the numerator in the Hilbert series
--of the associated monomial subalgebra
tHilbNumer = method();
tHilbNumer(TVariety,List) := RingElement => (X,L) -> (
    if L == {} then return toCharRing(X,1);
    R := monomialRing(L);
    I := matToIdeal transpose matrix L;
    f := value numerator hilbertSeries (map(R,ring I,gens R))I;
    toCharRing(X,f)
)





--a TKClass C has data of:
--C.tvar = a TVariety X that the T-equiv K-class C lives on
--C.hilb = a hash table whose keys are points of X and values are the Hilbert series at the point
TKClass = new Type of HashTable

--a TKClass is given by a TVariety X and a list L of Laurent polynomials for each torus-invariant
--points in X, (listed in the order of X.points)
tKClass = method();
tKClass(TVariety,List) := TKClass => (X,L) -> (
    K := X.points;
    L = L/(f -> toCharRing(X,f));
    new TKClass from {
	symbol tvar => X,
	symbol hilb => hashTable apply(#K, i -> (K_i,L_i))
    }
)

--the trivial TKClass (where X^T --> R is a constant 1 function) of a TVariety X
--in other words, the TKClass of the structure sheaf of X
trivialTKClass = method();
trivialTKClass(TVariety) := TKClass => X -> (
    L := apply(X.points, p -> 1);
    tKClass(X,L)
)


--multiplying two TKClasses
TKClass * TKClass := (C1,C2) -> (
    X1 := C1.tvar; X2 := C2.tvar;
    if not X1 === X2 then << "different varieties" << return error;
    L := apply(X1.points, p -> C1.hilb#p * C2.hilb#p);
    tKClass(X1,L)
)

TKClass ^ ZZ := (C,d) -> (
    if d > 0 then return product(d, i -> C)
    else if d == 0 then return trivialTKClass C.tvar
    else if d < 0 then (
	L := apply(C.tvar.points, p -> (C.hilb#p)^(-1));
	Cneg := tKClass(C.tvar,L);
	return product(-d, i -> Cneg)
	)
    )


--adding two TKClasses
TKClass + TKClass := (C1,C2) -> (
    X1 := C1.tvar; X2 := C2.tvar;
    if not X1 === X2 then << "different varieties" << return error;
    L := apply(X1.points, p -> C1.hilb#p + C2.hilb#p);
    tKClass(X1,L)
)

--if a TVariety X has a distinguished O(1) T-equivariant line bundle, then returns its TKClass
ampleTKClass = method();
ampleTKClass(TVariety) := TKClass => X -> (
    if X.cache.?ampleTKClass then return X.cache.ampleTKClass
    else << "no distinguished ample line bundle on this T-variety" << return error;
)


--a TMap f has data of:
--f.source: the source tVariety X, f.target: the target tVariety Y
--f.ptsMap: a hash table whose keys are X.points and values are the point Y.points that it maps to 
TMap = new Type of HashTable

--a TMap is given by providing the source X and target Y and a list L of pairs (X point, Y point)
tMap = method();
tMap(TVariety,TVariety,List) := TMap => (X,Y,L) -> (
    if not X.charRing === Y.charRing then << "character rings need be same" << return error;
    new TMap from {
	symbol source => X,
	symbol target => Y,
	symbol ptsMap => hashTable L,
	cache => new CacheTable	
    }
)

--pullback map of TKClasses given a TMap
pullback = method();
pullback(TMap) := FunctionClosure => phi -> (
    X := phi.source; Y := phi.target;
    if not phi.cache.?pullback then phi.cache.pullback = C -> (
	if not C.tvar === Y then << "t-K-Class not of target" << return error;
	L := apply(X.points, p -> C.hilb#((phi.ptsMap)#p));
	tKClass(X,L)
    );
    phi.cache.pullback
)


--default Macaulay2 can't do fraction fields for Laurent rings
--but by shifting so that we exit Laurent rings and then going back, we can.
--Here given two Laurent polynomials f,g and a ring S with Inverses=>false, outputs the ratio
--of f and g by multiplying both of them with a big enough monomial to make them polynomials
--the output is the ratio in the new ring S, and a function to go back to the ring of f.
toFraction = method();
toFraction(RingElement,RingElement,Ring) := (f,g,S) -> (
    R := ring f;
    if not R === ring g then << "the two polynomials live in different rings" << return error;
    if not (#gens R) == (#gens S) then << "the temporary ring is no good" << return error;
    Exps := -((transpose (exponents f | exponents g))/min);
    mult := product(#Exps, i -> R_i^(Exps_i));
    numer := sub(mult*f, apply(#(gens R), i -> R_i=>S_i));
    denom := sub(mult*g, apply(#(gens R), i -> R_i=>S_i));
    rat := numer / denom;
    goBack := ratFct -> (
	N := numerator ratFct;
	D := denominator ratFct;
	if not #(terms D) == 1 then << "denominator not a monomial" << return error;
	sub(N, apply(#(gens S), i -> S_i=>R_i))*(sub(D, apply(#(gens S), i -> S_i=>R_i)))^(-1)
	);
    {rat, goBack}
)


--given a TMap, computes the pushforward map as a function
pushforward = method();
pushforward(TMap) := FunctionClosure => phi -> (
    if phi.cache.?pushforward then return phi.cache.pushforward;
    X := phi.source; Y := phi.target; R := X.charRing;
    Ydenominators := hashTable apply(Y.points, q -> (q, tHilbDenom(Y,(Y.charts)#q)));
    Xdenominators := hashTable apply(X.points, p -> (p, tHilbDenom(X,(X.charts)#p)));
    x := symbol x; S := QQ[x_0..x_(#(gens R)-1)];
    pushforwardFct := C -> (
	if not C.tvar === X then << "TKClass not in source" << return error;
	L := apply(Y.points, q -> (
	    preimages := select(X.points, p -> (phi.ptsMap)#p === q);
	    if #preimages == 0 then return 0_R;
	    Ydenom := Ydenominators#q;
	    toSum := apply(preimages, p -> toFraction(Ydenom * C.hilb#p, Xdenominators#p, S));
	    val := (last first toSum) sum(toSum/first)
	    )
	);
	tKClass(Y,L)
    );
    phi.cache.pushforward = pushforwardFct
)


--given a TVariety X outputs the diagonal map X -> X x X
diagonalTMap = method();
diagonalTMap(TVariety) := TMap => X -> (
    Y := X ** X;
    ptPairs := apply(X.points, p -> (
	Q := first select(Y.points, q -> (first q) === (last q) and (first q) === p);
	(p,Q)
	)
    );
    tMap(X,Y,ptPairs)
)

--given two maps, takes Cartesian product of them
TMap ** TMap := TMap => (phi,psi) -> (
    X := phi.source ** psi.source;
    Y := phi.target ** psi.target;
    ptPairs := apply(X.points, p -> (
	Q := ((phi.ptsMap)#(first p),(psi.ptsMap)#(last p));
	(p,Q)
	)
    );
    tMap(X,Y,ptPairs)
)

--composition f o g of two TMaps g: X --> Y, f: Y --> Z;
compose(TMap,TMap) := TMap => (f,g) -> (
    if not f.source === g.target then << "did you mean two switch the order?" << return error;
    Y := f.target; X := g.source;
    ptPairs := apply(X.points, p -> (p, f.ptsMap#(g.ptsMap#p)));
    tMap(X,Y,ptPairs)
)


------------------< Fink-Speyer formuation of the K-theoretic Tutte polynomial >--------------------

--MOAR auxiliary functions:
--given an indicator vector for a flag of subsets converts into a list of sets
vecToFlag = v -> sort apply(max v, k -> set (positions(v, i -> i > k)))

--given a subset s of [n] = {0..(n-1)}, outputs the indicator vector for s
setIndicator = (s,n) -> apply(n, i -> if member(i,s) then 1 else 0)

--converts a list of subsets into the indicator vector
flagToVec = (F,n) -> sum(F/(s -> setIndicator(s,n)))

--swaps the (l_0,l_1)th places in the vector v
swap = (v,l) -> apply(#v, i -> 
    if i == first l then v_(last l)
    else if i == last l then v_(first l)
    else v_i
)

--converts a pair s into a length n list where s_0 has -1 and s_1 has 1 and 0 everywhere else
swapIndicator = (s,n) -> (
    apply(n, i -> 
	if i == first s then -1
	else if i == last s then 1
	else 0
    )
)


--Given a list K and a number n defining flag variety Fl(K;n), returns the TVariety where
--the action of T on the vector space kk^n is (t.v) = (t^-1v)
tFlagVariety = method();
tFlagVariety(List,ZZ,Ring) := TVariety => (K,n,R) -> (
    if not #(gens R) == n then << "check character ring" << return error;
    if not max K < n then << "check rank sequence K" << return error;
    E = set toList(0..(n-1));
    pts := (unique permutations sum(K/(k -> toList(k:1) | toList(n-k:0))))/vecToFlag;
    chrts := apply(pts, p -> (toList sum(p, l -> l**(E-l)))/(s -> swapIndicator(s,n)));
    X := tVariety(pts,chrts,R);
    L := apply(X.points, p -> (Exps := flagToVec(p,n); product(#Exps, i -> R_i^(Exps_i))));
    X.cache.ampleTKClass = tKClass(X,L);
    X
)


--if one does not want to make the charRing beforehand.
--NOT recommended for use...
tFlagVariety(List,ZZ) := TVariety => (K,n) -> tFlagVariety(K,n,makeCharRing n)



--Given two lists Kso, Kta (for rank sequences of source and target flag varieties) and n,
--creates the two TVarieties and makes a TMap between them (given the charRing R)
--corresponding to the "forgetting appropriate linear subspaces map"
tFlagMap = method();
tFlagMap(List,List,ZZ,Ring) := TMap => (Kso,Kta,n,R) -> (
    Flso := tFlagVariety(Kso,n,R);
    Flta := tFlagVariety(Kta,n,R);
    ptPairs := apply(Flso.points, p -> (p,first select(Flta.points, q -> isSubset(q,p))));
    tMap(Flso,Flta,ptPairs)
)

--if X and Y are already tFlagVarieties, then the following outputs the corresponding tFlagMap
-- X --> Y
tFlagMap(TVariety,TVariety) := TMap => (X,Y) -> (
    if not X.charRing === Y.charRing then << "character ring need be same" << return error;
    ptPairs := apply(X.points, p -> (p,first select(Y.points, q -> isSubset(q,p))));
    tMap(X,Y,ptPairs)
)


--given a flag matroid M, returns the TKClass of its 'torus-orbit' in the flag-variety X
tKClass(TVariety,FlagMatroid) := TKClass => (X,M) -> (
    E := M.groundSet;
    K := M.constituents/rank;
    R := X.charRing;
    if not (#(gens R) == #E and (first X.points)/(s -> #s) == K) then 
    	<< "wrong flag variety for the flag matroid" << return error;
    B := bases M;
    L := apply(X.points, p -> (
	if not member(p,B) then return 0_R;
	vp := flagToVec(p,#E);
	rays := select(X.charts#p, r -> member(vecToFlag(swap(vp,positions(r, i -> not i == 0))),B));
	nonrays := select(X.charts#p, r -> not member(r,rays));
	ConeP := tHilbNumer(X,rays);
	if #nonrays == 0 then ConeP else ConeP * tHilbDenom(X,nonrays)
	)
    );
    tKClass(X,L)    
)

--sets up the varieties involved in the Fourier-Mukai themed push-pull diagram for the
--flag geometric Tutte polynomial of a flag matroid.  The charRing R should be given.
--Given a list K and an integer n, sets up Fl(K;n) <-f- Fl(1,K,n-1;n) -g-> Fl(1;n) x Fl(n-1;n)
fourierMukai = method();
fourierMukai(List,ZZ,Ring) := List => (K,n,R) -> (
    FlK := tGeneralizedFlagVariety("A",n-1,K,R);
    Fl1Kn1 := tGeneralizedFlagVariety("A",n-1,unique ({1}|K|{n-1}), R);
    Fl1 := tGeneralizedFlagVariety("A",n-1,{1},R);
    Fln1 := tGeneralizedFlagVariety("A",n-1,{n-1},R);
    piK := tFlagMap(Fl1Kn1,FlK);
    f := tFlagMap(Fl1Kn1,Fl1); g := tFlagMap(Fl1Kn1,Fln1);
    pi1n1 := compose(f ** g, diagonalTMap Fl1Kn1);
    --<< "{{Fl(K;n), Fl(1,K,n-1;n),  Fl(1;n) x Fl(n-1;n)},{pi_K, pi_(1(n-1))}}" <<
    {{FlK, Fl1Kn1, pi1n1.target}, {piK, pi1n1}}
)

fourierMukai(List,ZZ) := List => (K,n) -> fourierMukai(K,n, makeCharRing n)


--these are auxiliary functions for converting T-equivariant class to in terms of alpha,beta

--takes in a tKClass C of Fl(1;n) x Fl(n-1;n) and outputs the matrix whose (i,j)th entry
--is the hilb value at the point {set{i}, [n]\set{j}}
toMatrix = C -> (
    R := C.tvar.charRing;
    n := #(gens R);
    E := set toList(0..(n-1));
    matrix apply(n, i -> apply(n, j -> (
        if i == j then return 0_R;
	C.hilb#(({set{i}}, {E - set{j}}))
	))
    )
)

aa = (i0,i,j,R) -> (R_i-R_(i0))*(R_i^(-1));
bb = (j0,i,j,R) -> (R_(j0) - R_j)*(R_(j0)^(-1));

equiToNonEquiStep = (i0,j0,X) -> (
    R := ring X_(0,0); n := #(gens R);
    t := symbol t; S := QQ[t_0..t_(n-1)]; 
    denom := sub(product(i0, k -> aa(k,i0,j0,R))*product(j0, k -> bb(k,i0,j0,R)),R);
    c := first toFraction(X_(i0,j0), denom,S);
    cVal := sub(c, apply(gens ring c, r -> r=>1));
    XX := matrix apply(n, i -> apply(n, j -> (
	numer := X_(i0,j0) * product(i0, k -> aa(k,i,j,R)) * product(j0, k -> bb(k,i,j,R));
	ratio := toFraction(numer,denom,S);
	(last ratio) (first toFraction(X_(i,j),1_R,S) - first ratio)
	))
    );
    {cVal,XX}
)

--Given a TKClass in P^(n-1) x P^(n-1) = Gr(n-1;n) x Gr(1;n),
--outputs the polynomial in representing its K-class
--where x,y are the structure sheaves of the two hyperplanes
toPolynomial = method();
toPolynomial(TKClass) := Matrix => C -> (
    T := toMatrix C;
    n := numcols T;
    TList := {T};
    M := matrix apply(n, i -> apply(n, j -> (
    	out := equiToNonEquiStep(i,j,TList_(-1));
	TList = append(TList, last out);
	first out
	))
    );
    x := symbol x; y := symbol y;
    S := ZZ[x,y];
    sum(n, i -> sum(n, j -> M_(i,j) * S_0^j * S_1^i))
)


--given a flag matroid M, outputs the flag-geometric Tutte polynomial of M
kTutte = method();
kTutte(FlagMatroid) := RingElement => M -> (
    if not M.cache.?kTutte then M.cache.kTutte = (
	n := #M.groundSet;
    	R := makeCharRing n;
    	K := M.constituents/rank;
    	FM := fourierMukai(K,n,R);
    	FlK := first first FM;
    	f := first last FM; g := last last FM;
    	yM := tKClass(FlK,M);
    	YM := (pushforward g)( (pullback f)(yM * (ampleTKClass FlK)) );
    	toPolynomial YM
    );
    M.cache.kTutte
)

--computes the bivariate polynomial representing the class of the pushforward of yM.O(d) in P x P
-- for a flag matroid M and any integer (can be negative) d.
kTutteEhrhartValue = method();
kTutteEhrhartValue(FlagMatroid,ZZ) := RingElement => (M,d) -> (
    if M.cache.?kTutteEhrhart then return (
	T := M.cache.kTutteEhrhart;
	S := ring T;
	sub(T, {S_0 => S_0, S_1 => S_1, S_2 => d})
	);
    n := #M.groundSet;
    R := makeCharRing n;
    K := M.constituents/rank;
    FM := fourierMukai(K,n,R);
    FlK := first first FM;
    f := first last FM; g := last last FM;
    yM := tKClass(FlK,M);
    YMd := (pushforward g)( (pullback f)(yM * (ampleTKClass FlK)^d) );
    toPolynomial YMd
)


--auxiliary code for the kTutteEhrhart function
vandermonde = d -> (
    sub(matrix apply(d, i -> apply(d, j -> i^j)),QQ)
)

--the i-th elementary symmetric function of a list L
elemSym = (L,i) -> sum(subsets(L,i), j -> product j)

--given a polynomial P and a variable x of P, outputs a polynomial P written in terms of the
--binomial polynomials in x.  The output is the vector (in column) of coordinates where the
--basis is { (x C d), (x C d-1), ... (x C 0) }
toProjective = (P,x) -> (
    R := ring P;
    d := degree(x,P);
    transMat := transpose matrix reverse apply(d+1, i -> apply(d+1, j -> 1/(i!) * elemSym(toList(-(i-1)..0),(j-d+i))));
    sub(inverse transMat,R) * (last coefficients(P, Variables => x))
)


--computes the (x,y,m)-multivariate k-Tutte valued Ehrhart polynomial of flag matroid M
--requires package "Polyhedra"
kTutteEhrhart = method(Options => {Projective => false});
kTutteEhrhart(FlagMatroid) := RingElement => opts -> M -> (
    if not M.cache.?kTutteEhrhart then (
    	n := #M.groundSet;
    	R := makeCharRing n;
    	K := M.constituents/rank;
    	FM := fourierMukai(K,n,R);
    	d := dim convexHull transpose matrix apply(bases M, b -> flagToVec(b,n));
	x := symbol x; y := symbol y; m := symbol m;
	S := QQ[x,y,m];
	Vmat := sub(vandermonde(d+1),S);
	coeff := Vmat^(-1) * transpose matrix {apply(d+1, i -> kTutteEhrhartValue(M,FM,i,S))};
	exps := matrix {apply(d+1, i -> S_2^i)};
	M.cache.kTutteEhrhart = (exps * coeff)_(0,0)
	);
    f := M.cache.kTutteEhrhart;
    if not opts.Projective then f else toProjective(f, (ring f)_2)
)


--often we don't want to initialize the fourierMukai multiple times (makes things slow!)
--one can give it as an input:
kTutte(FlagMatroid,List) := RingElement => (M,FM) -> (
    if not M.cache.?kTutte then M.cache.kTutte = (
    	FlK := first first FM;
    	f := first last FM; g := last last FM;
    	yM := tKClass(FlK,M);
    	YM := (pushforward g)( (pullback f)(yM * (ampleTKClass FlK)) );
    	toPolynomial YM
    );
    M.cache.kTutte
)

--same thing as above; don't want to initialize fourierMukai all the time
--actually, maybe I should cache the fourierMukai?? TODO later...
kTutteEhrhartValue(FlagMatroid,List,ZZ) := RingElement => (M,FM,d) -> (
    if M.cache.?kTutteEhrhart then return (
	T := M.cache.kTutteEhrhart;
	S := ring T;
	sub(T, {S_0 => S_0, S_1 => S_1, S_2 => d})
	);
    FlK := first first FM;
    f := first last FM; g := last last FM;
    yM := tKClass(FlK,M);
    YMd := (pushforward g)( (pullback f)(yM * (ampleTKClass FlK)^d) );
    toPolynomial YMd
)

--computes the kTutte polynomial and puts the element in a given ring R
kTutte(FlagMatroid,Ring) := RingElement => (M,R) -> (
    t := kTutte M; S := ring t;
    sub(t, {S_0 => R_0, S_1=>R_1})
)

kTutteEhrhartValue(FlagMatroid,ZZ,Ring) := RingElement => (M,d,R) -> (
    t := kTutteEhrhartValue(M,d); S := ring t;
    sub(t, {S_0 => R_0, S_1=>R_1})
)    

--computes the kTutte polynomial and puts the element in a given ring R
kTutte(FlagMatroid,List,Ring) := RingElement => (M,FM,R) -> (
    t := kTutte(M,FM); S := ring t;
    sub(t, {S_0 => R_0, S_1=>R_1})
)

kTutteEhrhartValue(FlagMatroid,List,ZZ,Ring) := RingElement => (M,FM,d,R) -> (
    t := kTutteEhrhartValue(M,FM,d); S := ring t;
    sub(t, {S_0 => R_0, S_1=>R_1})
)

--evaluates the kTutte of a flag matroid M at values (in a ring) x = a, y = b
kTutteEvaluate = method();
kTutteEvaluate(FlagMatroid,RingElement,RingElement) := ZZ => (M,a,b) -> (
    t := kTutte M; S := ring t;
    sub(t, {S_0 => a, S_1 => b})
)

kTutteEvaluate(FlagMatroid,ZZ,ZZ) := ZZ => (M,a,b) -> (
    t := kTutte M; S := ring t;
    sub(t, {S_0 => a, S_1 => b})
)

--gives the flag-geometric K-theoretic characteristic polynomial of a flag matroid M
kCharPol = method(Options => {Reduced => false});
kCharPol(FlagMatroid) := RingElement => opts -> M -> (
    T := kTutte M;
    R := ring T;
    P := (-1)^(rank last M.constituents) * sub(T, {R_0 => 1-R_0, R_1 => 0});
    if not opts.Reduced then return P;
    sub(P / (R_0 - 1),R)
)



------------------------< T-equivariant K-theoretic Tutte polynomials >----------------------------

--GENERAL WARNING: codes here can be buggy if the flag matroid has U(0,n) or U(n,n) constituents

--Let L be a list of list of integers representing monomials in a Laurant polynomial
-- ring S; the associated toric ring R is multigraded as it is a subalgebra of S.
-- Output below is the multigraded hilbert series of R. 
hilbertSeries(List) := opts -> L -> (
    A := monomialRing L;
    I := if L == {} then ideal(0_A) else matToIdeal transpose matrix L;
    hilbertSeries (map(A,ring I,gens A))I
)

--computes the T-equivariant version of the K-theoretic Tutte polynomial of a flag matroid
--informed by the "fundamental computation"
tKTutte = method(Options => {inUV => false});
tKTutte(FlagMatroid) := RingElement => opts -> M -> (
    if not M.cache.?tKTutte then M.cache.tKTutte = (
    	E := M.groundSet;
    	B := bases M;
    	dtop := rank last M.constituents;
    	t := symbol t;
    	S := QQ[t_0..t_(#E-1)];
    	u := symbol u; v := symbol v;
    	Suv := (frac S)[u,v];
    	f := sum(B, b -> (
	    	candRays := (toList sum(b, l -> l**(E-l)))/(s -> swapIndicator(s,#E));
	    	vp := flagToVec(b,#E);
	    	realRays := select(candRays, r -> 
		    member(vecToFlag(swap(vp,positions(r, i -> not i == 0))),B));
	    	hilbCone := (
		    if realRays == {} then 1_S else (
		    	H := hilbertSeries realRays;
		    	f := first toFraction(value numerator H, value denominator H, S)
		    	)
		    );
	    	O1part := product(drop(b,-1), i -> product(elements i, j -> S_j));
	    	uPart := sum apply(subsets last b, s -> 
		    sub(product(elements s, i -> S_i),Suv)*Suv_0^(dtop - #s));
	    	vPart := sum apply(subsets (E - first b), s -> 
		    sub(product(elements s, i -> S_i),Suv)*Suv_1^(#s));
	    	sub(hilbCone,Suv)*sub(O1part,Suv)*uPart*vPart
	    	)
	    )
    	);
    tk := M.cache.tKTutte;
    if opts.inUV then tk else sub(tk,QQ[gens ring tk][gens coefficientRing ring tk])
)

fastKTutte = method(Options => {inUV => false});
fastKTutte(FlagMatroid) := RingElement => opts -> M -> (
    if not M.cache.?tKTutte then M.cache.tKTutte = (
    	E := M.groundSet;
    	B := bases M;
    	dtop := rank last M.constituents;
    	t := symbol t;
    	S := QQ[t_0..t_(#E-1)];
	g := gens S;
	Ering := frac(QQ[e]);
	ee :=(gens Ering)#0;
    	u := symbol u; v := symbol v;
	Euv := Ering[u,v];
    	f := sum(B, b -> (
	    	candRays := (toList sum(b, l -> l**(E-l)))/(s -> swapIndicator(s,#E));
	    	vp := flagToVec(b,#E);
	    	realRays := select(candRays, r -> 
		    member(vecToFlag(swap(vp,positions(r, i -> not i == 0))),B));
	    	hilbCone := (
		    if realRays == {} then 1_S else (
		    	H := hilbertSeries realRays;
		    	f := first toFraction(value numerator H, value denominator H, S);
			g1 = gens ring f;
			subf := sub(f,apply(length g1,i->(g1#i=>1+2^i*ee)))
		    	)
		    );
	    	O1part := product(drop(b,-1), i -> product(elements i, j -> S_j));
		subO1part := sub(O1part,apply(length g,i->(g#i=>1+2^i*ee)));
		uPart := sum apply(subsets last b, s -> 
		    sub((product(elements s, i -> S_i))_S,apply(length g,i->(g#i=>1+2^i*ee)))*Euv_0^(dtop - #s));
		vPart := sum apply(subsets (E - first b), s -> 
		    sub((product(elements s, i -> S_i))_S,apply(length g,i->(g#i=>1+2^i*ee)))*Euv_1^(#s));
	    	sub(hilbCone,Euv)*subO1part*uPart*vPart
	    	)
	    );
	sub(f,{ee=>0})
    	);
    tk := M.cache.tKTutte;
    if opts.inUV then tk else sub(tk,QQ[gens ring tk][gens coefficientRing ring tk])
)

--converts tKTutte into x,y by setting u = x - 1 and v = y - 1
tKTutteTXY = method();
tKTutteTXY(FlagMatroid) := RingElement => M -> (
    f := tKTutte(M, inUV => true);
    R := ring f;
    S := (coefficientRing R)[symbol x, symbol y];
    sub(f, apply(#(gens R), i -> R_i => S_i - 1))   
)

--pushes down to non-equivariant K-theory (i.e. evaluates the t_i's to 1)
--and then converts into x,y by setting u = x - 1 and v = y - 1
tKTutteXY = method();
tKTutteXY(FlagMatroid) := RingElement => M -> (
    f := tKTutte(M);
    f = sub(f, apply(gens ring f, r -> r => 1));
    R := ring f;
    S := QQ[symbol x, symbol y];
    sub(f, apply(#(gens R), i -> R_i => S_i - 1))   
)


--K-theoretic approach to the Las Vergnas Tutte polynomial
--outputs the T-equivariant Las Vergnas Tutte polynomial of a 2-step flag matroid (M1, M2)
KLVTutte = method(Options => {inUVW => false});
KLVTutte(FlagMatroid) := RingElement => opts -> M -> (
    if #F.constituents > 2 then << "must be TWO-step flag matroid" << return error;
    M1 := first M. constituents;
    M2 := last M.constituents;
    if not F.cache.?KLVTutte then F.cache.KLVTutte = (
    	E := M.groundSet;
    	B := bases M;
    	r1 := rank M1;
    	r2 := rank M2;
    	t := symbol t;
    	S := QQ[t_0..t_(#E-1)];
    	u := symbol u; v := symbol v; w := symbol w;
    	Suvw := (frac S)[u,v,w];
    	f := sum(B, b -> (
	    	candRays := (toList sum(b, l -> l**(E-l)))/(s -> swapIndicator(s,#E));
	    	vp := flagToVec(b,#E);
	    	realRays := select(candRays, r -> 
		    member(vecToFlag(swap(vp,positions(r, i -> not i == 0))),B));
	    	hilbCone := (
		    if realRays == {} then 1_S else (
		    	H := hilbertSeries realRays;
		    	f := first toFraction(value numerator H, value denominator H, S)
		    	)
		    );
	    	uPart := sum apply(subsets first b, s -> 
		    sub(product(elements s, i -> S_i),Suvw)*Suvw_0^(r1 - #s));
	    	vPart := sum apply(subsets (E - last b), s -> 
		    sub(product(elements s, i -> S_i),Suvw)*Suvw_1^(#s));
	    	wPart := sum apply(subsets (last b - first b), s ->
		    sub(product(elements s, i -> S_i),Suvw)*Suvw_2^(r2 - r1 -#s));
	    	sub(hilbCone,Suvw)*uPart*vPart*wPart
	    	)
	    )
	);
    lvt := F.cache.KLVTutte;
    if opts.inUVW then lvt else sub(lvt,QQ[gens ring lvt][gens coefficientRing ring lvt])
)

KLVTutteTXYZ = method();
KLVTutteTXYZ(FlagMatroid) := RingElement => F -> (
    lvt := KLVTutte(F, inUVW => true);
    R := ring lvt;
    S := (coefficientRing R)[symbol x, symbol y, symbol z];
    sub(lvt, {R_0 => S_0 - 1, R_1 => S_1 -1, R_2 => S_2})
)

KLVTutteXYZ = method();
KLVTutteXYZ(FlagMatroid) := RingElement => F -> (
    lvt := KLVTutte F;
    lvt = sub(lvt, apply(gens ring lvt, r -> r => 1));
    R := ring lvt;
    S := (coefficientRing R)[symbol x, symbol y, symbol z];
    sub(lvt, {R_0 => S_0 - 1, R_1 => S_1 -1, R_2 => S_2})
)


--the "H-polynomial" (cf. section 6 of pub version Fink-Speyer)
--this one is through the "Las Vergnas path"
--NOTE: divide everything at the end by the product t_0t_1..t_(n-1) to get the true answer
kHpoly = method();
kHpoly(FlagMatroid) := RingElement => M -> (
    E := M.groundSet;
    B := bases M;
    r1 := rank first M.constituents;
    r2 := rank last M.constituents;
    t := symbol t;
    S := QQ[t_0..t_(#E-1)];
    u := symbol u; v := symbol v; w := symbol w;
    Suvw := (frac S)[u,v,w];
    f := sum(B, b -> (
	    candRays := (toList sum(b, l -> l**(E-l)))/(s -> swapIndicator(s,#E));
	    vp := flagToVec(b,#E);
	    realRays := select(candRays, r -> 
		member(vecToFlag(swap(vp,positions(r, i -> not i == 0))),B));
	    hilbCone := (
		if realRays == {} then 1_S else (
		    H := hilbertSeries realRays;
		    f := first toFraction(value numerator H, value denominator H, S)
		    )
		);
	    uPart := sum apply(subsets first b, s -> 
		    sub(product(elements s, i -> S_i),Suvw)*Suvw_0^(r1 - #s));
	    vPart := sum apply(subsets (E - last b), s -> 
		    sub(product(elements s, i -> S_i),Suvw)*Suvw_1^(#s));
	    wPart := sum apply(subsets (last b - first b), s ->
		    sub(product(elements s, i -> S_i),Suvw)*Suvw_2^(r2 - r1 -#s));
	    sub(hilbCone,Suvw)*uPart*vPart*product(elements (E - last b), i -> S_i)
	    )
	);
    sub(f,QQ[gens Suvw][gens S])
)

--the "H-polynomial" (cf. section 6 of pub version Fink-Speyer)
--this one is through the "flag-geometric path"
--NOTE: divide everything at the end by the product t_0t_1..t_(n-1) to get the true answer
kHpoly2 = method();
kHpoly2(FlagMatroid) := RingElement => M -> (
    E := M.groundSet;
    B := bases M;
    dtop := rank last M.constituents;
    t := symbol t;
    S := QQ[t_0..t_(#E-1)];
    u := symbol u; v := symbol v;
    Suv := (frac S)[u,v];
    f := sum(B, b -> (
	    candRays := (toList sum(b, l -> l**(E-l)))/(s -> swapIndicator(s,#E));
	    vp := flagToVec(b,#E);
	    realRays := select(candRays, r -> 
		member(vecToFlag(swap(vp,positions(r, i -> not i == 0))),B));
	    hilbCone := (
		if realRays == {} then 1_S else (
		    H := hilbertSeries realRays;
		    f := first toFraction(value numerator H, value denominator H, S)
		    )
		);
	    uPart := sum apply(subsets last b, s -> 
		    sub(product(elements s, i -> S_i),Suv)*Suv_0^(dtop - #s));
	    vPart := sum apply(subsets (E - first b), s -> 
		    sub(product(elements s, i -> S_i),Suv)*Suv_1^(#s));
	    sub(hilbCone,Suv)*uPart*vPart*product(elements (E - last b), i -> S_i)
	    )
	);
    sub(f,QQ[gens Suv][gens S])
    )

--the "H-polynomial" defined by pullback and pushing forward y(MM) to P x P
--(through the more natural Fl(1,r_1,..,r_k,n-1;n) construction
--NOTE: should be the same as kHpoly2, but faster (and we forget the T-equivariant-ness)
kHPolynomial = method(Options => {inUV => false});
kHPolynomial(FlagMatroid) := RingElement => opts -> M -> (
    if not M.cache.?kHPolynomial then M.cache.kHPolynomial = (
	n := #M.groundSet;
    	R := makeCharRing n;
    	K := M.constituents/rank;
    	FM := fourierMukai(K,n,R);
    	FlK := first first FM;
    	f := first last FM; g := last last FM;
    	yM := tKClass(FlK,M);
    	YM := (pushforward g)( (pullback f)(yM));
    	toPolynomial YM
    );
    KHP := M.cache.kHPolynomial;
    if opts.inUV then (
	S := QQ[symbol u, symbol v];
	return (map(S, ring KHP, apply(2, i -> (ring KHP)_i => S_i + 1)))KHP
	);
    M.cache.kHPolynomial
)

end

------------------------------------------------------------------------------------------------
--------------------------------------------- END ----------------------------------------------


-----------------------------< Tests for T-equivariant K-theory >--------------------------------

restart
load "kTutteMatMorph_upload.m2"

R = makeCharRing(3)

X = tGeneralizedFlagVariety("A",2,{1,2},R)
C = ampleTKClass X
C.tvar === X

P2a = tGeneralizedFlagVariety("A",2,{1},R);
P2b = tGeneralizedFlagVariety("A",2,{2},R);
P2a ** P2b

f = tFlagMap(X,P2a)
g = tFlagMap(X,P2b)
h = f ** g
d = diagonalTMap(X)
--compose(h,d)


---------------------------< Tests for the flag-geometric Tutte polynomial >--------------------------

restart
load "kTutteMatMorph_upload.m2"

--uniform(1,3) matroid
R = makeCharRing 3;
FM = fourierMukai({1},3,R);
F1 = first first FM
f = (FM_1)_0;
g = (FM_1)_1;

M = flagMatroid({uniformMatroid(1,3)})
yM = tKClass(F1,M)

Y = time (pushforward(g))( (pullback(f))(yM * (ampleTKClass F1)))
toMatrix Y
toPolynomial Y

kTutte M

--complete graph on 4 vertices example
M = flagMatroid {matroid completeGraph 4}
time FM = fourierMukai({3},6); -- 3.5 seconds (yikes!)
time kTutte(M,FM) -- 27 seconds (ugh...)


--a robust check with the matroids
restart
load "kTutteMatMorph_upload.m2"
ML = drop(drop(allMatroids 4,1),-1)
TML = apply(ML, m -> time {tuttePolynomial m, kTutte flagMatroid({m})}) -- 1 sec per matroid
all(TML, l -> (map(ring first l, ring last l, gens ring first l))(last l) == first l)

FM = time fourierMukai({3},5);
ML = select(allMatroids 5, m -> rank m == 3)
TML = apply(ML, m -> time {tuttePolynomial m, kTutte(flagMatroid {m},FM)}) -- 3 sec per matroid
all(TML, l -> (map(ring first l, ring last l, gens ring first l))(last l) == first l)


--the first example in CDMS18
R = makeCharRing 3;
FM = fourierMukai({1,2},3,R);
F12 = (first FM)_1
PP = last first FM
pi12 = last last FM
N = flagMatroid(matrix{{1,1,1},{1,0,0}},{1,2})
yN = tKClass(F12,N)

EhrN = m -> toPolynomial (pushforward pi12)(yN * product(m, i -> ampleTKClass F12))
EhrN 1
EhrN 2

kTutte N
tKTutteXY N

--the second example in CDMS18
N = flagMatroid({uniformMatroid(2,5), uniformMatroid(3,5)})
FM = time fourierMukai({2,3},5);
time kTutte(N,FM) -- 7 seconds
time tKTutteXY N -- 48 seconds

--A test for the faster algorithm ("K-theory to cohomology")
restart
load "kTutteMatMorph_upload.m2"
--M=flagMatroid({uniformMatroid(2,5), uniformMatroid(3,5)})
time f=fastKTutte(flagMatroid {uniformMatroid(2,5), uniformMatroid(3,5)} ) --1.23875 seconds
time f=tKTutte(flagMatroid {uniformMatroid(2,5), uniformMatroid(3,5)} ) --59.9602 seconds


--all Fl(1,2;3) examples
M1 = flagMatroid({uniformMatroid(1,3), uniformMatroid(2,3)})
M2 = flagMatroid({matroid({0,1,2}, {{0},{2}}), uniformMatroid(2,3)})
M3 = flagMatroid({uniformMatroid(1,3), matroid({0,1,2}, {{0,1},{0,2}})})
M4 = flagMatroid({matroid({0,1,2},{{0}}), matroid({0,1,2}, {{0,1},{0,2}})})
M5 = flagMatroid({matroid({0,1,2},{{1},{2}}), matroid({0,1,2}, {{0,1},{0,2}})})
M6 = flagMatroid({matroid({0,1,2},{{0},{2}}), matroid({0,1,2}, {{0,2}})})
M7 = flagMatroid({matroid({0,1,2},{{0}}), matroid({0,1,2},{{0,2}})})

L = {M1,M2,M3,M4,M5,M6,M7};
L/isWellDefined
L/kTutte




----------------------< K-theoretic Las Vergnas Tutte polynomial checks >------------------------------
restart
load "kTutteMatMorph_upload.m2"

ML = select(allMatroids 6, m -> rank m == 4 and (#loops m + #coloops m) == 0);
M = first random ML;
e = set (random elements M.groundSet)_{0,1}
F = flagMatroid{M/e,M\e}
(F.constituents)/(m -> {loops m, coloops m})
time KLVTutte(F, inUVW => true)
KLVTutteXYZ F
LVTutte F


--------------------------------< "H-polynomial" experiments >------------------------------
restart
load "kTutteMatMorph_upload.m2"

--kHpoly is from the Las-Vergnas diagram, whereas kHpoly2 is from the flag-geometric diagram

F = flagMatroid {U(1,3),U(2,3)}
kHpoly F
kHpoly2 F

F = flagMatroid {U(1,4),U(3,4)}
kHpoly F
kHpoly2 F

--for kHpoly (from the Las-Vergnas diagram), need not be a polynomial in uv
F = flagMatroid {U(2,4),U(3,4)}
kHpoly F
kHpoly2 F
