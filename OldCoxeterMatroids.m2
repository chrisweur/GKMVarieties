-*----------------------------------------------------------------------------------------------
Coxeter Matroids M2 codes

Author: Christopher Eur
Last update: 8/10/2018

-----------------------------------------------------------------------------------------------*-


needsPackage "Matroids"
needsPackage "Polyhedra"

--outputs the last coefficient of the characteristic polynomial (no sign)
--WARNING: not the beta invariant of the matroid
beta = M -> (f := characteristicPolynomial M; (-1)^(rank M) * sub(f, (ring f)_0 => 0))

--computes what the evaluations of the variables in the G-invariant should be for beta
gToBeta = method();
gToBeta(ZZ,ZZ) := Matrix => (r,n) -> (
    ML := select(allMatroids n, m -> rank m == r);
    rkSeqL := subsets(n,r)/(l -> apply(n, i -> if member(i,l) then 1 else 0));
    A := matrix apply(ML, m -> (g := gInvar m; apply(rkSeqL, l -> if member(l,keys g) then g#l else 0)));
    A = sub(A,QQ); B := sub(transpose matrix{ML/(m -> n! * beta m)},QQ);
    if not ker A == 0 then << "solution not unique";
    Sol := flatten entries transpose solve(A,B);
    hashTable apply(#rkSeqL, i -> (rkSeqL_i, Sol_i))
)

--matroid basis polytope with the lattice chosen to be e_1, e_2, ..., e_n with e_0 := sum e_i's
polytope(Matroid) := Polyhedron => M -> (
    n := #M.groundSet;
    convexHull( (id_(ZZ^(n-1)) | transpose matrix {toList(n-1:-1)}) * basisIndicatorMatrix M)
)

multiTutte = method();
multiTutte(Matroid) := RingElement => M -> (
    q := symbol q; v := symbol v;
    E := elements M.groundSet; Evar := E/(e -> v_e);
    R := frac(QQ[q])[Evar];
    sum(subsets E, s -> 1/((coefficientRing R)_0^(rank(M,set s)))*product(s, i -> R_i))
)



----------------------------------< ordinary flag matroids >--------------------------------------

FlagMatroid = new Type of HashTable

flagMatroid = method();

--Given a list L of concordant matroids (M_1, ... , M_k) returns the flag matroid type.
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
    F := apply(L, m -> flats m);
    all(k-1, i -> isSubset(F_i,F_(i+1)) ) 
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

--flats of a flag matroids are flag of flats of the constituents (i think??)
flats(FlagMatroid) := List => M -> (
    if not M.cache.?flats then M.cache.flats = (
	ML := M.constituents; k := #ML;
	FL := flatten apply(k, i -> (flats ML_i)/(f -> (i,f)) );
	rel := (a,b) -> first a < first b and isSubset(last a, last b);
	P := poset(FL, rel , AntisymmetryStrategy => "none");
	(maximalChains P)/(c -> c/last)
    );
    M.cache.flats
)

latticeOfFlats(FlagMatroid) := Poset => M -> (
    if not M.cache.?latticeOfFlats then M.cache.latticeOfFlats = (
	F := flats M; k := #M.constituents;
	rel := (a,b) -> all(k, i -> isSubset(a_i, b_i));
	poset(F, rel)
    );
    M.cache.latticeOfFlats
)

characteristicPolynomial(FlagMatroid) := RingElement => opts -> M -> characteristicPolynomial latticeOfFlats M;

rank(FlagMatroid,Set) := ZZ => (M,A) -> sum(M.constituents, m -> rank(m,A))

higgsLift = method();
higgsLift(Matroid,Matroid) := Matroid => (M,N) -> (
    E := M.groundSet; r := rank M;
    if not #E == #N.groundSet then << "ground sets are different" << return error;
    if not r < rank N then << "rank of first need be smaller than first" << return error;
    I := independentSets N;
    B := select(subsets(E,r+1), s -> member(s,I));
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

polytope(FlagMatroid) := Polyhedron => M -> sum(M.constituents/polytope);

--The multi-parameter version of the Potts model for flag matroids
--multiplies by global parameter q_1^rk_1...q_k^rk_k to make the exponents positive
multiTutte(FlagMatroid) := RingElement => M -> (
    q := symbol q; v := symbol v;
    E := elements M.groundSet; Evar := E/(e -> v_e);
    L := M.constituents; k := #L; param := apply(k, i -> q_i);
    R := ZZ[param][Evar];
    sum(subsets E, s -> product(k, i -> (coefficientRing R)_i^(rank(L_i) - rank(L_i,set s)))*product(s, i -> R_i))
)

LVTutte = method();
LVTutte(Matroid,Matroid) := RingElement => (M1,M2) -> (
    x := symbol x; y := symbol y; z := symbol z;
    E := M1.groundSet; r1 := rank M1; r2 := rank M2;
    R := QQ[x,y,z];
    sum(subsets E, s -> (R_0-1)^(r1-rank(M1,s))*(R_1-1)^(#s-rank(M2,s))*R_2^(r2-r1-rank(M2,s)+rank(M1,s)))
)

fVec = method();
fVec(Matroid,Matroid) := RingElement => (M,N) -> (
    E := M.groundSet;
    B := select(subsets E, s -> rank(N,s) == #s and rank(M,s) == rank M);
    x := symbol x; R := QQ[x];
    sum(B, b -> R_0^(#b))
)


------------------------------< G-invariant of a polymatroid >--------------------------------
rankSeq := (M,L) -> (
    E := elements M.groundSet;
    apply(#L, i -> rank(M, set apply(i+1, j -> E_(L_j)))- rank(M, set apply(i, j -> E_(L_j))))
)

gInvar = method(Options => {symbol Greedy => false});
gInvar(Matroid) := HashTable => opts -> M -> (
    if not opts.Greedy then return tally (permutations (#M.groundSet))/(s -> rankSeq(M,s));
    if opts.Greedy then (
	tally apply(permutations (#M.groundSet), s -> set s_(elements maxWeightBasis(M,-s)))
	)
)

gInvar(FlagMatroid) := HashTable => opts -> M -> (
    if not opts.Greedy then return tally (permutations (#M.groundSet))/(s -> rankSeq(M,s));
--TODO    if opts.Greedy then (
--	tally apply(permutations (#M.groundSet), s -> set s_(elements maxWeightBasis(M,-s)))
--	)
)






-*-----------------------------------------------------------------------------------------------
Symplectic (flag) matroids (Coxeter matroids of type B)
-----------------------------------------------------------------------------------------------*-


--------------------------------< Signed sets >---------------------------------

SignedSet = new Type of HashTable

net SignedSet := S -> (
	"sgn-set " | toString(flatten elements S.positives | (elements S.negatives)/(i -> toString(i) |"*") )
)


--Given finie disjoint sets S and T, creates a signed set with S cup T as the support and
--S being the positive set and T the negative.
signedSet = method();
signedSet(Set,Set) := SignedSet => (S,T) -> (
    if not S * T === set {} then << "the two sets must be disjoint" << return error;
    new SignedSet from {
	symbol positives => S,
	symbol negatives => T
    }
)

--Given two (disjoint) lists instead
signedSet(List,List) := SignedSet => (S,T) -> signedSet(set S, set T)

--signed set from a indicator vector consisting of 1, 0, or -1
signedSet(Matrix) := SignedSet => C -> (
    if numcols C > 1 then << "input need be a column vector" << return error;
    n := numrows C; C = flatten entries transpose C;
    signedSet(select(n, i -> C_i == 1), select(n, i -> C_i == -1))
)

--Gives the support of a signed set
support(SignedSet) := Set => S -> S.positives + S.negatives

--Tests if a signed set is a subset of another
isSubset(SignedSet, SignedSet) := Boolean => (S,T) -> (
    isSubset(S.positives,T.positives) and isSubset(S.negatives,T.negatives)
)

--star: the involution of the signed set (reverses all the signs)
setStar = method();
setStar(SignedSet) := SignedSet => S -> signedSet(S.negatives, S.positives)

--intersection of two signed (admissible) sets
SignedSet * SignedSet := SignedSet => (S,T) -> signedSet(S.positives * T.positives, S.negatives * T.negatives)

--"union" of two signed (admissible) sets
SignedSet + SignedSet := SignedSet => (S,T) -> (
    pos := S.positives + T.positives; neg := S.negatives + T.negatives;
    signedSet(pos - neg, neg - pos)
)

--member function for a signed set;
--outputs 1 if an element e is a positive member, -1 if negative, 0 otherwise
member(Thing, SignedSet) := ZZ => (e,S) -> (
    if member(e,S.positives) then 1
    else if member(e,S.negatives) then -1
    else 0
)


--outputs all signed subsets whose support is in {0,...,n-1}
signedSubsets = method();
signedSubsets(ZZ) := List => n -> (
    E := set (0..(n-1));
    flatten apply(subsets E, s -> subsets(E-s)/(t -> signedSet(s,t)) )
)

--outputs all signed subsets whose support is of size m and is contained in {0,...,n-1}
signedSubsets(ZZ,ZZ) := List => (n,m) -> (
    E := set(0..(n-1));
    flatten flatten apply(m+1, i -> apply(subsets(E,i), s -> subsets(E-s, m-i)/(t -> signedSet(s,t))))
)





-------------------------------< symplectic matroids >----------------------------------



SpMatroid = new Type of HashTable

net SpMatroid := M -> (
	"a symplectic matroid of rank " | toString(M.rank) | " on n=" | toString(#M.groundSet) | " elements"
)

--Given the size of the ground set E, and a list L of signed subsets
--returns a symplectic matroid
spMatroid = method();
spMatroid(Set,List) := SpMatroid => (E,L) -> (
    if any(L, l -> not isSubset(support l, E)) then << "there's a basis not in ground set" << return error;
    r := #(support first L); if any(L, l -> not # support l == r) then << "bases need be same size" << return error;
    M := new SpMatroid from {
	symbol groundSet => set (0..(#E-1)),
	symbol bases => L/(l -> signedSet(positions(elements E, e -> member(e,l.positives)), positions(elements E,e -> member(e,l.negatives)))),
	symbol rank => r,
	cache => new CacheTable,
    };
    M.cache.rankFunction = new MutableHashTable;
    M
)


bases(SpMatroid) := List => M -> M.bases;

rank(SpMatroid) := ZZ => M -> M.rank;

basisIndicatorMatrix(SpMatroid) := Matrix => M -> (
    E := elements M.groundSet; B := bases M;
    matrix apply(E, i -> apply(B, b -> member(i,b)))
)

polytope(SpMatroid) := Polytope => M -> (
    if not M.cache.?polytope then M.cache.polytope = convexHull basisIndicatorMatrix M;
    M.cache.polytope
)


isBroot := v -> (
    v = entries v; l = select(v, i -> not i == 0);
    #l < 3 and abs(first l) == abs(last l)
)

isWellDefined(SpMatroid) := Boolean => M -> (
    Q := polytope M; V:= vertices Q;
    E := faces(dim Q - 1,Q)/first/(l -> V_(first l) - V_(last l));
    all(E, v -> isBroot v)
)


rank(SpMatroid,SignedSet) := ZZ => (M,S) -> (
    if not isSubset(support S, M.groundSet) then << "need be a signed subset of ground set" << return error;
    if not M.cache.rankFunction#?S then M.cache.rankFunction#S = (
	v := matrix {apply(elements M.groundSet, i -> member(i,S))};
	max flatten entries (v * basisIndicatorMatrix M) 
    );
    M.cache.rankFunction#S
)

--Given a list L of admissble weight function on the 2*(ground set) of a symplectic matroid M,
--returns a basis that is maximal w/r/t to the weight function 
--the algorithm is the greedy algorithm described in BGW ch.3
maxWeightBasis(SpMatroid,List) := SignedSet => (M,L) -> (
    E := M.groundSet; n := #E;
    mE := set ((elements E)/(e -> e+n));
    Esort := reverse sort(elements (E+mE), i -> L_i);
    Etotal := apply(Esort, i -> if i < n then signedSet(set {i}, set {}) else signedSet(set {}, set{i-n}));
    B := signedSet(set{}, set{}); r := rank M; i := -1; P := {};
    while #(support B) < r do (
	i = first select(2*n, j -> j > i and #(support B) +1 == #(support (B+Etotal_j)) and any(bases M, b -> isSubset(B + Etotal_j, b)));
	B = B + Etotal_i;
	P = P | {i};
    );
    B, signedSet(set select(P, j -> j < n), set select(P, j -> j > n-1)/(j -> 2*n-j-1))
)

--given a set S outputs all signed sets with support S
allSigns := S -> apply(subsets S, p -> signedSet(p, S - p))

--outputs all signed permutations of n elements
adWeights := n -> (
    elts := X -> apply(elements support X, e -> member(e,X)*e);
    flatten ( (allSigns set(1..n) )/elts/permutations )
)

--max downward chains given a signed subset
maxChains :=  S -> (
    E := elements support S; n := #E;
    P := permutations E;
    apply(P, p -> apply(n, i -> (
		Ei := p_(toList(0..i));
		signedSet(select(Ei, e -> member(e,S.positives)), select(Ei, e -> member(e,S.negatives)))
	)))
)

--The proposed G-invariant of a symplectic matroid
gInvar(SpMatroid) := HashTable => opts -> M -> (
    H := tally apply(flatten ((allSigns M.groundSet)/maxChains), c -> c/(i -> rank(M,i)));
    if not opts.Greedy then return H;
    n := #M.groundSet;
    WL := (adWeights n)/(w-> w | -w);
    tally apply(WL, w -> last maxWeightBasis(M,-w))
)

--input: a matrix M whose row span is a isotropic subspace
--output: the associated symplectic matroid
spMatroid(Matrix) := SpMatroid => M -> (
    if numcols M % 2 == 1 then << "numcols need be even" << return error else n = (numcols M)//2;
    A := M_(apply(n, i -> i)); B := M_(apply(n, i -> i+n)); r = numrows M;
    if not A * transpose B == B * transpose A then << "row span not isotropic" << return error;
    L := select(subsets(2*n, r), l -> not any(l, i -> member(i+n,l)) and rank M_l == r);
    spMatroid(set toList(0..(n-1)), L/(l -> signedSet(select(l, i -> i < n), select(l, i -> not i < n)/(j -> j-n))))
)

--uniform symplectic matroid of rank r on ground set of size n
uniformSpMatroid = method();
uniformSpMatroid(ZZ,ZZ) := SpMatroid => (r,n) -> (
    E := set toList(0..(n-1));
    L := select(subsets(2*n, r), l -> not any(l, i -> member(i+n,l)) );
    L = L/(l -> signedSet(select(l, i -> i < n), select(l, i -> not i < n)/(j -> j-n)));
    spMatroid(E,L)
)

--random sp matroid in a VERY inefficient way... okay at least up to r = 4 and n = 4...
randomSpMatroid = method();
randomSpMatroid(ZZ,ZZ) := SpMatroid => (r,n) -> (
    E := set (0..(n-1));
    L := signedSubsets(n,r); b := #L;
    M := spMatroid(E, (random L)_(apply(random(1,b),i -> i)) );
    while not isWellDefined M do M = spMatroid(E, (random L)_(apply(random(1,b),i -> i)) );
    M
)

ring(SpMatroid) := Ring => M -> (
    n := #M.groundSet;
    x := symbol x; y := symbol y;
    S := QQ[x_0..x_(n-1),y_0..y_(n-1)];
    Ass := (bases M)/(b -> (elements b.positives)/(i -> S_i) | (elements b.negatives)/(i -> S_(i+n)))/ideal;
    S/(intersect Ass)
)


--Given a symplectic matroid, returns the symmetric matrix associated to the quadric form
--coming from the multivariate Tutte polynomial
HRtest = method();
HRtest(SpMatroid) := Matrix => M -> (
    n := #M.groundSet;
    q := symbol q; S := QQ[q, MonomialOrder => GLex, Inverses=>true];
    nToSS := A -> signedSet(set select(A, i -> i < n),set select(A, i -> i > n-1)/(i -> i-n));
    rkDiff := (i,j) -> -rank(M, nToSS {i,j}) + rank(M, nToSS {i}) + rank(M, nToSS {j});
    A := matrix apply(2*n, i -> apply(2*n, j -> 
	    if i == j or i == j+n or j == i+n then 0
	    else q^(rkDiff(i,j))
	)
    );
    A = A | transpose matrix{toList(2*n:1)};
    A || matrix{toList(2*n:1) | {n/(n-1)}}
)

--output the max and min matroid of a Lagrangian matroid
Mmaxmin = method();
Mmaxmin(SpMatroid) := List => M -> (
    n := #M.groundSet;
    B := basisIndicatorMatrix M;
    dotted := flatten entries (matrix{toList(n:1)}*B);
    small := min dotted; big := max dotted;
    Bmax := select(numcols B, i -> dotted_i == big)/(j -> B_j)/entries;
    Bmin := select(numcols B, i -> dotted_i == small)/(j -> B_j)/entries;
    Mmax := matroid(toList(0..(n-1)), Bmax/(b -> select(n, i -> b_i > 0)));
    Mmin := matroid(toList(0..(n-1)), Bmin/(b -> select(n, i -> b_i > 0)));
    {Mmax, Mmin}
)



------------------------< conormal fan: special Lagrangian matroids >-----------------------------

spMatroid(Matroid) := SpMatroid => M -> (
    E := M.groundSet;
    spMatroid(E, (bases M)/(b -> signedSet(b, E-b)))
)







----------------------------------< flag symplectic matroids >------------------------------------



FlagSpMatroid = new Type of HashTable

--Given a list L of concordant symplectic matroids (M_1, ... , M_k) returns the flag.
--Does not check concordance of matroids in the list, but checks same cardinality of ground set.
flagSpMatroid = method();
flagSpMatroid(List) := FlagSpMatroid => L -> (
    E := (first L).groundSet; n := #E;
    if any(L, m -> not #m.groundSet == n) then << "ground set not all same size" << return error;
    new FlagSpMatroid from {
	symbol groundSet => E,
	symbol constituents => L,
	cache => new CacheTable
    }
)

flagSpMatroid(Matrix,List) := SpMatroid => (A,L) -> (
    flagSpMatroid(apply(L, r -> spMatroid(A^apply(r, i -> i))))
)

truncate(SpMatroid) := SpMatroid => M -> (
    r := rank M; n := #M.groundSet; B := bases M;
    L := select(signedSubsets(n,r-1), s -> any(B, b -> isSubset(s,b)));
    spMatroid(M.groundSet, L)
)


higgsLift(SpMatroid) := FlagSpMatroid => M -> (
    r := rank M; n := #M.groundSet;
    L := {M};
    apply(r-1, i -> L = {truncate first L} | L);
    apply(n-r, i -> L = L | {uniformSpMatroid(r+i+1,n)});
    flagSpMatroid(L)
)

polytope(FlagSpMatroid) := Polyhedron => M -> sum(M.constituents/polytope);


end



----------------------------------------< TESTS >----------------------------------------
restart
load "CoxeterMatroids.m2"

isLogConc = P -> (
    L := reverse flatten entries last coefficients P;
    all(#L-2, i -> L_i*L_(i+2) < L_(i+1)^2)
)

--this is now for n=8
isStrongLogConc = P -> (
    L := reverse flatten entries last coefficients P;
    D := reverse flatten ((terms P)/degree);
    all(#L-2, i -> not L_i*L_(i+2)*(D_(i+1)+1)*(8-D_(i+1)+1)/((D_(i+1))*(8-D_(i+1))) > L_(i+1)^2)
)


M = uniformMatroid(2,5); N = uniformMatroid(3,5);
fVec(M,N)

ML = select(allMatroids 8, m -> rank m == 4);
NL = select(allMatroids 8, m -> rank m == 6);

isQuot = (m,n) -> isSubset(flats m, flats n);

L = select(apply(100, i -> (first random ML, first random NL)), s -> isQuot(first s, last s));
#L
FL = L/fVec;
select(FL, f -> not isStrongLogConc f)



F = flagMatroid({uniformMatroid(2,5),uniformMatroid(3,5)})
T = multiTutte F
C = sub(T,apply(gens ring T, i -> i=>-1))
C' = sub(C, apply(gens ring C, i -> i=> i+1))
factor oo

F = higgsLift matroid({0,1,2},{{0,1},{0,2}})

F = flagMatroid({uniformMatroid(1,3),matroid({0,1,2},{{0,1},{0,2}})})
T = multiTutte F
C = sub(T,apply(gens ring T, i -> i=>-1))

F = higgsLift(uniformMatroid(3,5))
T = multiTutte F
C = sub(T,apply(gens ring T, i -> i=>-1))

F = higgsLift(uniformMatroid(1,1)++uniformMatroid(2,3))
T = multiTutte F
C = sub(T,apply(gens ring T, i -> i=>-1))

F = higgsLift(matroid graph{{a,b},{b,c},{c,d},{d,a},{a,c}})
T = multiTutte F
C = sub(T,apply(gens ring T, i -> i=>-1))


isQuot = (m,n) -> isSubset(flats m, flats n);

M = randomSpMatroid(4,4)
L = Mmaxmin M
L/bases
isQuot(last L, first L)

select(100, i -> time (L := Mmaxmin randomSpMatroid(5,5); not isQuot(last L, first L)))


E = set toList(0..2)
BL = drop(subsets signedSubsets(3,3),1)/(b -> spMatroid(E,b))
#BL
time ML = select(BL, m -> isWellDefined m);

E = set toList(0..3)
BL = drop(subsets signedSubsets(4,4),1)/(b -> spMatroid(E,b));
#BL
time ML = select(BL, m -> isWellDefined m);

select(100, i -> if isWellDefined BL_(14*i+30000) then (L:= Mmaxmin BL_(14*i+30000); not isQuot(last L, first L)) else false)



M = uniformSpMatroid(3,3)
A = HRtest M
sub(A,{(ring A)_0=>1})
eigenvalues oo

M = spMatroid(matrix{{1,0,1,1},{0,1,1,0}})
bases M
isWellDefined M
A = HRtest M
sub(A,{(ring A)_0 =>1/2})
eigenvalues oo

M = uniformSpMatroid(3,3)
bases M
A = HRtest M
sub(A,{(ring A)_0 =>1})
eigenvalues oo

S = QQ[q, MonomialOrder=>GLex, Inverses=>true];
A = matrix{{0,q,0,q^(-1),1},{q,0,q^2,0,1},{0,q^2,0,q,1},{q^(-1),0,q,0,1},{1,1,1,1,2}}
determinant A
factor oo

signedSet({0,1},{1,2})
X = signedSet({0,1},{2,3})
star X
support X

signedSubsets 3
signedSubsets(3,3)

A = matrix{{1,0,1,0},{0,1,0,1}}
M = spMatroid A
B = bases M
B_0 + B_1
B_0 * B_1


M = uniformSpMatroid(2,2)
A = basisIndicatorMatrix M
L = apply(4, i -> A_i)
L/matrix/signedSet
P = polytope M
vertices P
faces(1,P)
rank(M, signedSet({1},{}))


M = uniformSpMatroid(3,6)
Q = polytope M
time L = faces(dim Q - 1,Q);
V = vertices Q;
E = apply(L/first, l -> V_(first l) - V_(last l));

M = randomSpMatroid(2,3);
gInvar M, gInvar(M, Greedy=>true)



---------------------< A bank of small symplectic matroids >----------------------
restart
load "CoxeterMatroids.m2"

E = set{0,1,2}
time L = signedSubsets(3,3);
time L2 = drop(subsets(L),1);
apply(20, i -> time isWellDefined spMatroid(E,first random L2))

apply(10, i -> time randomSpMatroid(4,4))


-------------------< The proposed G-invariant of symplectic matroid >---------------------
restart
load "CoxeterMatroids.m2"


--matroid subdivion of the uniform symplectic matroid of rank 2 on 2 elements
M = uniformSpMatroid(2,2)
E = M.groundSet; B = bases M
M1 = spMatroid(E, B_{0,1,3})
M2 = spMatroid(E, B_{0,2,3})
M12 = spMatroid(E, B_{0,3}) 
H = {M, M1, M2, M12}/gInvar
H_0 + H_3 === H_1 + H_2

--matroid subdivion of the uniform symplectic matroid of rank 3 on 3 elements
--it is a 3-cube, so we picked a subdivision into simplices coming from regular tetrahedron inside
M = uniformSpMatroid(3,3)
E = M.groundSet; B = bases M
M1 = spMatroid(E, B_{0,1,2,4})
bases M1
M2 = spMatroid(E, B_{3,2,1,7})
bases M2
M3 = spMatroid(E, B_{6,7,4,2})
bases M3
M4 = spMatroid(E, B_{5,4,7,1})
bases M4
M5 = spMatroid(E, B_{7,1,2,4})
bases M5
M15 = spMatroid(E, elements(set bases M1 * set bases M5))
bases M15
M25 = spMatroid(E, elements(set bases M2 * set bases M5))
M35 = spMatroid(E, elements(set bases M3 * set bases M5))
M45 = spMatroid(E, elements(set bases M4 * set bases M5))
H = {M,M1,M2,M3,M4,M5,M15,M25,M35,M45}/gInvar
H_0 + H_6 + H_7 + H_8 + H_9 === H_1 + H_2 + H_3 + H_4 + H_5


-----------------------< G-invariant specializations >--------------------------
restart
load "CoxeterMatroids.m2"

M = randomSpMatroid(3,4)
bases M
gInvar M
F = higgsLift M


M = randomSpMatroid(4,4)
P = polytope M
vertices P
isNormal P

apply(10, i -> (M := randomSpMatroid(2,4); (M, isNormal polytope M)))

--------------------------< random stuff >-------------------------------
restart
load "CoxeterMatroids.m2"


M = flagMatroid({uniformMatroid(2,5), uniformMatroid(3,5)})
multiTutte M

needsPackage "NormalToricVarieties"

V = {{-1,-1,-1},{1,0,0},{0,1,0},{0,0,1}};
S = drop(drop(subsets(4),1),-1)
L = poset(S,isSubset)
M = (maximalChains L)/(c -> apply(c, i -> position(S, j -> j == i)))
R = S/(s -> sum(V_s))



X = normalToricVariety(R,M)
isSmooth X
A = nefGenerators X
fromWDivToCl X



M = randomSpMatroid(3,3)
bases M
ring M
ideal oo
primaryDecomposition oo


M = uniformMatroid(2,3);
P = multiTutte M
sub(P,apply(gens ring P, i -> i => -1))

N = higgsLift M
P = multiTutte N
sub(P,apply(gens ring P, i -> i => -1))

M = uniformMatroid(1,1) ++ uniformMatroid(1,2)
N = higgsLift M
P = multiTutte N
sub(P,apply(gens ring P, i -> i => -1))

N = flagMatroid({uniformMatroid(2,5), uniformMatroid(3,5)})
isWellDefined N
P = multiTutte N
sub(P,apply(gens ring P, i -> i => -1))
numerator oo
oo/((ring oo)_0-1)



F = flagMatroid(matrix{{1,1,1},{1,0,0}},{1,2})
LVTutte(first F.constituents, last F.constituents)

R = QQ[x,y]
T = x^2*y^2+x^2*y+x*y^2+x^2+x*y
sub(T,{x=>x+1,y=>y+1})



------- M-permutohedron experiments
load "CoxeterMatroids.m2"




------finding that G-invar may not span in anything not minuscule type
restart
load "CoxeterMatroids.m2"

F = ZZ/3
A = random(F^5,F^5)
M = flagMatroid(A,{1,2,3,4})
gInvar M
volume polytope M
bases M

time L = apply(100, i -> ( m := flagMatroid(random(F^5,F^5),{1,3}); (m, gInvar m) ));
T = tally L/last;
K = select(keys T, i -> T#i > 1);
time apply(K, i -> (select(L, j -> (last j) === i))/first/polytope/volume)
