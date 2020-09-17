-*-----------------------------------------------------------------------

Tests for GKMManifolds.m2

-----------------------------------------------------------------------*-

--------------------------------
-- kTutte test
--------------------------------
TEST ///
ML = (random drop(drop(allMatroids 4,1),-1))_{0,1,2}
TML = apply(ML, m -> {tuttePolynomial m, kTutte flagMatroid({m})}) -- 4 seconds
assert all(TML, l -> (map(ring first l, ring last l, gens ring first l))(last l) == first l)
N = flagMatroid(matrix{{1,1,1},{1,0,0}},{1,2}) 
f = kTutte N --see Example 8.24 of [CDMS18]
x = (ring f)_0
y = (ring f)_1
assert (f == x^2*y^2 + x^2*y + x*y^2 + x^2 + x*y)
A = matrix{{1,1,1,1},{0,1,2,3}}
FMr = flagMatroid(A,{1,2})
FM = flagMatroid {uniformMatroid(1,4),uniformMatroid(2,4)}
assert (FMr === FM)
g = kTutte FM
assert (16 == sub(g, {(ring g)_0 => 1, (ring g)_1 => 1}))
///

--------------------------------
-- Type A tests (ordinary flag variety)
--------------------------------
TEST ///
X = tGeneralizedFlagVariety("A",3,{1,2})
A = matrix{{1,1,1,1},{0,1,2,3}}
U = flagMatroid {uniformMatroid(1,4),uniformMatroid(2,4)}
C = tKClass(X,U)
D = tOrbitClosure(X,A)
O1 = ampleTKClass X
assert (C === D)
assert (isWellDefined C)
assert (isWellDefined (C * O1))
Eu1 = tChi O1
Eu2 = tChi (C * O1)
assert (set exponents Eu1 === set exponents Eu2)
R = X.charRing
assert (sum(latticePts U, i -> R_i) == Eu2)
Y = tGeneralizedFlagVariety("A",3,{2},R)
assert (set Y.points === set apply(subsets(4,2), i -> {set i}))
Z = tGeneralizedFlagVariety("A",3,{1},R)
f = tFlagMap(X,Y)
g = tFlagMap(X,Z)
assert (O1 === (pullback f)(ampleTKClass Y) * (pullback g)(ampleTKClass Z))
dualTautSub = tKClass(Y,apply(Y.points, p -> sum(elements first p, i -> R_i)))
assert ((pushforward f)((pullback g)(ampleTKClass Z)) === dualTautSub)
assert (tChi ampleTKClass X == tChi (pushforward f)(ampleTKClass X))
///

----------------------------------------------------------------
-- Type C tests: Lagrangian Grassmannian LGr(2,4), and complete flag SpFl(1,2;4)
----------------------------------------------------------------
TEST ///
R = makeCharRing 2
X = tGeneralizedFlagVariety("C",2,{2},R)
assert (set X.points === set {{set{0,1}},{set{"0*",1}},{set{0,"1*"}},{set{"0*","1*"}}})
assert (6 == #(momentGraph X).edges)
O1 = ampleTKClass X
assert (tChi O1 == sum({{1,1},{-1,1},{1,-1},{-1,-1},{0,0}}, i -> R_i))
A1 = matrix{{1,0,1,0},{0,1,0,1}}
A2 = matrix{{1,0,1,2},{0,1,2,1}}
C1 = tOrbitClosure(X,A1)
C2 = tOrbitClosure(X,A2)
assert (tChi (C1*O1) == sum({{1,1},{-1,1},{1,-1},{-1,-1}}, i -> R_i))
assert (tChi (C2*O1) == sum({{1,1},{-1,1},{1,-1},{-1,-1},{0,0}}, i -> R_i))
Y = tGeneralizedFlagVariety("C",2,{1},R)
Z = tGeneralizedFlagVariety("C",2,{1,2},R)
f = tFlagMap(Z,Y)
g = tFlagMap(Z,X)
assert (tChi ampleTKClass Y == R_0 + R_1 + R_0^(-1) + R_1^(-1))
dualTautSub = tKClass(X,apply(X.points, p -> (
	    l := setIndicator(first p,2);
	    R_0^(l_0) + R_1^(l_1)
	    )
	))
assert (dualTautSub === (pushforward g)(pullback f)(ampleTKClass Y))
///


----------------------------------------------------------------
-- Type B tests: OG(3,7) and OGFl(2,3;7) and OGFl(1,2,3;7)
----------------------------------------------------------------
TEST ///
-- Checking points
R = makeCharRing 3
X = tGeneralizedFlagVariety("B",3,{3,3},R)
Y = tGeneralizedFlagVariety("B",3,{2,3,3},R)
Z = tGeneralizedFlagVariety("B",3,{1,2,3,3},R)
assert(set X.points === set {{set {0, 1, 2}}, {set {0, 1, "2*"}}, {set {0, "1*", 2}}, {set {0, "1*", "2*"}}, 
	{set {"0*", 1, 2}}, {set {"0*", 1, "2*"}}, {set {"0*", "1*", 2}}, {set {"0*", "1*", "2*"}}})
X' = tGeneralizedFlagVariety("C",3,{3,3},R)
Y' = tGeneralizedFlagVariety("C",3,{2,3,3},R)
Z' = tGeneralizedFlagVariety("C",3,{1,2,3,3},R)
assert(set X'.points === set X.points and set Y'.points === set Y.points and set Z'.points === set Z.points)

-- Checking charts
peek X;
HX = X.charts;
HX' = X'.charts;
missingCharX = apply(X.points, v -> toList(set HX'#v - set HX#v));
assert all(missingCharX, v -> #v == 3 and all(v, w -> max w ==2 or min w == -2))

HY = Y.charts;
HY' = Y'.charts;
missingCharY = apply(Y.points, v -> toList(set HY'#v - set HY#v));
missingCharY#0
assert all(missingCharY, v -> #v == 3 and all(v, w -> max w ==2 or min w == -2))

HZ = Z.charts;
HZ' = Z'.charts;
missingCharZ = apply(Z.points, v -> toList(set HZ'#v - set HZ#v));
assert all(missingCharZ, v -> #v == 3 and all(v, w -> max w ==2 or min w == -2))

-- Checking maps
f = tFlagMap(Z,Y)
g = tFlagMap(Y,X)
h = tFlagMap(Z,X)
assert(h === compose(g,f))

-- Checking TOrbClosure
M = matrix{{-2,0,0,1,0,0,2},{0,-2,0,2,1,0,2},{0,0,-2,2,2,1,2}}
-- Verifying M is isotropic
A = matrix{{-2,0,0},{0,-2,0},{0,0,-2}};
B = matrix{{1,0,0},{2,1,0},{2,2,1}};
D = matrix{{2},{2},{2}};
assert(A* transpose(B)  + B *transpose(A) + D*transpose(D) == 0)


time C1 = tOrbitClosure(X,M);
time C1'= tOrbitClosure(X,M, RREFMethod => true);
assert isWellDefined C1
assert(C1 === C1')

time C2 = tOrbitClosure(Y,M)
time C2'= tOrbitClosure(Y,M, RREFMethod => true)
assert isWellDefined C2
assert(C2 === C2')

time C3 = tOrbitClosure(Z,M)
time C3'= tOrbitClosure(Z,M, RREFMethod => true)
assert isWellDefined C3
assert(C3 === C3')
///



----------------------------------------
-- Type D tests: OG(4,8) and OGFl(1,4;8)
----------------------------------------
TEST ///
-- Checking points
R = makeCharRing 4
X1 = tGeneralizedFlagVariety("D",4,{4,4},R)
X2 = tGeneralizedFlagVariety("D",4,{3,3},R)
Z = tGeneralizedFlagVariety("C",4,{4},R)

PointSet = set{
    {set{0,1,2,3}}, {set{0,1,2,"3*"}}, {set{0,1,"2*",3}}, {set{0,1,"2*","3*"}},
    {set{"0*",1,2,3}}, {set{"0*",1,2,"3*"}}, {set{"0*",1,"2*",3}}, {set{"0*",1,"2*","3*"}},
    {set{0,"1*",2,3}}, {set{0,"1*",2,"3*"}}, {set{0,"1*","2*",3}}, {set{0,"1*","2*","3*"}},
    {set{"0*","1*",2,3}}, {set{"0*","1*",2,"3*"}}, {set{"0*","1*","2*",3}}, {set{"0*","1*","2*","3*"}}}
assert(set X1.points + set X2.points === PointSet)

assert(set X1.points + set X2.points === set Z.points)


-- Checking charts
HX1 = X1.charts;
HX2 = X2.charts;
HZ = Z.charts;
missingCharX = (apply(X1.points, v -> toList(set HZ#v - set HX1#v)) |apply(X2.points, v -> toList(set HZ#v - set HX2#v)))
assert all(missingCharX, v -> #v == 4 and all(v, w -> max w ==2 or min w == -2))


-- Checking TOrbClosure
-- Example that is 0 on X2
A = matrix{{1,3,-2,-1/4},{-1,-1,19,-61/4},{0,1,19,-73/4},{2,0,22,-89/4}}
B = matrix(QQ,{{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}})
M = A | B
-- Verifying M is isotropic
assert(A* transpose(B)  + B *transpose(A) == 0)

time C1 = tOrbitClosure(X1,M);
time C2 = tOrbitClosure(X2,M);
assert(isWellDefined C1 and isWellDefined C2)
assert all(X2.points, v -> C2.hilb#v == 0)

-- Example that is 0 on X1
A = matrix{{1,0,0,1},{0,1,0,1},{0,0,1,1},{0,0,0,0}}
B = matrix{{0,-1/9,-6/5,0},{1/9,0,-6/5,0},{6/5,6/5,0,0},{-1,-1,-1,1}}
N = A|B
-- Verifying M is isotropic
assert(A* transpose(B)  + B *transpose(A) == 0)

time D1 = tOrbitClosure(X1,N);
time D2 = tOrbitClosure(X2,N);
assert(isWellDefined D1 and isWellDefined D2)
assert all(X1.points, v -> D1.hilb#v == 0)



--- Same thing happens with the flag variety
R = makeCharRing 4
Y1 = tGeneralizedFlagVariety("D",4,{3,4,3,3},R)
Y2 = tGeneralizedFlagVariety("D",4,{3,4,4,4},R)
Z = tGeneralizedFlagVariety("D",4,{3,4},R)
assert(set ((X.points)/first) === set ((Z.points)/first))
assert(set ((Y.points)/first) === set ((Z.points)/first))

time C1 = tOrbitClosure(Y1,M);
time C2 = tOrbitClosure(Y2,M);
assert all(Y1.points, v -> C1.hilb#v == 0)
assert isWellDefined C2

f = tFlagMap(X,Z)
assert (trivialTKClass Z === (pushforward f)(trivialTKClass X) )
///



-------------------------
-- Low degree isogenies
-------------------------
TEST ///
-- OG(1,3), SGr(1,2) and P1^* = Gr(1,2)
X1 = tGeneralizedFlagVariety("B",1,{1,1})
X2 = tGeneralizedFlagVariety("C",1,{1})
X3 = tGeneralizedFlagVariety("A",1,{1})
assert(#X1.points === #X2.points)
assert(apply(X1.points, v-> #(X1.charts)#v) === apply(X2.points, v-> #(X2.charts)#v))
assert(apply(X1.points, v-> #(X1.charts)#v) === apply(X3.points, v-> #(X3.charts)#v))


-- SOGr(3,6), SOGr(2,6), and P3^* = Gr(3,4)
X1 = tGeneralizedFlagVariety("D",3,{3,3})
X2 = tGeneralizedFlagVariety("D",3,{2,2})
X3 = tGeneralizedFlagVariety("A",3,{3})
assert(#X1.points === #X2.points and #X2.points === #X3.points)
assert(apply(X1.points, v-> #(X1.charts)#v) == apply(X2.points, v-> #(X2.charts)#v) and
    apply(X1.points, v-> #(X1.charts)#v) === apply(X3.points, v-> #(X3.charts)#v))

--SOGr(1,6) and Gr(2,4)
X1 = tGeneralizedFlagVariety("D",3,{1})
X2 = tGeneralizedFlagVariety("A",3,{2})
assert(#X1.points === #X2.points)
assert(apply(X1.points, v-> #(X1.charts)#v) === apply(X2.points, v-> #(X2.charts)#v))

///

-------------------------------------
-- Boundary cases for TOrbitClosure
-------------------------------------
TEST ///
M = matrix(QQ,{{1,0,0,0},{0,1,0,0}})
R = makeCharRing 4
X = tGeneralizedFlagVariety("A",3,{2},R)
C = tOrbitClosure(X,M)
peek C
assert all(keys C.hilb, v -> (
	if not v === {set{0,1}} then C.hilb#v == 0 else 
	sub(C.hilb#v,R)  === product apply({0,1} ** {2,3}, v -> 1- R_(v_0)^(-1)*R_(v_1)) 
	)
    )
///

-----------------------
-- A toric variety test
-----------------------
TEST ///
X = kleinschmidt(2,{2})
Y = tVariety X
assert (normalToricVariety Y === X)
R = Y.charRing
assert (2 == numgens R)
G = momentGraph Y
assert all(keys G.edges, e -> 0 == (matrix{(rays X)_(first elements (set e_0 * set e_1))} * transpose matrix{G.edges#e}))
antiK = - toricDivisor X
TantiK = tKClass(Y,antiK)
ltsPts = {{1,1},{0,-1},{0,0},{0,1},{-1,-3},{-1,-2},{-1,-1},{-1,0},{-1,1}}
assert (sum(ltsPts, i -> R_i) == tChi TantiK)
///

-*---- old tests during developing

----------< tGeneralizedFlagVariety tests >-------------
restart
needsPackage "GKMManifolds"

TEST ///
--the flag variety of 1-dim in 2-dim in ambient 3-dim space, as proj var embedded in P7
--(here, dim means linear space dim, not projective dim)
X = tGeneralizedFlagVariety("A",2,{1,2})
peek X
H = X.charts
G = momentGraph X
peek G
G.edges
graph G
C = ampleTKClass X
peek C
isWellDefined C
lieType X
P = bruhatOrder X
displayPoset(P, SuppressLabels => false)

--if a character ring is pre-determined
R = makeCharRing 3
X = tGeneralizedFlagVariety("A",2,{1,2},R)
C = ampleTKClass X
ring first values C.hilb === R
peek C
isWellDefined C

--ordinary Grassmannian(2,4)
X = tGeneralizedFlagVariety("A",3,{2})
peek X
C = ampleTKClass X
tChi C
Y = tGeneralizedSchubertVariety(X,{set{0,1}}) --basically same as X
peek Y
Z = tGeneralizedSchubertVariety(X,{set{0,2}})
peek Z
W = tGeneralizedSchubertVariety(X,{set{1,2}})
peek W
displayPoset cellOrder momentGraph X
displayPoset cellOrder momentGraph Z


--same variety but embedded differently in P14
X = tGeneralizedFlagVariety("A",2,{1,1,2})
peek X
H = X.charts
G = momentGraph X
peek G
G.edges
C = ampleTKClass X
peek C

--Lagrangian Grassmannian(2,4)
X = tGeneralizedFlagVariety("C",2,{2})
peek X
H = X.charts
G = momentGraph X
peek G
G.edges
graph G
C = ampleTKClass X
peek C
isWellDefined C
lieType X
P = bruhatOrder X
displayPoset(P, SuppressLabels => false)

--Lagrangian Grassmannian(3,6)
X = tGeneralizedFlagVariety("C",3,{3})
peek X
H = X.charts
netList {(first keys H), H#(first keys H)}
G = momentGraph X
peek G
G.edges
C = ampleTKClass X
peek C

--isotropic Grassmannian(1,2,3;6) (i.e. SpGr(1,2,3;6))
X = tGeneralizedFlagVariety("C",3,{1,2,3})
peek X
H = X.charts
netList {(first keys H), H#(first keys H)}
G = momentGraph X
peek G
G.edges
C = ampleTKClass X
peek C
isWellDefined C
P = bruhatOrder X
displayPoset P

X = tGeneralizedFlagVariety("B",3,{2})
peek X
H = X.charts
G = momentGraph X
peek G
G.edges
C = ampleTKClass X
peek C


--spin groups are not implemented yet
X = tGeneralizedFlagVariety("B",3,{3})





--even orthogonal Grassmannians are very tricky!
--All fine as long as k < n-1, where we consider k-dim'l isotropic subspace of a 2n-dim'l space.
--When k = n-1, we need consider the sum of the two fundamental weights
--When k = n, it seems like we can one conn component per fundamental weight (? check)


--orthogonal Grassmannian(2,8)
X = tGeneralizedFlagVariety("D",4,{2})
peek X
H = X.charts
netList {(first keys H), H#(first keys H)}
G = momentGraph X
graph G
peek G
G.edges
C = ampleTKClass X
peek C


--orthogonal Grassmannian(3,8) (i.e. SOGr(3,8))
--I think? This needs to be checked
X = tGeneralizedFlagVariety("D",4,{3,4})
X.points
peek X
H = X.charts
H#{set{0,1,2}}
G = momentGraph X
graph G
peek G
G.edges
C = ampleTKClass X
peek C

--again, spin groups not implemented
X = tGeneralizedFlagVariety("D",4,{4})
X = tGeneralizedFlagVariety("D",4,{3})


--one of two connected components (I think?) 
--of orthogonal Grassmannian(4,8) (i.e. SOGr(4,8))
X = tGeneralizedFlagVariety("D",4,{4,4})
peek X
H = X.charts
netList {(first keys H), H#(first keys H)}
G = momentGraph X
graph G
peek G
G.edges
C = ampleTKClass X
peek C

--the other half of orthogonal Grassmannian(4,8) (i.e. SOGr(4,8))
X = tGeneralizedFlagVariety("D",3,{3,3})
peek X
H = X.charts
netList {(first keys H), H#(first keys H)}
G = momentGraph X
graph G
peek G
G.edges
C = ampleTKClass X
peek C


-------------------------------------------< speed tests >-------------------------------------
--Seems like for Grassmannians the functions "tGeneralizedFlagVariety" and "tFlagVariety"
--are about the same speed because the limitation is in the function "permutations"
--However, for partial flag varities, the "tGeneralizedFlagVariety" is considerably slower
--Similarly in other Lie types than "A", Grassmannians are pretty fast, limited by the
--function "signedPermutations", but much slower for partial flag varieties.


restart
uninstallPackage"GKMManifolds"
installPackage "GKMManifolds"


--Grassmannians
time X = tGeneralizedFlagVariety("A", 8, {3}) -- 3 seconds
time X = tGeneralizedFlagVariety("A", 8, {4}) -- 3 seconds
time X = tGeneralizedFlagVariety("A", 9, {4}) -- 29 seconds
time permutations {1,2,3,4,5,6,7,8,9,10}; -- 27 seconds
--so Gr(4,10) above is slow due to permutations being slow
time X = tFlagVariety({4},9) -- 4 seconds
time X = tFlagVariety({3},9) -- 4 seconds

--full flag varieties
time X = tGeneralizedFlagVariety("A", 4, {1,2,3,4}) -- 0.3 seconds
time X = tGeneralizedFlagVariety("A", 5, {1,2,3,4,5}) -- 6 seconds
--time X = tGeneralizedFlagVariety("A", 6, {1,2,3,4,5,6}) -- don't try on laptop
time X = tFlagVariety({1,2,3,4,5,6},7) -- 3 seconds!!

--isotropic Grassmannians
time X = tGeneralizedFlagVariety("B", 7, {4}) -- 3 seconds
time X = tGeneralizedFlagVariety("B", 7, {6}) -- 3 seconds
time X = tGeneralizedFlagVariety("B", 7, {7,7}) -- 2 seconds
time X = tGeneralizedFlagVariety("B", 8, {8,8}) -- 34 seconds
time X = tGeneralizedFlagVariety("B", 8, {4}) -- 34 seconds
time signedPermutations {1,2,3,4,5,6,7,8}; -- 28 seconds
--again, the slowness for isotropic Grassmannians is largely due to signPermutations

--flag isotropic
time X = tGeneralizedFlagVariety("B", 4, {1,2,3,4,4}) -- 3 seconds
time X = tGeneralizedFlagVariety("B", 5, {1,2,3,4}) -- 40 seconds

--few other types
time X = tGeneralizedFlagVariety("C", 7, {4}) -- 3 seconds
time X = tGeneralizedFlagVariety("D", 7, {6,7}) -- 3 seconds
time X = tGeneralizedFlagVariety("C", 5, {1,2,4}) -- 7 seconds
time X = tGeneralizedFlagVariety("D", 5, {2,3,4,5}) -- 5 seconds

--a test for normal toric varieties: Hirezebruch surface H_2
restart
needsPackage "GKMManifolds"
X = kleinschmidt(2,{2})
peek X
Y = tVariety X
peek Y
peek momentGraph Y
normalToricVariety Y
antiK = - toricDivisor X
TantiK = tKClass(Y,antiK)
peek TantiK
tChi TantiK

--a sanity check with projective space P^3
P3 = tProjectiveSpace 3
R = P3.charRing;
describe R
charts P3
G = momentGraph P3
G.edges
peek G
graph G


assert isWellDefined C
badTKC = tKClass(P3,{R_0,R_0,R_0,R_1})
peek badTKC
assert not isWellDefined badTKC


--product of T-varieties check with P^1
P1 = tProjectiveSpace 1
P1xP1 = P1 ** P1
peek P1xP1
G = P1.momentGraph
GxG = G ** G
peek GxG
graph GxG

--test of affineToricRing
B = {{1,0},{0,1},{1,1},{1,-1}}
affineToricRing B
hilbertSeries(oo, Reduce => true)


----------< TOrbClosure tests >-------------

restart
needsPackage "GKMManifolds"

-- Easy test (can compre with Speyer/Fink pg 17)
M = matrix(QQ,{{1,1,0,1},{0,0,1,1}})
N = matrix(QQ,{{1,1,1,3},{1,1,2,1}})
X = tGeneralizedFlagVariety("A",3,{2})
time C1 = TOrbClosure(X,{M})
time C2 = TOrbClosure(X,{M})
peek C1
peek C2
C = tKClass(X,flagMatroid(M,{2})); peek C


--Lagrangian Grassmannian test
M = matrix(QQ,{{1,0,1,2},{0,1,2,1}})
X = tGeneralizedFlagVariety("C",2,{2})
C = TOrbClosure(X,{M})
peek C
isWellDefined C


-- Type "A"
M = matrix(QQ,{{1,1,0,1},{0,0,1,1}})
N = matrix(QQ,{{1,0,0,0},{0,1,0,0}})
X = tGeneralizedFlagVariety("A",3,{2})
time C = TOrbClosure(X,{M});
peek C
isWellDefined C
time D = TOrbClosure(X,{N}); --example where the matroid has only one basis
peek D
D.hilb#{set{0,1}}
factor oo

Y = tGeneralizedFlagVariety("A",2,{2,1})
M = matrix(QQ,{{1,0,0},{0,1,1}})
N = matrix(QQ,{{1,1,1}})
time C = TOrbClosure(Y,{N,M}); peek C



-- Type "C"
M = matrix(QQ,{{1,0,1,3},{0,1,3,1}})
N = matrix(QQ,{{1,0,1,0},{0,1,0,1}})
X = tGeneralizedFlagVariety("C",2,{2})
time C = TOrbClosure(X,{M}); D = TOrbClosure(X,{N});
peek C
peek D
{C,D}/isWellDefined



M = matrix{{1,1,1,0,0,0},{0,0,0,1,1,-2}}
X = tGeneralizedFlagVariety("C",3,{2})
time C = TOrbClosure(X,{M})

N = matrix{{1,1,1,1,1,-2}}
Y = tGeneralizedFlagVariety("C",3,{2,1})
time D = TOrbClosure(Y,{N,M});
peek D






-- Error testing
M = matrix(QQ,{{1,0,1,3},{0,1,3,1}})
X = tProjectiveSpace 3
C = tOrbitClosure(X,M)

X = tGeneralizedFlagVariety("A",3,{3})
C = TOrbClosure(X,{M})




-- Sanity check
-- The closure of the following is just a point




A = matrix(QQ,{{1,0,2,1},{0,1,1,2}})
X = tGeneralizedFlagVariety("C",2,{2})
C = tOrbitClosure(X,A); peek C



M = matrix(QQ,{{1,2,0,0},{1,2,0,0}})
X = tGeneralizedFlagVariety("A",3,{1})
C = tOrbitClosure(X,M); peek C
M = matrix(QQ,{{1,2,0,0}})



----------------------------------------------------------*-
