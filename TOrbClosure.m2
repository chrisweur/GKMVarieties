restart
uninstallPackage "GKMManifolds"
installPackage "GKMManifolds"
needsPackage "GKMManifolds"

--------------------------------------------------------------------------------------------------
--------------------------------< T-orbit closures of points >------------------------------------

---------------------------< auxiliary functions for TOrbClosure >--------------------------------
convertToNum = (n,L) -> (
    return apply(toList L, v -> if v === unastrsk(v) then v else n + unastrsk v)
    )

revMat = M -> return matrix apply( reverse entries M, v-> reverse v)

-- Takes in a mutable matrix and outputs the RREF with the identitiy block in the beginning
rowRed = M -> (
    Mat := matrix revMat(M) ;
    return mutableMatrix revMat transpose gens gb image transpose Mat;
    )



--the T-equivariant K-class of a torus orbit closure of a point in a generalized flag variety
--input:  X a tGeneralizedFlagVariety and MatLst a list of matrices, {M1,...,Mn}, that defines a point in X.
--    	  For convenience we assume that the ranks of the M_i are distinct.    	
--output: The tKClass of the closure of the orbit of the point corresponding to {M1,...,Mn}
TOrbClosure = method()
TOrbClosure(TVariety,List) := TKClass => (X,MatLst) -> ( 
    if X.cache.?lieType then Typ := X.cache.lieType else (
	 << "T-orbit closures are only implemented for Lie types" << 
	 return error
	 );
    R := X.charRing;
    m := numgens X.charRing;
    x := symbol x;
    S := QQ[x_0..x_(m-1)];
    QS := frac(QQ[x_0..x_(m-1), Degrees => apply(m, i -> setIndicator(set {i},m))]);
    rks :=  apply(first X.points, v -> #(elements v));
    if not all(MatLst, v -> numcols v == (
	    if Typ === "A" then m
	    else if Typ === "B" then 2*m+1
	    else 2*m )
	) then << "the column size of the matrices are incorrect" << return error;
    if not rks === apply(MatLst, v -> rank v) then (
	<< " the rank of the matrices are incorrect" << 
	return error
	);
    (if Typ === "A" then Gens := apply(gens QS, v -> v^(-1))
	else if Typ  === "C" or Typ === "D" then Gens = apply(gens QS, v -> v^(-1)) | gens QS
	else if Typ === "B" then Gens = apply(gens QS, v -> v^(-1)) | gens QS |{1}
	);
    Lst := apply(X.points, pt -> (
	    degList := {};
	    convertPt := apply(pt, v -> convertToNum(m,v));
	    if any(#rks, v -> determinant (rowRed MatLst_v)_(convertPt_v) == 0) then 0_R else (
		for j in toList(0..#rks-1) do (
		    L := (reverse convertPt)#j;
		    M := mutableMatrix sub((reverse MatLst)#j,QS);
		    for i in toList(0..#L-1) do M = columnSwap(M,i,L#i);
		    M = rowRed M;
		    for i in reverse toList(0..#L-1) do M = columnSwap(M,i,L#i);
		    for i in L do M = columnMult(M,i,0);
		    for v in toList(0..#L-1) do M = rowMult(M,v,sub(Gens_(L_v)^(-1),ring M));
		    for v in toList(0..#Gens-1) do M = columnMult(M,v,sub(Gens_v, ring M));
		    degList = degList | apply(flatten entries M, v-> degree v);
		    );
		degs := - unique delete(-infinity, degList);
		tHilbSer := hilbertSeries affineToricRing degs;
		numer := toCharRing(X,value numerator tHilbSer);
		denom := toCharRing(X,value denominator tHilbSer);
		fracVal := toFraction(numer * product(X.charts#pt, l -> (1-R_l)), denom, S);
		(last fracVal)(first fracVal)
		-* degs *-
	    	)
	    )
    	);
    tKClass(X,Lst)
    -* hashTable apply(#X.points, i -> (X.points_i, Lst_i)) --and comment this *-
    )





------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------

restart
needsPackage "GKMManifolds"


-- Type "A"
M = matrix(QQ,{{1,1,0,1},{0,0,1,1}})
X = tGeneralizedFlagVariety("A",3,{2})
C = tOrbitClosure(X,M)
peek C
D = tOrbitClosure(X,M,RREFMethod => true)
peek D
C === D

M = random(QQ^2,QQ^4)
X = tGeneralizedFlagVariety("A",3,{1,2})
time C = tOrbitClosure(X,M); peek C
time D = tOrbitClosure(X,M,RREFMethod => true); peek D
C === D


Y = tGeneralizedFlagVariety("A",2,{2,1})
M = matrix(QQ,{{1,0,0},{0,1,1}})
N = matrix(QQ,{{1,1,1}})
time C = TOrbClosure(Y,{N,M}); peek C



-- Type "C"
M = matrix(QQ,{{1,0,1,3},{0,1,3,1}})
N = matrix(QQ,{{1,0,1,0},{0,1,0,1}})
X = tGeneralizedFlagVariety("C",2,{2})
time C = tOrbitClosure(X,M); D = tOrbitClosure(X,N);
peek C
peek D
time C1 = tOrbitClosure(X,M,RREFMethod => true); D1 = tOrbitClosure(X,N,RREFMethod => true);
C1 === C, D1 === D

M = matrix{{1,1,1,0,0,0},{0,0,0,1,1,-2}}
X = tGeneralizedFlagVariety("C",3,{2})
time C = TOrbClosure(X,{M})

N = matrix{{1,1,1,1,1,-2}}
Y = tGeneralizedFlagVariety("C",3,{2,1})
time D = TOrbClosure(Y,{N,M});
peek D



-- Type "B"
X = tGeneralizedFlagVariety("B",3,{1})
peek X
M = matrix(QQ,{{-2,0,0,1,0,0,2}})
--matrix{{-2,0,0}} * transpose matrix{{1,0,0}}  + ( matrix{{1,0,0}} * transpose matrix{{-2,0,0}} ) + 4
C = TOrbClosure(X,{M})
peek C



-- Error testing
M = matrix(QQ,{{1,0,1,3},{0,1,3,1}})
X = tProjectiveSpace 3
C = TOrbClosure(X,{M})

X = tGeneralizedFlagVariety("A",3,{3})
C = TOrbClosure(X,{M})




-- Sanity check
-- The closure of the following is just a point
M = matrix(QQ,{{1,0,0,0},{0,1,0,0}})
X = tGeneralizedFlagVariety("A",3,{2})
C = TOrbClosure(X,{M}); peek C

Y = tGeneralizedFlagVariety("A",3,{2,2})
C = TOrbClosure(Y,{M,M}); peek C

Z = tGeneralizedFlagVariety("C",2,{2})
C = TOrbClosure(Z,{M}); peek C



M = matrix(QQ,{{1,2,0,0},{1,2,0,0}})
X = tGeneralizedFlagVariety("A",3,{1})
C = TOrbClosure(X,{M}); peek C









--- 
-- Needs a better name? Projection maps?
tGeneralizedMap = method()
tGeneralizedMap(TVariety,TVariety) := TMap => (X,Y) -> (
    if not (X.cache.?lieType and Y.cache.?lieType) then << "projection maps are only implemented for Lie types" << return error;
    if not numgens X.charRing === numgens Y.charRing then << "character ring need be same" << return error;
    TypX := X.cache.lieType;
    TypY := Y.cache.lieType;
    if not TypX === TypY then << "the varieties are of different Lie types" << return error;
    LX := apply(first X.points, v -> #(elements v));
    LY := apply(first Y.points, v -> #(elements v));
    if not isSubset(LY,LX) then << "there is no projection map" << return error;
    Typ := X.cache.lieType;
    m := numgens X.charRing;
    rk := if Typ === "A" then m-1 else m;
    T := symbol T;
    R := makeCharRing(m);
    X' := tGeneralizedFlagVariety(TypX,rk,LX,R);
    Y' := tGeneralizedFlagVariety(TypY,rk,LY,R);
    ptPairs := apply(X'.points, p -> (p,first select(Y'.points, q -> isSubset(q,p))));
    tMap(X',Y',ptPairs)
)

-- Tests
X = tGeneralizedFlagVariety("A",3,{1,2})
Y = tGeneralizedFlagVariety("A",3,{2})
peek tGeneralizedMap(X,Y)
tGeneralizedMap(Y,X)


X = tGeneralizedFlagVariety("B",3,{2})
Y = tGeneralizedFlagVariety("A",2,{2})
tGeneralizedMap(X,Y)


X = tGeneralizedFlagVariety("B",3,{1,2})
Y = tGeneralizedFlagVariety("B",3,{1})
peek tGeneralizedMap(X,Y)

