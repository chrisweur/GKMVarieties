
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
--input: X a tGeneralizedFlagVariety
--output:
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
		    for i in toList(0..#L-1) do M = columnSwap(M,i,L#i);
		    for i in L do M = columnMult(M,i,0);
		    for v in toList(0..#L-1) do M = rowMult(M,v,sub(Gens_(L_v)^(-1),ring M));
		    for v in toList(0..#Gens-1) do M = columnMult(M,v,sub(Gens_v, ring M));
		    degList = degList | apply(flatten entries M, v-> degree v);
		    );
		degs := - unique delete(-infinity, degList);
		-*-- uncomment the next five lines once bug for "degs" is fixed
		tHilbSer := hilbertSeries affineToricRing degs;
		numer := toCharRing(X,value numerator tHilbSer);
		denom := toCharRing(X,value denominator tHilbSer);
		fracVal := toFraction(numer * product(X.charts#pt, l -> (1-R_l)), denom, S);
		(last fracVal)(first fracVal)
		--*-
		degs --comment out this line then
	    	)
	    )
    	);
    -*-- uncomment this too then
    tKClass(X,Lst)
    --*-
    hashTable apply(#X.points, i -> (X.points_i, Lst_i)) --and comment this
    )


end

------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------


restart
load "TOrbClosure.m2"


-- Easy test (can compre with Speyer/Fink pg 17)
M = matrix(QQ,{{1,1,0,1},{0,0,1,1}})
X = tGeneralizedFlagVariety("A",3,{2})
time C = TOrbClosure(X,{M})
peek C --the degrees at {set{1,2}} is incorrect: {1,0,-1,0} must be replaced with {1,-1,0,0}


Creal = tKClass(X,flagMatroid(M,{2}))
peek Creal
isWellDefined Creal


--matrix with same underlying matroid as M above
N = matrix(QQ,{{1,1,1,3},{1,1,2,1}})
CN = TOrbClosure(X,{N})
C === CN





--A Lagrangian Grassmannian example
M = matrix(QQ,{{1,0,1,2},{0,1,2,1}})
N = matrix(QQ,{{1,0,1,0},{0,1,0,1}})
X = tGeneralizedFlagVariety("C",2,{2})
CM = TOrbClosure(X,{M})
CN = TOrbClosure(X,{N})




--- Basic functionality tests (not checked for math yet)
M = matrix(QQ,{{1,1,0,1},{0,0,1,1}})
N = matrix(QQ,{{1,1,0,1}})
X = tGeneralizedFlagVariety("A",3,{1,2})
peek X
time C = TOrbClosure(X,{N,M})
peek C


M = matrix(QQ,{{1,1,0,1},{1,1,0,1}})
M = matrix(QQ,{{1,1,0,1}})
X = tGeneralizedFlagVariety("A",3,{1})
peek X
time C = TOrbClosure(X,{M})
peek C


--
M = random(QQ^2,QQ^4)
X = tGeneralizedFlagVariety("A",3,{2})
TOrbClosure(X,{M})
peek C




-- 
X = tProjectiveSpace 3
C = TOrbClosure(X,{M})


--
X = tGeneralizedFlagVariety("B",3,{1})
peek X
M = matrix(QQ,{{-2,0,0,1,0,0,2}})
--matrix{{-2,0,0}} * transpose matrix{{1,0,0}}  + ( matrix{{1,0,0}} * transpose matrix{{-2,0,0}} ) + 4
C = TOrbClosure(X,{M})
peek C






-- Sanity check
M = matrix(QQ,{{1,0,0,0},{0,1,0,0}})
X = tGeneralizedFlagVariety("A",3,{2})
C = TOrbClosure(X,{M}); peek C
X = tGeneralizedFlagVariety("C",2,{2})
C = TOrbClosure(X,{M}); peek C


