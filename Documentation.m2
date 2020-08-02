newPackage("EquivariantLocalizations",
	Version => "0.1",
	Date => "August xx, 2020",
	Authors => {
	    {Name => "Chris Eur",
       	     Email => "ceur@math.berkeley.edu",
       	     HomePage => "https://math.berkeley.edu/~ceur"},
	    {Name => "Ritvik Ramkumar",
	     Email => "ritvik@math.berkeley.edu",
	     HomePage => "https://math.berkeley.edu/~ritvik"}
	    },
	Headline => "a package for computations with flag matroids and equivariant localizations",
	HomePage => "https://github.com/chrisweur/EquivariantLocalizations",
	PackageImports => {"Matroids"},
	PackageExports => {"Matroids"},	
	DebuggingMode => true
)
export {
	"TVariety",
	"tVariety", -- Do I need to export product of T-varieties?
    	"TKClass",
	"tKClass", -- Again, do I need to export product/addition of K classes?
	"ampleTKClass",
	"TMap",
	"tMap",
	"pullback",
	"pushforward",-- Again, product of T-maps?	
	"diagonalTMap",
	--"compose",
	"FlagMatroid",
	"flagMatroid",
	--"isWellDefined",
	--"rank",
	"tFlagVariety",
	"tFlagMap",
	"fourierMukai",
	"kTutte",
	"kCharPol"
}



beginDocumentation()

-- Documentation --
-- <<docTemplate
doc ///
	Key
		EquivariantLocalizations
	Headline
		a package for computations with flag matroids and equivariant localization
	Description
		Text
			Define in words T-Varieties and equivariant localization......

			This Macaulay2 package is designed to .....

			References.....
			

///


doc ///
	Key
		TVariety
	Headline
		the class of all T-varieties
	Description
		Text
			To see how to specify a T-variety, see @TO tVariety@.

			Describe basic functionality of the package...
		
		Example
			Here....

///


doc ///
	Key
		tVariety
		(tVariety, List, HashTable, List, Ring)
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
			
		Example
			Here...

	Caveat
		This function does not check if X defines a T-variety - see 
		@TO2{(isWellDefined, TVariety), "isWellDefined"}@.
	
	SeeAlso
		(isWellDefined, TVariety)
		(symbol **, TVariety, TVariety)
		tFlagVariety
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
			
		Example
			Here...

	Caveat
		This function does not check if X defines a T-variety - see 
		@TO2{(isWellDefined, TVariety), "isWellDefined"}@.
	
	SeeAlso
		(isWellDefined, TVariety)
///


doc ///
	Key
		TKClass
	Headline
		a type of HashTable
	Description
		Text
			Describe it...
			
	SeeAlso
	    

///


doc ///
	Key
		tKClass
		(tKClass, TVariety, List)
	Headline
		constructs a equivariant K-class
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
			Here..
			
		Example
			Here...

	Caveat
		This function does not check if X defines a T-variety - see 
		@TO2{(isWellDefined, TKClass), "isWellDefined"}@.
	
	SeeAlso
		(isWellDefined, TKClass)
		(symbol *, TKClass, TKClass)
		(symbol +, TKClass, TKClass)
		assignAmpleTKClass
		pullback
		pushforward
///





end


