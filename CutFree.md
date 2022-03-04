# IUB_Codes

## Description:
Dictionary consisting of IUB keys and associated nucleotide values.

## Examples:
```julia-repl
julia> println(IUB_CODES["R"][2])
G
```

# get_blocking_codes(code, names)

## Description:
Find IUB codes that do not contain a specified nucleotide code.

If 'names' is true, a vector of the possible IUB codes is returned. 

If 'names' is false, a vector of the possible nucleotide groups is returned.

## Examples:
```julia-repl
julia> print(get_blocking_codes("A", true))
["T", "C", "G", "Y", "B", "S", "K"]
julia> print(get_blocking_codes("A", false))
[["T"], ["C"], ["G"], ["C", "T"], ["C", "G", "T"], ["C", "G"], ["G", "T"]]
```

# str_to_vector(str)

## Description:
Transform a string into a vector of single characters.

## Examples:
```julia-repl
julia> println(str_to_vector("CGTACCCGG"))
SubString{String}["C", "G", "T", "A", "C", "C", "C", "G", "G"]
```

# make_oligo_block(oligo)

## Description:
Creates a text block showing which bases are allowed by the oligo. 

Returns four strings denoting the presence of each base.

## Examples:
```julia-repl
julia> print(make_oligo_block("ACGC"))
["A---", "----", "-C-C", "--G-"]
```

# print_oligo_block(oligo)

## Description:
Prints the oligo and the output of make_oligo_block.

## Examples:
```julia-repl
julia> print_oligo_block("NAGN")
NAGN
AA-A
T--T
C--C
G-GG
```

# xiny(x,y)

## Description:
Returns a boolean value determining if x is contained within a vector y

## Examples:
```julia-repl
julia> print(xiny("A", IUB_CODES["N"]))
true
julia> print(xiny("A", IUB_CODES["C"]))
false
```

# subcodes(code)

## Description:
Return array of IUB codes that are subsets of provided IUB code.

## Examples:
```julia-repl
julia> print(subcodes("H"))
SubString{String}["A", "W", "T", "C", "Y", "M", "H"]
```

# complement(oligo)

## Description:
Return array the complements for any given IUB codes.

## Examples:
```julia-repl
julia> print(complement("HA"))
SubString{String}["D", "T"]
```