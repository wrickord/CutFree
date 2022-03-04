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
julia> println(get_blocking_codes("A", true))
Any["T", "C", "G", "Y", "B", "S", "K"]
julia> print(get_blocking_codes("A", false))
Any[["T"], ["C"], ["G"], ["C", "T"], ["C", "G", "T"], ["C", "G"], ["G", "T"]]
```

# vector_to_str(vctr)

## Description:
Transform a string into a vector of single strings.

## Examples:
```julia-repl
julia> println(str_to_vector("CGTACCCGG"))
["C", "G", "T", "A", "C", "C", "C", "G", "G"]
```

# str_to_vector(str)

## Description:
Transform a vector of strings into a single string.

## Examples:
```julia-repl
julia> println(vector_to_str(["C", "G", "T", "A", "C", "C", "C", "G", "G"]))
CGTACCCGG
```

# make_oligo_block(oligo)

## Description:
Creates a text block showing which bases are allowed by the oligo. 

Returns four strings denoting the presence of each base.

## Examples:
```julia-repl
julia> println(make_oligo_block("ACGC"))
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
julia> println(xiny("A", IUB_CODES["N"]))
true
julia> println(xiny("A", IUB_CODES["C"]))
false
```

# subcodes(code)

## Description:
Return array of IUB codes that are subsets of provided IUB code.

## Examples:
```julia-repl
julia> println(subcodes("H"))
["A", "W", "T", "C", "Y", "M", "H"]
```

# complement(oligo)

## Description:
Return array the complements for any given IUB codes.

## Examples:
```julia-repl
julia> println(complement("HA"))
["D", "T"]
```

# reverse_complement(oligo)

## Description:
Return the reverse of an array of the complements for any given IUB codes.

## Examples:
```julia-repl
julia> println(reverse_complement("HDN"))
["N", "H", "D"]
```

# degeneracy(oligo)

## Description:
Return the number of sequences in a degenerate oligo.

## Examples:
```julia-repl
julia> println(degeneracy("NNN"))
64
```