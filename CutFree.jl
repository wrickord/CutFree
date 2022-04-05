import Pkg;

Pkg.add("Clp")
Pkg.add("Gurobi")
Pkg.add("JuMP")

Pkg.update("Clp")
Pkg.update("Gurobi")
Pkg.update("JuMP")

using JuMP, Gurobi

IUB_CODES = Dict(
    "A" => ["A"], 
    "C" => ["C"], 
    "G" => ["G"], 
    "T" => ["T"], 
    "R" => ["A", "G"], 
    "Y" => ["C", "T"], 
    "M" => ["A", "C"], 
    "K" => ["G", "T"], 
    "S" => ["C", "G"],
    "W" => ["A", "T"], 
    "H" => ["A", "C", "T"], 
    "B" => ["C", "G", "T"], 
    "V" => ["A", "C", "G"], 
    "D" => ["A", "G", "T"], 
    "N" => ["A", "C", "G", "T"]
)

"""
get_blocking_codes(code, names)

Description:
    Find IUB codes that do not contain a specified nucleotide code.
    If 'names' is true, a vector of the possible IUB codes is returned. 
    If 'names' is false, a vector of the possible nucleotide groups is returned.

Examples:
    julia> println(get_blocking_codes("A", true))
    Any["T", "C", "G", "Y", "B", "S", "K"]
    julia> print(get_blocking_codes("A", false))
    Any[["T"], ["C"], ["G"], ["C", "T"], ["C", "G", "T"], ["C", "G"], ["G", "T"]]
"""
function get_blocking_codes(code, names)
    blocking_codes = Vector()

    for (key, value) in IUB_CODES
        if filter(x -> x != code, value) == IUB_CODES[key] && !names
            push!(blocking_codes, value)
        end

        if filter(x -> x != code, value) == IUB_CODES[key] && names
            push!(blocking_codes, key)
        end
    end

    return blocking_codes
end

"""
str_to_vector(str)

Description:
Transform a string into a vector of single strings.

Example:
julia> println(str_to_vector("CGTACCCGG"))
["C", "G", "T", "A", "C", "C", "C", "G", "G"]
"""
function str_to_vector(str)
    return string.(collect(str))
end

"""
vector_to_str(vctr)

Description:
Transform a vector of strings into a single string.

Example:
julia> println(vector_to_str(["C", "G", "T", "A", "C", "C", "C", "G", "G"]))
CGTACCCGG
"""
function vector_to_str(vctr)
    return join(vctr)
end

"""
make_oligo_block(oligo)

Description:
Creates a text block showing which bases are allowed by the oligo. 

Returns four strings denoting the presence of each base.

Example:
julia> println(make_oligo_block("ACGC"))
["A---", "----", "-C-C", "--G-"]
"""
function make_oligo_block(oligo)
    oligo_bases = str_to_vector(oligo)
    strings = Dict("A" => "", "C" => "", "G" => "", "T" => "")

    i = 0
    for code in oligo_bases
        for base in IUB_CODES[code]
            for (key, value) in strings
                if (key == base)
                    strings[key] *= key
                end
            end
        end

        i += 1
        for (key, value) in strings
            if (sizeof(value) != i)
                strings[key] *= "-"
            end
        end
    end

    return collect(values(strings))
end

"""
print_oligo_block(oligo)

Description:
Prints the oligo and the output of make_oligo_block.

Example:
julia> print_oligo_block("NAGN")
NAGN
AA-A
T--T
C--C
G-GG
"""
function print_oligo_block(oligo)
    return print(oligo * "\n" * make_oligo_block(oligo)[1] * "\n" * make_oligo_block(oligo)[2]
    * "\n" * make_oligo_block(oligo)[3] * "\n" * make_oligo_block(oligo)[4])
end

"""
xiny(x,y)

Description:
Returns a boolean value determining if x is contained within a vector y

Examples:
julia> println(xiny("A", IUB_CODES["N"]))
true
julia> println(xiny("A", IUB_CODES["C"]))
false
"""
function xiny(x, y)
    return (x ∈ y)
end

"""
subcodes(code)

Description:
Return array of IUB codes that are subsets of provided IUB code.

Example:
julia> println(subcodes("H"))
["A", "W", "T", "C", "Y", "M", "H"]
"""
function subcodes(code)
    codes = ""

    for (key, value) in IUB_CODES
        if (issubset(value, IUB_CODES[code]))
            codes *= key
        end
    end
    
    return str_to_vector(codes)
end

"""
complement(oligo)

Description:
Return array the complements for any given IUB codes.

Example:
julia> println(complement("HA"))
["D", "T"]
"""
function complement(oligo)
    codes = ""
    oligo_bases = str_to_vector(oligo)
    base_comps = Dict(
        "A" => "T", 
        "C" => "G", 
        "G" => "C", 
        "T" => "A",
    )

    for code in oligo_bases
        temp = ""
        for base in IUB_CODES[code]
            temp *= base_comps[base]
        end

        for (key, value) in IUB_CODES
            if (issetequal(str_to_vector(temp), value))
                codes *= key
            end
        end
    end

    return str_to_vector(codes)
end

"""
reverse_complement(oligo)

Description:
Return the reverse of an array of the complements for any given IUB codes.

Example:
julia> println(reverse_complement("HDN"))
["N", "H", "D"]
"""
function reverse_complement(oligo)
    return reverse(complement(oligo))
end

"""
expand_asymmetric(oligo)

Description:
Return both strands if reverse complement is different than oligo.

Example:
julia> println(expand_asymmetric("HDN"))
["N", "H", "D"]
"""
function expand_asymetric(oligo)
    if (isequal(str_to_vector(oligo), reverse_complement(oligo)))
        return [oligo]
    else
        return [oligo, vector_to_str(reverse_complement(oligo))]
    end
end

"""
degeneracy(oligo)

Description:
Return the number of sequences in a degenerate oligo.

Example:
julia> println(degeneracy("NNN"))
64
"""
function degeneracy(oligo)
    value = 1
    oligo_bases = str_to_vector(oligo)

    for code in oligo_bases
        value *= length(IUB_CODES[code])
    end

    return value
end

"""
cutfree()

Description:
Return a site-free random of maximum diversity (of nucleotides). This single degenerate sequence can be used to synthesize a pool of oligos with mixed bases.

Example:
julia>

"""
function cutfree(;
            starting_oligo = "NNNNNNNNNNNNNNNNNNNN",
            restriction_sites = ["GGTCTC", "GGCCGG"],
            min_blocks = 1,
            increase_diversity = true,
            )

    m = length(starting_oligo)
    k = []
    for i in restriction_sites
        push!(k, length(i))
    end

    starting_oligo = str_to_vector(starting_oligo)
    number_binary = 15 * m
    allowed_codes = []

    for i in 1:m
        subs = subcodes(starting_oligo[i])
        blocking = get_blocking_codes(starting_oligo[i], false)
        for b in subs
            
        end
        
    end

    # Find the maximum acheivable diversity and optimal objective value

    # Run CutFree with a constraint addeed requiring the solution to be equal to the optimal objective value

    # Change objective function to the sum of all variables multiplied by random coefficients (randomize)
    
    #=  model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "TimeLimit", maxtime)
        set_optimizer_attribute(model, "Presolve", 0)
    =#

    return nothing
end

cutfree()
