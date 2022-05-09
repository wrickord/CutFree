using Clp, Gurobi, GLPK, JuMP, DataStructures, NamedArrays

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
    blocking_codes = []
    for (key, value) in IUB_CODES
        if findall(x -> x in value, IUB_CODES[code]) == Int64[] && !names
            push!(blocking_codes, value)
        end

        if findall(x -> x in value, IUB_CODES[code]) == Int64[] && names
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
    strings = Dict("A" => "", "C" => "", "G" => "", "T" => "")

    i = 0
    for code in str_to_vector(oligo)
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
    println(oligo * "\n" * make_oligo_block(oligo)[1] * "\n" * make_oligo_block(oligo)[2]
    * "\n" * make_oligo_block(oligo)[3] * "\n" * make_oligo_block(oligo)[4])
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
    Return array of the complements for any given IUB codes.

    Example:
    julia> println(complement("HA"))

    ["D", "T"]
"""
function complement(oligo)
    base_comps = Dict("A" => "T", "C" => "G", "G" => "C", "T" => "A")

    codes = ""
    for code in str_to_vector(oligo)
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

    ["HDN", "NHD"]
"""
function expand_asymmetric(oligo)
    return [oligo, vector_to_str(reverse_complement(oligo))]
end

"""
degeneracy(oligo)

    Description:
    Return the log of the number of sequences in a degenerate oligo.

    Example:
    julia> println(degeneracy("NNN"))

    4.1588830833596715
"""
function degeneracy(oligo)
    value = 0
    for code in str_to_vector(oligo)
        value += log(length(IUB_CODES[code]))
    end

    return value
end

"""
cutfree()

    Description:
    Return a site-free random of maximum diversity (of nucleotides). This single degenerate sequence can be used to synthesize a pool of oligos with mixed bases.

    Example:
    julia> cutfree()

    NNNNNWDNHDNHDNHDNDNN
    AAAAAAAAAAAAAAAAAAAA
    TTTTTTTTTTTTTTTTTTTT
    CCCCC--CC-CC-CC-C-CC
    GGGGG-GG-GG-GG-GGGGG
    24.731283462223622
    "NNNNNWDNHDNHDNHDNDNN"

    julia> cutfree(starting_oligo = "GNGNNNYBDKVDNGCTNNNNN", restriction_sites = ["GGTCTC", "GGCCGG"], min_blocks = 1, increase_diversity = true)

    GNGNDNYKWKVDNGCTNNNNN
    -A-AAA--A-AAA---AAAAA
    -T-TTTTTTT-TT--TTTTTT
    -C-C-CC---C-C-C-CCCCC
    GGGGGG-G-GGGGG--GGGGG
    18.545074838323124
    "GNGNDNYKWKVDNGCTNNNNN"
"""
function cutfree(;
            starting_oligo = "NNNNNNNNNNNNNNNNNNNN",
            restriction_sites = ["GGTCTC", "GGCCGG"],
            min_blocks = 1,
            increase_diversity = true,
            )

    m = length(starting_oligo)

    sites = []
    for i in 1:size(restriction_sites)[1]
        push!(sites, expand_asymmetric(restriction_sites[i])[1])
        push!(sites, expand_asymmetric(restriction_sites[i])[2])
    end

    A = NamedArray(zeros(15, m), (["A", "C", "G", "T", "R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N"],
        1:m), ("IUB_CODES", "Position"))
    for i in 1:m
        subs = subcodes(string(starting_oligo[i]))
        for b in subs
            A[b, i] = 1
        end
    end

    B = []
    for rs in sites
        B_i = NamedArray(zeros(15, length(rs)), (["A_blocked", "C_blocked", "G_blocked", "T_blocked", "R_blocked", "Y_blocked",
            "M_blocked", "K_blocked", "S_blocked", "W_blocked", "H_blocked", "B_blocked", "V_blocked", 
            "D_blocked", "N_blocked"], 1:length(rs)), ("IUB_CODES", "Position"))
        for j in 1:length(rs)
            blocked = get_blocking_codes(string(rs[j]), true)
            for b in blocked
                B_i[string(b, "_blocked"), j] = 1
            end
        end

        push!(B, B_i)
    end
    
    model = Model(Gurobi.Optimizer) # Gurobi can be changed to any open source solver you wish to use

    C = zeros(15, m) # Degeneracy matrix
    for i in 1:15, j in 1:m
        C[i, j] = log(length(IUB_CODES[names(A,1)[i]])) # Assign degeneracy values to matrix
    end

    @variable(model, output[1:15, 1:m], Bin) # Create output variable
    @constraint(model, sum(output[j, i] for i in 1:m, j in findall(x->x==0, A[:, i])) .== 0) # Enforce selected codes exist within subcodes(starting_oligo)
    @constraint(model, sum(output[i, 1:m] for i=1:15) .<= 1) # Enforce one code selected for each position
    
    for rs in sites, i in 1:length(sites), j in 1:m-length(rs)+1
        @constraint(model, sum(output[:, j:j+length(rs)-1] .* B[i]) .>= min_blocks) # Enforce presence of blocking codes
    end
    
    @objective(model, Max, sum(C .* output)) # Maximize degeneracy   
    optimize!(model)

    # Increase diversity of solution
    if (increase_diversity == true)
        # Run CutFree with a constraint addeed requiring the solution to be equal to the optimal objective value
        @constraint(model, sum(C .* output) == objective_value(model))

        # Change objective function to the sum of all variables multiplied by random coefficients (randomize)
        D = rand(15, m)
        @objective(model, Max, sum(D .* output))
        optimize!(model)
    end

    # Find oligo from output matrix, return empty if no possible output
    oligo = ""
    for i in 1:m, j in 1:15
        try value.(output[j, i])
            if value.(output[j, i]) == 1
                oligo *= names(A, 1)[j]
            end
        catch
            return "No possible output with given arguments."
        end
    end

    print_oligo_block(oligo)
    println(degeneracy(oligo))

    return oligo
end

cutfree()